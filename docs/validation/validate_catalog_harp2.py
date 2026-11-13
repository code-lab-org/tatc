"""
Validates the HARP2 catalog entry against the PACE HARP2 Level 1B product, using the
exact values of the entry through public TAT-C methods. HARP2 is a wide-angle
imaging polarimeter whose frames give 90 views of each place at angles along track;
the entry models the near-nadir view (index 39, at 664 nm, 0.6 deg along track),
whose positions and times are read for each line. Its angles are fitted to the
positions of the samples on one day (run with `--fit`), from the spacecraft
positions the product records, so the fit needs no element set. Running on a day at
least 30 days from it tests the fitted values on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
near-nadir view's positions, times, and spacecraft positions are read, with range
requests:

    python docs/validation/validate_catalog_harp2.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_harp2.py <data directory> 2026-11-07

A day is reduced to every fourth line of the view (about 1.8 s) and every eighth
sample (and the last). The product covers the sunlit part of each orbit, whose dark
part is a gap.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
from datetime import datetime, timezone

import h5py
import numpy as np
import pandas as pd

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import (
    Product,
    Reference,
    Setup,
    fit,
    fit_yaw,
    load_day,
    validate_swath,
)
from validation_io import RemoteFile

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
PACE_ELEMENTS = (
    "1 58928U 24025A   26313.27483190  .00000501  00000+0  10091-3 0  9993",
    "2 58928  98.1040 242.4857 0002666 129.8006 230.3433 14.65644467147180",
)
ELEMENTS = {"2026-11-07": PACE_ELEMENTS, "2026-11-08": PACE_ELEMENTS}
VIEW = 39
# the samples with positions in the view's lines
FIRST_SAMPLE, LAST_SAMPLE = 3, 638
LINE_STRIDE, SAMPLE_STRIDE = 4, 8
COLUMNS = np.r_[FIRST_SAMPLE : LAST_SAMPLE + 1 : SAMPLE_STRIDE, LAST_SAMPLE]


def read_harp2(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the near-nadir view's lines of a granule with range requests: their times
    (seconds from the start of the day), the positions of their samples, and the
    spacecraft's position, keeping the lines with positions.
    """
    with h5py.File(RemoteFile(url, session, block_size=2**20)) as f:
        start = pd.Timestamp(f.attrs["time_coverage_start"].decode())
        lat = f["geolocation_data/latitude"][VIEW]
        lon = f["geolocation_data/longitude"][VIEW]
        seconds = f["navigation_data/seconds_of_day"][VIEW].astype(float)
        sc_lat = f["navigation_data/orbit_latitude"][VIEW]
        sc_lon = f["navigation_data/orbit_longitude"][VIEW]
        sc_alt = f["navigation_data/orbit_altitude"][VIEW]
    # the times are seconds of the day on which the granule starts
    seconds = seconds + (start.normalize() - pd.Timestamp(day)).total_seconds()
    seconds = np.where(seconds < seconds[0] - 43200, seconds + 86400, seconds)
    lines = np.flatnonzero(np.all(np.abs(lat[:, COLUMNS]) <= 90, axis=1))
    lines = lines[::LINE_STRIDE]
    return {
        "time": seconds[lines],
        "lat": lat[lines][:, COLUMNS],
        "lon": lon[lines][:, COLUMNS],
        "sc_lat": sc_lat[lines],
        "sc_lon": sc_lon[lines],
        "sc_alt": sc_alt[lines],
    }


PRODUCT = Product(
    short_name="PACE_HARP2_L1B_SCI",
    version="4",
    read=read_harp2,
    host="https://obdaac-tea",
)
SETUP = Setup(
    entry_id="harp2",
    platform="pace",
    swath="swath",
    beams="",
    size_instrument="HARP2",
    published_footprints=(),
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "harp2.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "harp2.npz",
    source="docs/validation/validate_catalog_harp2.py",
    reference="PACE_HARP2_L1B_SCI version 4 (positions and times of the near-nadir "
    "view)",
    published_reference="",
    fit_day=FIT_DAY,
    time_steps=(10, 30, 120),
    pass_radius=20,
    point_spacing=200e3,
)


def load_reference(day: datetime, directory: pathlib.Path) -> Reference:
    """Loads a reduced day of the near-nadir view."""
    data = load_day(PRODUCT, day, directory)
    data.update(pixel=COLUMNS - FIRST_SAMPLE, pixels=LAST_SAMPLE - FIRST_SAMPLE + 1)
    # each line of the view is about 0.44 s after the last
    period = float(np.median(np.diff(data["time"])))
    data["line_stride"] = LINE_STRIDE
    return Reference(day, data, period)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    reference = load_reference(day, arguments.directory)
    print(
        f"{arguments.day}: {len(reference.time)} lines from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("earth_fixed", "inertial"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate_swath(
        SETUP, day, tle, get_entry("harp2"), orbit, reference, arguments.write
    )
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
