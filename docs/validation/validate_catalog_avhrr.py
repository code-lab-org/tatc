"""
Validates the AVHRR catalog entry against the full-resolution AVHRR swath of Metop-C
in NOAA STAR's ACSPO sea surface temperature product (GHRSST Level 2P), which gives
the position of every sample, using the exact values of the entry through public
TAT-C methods. The entry's scan angles are fitted to the positions of the samples on
one day (run with `--fit`): running on that day checks the fit, while running on a
day at least 30 days from it tests the fitted angles on data held out from the fit.
Element sets for other days are read from a file of the history of Metop's element
sets.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
positions and times are read, with range requests (about 2 GB a day):

    python docs/validation/validate_catalog_avhrr.py <data directory> 2026-11-05 --fit
    python docs/validation/validate_catalog_avhrr.py <data directory> 2026-11-05
    python docs/validation/validate_catalog_avhrr.py <data directory> 2026-10-05

A day is reduced to one scan line a second and every 16th sample (and the last). The
product records no spacecraft position, which the fit takes from the element set.

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with the footprint published
by the NWP SAF, and the projected pixels with the samples' positions.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results. A record of the day of the fit has no
holdout, which the entry's validation reports as a gap.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
from datetime import datetime, timezone

import h5py
import numpy as np
import pandas as pd
from skyfield.api import wgs84

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import (
    Product,
    Reference,
    Setup,
    fit,
    load_day,
    validate_beams,
    validate_swath,
)
from validation_io import RemoteFile, read_element_set_file

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
NORAD_ID = 43689
FIT_DAY = datetime(2026, 11, 5, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-05": (
        "1 43689U 18087A   26312.29310542  .00000136  00000+0  82226-4 0  9997",
        "2 43689  98.7256   9.6406 0001524 126.9326 233.1991 14.21484375415305",
    ),
}
ELEMENT_SET_FILE = REPOSITORY / "docs" / "validation" / "metop_2026-10-05.3le"
GHRSST_EPOCH = pd.Timestamp("1981-01-01", tz="UTC")
SAMPLES = 2048
LINE_RATE = 6  # scan lines a second
LINE_STRIDE = LINE_RATE
COLUMNS = np.r_[0:SAMPLES:16, SAMPLES - 1]
# the middle sample, whose time is fitted
TIME_SAMPLE = SAMPLES // 2


def read_avhrr(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the positions and times of a granule's samples with range requests, and
    reduces them to one scan line a second and some samples. The product gives times
    in quarter seconds, so each line's time is that of a line fitted to the times by
    line, which is checked against them.
    """
    with h5py.File(RemoteFile(url, session, block_size=2**20)) as f:
        lat = f["lat"][:]
        lon = f["lon"][:]
        base = float(f["time"][0])
        dtime = f["sst_dtime"][0, :, TIME_SAMPLE].astype(float)
        fill = f["sst_dtime"].attrs["_FillValue"][0]
        scale = float(f["sst_dtime"].attrs["scale_factor"][0])
    rows = np.arange(len(dtime))
    valid = dtime != fill
    slope, intercept = np.polyfit(rows[valid], dtime[valid] * scale, 1)
    residual = np.abs(intercept + slope * rows[valid] - dtime[valid] * scale)
    if abs(slope * LINE_RATE - 1) > 1e-3 or residual.max() > 0.3:
        raise ValueError(
            f"{url}: lines do not follow a constant rate (slope {slope:.6f} s, "
            f"largest residual {residual.max():.3f} s)"
        )
    keep = rows[::LINE_STRIDE]
    seconds = (
        base
        - (pd.Timestamp(day) - GHRSST_EPOCH).total_seconds()
        + intercept
        + slope * keep
    )
    return {
        "time": seconds,
        "lat": lat[keep][:, COLUMNS],
        "lon": lon[keep][:, COLUMNS],
    }


PRODUCT = Product(
    short_name="AVHRRF_MC-STAR-L2P-v2.80",
    version="2.80",
    read=read_avhrr,
    host="https://archive.podaac",
)
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the NWP SAF's description of
# AVHRR/3, which gives a square 1.1 km at nadir and states no altitude, at Metop's
# nominal altitude
PUBLISHED_FOOTPRINTS = (("AVHRR", 817, 98.7, 0, 1.1, 1.1),)
SETUP = Setup(
    entry_id="avhrr",
    platform="metop-c",
    swath="swath",
    beams="pixels",
    size_instrument="AVHRR",
    published_footprints=PUBLISHED_FOOTPRINTS,
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "avhrr.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "avhrr.npz",
    source="docs/validation/validate_catalog_avhrr.py",
    reference="AVHRRF_MC-STAR-L2P version 2.80 (sample positions and times)",
    published_reference="Footprint published by the NWP SAF; sample positions of "
    "AVHRRF_MC-STAR-L2P version 2.80",
    fit_day=FIT_DAY,
    time_steps=(10, 30, 120),
    pass_radius=40,
    point_spacing=200e3,
    pixel_stride=100,
    spacing_stride=150,
)


def load_reference(day: datetime, directory: pathlib.Path, tle: list[str]):
    """Loads a reduced day, with the spacecraft's position from the element set."""
    data = load_day(PRODUCT, day, directory)
    times = list(
        (pd.Timestamp(day) + pd.to_timedelta(data["time"], unit="s")).to_pydatetime()
    )
    position = wgs84.geographic_position_of(
        GeneralPerturbationsOrbit.from_tle(tle).get_orbit_track(times)
    )
    data.update(
        sc_lat=position.latitude.degrees,
        sc_lon=position.longitude.degrees,
        sc_alt=position.elevation.m,
        pixel=COLUMNS,
        pixels=SAMPLES,
        line_stride=LINE_STRIDE,
    )
    return Reference(day, data, LINE_STRIDE / LINE_RATE)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--sample-sizes-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    tle = list(
        ELEMENTS[arguments.day]
        if arguments.day in ELEMENTS
        else read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    )
    print(f"element set: {tle[0][:32]}", flush=True)
    reference = load_reference(day, arguments.directory, tle)
    print(
        f"{arguments.day}: {len(reference.time)} scans from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("earth_fixed", "inertial"):
            fit(reference, frame)
        return
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("avhrr")
    validate_beams(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
