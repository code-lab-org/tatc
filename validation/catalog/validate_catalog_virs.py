"""
Validates the VIRS catalog entry against the GPM-format Level 1B VIRS product of TRMM,
using the exact values of the entry through public TAT-C methods. The product gives
no spacecraft positions, so they are interpolated to its scans' times from those the
GPM-format Level 2 PR product records for the same spacecraft (every 0.6 s), and the
fit needs no element set. The entry's scan angles are fitted to the positions of the
samples on one day (run with `--fit`): running on that day checks the fit, while
running on a day at least 30 days from it tests the fitted angles on data held out
from the fit. Both days are in TRMM's forward orientation (the PR product's
`SCorientation` of 0) throughout. Element sets are read from a file of the history of
TRMM's element sets.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. The positions
of samples are read whole, and reduced to every third scan and some samples:

    python validation/catalog/validate_catalog_virs.py <data directory> 2013-06-01 --fit
    python validation/catalog/validate_catalog_virs.py <data directory> 2013-06-01
    python validation/catalog/validate_catalog_virs.py <data directory> 2013-07-10

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with the footprint published in
the TRMM Version 7 README of GES DISC, and the projected pixels with the samples'
positions.

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
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    ecef,
    read_element_set_file,
    RECORDS,
    FIXTURES,
    source_of,
)
from validate_catalog_trmm_pr import ELEMENT_SET_FILE, NORAD_ID
from validate_catalog_trmm_pr import PRODUCT as PR_PRODUCT

SAMPLES = 261
LINE_STRIDE = 3
COLUMNS = np.r_[0:SAMPLES:4, SAMPLES - 1]


def read_virs(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the times and positions of a granule's samples, and reduces them to every
    third scan and some samples. The spacecraft's position is added afterward.
    """
    with h5py.File(RemoteFile(url, session, block_size=2**20)) as f:
        g = f["Swath"]
        date = pd.to_datetime(
            pd.DataFrame(
                {
                    "year": g["ScanTime/Year"][:],
                    "month": g["ScanTime/Month"][:],
                    "day": g["ScanTime/DayOfMonth"][:],
                }
            ),
            utc=True,
        )
        seconds = (date - pd.Timestamp(day)).dt.total_seconds().to_numpy() + g[
            "ScanTime/SecondOfDay"
        ][:].astype(float)
        lat = g["Latitude"][:]
        lon = g["Longitude"][:]
    keep = np.arange(0, len(seconds), LINE_STRIDE)
    return {
        "time": seconds[keep],
        "lat": lat[keep][:, COLUMNS],
        "lon": lon[keep][:, COLUMNS],
    }


PRODUCT = Product(short_name="GPM_1BVIRS", version="08", read=read_virs)
# the scan period is measured on the day of the fit
SCAN_PERIOD = 0.30467
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the TRMM Version 7 README of GES
# DISC, which gives a horizontal resolution of 2.4 km after TRMM's orbit boost of
# August 2001, identical in both directions, at nadir
PUBLISHED_FOOTPRINTS = (("VIRS", 403, 35, 0, 2.4, 2.4),)
SETUP = Setup(
    entry_id="virs",
    platform="trmm",
    swath="swath",
    beams="pixels",
    size_instrument="VIRS",
    published_footprints=PUBLISHED_FOOTPRINTS,
    record_path=RECORDS / "virs.json",
    fixture_path=FIXTURES / "virs.npz",
    source=source_of(__file__),
    reference="GPM_1BVIRS version 08 (sample positions), with spacecraft positions of "
    "GPM_2APR version 07",
    published_reference="Footprint of the TRMM Version 7 README of GES DISC; sample "
    "positions of GPM_1BVIRS version 08",
    time_steps=(10, 30, 120),
    pass_radius=40,
    point_spacing=200e3,
    pixel_stride=100,
    spacing_stride=150,
)


def load_reference(day: datetime, directory: pathlib.Path) -> Reference:
    """
    Loads a reduced day, with the spacecraft's position interpolated, in Earth-fixed
    coordinates, from the positions the PR product records.
    """
    data = load_day(PRODUCT, day, directory)
    radar = load_day(PR_PRODUCT, day, directory)
    altitude = radar["sc_alt"].astype(float)
    altitude = altitude * 1e3 if np.median(altitude) < 1e4 else altitude
    position = ecef(
        radar["sc_lon"].astype(float), radar["sc_lat"].astype(float), altitude
    )
    order = np.argsort(radar["time"])
    within = (data["time"] >= radar["time"][order[0]]) & (
        data["time"] <= radar["time"][order[-1]]
    )
    data = {k: v[within] for k, v in data.items() if k != "granules"}
    interpolated = np.stack(
        [
            np.interp(data["time"], radar["time"][order], position[order, axis])
            for axis in range(3)
        ],
        axis=-1,
    )
    lon, lat, height = TO_GEODETIC.transform(*interpolated.T)
    data.update(
        sc_lat=lat,
        sc_lon=lon,
        sc_alt=height,
        pixel=COLUMNS,
        pixels=SAMPLES,
        line_stride=LINE_STRIDE,
    )
    reference = Reference(day, data, LINE_STRIDE * SCAN_PERIOD)
    # scans in gaps of the PR's positions longer than a few of its scans are dropped
    gap = np.diff(radar["time"][order]) > 5
    if np.any(gap):
        starts, ends = radar["time"][order][:-1][gap], radar["time"][order][1:][gap]
        print(
            f"{len(starts)} gaps of PR positions over 5 s "
            f"(longest {np.max(ends - starts):.1f} s)",
            flush=True,
        )
    return reference


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
    reference = load_reference(day, arguments.directory)
    print(
        f"{arguments.day}: {len(reference.time)} scans (every {LINE_STRIDE}), "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
        print(
            f"scan period of the product: {np.median(np.diff(reference.seconds)) / LINE_STRIDE:.6f} s"
        )
        return
    tle = read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("virs")
    validate_beams(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
