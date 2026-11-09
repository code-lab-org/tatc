"""
Validates the AIRS catalog entry against the AIRS Level 1B infrared product of Aqua,
using the exact values of the entry through public TAT-C methods. The entry's scan
angles are fitted to the positions of the footprints on one day (run with `--fit`),
which needs no element set: running on a day at least 30 days from it tests the
entry on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python docs/validation/validate_catalog_airs.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_airs.py <data directory> 2026-11-05
    python docs/validation/validate_catalog_airs.py <data directory> 2026-11-04

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with footprints published by
the AIRS project, and the projected pixels with the footprints' positions.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
from datetime import datetime, timezone

import h5py
import numpy as np

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import (
    Product,
    Reference,
    Setup,
    fit,
    load_day,
    tai93_to_seconds,
    validate_beams,
    validate_swath,
)
from validation_io import RemoteFile

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
# the public element set of `validate_catalog_modis.py` for Aqua, with an epoch after
# the days of data
AQUA_ELEMENTS = (
    "1 27424U 02022A   26311.28865478  .00000831  00000+0  17465-3 0  9993",
    "2 27424  98.4481 283.5818 0000983 104.9545   2.1913 14.62267061304185",
)
ELEMENTS = {"2026-11-05": AQUA_ELEMENTS, "2026-11-04": AQUA_ELEMENTS}
SCAN_PERIOD = 8 / 3


def read_airs(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the fields of a granule with range requests: the middle time of each scan's
    Earth view and each footprint's time (seconds from the start of the day), the
    footprints' positions on the geoid, the subsatellite point and altitude, and
    whether any footprint of a scan is not usable for science.
    """
    with h5py.File(RemoteFile(url, session)) as f:
        return {
            "time": tai93_to_seconds(f["scan_mid_time"][:], day),
            "beam_time": tai93_to_seconds(f["obs_time_tai93"][:], day),
            "lat": f["lat_geoid"][:],
            "lon": f["lon_geoid"][:],
            "sc_lat": f["subsat_lat"][:],
            "sc_lon": f["subsat_lon"][:],
            "sc_alt": f["sat_alt"][:],
            "quality": np.max(f["instrument_state"][:], axis=1),
        }


PRODUCT = Product(short_name="AIRSAQIRL1B", version="8.0", read=read_airs)
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the AIRS project at Aqua's
# nominal altitude: 13.5 km at nadir, and 41 x 21.4 km at the scan's extremes, the
# larger of which is across track, at the outermost footprint's angle
PUBLISHED_FOOTPRINTS = (
    ("AIRS", 705, 98.2, 0, 13.5, 13.5),
    ("AIRS", 705, 98.2, 48.95, 21.4, 41.0),
)
SETUP = Setup(
    entry_id="airs",
    platform="aqua",
    swath="swath",
    beams="footprints",
    size_instrument="AIRS",
    published_footprints=PUBLISHED_FOOTPRINTS,
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "airs.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "airs.npz",
    source="docs/validation/validate_catalog_airs.py",
    reference="AIRSAQIRL1B version 8.0 (footprint positions and times)",
    published_reference="Footprints published by the AIRS project; footprint "
    "positions of AIRSAQIRL1B version 8.0",
    fit_day=FIT_DAY,
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=30,
    point_spacing=200e3,
    pixel_stride=100,
    spacing_stride=150,
)


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
    reference = Reference(day, load_day(PRODUCT, day, arguments.directory), SCAN_PERIOD)
    print(
        f"{arguments.day}: {len(reference.time)} scans from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
        offsets = reference.beam_seconds - reference.seconds[:, np.newaxis]
        print(f"footprint times from the scan's middle: {np.median(offsets, axis=0)}")
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("airs")
    validate_beams(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
