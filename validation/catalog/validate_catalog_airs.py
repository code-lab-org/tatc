"""
Validates the AIRS catalog entry against the AIRS Level 1B infrared product of Aqua,
using the exact values of the entry through public TAT-C methods. The entry's scan
angles are fitted to the positions of the footprints on one day (run with `--fit`),
which needs no element set: running on a day at least 30 days from it tests the
entry on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python validation/catalog/validate_catalog_airs.py <data directory> 2026-10-04 --fit
    python validation/catalog/validate_catalog_airs.py <data directory> 2026-11-05
    python validation/catalog/validate_catalog_airs.py <data directory> 2026-11-04

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
Aqua's element sets are those registered in `elements/element_sets.yml`, or those
of a file given with `--elements`.

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with footprints published by
the AIRS project, and the projected pixels with the footprints' positions.
`--configuration footprints` (or `swath`) tests one configuration.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime

import h5py
import numpy as np

from tatc.schemas import GeneralPerturbationsOrbit
from scanner import (
    Product,
    Reference,
    Setup,
    fit,
    load_day,
    tai93_to_seconds,
    validate_beams,
    validate_swath,
)
from validation_io import FIXTURES, RemoteFile, parse_run

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
    swath="swath",
    beams="footprints",
    size_instrument="AIRS",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "airs.npz",
    reference="AIRSAQIRL1B version 8.0 (footprint positions and times)",
    published_reference="Footprints published by the AIRS project; footprint "
    "positions of AIRSAQIRL1B version 8.0",
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=30,
    point_spacing=200e3,
    pixel_stride=100,
    spacing_stride=150,
)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("airs", __doc__)
    reference = Reference(
        run.day, load_day(PRODUCT, run.day, run.directory), SCAN_PERIOD
    )
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
        offsets = reference.beam_seconds - reference.seconds[:, np.newaxis]
        print(f"footprint times from the scan's middle: {np.median(offsets, axis=0)}")
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    if SETUP.beams in run.configurations:
        validate_beams(SETUP, run, reference, orbit, tle)
    if SETUP.swath in run.configurations:
        validate_swath(SETUP, run, reference, orbit, tle)


if __name__ == "__main__":
    main()
