"""
Validates the MHS catalog entry against the GPM Level 1C MHS product of Metop-C,
using the exact values of the entry through public TAT-C methods. The entry's scan
angles are fitted to the positions of the beams on one day (run with `--fit`):
running on that day checks the fit, while running on a day at least 30 days from it
tests the fitted angles on data held out from the fit. Element sets for other days
are read from a file of the history of Metop's element sets.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_mhs.py <data directory> 2026-11-05 --fit
    python docs/validation/validate_catalog_mhs.py <data directory> 2026-11-05
    python docs/validation/validate_catalog_mhs.py <data directory> 2026-10-05

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with footprints published in
the TC PRIMED documentation, and the projected pixels with the beams' positions.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results. A record of the day of the fit has no
holdout, which the entry's validation reports as a gap.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
from datetime import datetime, timezone

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
from validation_io import read_element_set_file

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
PRODUCT = Product(
    short_name="GPM_1CMETOPCMHS",
    version="08",
    group="S1",
    spacecraft=("SCstatus/SClatitude", "SCstatus/SClongitude", "SCstatus/SCaltitude"),
)
SCAN_PERIOD = 8 / 3
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the TC PRIMED documentation,
# which states no altitude, at Metop's nominal altitude and the outermost beam
PUBLISHED_FOOTPRINTS = (
    ("MHS", 817, 98.7, 0, 15.9, 15.9),
    ("MHS", 817, 98.7, 49.444, 27.1, 52.8),
)
SETUP = Setup(
    entry_id="mhs",
    platform="metop-c",
    swath="swath",
    beams="beams",
    size_instrument="MHS",
    published_footprints=PUBLISHED_FOOTPRINTS,
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "mhs.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "mhs.npz",
    source="docs/validation/validate_catalog_mhs.py",
    reference="GPM_1CMETOPCMHS version 08 (beam positions)",
    published_reference="Footprints of the TC PRIMED documentation; beam positions "
    "of GPM_1CMETOPCMHS version 08",
    fit_day=FIT_DAY,
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=60,
    point_spacing=500e3,
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
        f"{arguments.day}: {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("earth_fixed", "inertial"):
            fit(reference, frame)
        return
    tle = list(
        ELEMENTS[arguments.day]
        if arguments.day in ELEMENTS
        else read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    )
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("mhs")
    validate_beams(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
