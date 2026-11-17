"""
Validates the TRMM PR catalog entry against the GPM-format Level 2 PR product of TRMM,
using the exact values of the entry through public TAT-C methods. The entry's scan
angles are fitted to the positions of the beams on one day (run with `--fit`), from
the spacecraft positions the product records, so the fit needs no element set:
running on that day checks the fit, while running on a day at least 30 days from it
tests the fitted angles on data held out from the fit. Both days are in TRMM's
forward orientation (the products' `SCorientation` of 0) throughout. Element sets
are read from a file of the history of TRMM's element sets.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python validation/catalog/validate_catalog_trmm_pr.py <data directory> 2013-06-01 --fit
    python validation/catalog/validate_catalog_trmm_pr.py <data directory> 2013-06-01
    python validation/catalog/validate_catalog_trmm_pr.py <data directory> 2013-07-10

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
from validation_io import (
    read_element_set_file,
    FIXTURES,
    source_of,
    ELEMENT_SET_FOLDER,
)

NORAD_ID = 25063
ELEMENT_SET_FILE = ELEMENT_SET_FOLDER / "trmm_2013-05-20.3le"
PRODUCT = Product(
    short_name="GPM_2APR",
    version="07",
    group="FS",
    spacecraft=("navigation/scLat", "navigation/scLon", "navigation/scAlt"),
    quality="scanStatus/dataQuality",
)
# the scan period is measured on the day of the fit
SCAN_PERIOD = 0.5997
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the TC PRIMED documentation,
# after TRMM's orbit boost of August 2001; the edge of the scan is at its nominal
# 17 deg
PUBLISHED_FOOTPRINTS = (
    ("PR", 403, 35, 0, 5.0, 5.0),
    ("PR", 403, 35, 17.0, 5.0, 5.0),
)
SETUP = Setup(
    entry_id="trmm_pr",
    platform="trmm",
    swath="swath",
    beams="beams",
    size_instrument="PR",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "trmm_pr.npz",
    source=source_of(__file__),
    reference="GPM_2APR version 07 (beam positions of the full swath)",
    published_reference="Footprints of the TC PRIMED documentation; beam positions "
    "of GPM_2APR version 07",
    time_steps=(10, 30, 120),
    pass_radius=10,
    point_spacing=100e3,
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
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
        return
    tle = read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("trmm_pr")
    validate_beams(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
