"""
Validates the DPR catalog entry against the GPM Level 2 DPR product, using the exact
values of the entry through public TAT-C methods. The entry's scan angles are fitted
to the positions of the beams on one day (run with `--fit`): running on that day
checks the fit, while running on a day at least 30 days from it tests the fitted
angles on data held out from the fit. Element sets for other days are read from a
file of the history of GPM's element sets.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python docs/validation/validate_catalog_dpr.py <data directory> 2026-10-03 --fit
    python docs/validation/validate_catalog_dpr.py <data directory> 2026-10-03
    python docs/validation/validate_catalog_dpr.py <data directory> 2026-07-20
    python docs/validation/validate_catalog_dpr.py <data directory> 2026-11-10

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
NORAD_ID = 39574
FIT_DAY = datetime(2026, 10, 3, tzinfo=timezone.utc)
# the element set of `ValidateConicalGMI` for the day of the fit, and a public element
# set from CelesTrak, with an epoch after the day of data, for a later held-out day
ELEMENTS = {
    "2026-10-03": (
        "1 39574U 14009C   26277.54227627  .00004528  00000+0  10004-3 0  9997",
        "2 39574  64.9701 301.6724 0010793 268.6572  91.3332 15.45217995714704",
    ),
    "2026-11-10": (
        "1 39574U 14009C   26315.61017546  .00006513  00000+0  14350-3 0  9999",
        "2 39574  64.9689 174.1431 0011469 256.2485 103.7379 15.44776426720584",
    ),
}
ELEMENT_SET_FILE = REPOSITORY / "docs" / "validation" / "gpm_2026-07-15.3le"
PRODUCT = Product(
    short_name="GPM_2ADPR",
    version="08",
    group="FS",
    spacecraft=("navigation/scLat", "navigation/scLon", "navigation/scAlt"),
    quality="scanStatus/dataQuality",
)
# the scan period is measured on the day of the fit
SCAN_PERIOD = 0.666
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the TC PRIMED documentation,
# before and after GPM's orbit boost of November 2023; the edge of the scan is at
# its nominal 17 deg
PUBLISHED_FOOTPRINTS = (
    ("DPR", 400, 65, 0, 5.04, 5.04),
    ("DPR", 400, 65, 17.0, 5.04, 5.57),
    ("DPR", 435, 65, 0, 5.48, 5.48),
    ("DPR", 435, 65, 17.0, 5.48, 6.05),
)
SETUP = Setup(
    entry_id="dpr",
    platform="gpm",
    swath="swath",
    beams="beams",
    size_instrument="DPR",
    published_footprints=PUBLISHED_FOOTPRINTS,
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "dpr.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "dpr.npz",
    source="docs/validation/validate_catalog_dpr.py",
    reference="GPM_2ADPR version 08 (beam positions of the full swath)",
    published_reference="Footprints of the TC PRIMED documentation; beam positions "
    "of GPM_2ADPR version 08",
    fit_day=FIT_DAY,
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
    tle = list(
        ELEMENTS[arguments.day]
        if arguments.day in ELEMENTS
        else read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    )
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("dpr")
    validate_beams(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
