"""
Validates the TRMM PR catalog entry against the GPM-format Level 2 PR product of TRMM,
using the exact values of the entry through public TAT-C methods. The entry's scan
angles are fitted to the positions of the beams on one day (run with `--fit`), from
the spacecraft positions the product records, so the fit needs no element set:
running on that day checks the fit, while running on a day at least 30 days from it
tests the fitted angles on data held out from the fit. Both days are in TRMM's
forward orientation (the products' `SCorientation` of 0) throughout. Element sets
are read from the file of the history of TRMM's element sets in `elements/` (or
from a file given with `--elements`).

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python validation/catalog/validate_catalog_trmm_pr.py <data directory> 2013-06-01 --fit
    python validation/catalog/validate_catalog_trmm_pr.py <data directory> 2013-06-01
    python validation/catalog/validate_catalog_trmm_pr.py <data directory> 2013-07-10

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with footprints published in
the TC PRIMED documentation, and the projected pixels with the beams' positions.
`--configuration beams` (or `swath`) tests one configuration.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results. A record of the day of the fit has no
holdout, which the entry's validation reports as a gap.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from tatc.schemas import GeneralPerturbationsOrbit
from scanner import (
    Product,
    Reference,
    Setup,
    fit,
    load_day,
    validate_beams,
    validate_swath,
)
from validation_io import FIXTURES, parse_run

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
    swath="swath",
    beams="beams",
    size_instrument="PR",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "trmm_pr.npz",
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
    run = parse_run("trmm_pr", __doc__)
    reference = Reference(
        run.day, load_day(PRODUCT, run.day, run.directory), SCAN_PERIOD
    )
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
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
