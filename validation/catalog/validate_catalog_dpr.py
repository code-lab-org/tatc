"""
Validates the DPR catalog entry against the GPM Level 2 DPR product, using the exact
values of the entry through public TAT-C methods. The entry's scan angles are fitted
to the positions of the beams on one day (run with `--fit`): running on that day
checks the fit, while running on a day at least 30 days from it tests the fitted
angles on data held out from the fit. GPM's element sets are those registered in
`elements/element_sets.yml`, or else read from the file of the history of GPM's
element sets in `elements/` (or from a file given with `--elements`).

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python validation/catalog/validate_catalog_dpr.py <data directory> 2026-10-03 --fit
    python validation/catalog/validate_catalog_dpr.py <data directory> 2026-10-03
    python validation/catalog/validate_catalog_dpr.py <data directory> 2026-07-20
    python validation/catalog/validate_catalog_dpr.py <data directory> 2026-11-10

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
    swath="swath",
    beams="beams",
    size_instrument="DPR",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "dpr.npz",
    reference="GPM_2ADPR version 08 (beam positions of the full swath)",
    published_reference="Footprints of the TC PRIMED documentation; beam positions "
    "of GPM_2ADPR version 08",
    time_steps=(10, 30, 120),
    pass_radius=10,
    point_spacing=100e3,
    pixel_stride=100,
    spacing_stride=150,
)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("dpr", __doc__)
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
