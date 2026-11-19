"""
Validates the MHS catalog entry against the GPM Level 1C MHS product of Metop-C,
using the exact values of the entry through public TAT-C methods. The entry's scan
angles are fitted to the positions of the beams on one day (run with `--fit`):
running on that day checks the fit, while running on a day at least 30 days from it
tests the fitted angles on data held out from the fit. Element sets are those
registered in `elements/element_sets.yml`, or else read from the file of the history
of Metop's element sets in `elements/` (or from a file given with `--elements`).

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python validation/catalog/validate_catalog_mhs.py <data directory> 2026-11-05 --fit
    python validation/catalog/validate_catalog_mhs.py <data directory> 2026-11-05
    python validation/catalog/validate_catalog_mhs.py <data directory> 2026-10-05

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metrics compare sizes with footprints published in
the TC PRIMED documentation, and the projected pixels with the beams' positions.
`--configuration beams` (or `swath`) tests one configuration.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results. A record of the day of the fit has no
holdout, which the entry's validation reports as a gap.

With `--platform metop-b`, it tests the same values on a sibling satellite, and
adds the results to the record as a transfer, which does not count toward the tier.

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

SPACECRAFT = ("SCstatus/SClatitude", "SCstatus/SClongitude", "SCstatus/SCaltitude")
# the product of each platform: Metop-C's, on which the entry is fitted and tested,
# and a sibling's, on which results are reported as a transfer
PRODUCTS = {
    "metop-c": Product(
        short_name="GPM_1CMETOPCMHS", version="08", group="S1", spacecraft=SPACECRAFT
    ),
    "metop-b": Product(
        short_name="GPM_1CMETOPBMHS", version="08", group="S1", spacecraft=SPACECRAFT
    ),
}
SCAN_PERIOD = 8 / 3
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the TC PRIMED documentation,
# which states no altitude, at Metop's nominal altitude and the outermost beam
PUBLISHED_FOOTPRINTS = (
    ("MHS", 817, 98.7, 0, 15.9, 15.9),
    ("MHS", 817, 98.7, 49.444, 27.1, 52.8),
)
SETUP = Setup(
    swath="swath",
    beams="beams",
    size_instrument="MHS",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "mhs.npz",
    reference="GPM_1CMETOPCMHS version 08 (beam positions)",
    published_reference="Footprints of the TC PRIMED documentation; beam positions "
    "of GPM_1CMETOPCMHS version 08",
    time_steps=(SCAN_PERIOD, 32, 128),
    pass_radius=60,
    point_spacing=500e3,
    pixel_stride=100,
    spacing_stride=150,
)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("mhs", __doc__)
    reference = Reference(
        run.day, load_day(PRODUCTS[run.platform], run.day, run.directory), SCAN_PERIOD
    )
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        for frame in ("earth_fixed", "inertial"):
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
