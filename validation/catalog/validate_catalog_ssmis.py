"""
Validates the SSMIS catalog entry against the GPM Level 1C SSMIS product of DMSP F18,
using the exact values of the entry through public TAT-C methods. The entry's scan is
fitted to the positions of the samples of the low-frequency group (S1, 19 and 22 GHz)
on one day (run with `--fit`), from the spacecraft positions the product records, so
the fit needs no element set; running on a day at least 30 days from it tests the
fitted scan on data held out from the fit. The tests follow those of AMSR2:

* A point is observed on a pass if it lies in a cell of the lattice of S1 samples, at
  the time it crosses the cone, from the recorded spacecraft positions.
* The swath width is compared with the distance between the first and last samples.
* The swath area over a span of the day (from noon) is compared with the union of
  the cells.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python validation/catalog/validate_catalog_ssmis.py <data directory> 2026-10-04 --fit
    python validation/catalog/validate_catalog_ssmis.py <data directory> 2026-11-03

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
Element sets for days not registered in `elements/element_sets.yml` can be given in
a file of DMSP F18's element sets (`--elements`).

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import numpy as np

from tatc.schemas import GeneralPerturbationsOrbit
from conical import fit_cone_scan, validate_cone_swath
from scanner import Product, Reference, load_day
from validation_io import parse_run

PRODUCT = Product(
    short_name="GPM_1CF18SSMIS",
    version="08",
    group="S1",
    spacecraft=("SCstatus/SClatitude", "SCstatus/SClongitude", "SCstatus/SCaltitude"),
)
SCAN_PERIOD = 1.899


def main() -> None:
    """Validates the entry on a day of data, or fits its scan."""
    run = parse_run("ssmis", __doc__)
    reference = Reference(
        run.day, load_day(PRODUCT, run.day, run.directory), SCAN_PERIOD
    )
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans of {reference.beams} "
        f"samples, {len(reference.gaps)} gaps over a minute, median altitude "
        f"{np.median(reference.altitude) / 1e3:.1f} km",
        flush=True,
    )
    if run.fit:
        fit_cone_scan(reference, range(-90, 90, 30))
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    validate_cone_swath(
        run,
        {"S1": reference},
        GeneralPerturbationsOrbit.from_tle(tle),
        tle,
        key="swath",
        reference="GPM_1CF18SSMIS version 08 (positions of the S1 samples)",
        scan_period=SCAN_PERIOD,
        pass_radius=60,
        point_spacing=300e3,
        fixture_passes=80,
    )


if __name__ == "__main__":
    main()
