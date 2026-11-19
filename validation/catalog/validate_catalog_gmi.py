"""
Validates the GMI catalog entry against the GPM Level 1C GMI product, using the
exact values of the entry through public TAT-C methods. The test of swath width
follows `ValidateConicalGMI`, on whose day of data the entry's angles were fitted;
the tests of passes, epochs, and swath area are this script's own. Running on that
day checks the shared test against the notebook, while running on a day at least 30
days from it tests the fitted angles on data held out from the fit. The day must be
one on which GPM flew forward, as it does for about 40 days at a time, because the
configuration represents forward flight. Element sets for days not registered in
`elements/element_sets.yml` are taken from the file of the history of GPM's element
sets in `elements/` (or one given with `--elements`), of which every set within 1.5
days of the middle of the day is used.

Needs the `validation` extra and a NASA Earthdata login:

    python validation/catalog/validate_catalog_gmi.py <data directory> 2026-10-03 --fit
    python validation/catalog/validate_catalog_gmi.py <data directory> 2026-10-03
    python validation/catalog/validate_catalog_gmi.py <data directory> 2026-07-20

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).

With `--fit`, it measures the angles of the entry's configurations instead, on the
day of the fit, as `ValidateConicalGMI` did: for each feed, the cone angle of the
lines of sight from geodetic nadir and the center and half width of the scan sector
in azimuth from the inertial velocity (medians over every 20th scan with contiguous
neighbors). It prints them, without writing.

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metric compares sizes with footprints published
for two altitudes, which are independent of the fit of the angles. The spacing of
samples on the day is reported, as the standard does not gate it. To run only
these tests, use `--configuration sample_size`.

With `--write`, it writes the entry's validation record and the reduced fixture
from which unit tests recompute the results. A record of the day of the fit has
no holdout, which the entry's validation reports as a gap.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
from datetime import timedelta

import numpy as np
import pandas as pd

from tatc.analysis import collect_observations
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit, Point
from conical import (
    GATED,
    Product,
    Reference,
    SampleSizes,
    load_day,
    orbit_period,
    pick_cone_passes,
    score_swath,
    validate_sample_sizes,
)
from validation_io import fit_cone, parse_run, write_results

PRODUCT = Product(
    short_name="GPM_1CGPMGMI",
    cache_prefix="gmi_l1c",
    groups=("S1", "S2"),
    scan_period=1.875,
    contiguous_gap=3,
    pass_radius=40,
    orientation=True,
)
FEEDS = PRODUCT.groups
FIXTURE_PASSES = 40
# published footprints: instrument, altitude (km), along and across track (km), from
# the GMI L1B ATBD, Table 1-1 (407 km) and TC PRIMED after the orbit boost (435 km);
# an instrument of each feed, which share their spacing
SAMPLE_SIZES = SampleSizes(
    key="sample_size",
    footprints=(
        ("GMI 10.65 GHz", 407, 32.1, 19.4),
        ("GMI 18.7 GHz", 407, 18.1, 10.9),
        ("GMI 23.8 GHz", 407, 16.0, 9.7),
        ("GMI 36.64 GHz", 407, 15.6, 9.4),
        ("GMI 89 GHz", 407, 7.2, 4.4),
        ("GMI 166 GHz", 407, 6.3, 4.1),
        ("GMI 183.31 GHz", 407, 5.8, 3.8),
        ("GMI 10.65 GHz", 435, 35.0, 21.0),
        ("GMI 36.64 GHz", 435, 16.7, 10.0),
        ("GMI 89 GHz", 435, 7.7, 4.6),
        ("GMI 166 GHz", 435, 6.7, 4.3),
    ),
    inclination=65.0,
    spacing_instruments={"S1": "GMI 89 GHz", "S2": "GMI 166 GHz"},
    reference="Footprints of the GMI L1B ATBD and TC PRIMED; spacing of "
    "GPM_1CGPMGMI version 08",
)


def load_reference(run) -> Reference:
    """Loads the day of a run, which must be one on which GPM flew forward."""
    data = load_day(PRODUCT, run.day, run.directory)
    orientations = np.unique(data["orientation"])
    if not np.array_equal(orientations, [0]):
        raise ValueError(
            f"GPM did not fly forward all of {run.day:%Y-%m-%d}: {orientations}; "
            "the configuration represents forward flight"
        )
    return Reference(PRODUCT, run.day, data)


def recompute_passes(satellite, day, passes, feeds) -> dict[str, float]:
    """
    Recomputes the observed passes of the fixture, each by its feed, as the unit
    tests do: the percent found, and the 95th percentile of the earliest epoch's
    difference (not the nearest, as `validation_io.recompute_passes` takes).
    """
    differences = []
    for (lat, lon, seconds), index in zip(passes, feeds):
        time = day + timedelta(seconds=float(seconds))
        found = collect_observations(
            Point(latitude=float(lat), longitude=float(lon)),
            satellite,
            time - timedelta(minutes=10),
            time + timedelta(minutes=10),
            int(index),
        )
        differences.append(
            min((e - pd.Timestamp(time)).total_seconds() for e in found.epoch)
            if not found.empty
            else np.nan
        )
    differences = np.array(differences, dtype=float)
    return {
        "passes_found": 100 * float(np.isfinite(differences).mean()),
        "epoch_p95": summarize(differences)["p95"],
    }


def fit(reference: Reference) -> None:
    """
    Prints the cone and scan sector of each feed, measured on every 20th scan with
    contiguous neighbors, as `ValidateConicalGMI` measured them.
    """
    scans = np.nonzero(reference.contiguous)[0][::20]
    print(f"{len(reference.time)} scans, {len(scans)} measured")
    for feed in FEEDS:
        lon, lat = reference.foot[feed]
        values = fit_cone(
            reference.position[scans],
            reference.velocity[scans],
            lon[scans],
            lat[scans],
        )
        print(
            f"{feed}: cone_angle {values['cone_angle']:.3f}, "
            f"scan_center_azimuth {values['scan_center_azimuth']:.3f}, "
            f"scan_half_width {values['scan_half_width']:.3f} "
            f"(azimuth step {values['scan_sample_spacing']:.4f}; center from the "
            f"Earth-fixed velocity {values['scan_center_azimuth_earth_fixed']:.3f})"
        )


def validate_feeds(run, orbit, tle, reference: Reference) -> None:
    """
    Scores each feed of the swath configuration on its samples, and writes the
    results if the run writes.
    """
    entry = run.entry
    specifications = entry.configurations["feeds"].instruments

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, time_step=time_step, platform=run.platform)

    # a lattice with a typical spacing of 500 km within the latitudes GMI observes
    points = generate_points_fibonacci_lattice(500e3)
    points = points[np.abs(points.geometry.y) < 72].reset_index(drop=True)
    # one orbit from 12:00 UTC
    scans = np.nonzero(
        (reference.time >= 12 * 3600)
        & (reference.time < 12 * 3600 + orbit_period(orbit))
    )[0]
    results, observed = {}, {}
    for index, feed in enumerate(FEEDS):
        cone_angle = specifications[index].parameters["cone_angle"]
        results[feed], observed[feed] = score_swath(
            PRODUCT, reference, satellite_for, points, feed, index, cone_angle, scans
        )
        print(feed, json.dumps(results[feed]))
    if not run.write:
        return
    rng = np.random.default_rng(1)
    passes = np.concatenate(
        [
            pick_cone_passes(points, observed[feed], FIXTURE_PASSES, rng)
            for feed in FEEDS
        ]
    )
    feeds = np.concatenate(
        [np.full(FIXTURE_PASSES, index) for index in range(len(FEEDS))]
    )
    write_results(
        run,
        "feeds",
        reference="GPM_1CGPMGMI version 08",
        metrics=[
            {"name": name, "value": float(results[feed][name]), "scope": feed}
            for name in GATED
            for feed in FEEDS
        ],
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "passes": passes,
            "feeds": feeds,
        },
        fixture_metrics=recompute_passes(
            satellite_for(PRODUCT.scan_period), run.day, passes, feeds
        ),
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("gmi", __doc__)
    reference = load_reference(run)
    if run.fit:
        fit(reference)
        return
    tle = run.element_set(within=timedelta(days=1.5))
    print(f"{len(tle) // 2} element sets")
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute"
    )
    if SAMPLE_SIZES.key in run.configurations:
        validate_sample_sizes(SAMPLE_SIZES, PRODUCT, run, orbit, tle, reference)
    if "feeds" in run.configurations:
        validate_feeds(run, orbit, tle, reference)


if __name__ == "__main__":
    main()
