"""
Validates the AMSR2 catalog entry against the GPM Level 1C AMSR2 product, using
the exact values of the entry through public TAT-C methods. The tests of passes,
epochs, and swath width follow `ValidateConicalAMSR2`, on whose day of data the
entry's angles were fitted; the test of swath area is this script's own. Running on
that day checks the shared tests against the notebook, while running on a day at
least 30 days from it tests the fitted angles on data held out from the fit. Element
sets for days not registered in `elements/element_sets.yml` can be given in a file of
the history of GCOM-W's element sets (`--elements`), of which every set within 1.5
days of the middle of the day is used.

Needs the `validation` extra and a NASA Earthdata login:

    python validation/catalog/validate_catalog_amsr2.py <data directory> 2026-09-26 --fit
    python validation/catalog/validate_catalog_amsr2.py <data directory> 2026-09-26
    python validation/catalog/validate_catalog_amsr2.py <data directory> 2026-11-01

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).

With `--fit`, it measures the angles of the entry's configurations instead, on the
day of the fit, as `ValidateConicalAMSR2` did: the cone angle of the lines of sight
from geodetic nadir and the center and half width of the scan sector in azimuth
from the inertial velocity (medians over every 20th scan with contiguous
neighbors), of the 10.65 GHz group (S1), with the azimuth step between its samples.
It prints them, without writing.

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metric compares sizes with footprints published
by JAXA, which are independent of the fit of the angles. The spacing of samples on
the day is reported, as the standard does not gate it. To run only these tests, use
`--configuration sample_size`.

With `--write`, it writes the entry's validation record and the reduced fixture
from which unit tests recompute the results. A record of the day of the fit has
no holdout, which the entry's validation reports as a gap.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
from datetime import timedelta

import numpy as np

from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from conical import (
    Product,
    Reference,
    SampleSizes,
    load_day,
    orbit_period,
    pick_cone_passes,
    score_swath,
    validate_sample_sizes,
)
from validation_io import fit_cone, parse_run, recompute_passes, write_results

# the group whose samples the scan was fitted to, and that of the 89 GHz channels
PRODUCT = Product(
    short_name="GPM_1CGCOMW1AMSR2",
    cache_prefix="amsr2_l1c",
    groups=("S1", "S5"),
    scan_period=1.5,
    contiguous_gap=2,
    pass_radius=60,
)
FIXTURE_PASSES = 80
# published footprints: instrument, altitude (km), along and across track (km), from
# JAXA's AMSR2 channel set (ground resolution at the scan center); an instrument of
# each group, whose spacing is measured
SAMPLE_SIZES = SampleSizes(
    key="sample_size",
    footprints=(
        ("AMSR2 6.925 and 7.3 GHz", 700, 62, 35),
        ("AMSR2 10.65 GHz", 700, 42, 24),
        ("AMSR2 18.7 GHz", 700, 22, 14),
        ("AMSR2 23.8 GHz", 700, 26, 15),
        ("AMSR2 36.5 GHz", 700, 12, 7),
        ("AMSR2 89.0 GHz", 700, 5, 3),
    ),
    inclination=98.2,
    spacing_instruments={"S1": "AMSR2 10.65 GHz", "S5": "AMSR2 89.0 GHz"},
    reference="Footprints of JAXA's AMSR2 channel set; spacing of "
    "GPM_1CGCOMW1AMSR2 version 08",
)


def get_orbit_scans(reference: Reference, period: float) -> np.ndarray:
    """
    Gets the scans of the first orbit from 12:00 UTC without a gap in the data, as
    an earlier version of `ValidateConicalAMSR2` used an orbit after the gap of its day.
    """
    time = reference.time
    start = 12 * 3600
    while start + period < time[-1]:
        scans = np.nonzero((time >= start) & (time < start + period))[0]
        if len(scans) < 2:
            start += 600
            continue
        breaks = np.nonzero(np.diff(time[scans]) >= 2)[0]
        if len(breaks):
            start = time[scans[breaks[-1] + 1]]
        elif time[scans[0]] - start >= 2 or time[scans[-1]] < start + period - 2:
            start += 600
        else:
            return scans
    raise ValueError("no orbit without a gap after 12:00 UTC")


def fit(reference: Reference) -> None:
    """
    Prints the cone and scan sector of the 10.65 GHz group (S1), measured on every
    20th scan with contiguous neighbors, as `ValidateConicalAMSR2` measured them.
    """
    scans = np.nonzero(reference.contiguous)[0][::20]
    lon, lat = reference.foot["S1"]
    values = fit_cone(
        reference.position[scans], reference.velocity[scans], lon[scans], lat[scans]
    )
    print(f"{len(reference.time)} scans, {len(scans)} measured")
    for key in ("cone_angle", "scan_center_azimuth", "scan_half_width"):
        print(f"{key}: {values[key]:.3f}")
    print(f"scan_sample_spacing: {values['scan_sample_spacing']:.4f}")
    print(
        f"azimuths of the first and last samples: {values['first_azimuth']:.2f}, "
        f"{values['last_azimuth']:.2f}; center from the Earth-fixed velocity: "
        f"{values['scan_center_azimuth_earth_fixed']:.3f}"
    )


def validate_swath(run, orbit, tle, reference: Reference) -> None:
    """Scores the swath configuration, and writes its results if the run writes."""
    entry = run.entry

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, time_step=time_step, platform=run.platform)

    # a lattice with a typical spacing of 500 km, as in ValidateConicalAMSR2
    points = generate_points_fibonacci_lattice(500e3)
    cone_angle = entry.configurations["swath"].instruments[0].parameters["cone_angle"]
    scans = get_orbit_scans(reference, orbit_period(orbit))
    first, last = reference.time[scans[0]], reference.time[scans[-1]]
    print(f"swath orbit from {first / 3600:.2f} to {last / 3600:.2f} h UTC")
    results, observed = score_swath(
        PRODUCT, reference, satellite_for, points, "S1", 0, cone_angle, scans
    )
    print(json.dumps(results))
    if not run.write:
        return
    passes = pick_cone_passes(points, observed, FIXTURE_PASSES)
    write_results(
        run,
        "swath",
        reference="GPM_1CGCOMW1AMSR2 version 08",
        metrics=[
            {"name": name, "value": float(results[name])}
            for name in (
                "pass_recall",
                "pass_precision",
                "epoch_p95",
                "swath_width_median",
                "footprint_iou_median",
            )
        ],
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "passes": passes,
        },
        fixture_metrics=recompute_passes(
            satellite_for(PRODUCT.scan_period), 0, run.day, passes
        ),
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("amsr2", __doc__)
    reference = Reference(PRODUCT, run.day, load_day(PRODUCT, run.day, run.directory))
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
    if "swath" in run.configurations:
        validate_swath(run, orbit, tle, reference)


if __name__ == "__main__":
    main()
