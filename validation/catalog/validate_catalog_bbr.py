"""
Validates the BBR catalog entry against EarthCARE's Level 1B BBR product of single
pixels (BBR_SNG_1B), using the exact values of the entry through public TAT-C methods.
BBR has three fixed telescopes, viewing about 50 deg behind (aft), at, and about 50
deg ahead of (fore) nadir along track, each imaging a line of 30 pixels. The product
gives each pixel's position, each line's time, and the platform's position at it, by
view, so the fit needs no element set. The views' angles are fitted to one day (run
with `--fit`, then `--fit-view` to center each view in TAT-C's frame with its pitch):
running on a day at least 30 days from it tests the fitted angles on data held out
from the fit. Each view's swath is scored as for other pushbroom swaths: passes of
points within the lattice of its pixels, the swath width, and the swath area over a
span of the day. The views break for about 4.3 s about every 93 s (about 5% of the
day); these breaks, like longer ones, are periods without data, outside which passes
are scored and swath areas compared.

Needs the `validation` extra and a MAAP offline token in `ESA_MAAP_KEY`. Only the
pixels' positions and the lines' times are read (of the shortwave band, every second
line), with range requests:

    python validation/catalog/validate_catalog_bbr.py <data directory> 2026-09-27 --fit
    python validation/catalog/validate_catalog_bbr.py <data directory> 2026-09-27 --fit-view
    python validation/catalog/validate_catalog_bbr.py <data directory> 2026-11-07

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit;
`--fit-view` needs its day). The element set is the one nearest the middle of the day
in `elements/earthcare_2026-09-27.3le` (or in a file given with `--elements`). With
`--write`, it writes the entry's validation record and the reduced fixture from which
unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pathlib
from datetime import datetime, timedelta

import numpy as np

from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from earthcare_io import check_time_scale, load_day, to_seconds
from scanner import (
    GATED,
    Reference,
    fit,
    fit_view,
    fit_yaw,
    observe_orbit,
    score_swath_width,
)
from validation_io import (
    build_grid,
    parse_run,
    pick_passes,
    recompute_passes,
    score_passes,
    score_swath_area,
    write_results,
)

CONFIGURATION = "views"
PRODUCT = "BBR_SNG_1B"
# the product's order of views, and the nominal pitch (deg) from which each is fitted
VIEWS = {"AFT": -50.0, "NADIR": 0.0, "FORE": 50.0}
BAND = 0
LINE_STRIDE = 2
TIME_STEPS = (10, 30, 120)
POINT_SPACING = 100e3
PASS_RADIUS = 15  # km
# the product's views break for about 4.3 s about every 93 s, periods without data
# like those over a minute in other products
GAP_THRESHOLD = 1.0  # s
FIXTURE_PASSES = 60
WIDTH_STRIDE = 50


def read_bbr(f, day: datetime) -> dict:
    """Reads each view's pixel positions, line times, and platform positions."""
    g = f["ScienceData"]
    lines = np.arange(0, g["time"].shape[-1], LINE_STRIDE)
    time = to_seconds(g["time"][:, BAND, :], day)
    part = {"time": time[1, lines]}
    for v, view in enumerate(VIEWS):
        part[f"time_{view}"] = time[v, lines]
        part[f"lat_{view}"] = g["latitude"][v, BAND][lines]
        part[f"lon_{view}"] = g["longitude"][v, BAND][lines]
        part[f"sc_lat_{view}"] = g["platform_latitude"][v, BAND][lines]
        part[f"sc_lon_{view}"] = g["platform_longitude"][v, BAND][lines]
        part[f"sc_alt_{view}"] = g["platform_altitude"][v, BAND][lines]
    return part


def load_references(day: datetime, directory: pathlib.Path) -> dict[str, Reference]:
    """Loads a day, as a reference for each view, in the order of its own times."""
    data = load_day(PRODUCT, day, directory, read_bbr)
    references = {}
    for view in VIEWS:
        seconds = data[f"time_{view}"]
        order = np.argsort(seconds)
        keep = order[(seconds[order] >= 0) & (seconds[order] < 86400)]
        _, unique = np.unique(np.round(seconds[keep], 4), return_index=True)
        keep = keep[unique]
        part = {
            "time": seconds[keep],
            "lat": data[f"lat_{view}"][keep],
            "lon": data[f"lon_{view}"][keep],
            "sc_lat": data[f"sc_lat_{view}"][keep],
            "sc_lon": data[f"sc_lon_{view}"][keep],
            "sc_alt": data[f"sc_alt_{view}"][keep],
        }
        period = float(np.median(np.diff(part["time"])))
        references[view] = Reference(day, part, period, GAP_THRESHOLD)
    return references


def validate(run, orbit, tle, references) -> None:
    """
    Scores each view of the configuration on its pixels, and writes the results if
    the run writes, with the fixture of the nadir view's passes.
    """

    def satellite_for(time_step: float):
        return run.entry.to_satellite(
            orbit,
            configuration=CONFIGURATION,
            time_step=time_step,
            platform=run.platform,
        )

    satellite = satellite_for(TIME_STEPS[0])
    points = generate_points_fibonacci_lattice(POINT_SPACING)
    grid = build_grid()
    metrics, nadir_observed = [], None
    for index, view in enumerate(VIEWS):
        reference = references[view]
        observed = reference.observe_passes(points, PASS_RADIUS)
        start = reference.time[0].to_pydatetime()
        end = reference.time[-1].to_pydatetime() + timedelta(
            seconds=reference.scan_period
        )
        results = score_passes(
            points, satellite, index, observed, start, end, reference.gaps
        )
        results["swath_width_median"] = score_swath_width(
            reference,
            satellite.instruments[index],
            orbit,
            WIDTH_STRIDE,
            satellite.attitude,
        )
        first, last, compared, covered = observe_orbit(
            reference, orbit, grid, 2 * PASS_RADIUS
        )
        results["footprint_iou_median"] = score_swath_area(
            satellite_for, index, first, last, TIME_STEPS, grid, compared, covered
        )
        print(view, json.dumps(results), flush=True)
        metrics += [
            {"name": name, "value": float(results[name]), "scope": view}
            for name in GATED
        ]
        if view == "NADIR":
            nadir_observed = observed
    if not run.write:
        return
    passes = pick_passes(points, nadir_observed, run.day, FIXTURE_PASSES)
    nadir = list(VIEWS).index("NADIR")
    write_results(
        run,
        CONFIGURATION,
        reference=f"{PRODUCT} (pixel positions of the shortwave band, every "
        "second line)",
        metrics=metrics,
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "passes": passes,
        },
        fixture_metrics=recompute_passes(satellite, nadir, run.day, passes),
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its views."""
    run = parse_run(
        "bbr",
        __doc__,
        add_arguments=lambda parser: parser.add_argument(
            "--fit-view",
            action="store_true",
            help="center each view in TAT-C's frame, with the instrument's pitch",
        ),
    )
    references = load_references(run.day, run.directory)
    for view, reference in references.items():
        print(
            f"{run.day:%Y-%m-%d} {view}: {len(reference.time)} lines (every "
            f"{LINE_STRIDE}), {len(reference.gaps)} gaps over {GAP_THRESHOLD:g} s, "
            f"line period {reference.scan_period:.4f} s",
            flush=True,
        )
    if run.fit:
        check_time_scale(PRODUCT, run.day, "ScienceData/time")
        for view, reference in references.items():
            print(f"== {view}")
            for frame in ("inertial", "earth_fixed"):
                fit(reference, frame)
                fit_yaw(reference, frame)
        return
    if run.fit_view:
        for view, reference in references.items():
            for frame in ("inertial", "earth_fixed"):
                print(f"== {view}")
                fit_view(reference, frame, 0.0, VIEWS[view], instrument_pitch=True)
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate(run, orbit, tle, references)


if __name__ == "__main__":
    main()
