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
* The swath area over an orbit is compared with the union of the cells.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python docs/validation/validate_catalog_ssmis.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_ssmis.py <data directory> 2026-11-03

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import json
import pathlib
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.spatial import cKDTree

import tatc
from tatc.catalog import get_entry
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import Product, Reference, load_day, observe_orbit
from validation_io import (
    TO_GEODETIC,
    build_grid,
    ecef,
    get_holdout_gap,
    normalize,
    recompute_passes,
    score_passes,
    score_swath_area,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "ssmis.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "ssmis.npz"
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
PLATFORM = "dmsp-f18"
CONFIGURATION = "swath"
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-03": (
        "1 35951U 09057A   26312.25663183  .00000227  00000+0  14044-3 0  9998",
        "2 35951  98.9223 292.8431 0010542 197.2322 162.8494 14.14962193879968",
    ),
}
PRODUCT = Product(
    short_name="GPM_1CF18SSMIS",
    version="08",
    group="S1",
    spacecraft=("SCstatus/SClatitude", "SCstatus/SClongitude", "SCstatus/SCaltitude"),
)
SCAN_PERIOD = 1.899
TIME_STEPS = (SCAN_PERIOD, 10, 30)
POINT_SPACING = 300e3
PASS_RADIUS = 60  # km
FIXTURE_PASSES = 80
GEOD = Geod(ellps="WGS84")


def fit(reference: Reference) -> None:
    """
    Prints the fitted scan: the median angle of the samples from nadir (the cone
    angle), and the azimuths of the first and last samples from the velocity, in the
    inertial and Earth-fixed frames, with their center by latitude.
    """
    lon, lat, _ = TO_GEODETIC.transform(*reference.position.T)
    nadir = -up(lon, lat)
    sight = normalize(
        ecef(reference.lon, reference.lat) - reference.position[:, np.newaxis]
    )
    cone = np.degrees(np.arccos(np.sum(sight * nadir[:, np.newaxis], axis=-1)))
    print(
        f"cone angle: median {np.median(cone):.4f} deg, p5 "
        f"{np.percentile(cone, 5):.4f}, p95 {np.percentile(cone, 95):.4f}"
    )
    for frame in ("inertial", "earth_fixed"):
        motion = reference.velocity(frame)
        forward = normalize(motion - np.sum(motion * nadir, axis=-1)[:, None] * nadir)
        # to the left of the motion
        left = normalize(np.cross(forward, nadir))
        azimuth = np.degrees(
            np.arctan2(
                np.sum(sight * left[:, None], axis=-1),
                np.sum(sight * forward[:, None], axis=-1),
            )
        )
        first, last = np.median(azimuth[:, 0]), np.median(azimuth[:, -1])
        center = (azimuth[:, 0] + azimuth[:, -1]) / 2
        bands = ", ".join(
            f"{low:+d}..{low + 30:+d}: "
            f"{np.median(center[(np.asarray(lat) >= low) & (np.asarray(lat) < low + 30)]):+.3f}"
            for low in range(-90, 90, 30)
            if np.any((np.asarray(lat) >= low) & (np.asarray(lat) < low + 30))
        )
        print(
            f"{frame}: first sample {first:+.4f} deg, last {last:+.4f} deg; center "
            f"{np.median(center):+.4f} (p5 {np.percentile(center, 5):+.4f}, p95 "
            f"{np.percentile(center, 95):+.4f}), half width "
            f"{np.median(np.abs(azimuth[:, -1] - azimuth[:, 0]) / 2):.4f} deg; "
            f"center by latitude: {bands}"
        )
    print(f"scan period: median {np.median(np.diff(reference.seconds)):.6f} s")


def cone_crossing(reference: Reference, target, k: int, cone_angle: float) -> float:
    """Gets the time (s) a target crosses the cone near scan k, by bisection."""

    def angle(t: float) -> float:
        seconds = reference.seconds
        i = int(np.clip(np.searchsorted(seconds, t) - 1, 0, len(seconds) - 2))
        f = (t - seconds[i]) / (seconds[i + 1] - seconds[i])
        position = reference.position[i] + f * (
            reference.position[i + 1] - reference.position[i]
        )
        lon, lat, _ = TO_GEODETIC.transform(*position)
        return (
            np.degrees(np.arccos(np.dot(normalize(target - position), -up(lon, lat))))
            - cone_angle
        )

    lower, upper = reference.seconds[k] - 8, reference.seconds[k] + 8
    if angle(lower) * angle(upper) > 0:
        return np.nan
    for _ in range(30):
        middle = (lower + upper) / 2
        lower, upper = (middle, upper) if angle(middle) > 0 else (lower, middle)
    return (lower + upper) / 2


def observe_passes(reference: Reference, points, cone_angle: float) -> pd.DataFrame:
    """
    Finds the passes on which the samples observed points (in a cell of the lattice of
    samples), and their times (when they cross the cone).
    """
    point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
    flat = reference.unit.reshape(-1, 3)
    scan_of = np.repeat(np.arange(len(reference.seconds)), reference.beams)
    rows = []
    for p, hits in enumerate(
        cKDTree(flat).query_ball_point(
            point_unit, r=2 * np.sin(PASS_RADIUS / 6371 / 2), workers=-1
        )
    ):
        if not hits:
            continue
        hits = np.array(hits)
        hits = hits[np.argsort(reference.seconds[scan_of[hits]])]
        times = reference.seconds[scan_of[hits]]
        for group in np.split(hits, np.nonzero(np.diff(times) > 600)[0] + 1):
            k = int(
                scan_of[
                    group[
                        np.argmin(np.linalg.norm(flat[group] - point_unit[p], axis=-1))
                    ]
                ]
            )
            first, last = max(k - 2, 0), min(k + 2, len(reference.seconds) - 1)
            if not reference.in_cells(point_unit[[p]], first, last)[0]:
                continue
            crossing = cone_crossing(
                reference,
                ecef(points.geometry.x.iloc[p], points.geometry.y.iloc[p]),
                k,
                cone_angle,
            )
            if np.isfinite(crossing):
                rows.append({"point_id": p, "seconds": float(crossing)})
    observed = pd.DataFrame(rows).astype({"point_id": "int64"})
    observed["time"] = pd.Timestamp(reference.day) + pd.to_timedelta(
        observed.seconds, unit="s"
    )
    return observed


def score_swath_width(reference: Reference, instrument) -> float:
    """
    Gets the median absolute difference (km) of the swath width from the distance
    along the surface between the first and last samples of scans.
    """
    scans = np.arange(0, len(reference.seconds), 20)
    distance = GEOD.inv(
        reference.lon[scans, 0],
        reference.lat[scans, 0],
        reference.lon[scans, -1],
        reference.lat[scans, -1],
    )[2]
    width = np.vectorize(instrument.get_swath_width)(reference.altitude[scans])
    return float(np.median(np.abs(width - distance))) / 1e3


def main() -> None:
    """Validates the entry on a day of data, or fits its scan."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    reference = Reference(day, load_day(PRODUCT, day, arguments.directory), SCAN_PERIOD)
    print(
        f"{arguments.day}: {len(reference.time)} scans of {reference.beams} samples, "
        f"{len(reference.gaps)} gaps over a minute, median altitude "
        f"{np.median(reference.altitude) / 1e3:.1f} km",
        flush=True,
    )
    if arguments.fit:
        fit(reference)
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("ssmis")

    def satellite_for(time_step: float):
        return entry.to_satellite(
            orbit, CONFIGURATION, time_step=time_step, platform=PLATFORM
        )

    satellite = satellite_for(SCAN_PERIOD)
    cone_angle = (
        entry.configurations[CONFIGURATION].instruments[0].parameters["cone_angle"]
    )
    points = generate_points_fibonacci_lattice(POINT_SPACING)
    observed = observe_passes(reference, points, cone_angle)
    start = reference.time[0].to_pydatetime()
    end = reference.time[-1].to_pydatetime() + timedelta(seconds=SCAN_PERIOD)
    results = score_passes(points, satellite, 0, observed, start, end, reference.gaps)
    results["swath_width_median"] = score_swath_width(
        reference, satellite.instruments[0]
    )
    print(json.dumps(results), flush=True)
    grid = build_grid()
    first, last, compared, covered = observe_orbit(
        reference, orbit, grid, 2 * PASS_RADIUS
    )
    results["footprint_iou_median"] = score_swath_area(
        satellite_for, 0, first, last, TIME_STEPS, grid, compared, covered
    )
    print(json.dumps(results), flush=True)
    if not arguments.write:
        return
    rng = np.random.default_rng(1)
    picks = observed.iloc[rng.choice(len(observed), FIXTURE_PASSES, replace=False)]
    passes = np.stack(
        [
            points.geometry.y.to_numpy()[picks.point_id],
            points.geometry.x.to_numpy()[picks.point_id],
            picks.seconds.to_numpy(),
        ],
        axis=-1,
    )
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=passes
    )
    update_record(
        RECORD_PATH,
        "ssmis",
        CONFIGURATION,
        {
            "digest": entry.configurations[CONFIGURATION].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_ssmis.py",
            "reference": "GPM_1CF18SSMIS version 08 (positions of the S1 samples)",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": [
                {"name": name, "value": float(results[name])}
                for name in (
                    "pass_recall",
                    "pass_precision",
                    "epoch_p95",
                    "swath_width_median",
                    "footprint_iou_median",
                )
            ],
            "fixture_metrics": recompute_passes(satellite, 0, day, passes),
        },
    )
    print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()
