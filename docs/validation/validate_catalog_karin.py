"""
Validates the KaRIn catalog entry against SWOT's Level 2 KaRIn Low Rate Sea Surface
Height product (Basic, version D), using the exact values of the entry through public
TAT-C methods. Each file is a half-orbit pass on a grid of 69 pixels, 2 km apart across
track, with the time of each line. Pixels with good-quality heights (quality flag 0)
show where KaRIn observes over the ocean, as in `ValidateSwathSWOT`:

* The swath edges are the boundaries between good-quality pixels and those beside
  them, on the sides of lines over the open ocean (at least 20 km from a coast), where
  the boundaries are those of the swath rather than of the ocean. They are compared
  with the rays at the edges of the views at the times of the lines.
* The swath area is compared, pass by pass, on the grids of lines over the open ocean
  on both sides, at several time steps.
* Passes are those of lattice points in good-quality pixels, at the time of the line.

The left and right swaths' roll angles and fields of view come from `--fit`, which
measures the look angles of the edges on a day of data. Running on a day at least 30
days from it tests them on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_karin.py <data directory> 2026-10-01 --fit
    python docs/validation/validate_catalog_karin.py <data directory> 2026-11-04

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
from datetime import datetime, timedelta, timezone

import earthaccess
import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Geod
from scipy.spatial import cKDTree
from shapely import contains_xy, prepare
from shapely.ops import unary_union
from skyfield.framelib import itrs

import tatc
from tatc.analysis import compute_ground_track
from tatc.catalog import get_entry
from tatc.catalog.metrics import compute_mask_iou
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from tatc.utils import compute_projected_ray_position
from validation_io import (
    TO_GEODETIC,
    ecef,
    get_holdout_gap,
    normalize,
    predict_passes,
    match_passes,
    recompute_passes,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "karin.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "karin.npz"

FIT_DAY = datetime(2026, 10, 1, tzinfo=timezone.utc)
PLATFORM = "swot"
CONFIGURATION = "swaths"
# public element sets from CelesTrak, with epochs after each day of data (that of the
# fit day is the one of ValidateSwathSWOT)
ELEMENT_SETS = {
    "2026-10-01": [
        "1 54754U 22173A   26277.59902338  .00000094  00000+0  65386-4 0  9997",
        "2 54754  77.6084 247.6064 0000264 143.3265 216.7905 14.00173063194463",
    ],
    "2026-11-04": [
        "1 54754U 22173A   26312.32579237  .00000026  00000+0  19087-4 0  9998",
        "2 54754  77.6110 200.7326 0000419  61.8188 298.3006 14.00168167199321",
    ],
}
NADIR = 34  # grid column on the reference ground track
OPEN_OCEAN = 20e3  # m from a coast
MIN_GOOD = 0.9  # fraction of good pixels within the swath on a side of a line
INNER_CHECK, OUTER_CHECK = 14, 56  # km, the part of a swath that must be good
TIME_STEPS = (10, 30, 120)
POINT_SPACING = 100e3
SEARCH_WINDOW = timedelta(minutes=10)
EDGE_STRIDE = 10  # lines
PIXEL_CHORD = 2 * np.sin(1.45 / 6371 / 2)  # half the diagonal of a pixel
FIXTURE_PASSES = 150
EDGES = ("left inner", "left outer", "right inner", "right outer")
GEOD = Geod(ellps="WGS84")


def load_day(day: datetime, directory: pathlib.Path) -> list[dict]:
    """Downloads and reduces the passes that start on a day, or loads them."""
    path = directory / f"swot_karin_{day:%Y%m%d}.npz"
    if not path.exists():
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        results = earthaccess.search_data(
            short_name="SWOT_L2_LR_SSH_BASIC_D",
            temporal=(f"{day:%Y-%m-%d}T00:00:00Z", f"{day:%Y-%m-%d}T23:59:59Z"),
        )
        files = sorted(
            earthaccess.download(results, str(directory / f"swot_{day:%Y%m%d}"))
        )
        passes = []
        for file in files:
            with xr.open_dataset(file) as ds:
                seconds = (
                    pd.to_datetime(ds.time.values).tz_localize("UTC")
                    - pd.Timestamp(day)
                ).total_seconds()
                if seconds[0] < 0:
                    continue
                passes.append(
                    {
                        "time": np.asarray(seconds),
                        "lat": ds.latitude.values.astype(np.float32),
                        "lon": ((ds.longitude.values + 180) % 360 - 180).astype(
                            np.float32
                        ),
                        "good": ds.ssh_karin_qual.values == 0,
                        "coast": np.nan_to_num(ds.distance_to_coast.values)
                        >= OPEN_OCEAN,
                    }
                )
        np.savez_compressed(
            path,
            lines=[len(p["time"]) for p in passes],
            **{k: np.concatenate([p[k] for p in passes]) for k in passes[0]},
        )
    with np.load(path, allow_pickle=False) as data:
        cache = {key: data[key] for key in data.files}
    cuts = np.cumsum(cache["lines"])[:-1]
    keys = ("time", "lat", "lon", "good", "coast")
    passes = [
        dict(zip(keys, values))
        for values in zip(*(np.split(cache[k], cuts) for k in keys))
    ]
    for item in passes:
        add_distances(item)
    print(f"{len(passes)} passes, {int(cache['lines'].sum())} lines", flush=True)
    return passes


def add_distances(item: dict) -> None:
    """Adds the distance (km) of each pixel from the ground track, positive left."""
    lat, lon = item["lat"].astype(float), item["lon"].astype(float)
    following = np.r_[np.arange(1, len(lat)), len(lat) - 2]
    heading = GEOD.inv(
        lon[:, NADIR], lat[:, NADIR], lon[following, NADIR], lat[following, NADIR]
    )[0]
    heading[-1] = (heading[-1] + 180) % 360  # the last line looks back
    azimuth, _, distance = GEOD.inv(
        np.repeat(lon[:, [NADIR]], lat.shape[1], axis=1),
        np.repeat(lat[:, [NADIR]], lat.shape[1], axis=1),
        lon,
        lat,
    )
    item["distance"] = (
        -np.sign(np.sin(np.radians(azimuth - heading[:, None]))) * distance / 1e3
    )
    item["distance"][:, NADIR] = 0


def open_sides(item: dict) -> np.ndarray:
    """Gets whether each side (left, right) of each line is over the open ocean."""
    result = []
    for sign in (1, -1):
        swath = (sign * item["distance"] >= INNER_CHECK) & (
            sign * item["distance"] <= OUTER_CHECK
        )
        side = sign * item["distance"] > 0
        coast = np.all(item["coast"] | ~side, axis=1)
        good = np.sum(item["good"] & swath, axis=1) >= MIN_GOOD * np.sum(swath, axis=1)
        result.append(coast & good & np.isfinite(item["time"]))
    return np.stack(result, axis=1)


def find_edges(item: dict, lines: np.ndarray, side: int) -> dict[str, np.ndarray]:
    """
    Gets the inner and outer edges of a side's swath on lines, as the midpoints between
    the good-quality pixels at the ends of the swath and the pixels beside them.
    """
    sign = 1 if side == 0 else -1
    edges = {"inner": [], "outer": [], "line": []}
    for line in lines:
        # the columns of the side, outward from the ground track
        columns = np.nonzero(sign * item["distance"][line] > 0)[0]
        columns = columns[np.argsort(np.abs(item["distance"][line, columns]))]
        good = item["good"][line, columns]
        middle = np.searchsorted(np.abs(item["distance"][line, columns]), 35)
        if not good[middle]:
            continue
        # the run of good pixels around the middle of the swath
        first = middle
        while first > 0 and good[first - 1]:
            first -= 1
        last = middle
        while last < len(columns) - 1 and good[last + 1]:
            last += 1
        if first == 0 or last == len(columns) - 1:
            continue
        pairs = {
            "inner": (columns[first], columns[first - 1]),
            "outer": (columns[last], columns[last + 1]),
        }
        for edge, (a, b) in pairs.items():
            edges[edge].append(
                ecef(
                    (item["lon"][line, a] + item["lon"][line, b]) / 2,
                    (item["lat"][line, a] + item["lat"][line, b]) / 2,
                )
            )
        edges["line"].append(line)
    return {k: np.array(v) for k, v in edges.items()}


def to_utc(day: datetime, seconds) -> list[datetime]:
    """Converts seconds from the start of a day to times."""
    return [day + timedelta(seconds=float(s)) for s in np.atleast_1d(seconds)]


def look_angles(orbit, day, seconds, points) -> np.ndarray:
    """Gets the look angles (deg, positive left) of points from the propagated orbit."""
    track = orbit.to_gp_orbit().get_orbit_track(to_utc(day, seconds))
    position, velocity = track.frame_xyz_and_velocity(itrs)
    position, velocity = np.array(position.m).T, np.array(velocity.m_per_s).T
    lon, lat, _ = TO_GEODETIC.transform(*position.T)
    nadir = -up(lon, lat)
    left = normalize(np.cross(velocity, nadir))
    los = normalize(points - position)
    return np.degrees(
        np.arctan2(np.sum(los * left, axis=-1), np.sum(los * nadir, axis=-1))
    )


def measure_edges(passes, orbit, day) -> dict[str, np.ndarray]:
    """Gets the look angles of each edge on the open-ocean sides of lines."""
    angles = {edge: [] for edge in EDGES}
    for item in passes:
        sides = open_sides(item)
        for side, label in enumerate(("left", "right")):
            lines = np.nonzero(sides[:, side])[0][::EDGE_STRIDE]
            edges = find_edges(item, lines, side)
            if not len(edges["line"]):
                continue
            seconds = item["time"][edges["line"]]
            for edge in ("inner", "outer"):
                angles[f"{label} {edge}"].append(
                    look_angles(orbit, day, seconds, edges[edge])
                )
    return {k: np.concatenate(v) for k, v in angles.items()}


def fit(passes, orbit, day) -> None:
    """Measures the look angles of the edges of each swath."""
    angles = measure_edges(passes, orbit, day)
    for edge, values in angles.items():
        print(
            f"{edge}: {len(values)} lines; look angle median {np.median(values):.4f} "
            f"deg, p5 {np.percentile(values, 5):.4f}, p95 {np.percentile(values, 95):.4f}"
        )
    for label in ("left", "right"):
        inner = float(np.median(angles[f"{label} inner"]))
        outer = float(np.median(angles[f"{label} outer"]))
        print(
            f"{label}: roll {(inner + outer) / 2:.4f} deg, field of view "
            f"{abs(outer - inner):.4f} deg"
        )


def score_edges(passes, orbit, day, satellite) -> dict[str, np.ndarray]:
    """
    Gets the horizontal distances (km) between the measured edges and the rays at the
    edges of the views at the times of the lines.
    """
    distances = {edge: [] for edge in EDGES}
    for item in passes:
        sides = open_sides(item)
        for side, label in enumerate(("left", "right")):
            lines = np.nonzero(sides[:, side])[0][::EDGE_STRIDE]
            edges = find_edges(item, lines, side)
            if not len(edges["line"]):
                continue
            instrument = satellite.instruments[side]
            track = orbit.to_gp_orbit().get_orbit_track(
                to_utc(day, item["time"][edges["line"]])
            )
            half = instrument.cross_track_field_of_view / 2
            sign = 1 if instrument.roll_angle > 0 else -1
            for edge, offset in (("inner", -sign * half), ("outer", sign * half)):
                ray = compute_projected_ray_position(
                    track, 0, 0, instrument.roll_angle + offset
                )
                point = edges[edge]
                vertical = normalize(point)
                difference = ecef(ray.longitude.degrees, ray.latitude.degrees) - point
                horizontal = (
                    difference
                    - np.sum(difference * vertical, axis=-1)[:, None] * vertical
                )
                distances[f"{label} {edge}"].append(
                    np.linalg.norm(horizontal, axis=-1) / 1e3
                )
    return {k: np.concatenate(v) for k, v in distances.items()}


def score_pass_area(item, day, satellite_for) -> dict[str, float]:
    """Gets the intersection over union of a pass's open-ocean lines, by time step."""
    timed = np.isfinite(item["time"])
    lines = np.all(open_sides(item), axis=1) & timed
    if lines.sum() < 100:
        return {}
    # the grid of the open-ocean lines
    lon, lat = item["lon"][lines].ravel(), item["lat"][lines].ravel()
    observed = item["good"][lines].ravel()
    first, last = item["time"][timed][0], item["time"][timed][-1]
    results = {"lines": int(lines.sum())}
    for time_step in TIME_STEPS:
        n = max(1, int(round((last - first) / time_step)))
        step = (last - first) / n
        times = to_utc(day, first + step * (np.arange(n) + 0.5))
        satellite = satellite_for(step)
        geometry = unary_union(
            [
                compute_ground_track(
                    satellite,
                    times,
                    instrument_index=i,
                    time_step=timedelta(seconds=float(step)),
                ).geometry.iloc[0]
                for i in range(len(satellite.instruments))
            ]
        )
        prepare(geometry)
        predicted = contains_xy(geometry, lon, lat)
        results[f"iou_{time_step}"] = compute_mask_iou(
            observed, predicted, np.ones(len(lon))
        )
    return results


def observe_passes(points, passes, day) -> pd.DataFrame:
    """Finds the passes of lattice points in good-quality pixels."""
    unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
    rows = []
    for k, item in enumerate(passes):
        pixels = up(item["lon"].astype(float), item["lat"].astype(float)).reshape(-1, 3)
        distance, nearest = cKDTree(pixels).query(
            unit, distance_upper_bound=PIXEL_CHORD
        )
        found = np.isfinite(distance)
        found[found] &= item["good"].ravel()[nearest[found]]
        found[found] &= np.isfinite(item["time"][nearest[found] // 69])
        for point, pixel in zip(np.nonzero(found)[0], nearest[found]):
            rows.append(
                {
                    "point_id": int(point),
                    "pass": k,
                    "time": pd.Timestamp(day)
                    + pd.Timedelta(seconds=float(item["time"][pixel // 69])),
                }
            )
    return pd.DataFrame(rows)


def validate(day, tle, entry, orbit, passes, write) -> None:
    """Scores the configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, time_step=time_step, platform=PLATFORM)

    satellite = satellite_for(TIME_STEPS[0])
    edges = score_edges(passes, orbit, day, satellite)
    for edge, values in edges.items():
        print(
            f"swath edge, {edge}: {len(values)} lines, median {np.median(values):.2f} "
            f"km, p95 {np.percentile(values, 95):.2f} km",
            flush=True,
        )
    lattice = generate_points_fibonacci_lattice(POINT_SPACING)
    observed = observe_passes(lattice, passes, day)
    ids = np.sort(observed.point_id.unique())
    subset = lattice.iloc[ids].reset_index(drop=True)
    observed["point_id"] = np.searchsorted(ids, observed.point_id)
    start = to_utc(day, np.nanmin(passes[0]["time"]))[0] - SEARCH_WINDOW
    end = to_utc(day, np.nanmax(passes[-1]["time"]))[0] + SEARCH_WINDOW
    print(f"predicting passes of {len(subset)} points; {len(observed)} observed")
    predicted = predict_passes(subset, satellite, None, start, end, [])
    results = match_passes(predicted, observed)
    print(json.dumps(results), flush=True)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        area = pd.DataFrame(
            [
                r
                for r in pool.map(
                    lambda item: score_pass_area(item, day, satellite_for), passes
                )
                if r
            ]
        )
    for time_step in TIME_STEPS:
        print(
            f"swath IoU at {time_step} s over {len(area)} passes: median "
            f"{area[f'iou_{time_step}'].median():.5f}, p5 "
            f"{area[f'iou_{time_step}'].quantile(0.05):.5f}",
            flush=True,
        )
    metrics = [
        {"name": "pass_recall", "value": float(results["pass_recall"])},
        *[
            {
                "name": "swath_edge_p95",
                "value": float(np.percentile(values, 95)),
                "scope": edge,
            }
            for edge, values in edges.items()
        ],
        *[
            {"name": name, "value": float(value), "scope": f"time step {time_step} s"}
            for time_step in TIME_STEPS
            for name, value in (
                ("footprint_iou_median", area[f"iou_{time_step}"].median()),
                ("footprint_iou_p5", area[f"iou_{time_step}"].quantile(0.05)),
            )
        ],
        {"name": "epoch_p95", "value": float(results["epoch_p95"])},
        {"name": "epoch_median", "value": float(results["epoch_median"])},
    ]
    print(json.dumps(metrics), flush=True)
    if not write:
        return
    rng = np.random.default_rng(1)
    picks = observed.iloc[rng.choice(len(observed), FIXTURE_PASSES, replace=False)]
    fixture = np.column_stack(
        [
            subset.geometry.y.to_numpy()[picks.point_id],
            subset.geometry.x.to_numpy()[picks.point_id],
            (picks.time - pd.Timestamp(day)).dt.total_seconds().to_numpy(),
        ]
    )
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=fixture
    )
    update_record(
        RECORD_PATH,
        "karin",
        CONFIGURATION,
        {
            "digest": entry.configurations[CONFIGURATION].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_karin.py",
            "reference": "SWOT_L2_LR_SSH_BASIC_D (good-quality pixels and line times)",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": metrics,
            "fixture_metrics": recompute_passes(satellite, None, day, fixture),
        },
    )
    print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


def main() -> None:
    """Fits or validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--download-only", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    passes = load_day(day, arguments.directory)
    if arguments.download_only:
        return
    tle = ELEMENT_SETS[arguments.day]
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    if arguments.fit:
        fit(passes, orbit, day)
        return
    validate(day, tle, get_entry("karin"), orbit, passes, arguments.write)


if __name__ == "__main__":
    main()
