"""
What the validations of the side-looking radars (C-SAR and L-SAR) share: the products'
geolocation grids, with the zero-Doppler time of each line and the restituted orbit,
from which the roll angles of the swath edges are measured, and against which the
swath edges, passes, and epochs are scored.

A product is a dictionary of its grid's `time` (s from the start of the day, by line),
`lon` and `lat` (lines by samples, near range first), its `datatake`, and its orbit's
`orbit_time`, `orbit_position`, and `orbit_velocity` (Earth-fixed).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from tatc.analysis import collect_observations
from tatc.catalog.metrics import summarize
from tatc.schemas import Point
from tatc.utils import compute_argument_of_latitude, compute_projected_ray_position
from validation_io import (
    TO_GEODETIC,
    ecef,
    in_lattice_cells,
    match_passes,
    normalize,
    predict_passes,
    up,
)

PROFILE_BIN = 10  # deg of argument of latitude
TIME_STEPS = (10, 30, 120)
CHUNK_PRODUCTS = 8
MIN_CHUNK_PRODUCTS = 3
MAX_PRODUCT_GAP = 5  # s, between the last line of a product and the first of the next
EPOCH_PRODUCTS = 40
EPOCH_LINES = 3
EPOCH_SAMPLES = np.linspace(0.1, 0.9, 5)
POINT_SPACING = 100e3
PRODUCT_RADIUS = 2 * np.sin(350 / 6371 / 2)  # chord, on the unit sphere
SEARCH_WINDOW = timedelta(minutes=10)
FIXTURE_PASSES = 150


def to_utc(day: datetime, seconds) -> list[datetime]:
    """Converts seconds from the start of a day to times."""
    return [day + timedelta(seconds=float(s)) for s in np.atleast_1d(seconds)]


def interpolate_orbit(product: dict, seconds) -> tuple[np.ndarray, np.ndarray]:
    """Gets the restituted position and velocity at times (cubic Hermite)."""
    t, p, v = (
        product["orbit_time"],
        product["orbit_position"],
        product["orbit_velocity"],
    )
    seconds = np.atleast_1d(seconds)
    k = np.clip(np.searchsorted(t, seconds) - 1, 0, len(t) - 2)
    h = (t[k + 1] - t[k])[:, None]
    s = ((seconds - t[k]) / (t[k + 1] - t[k]))[:, None]
    position = (
        (2 * s**3 - 3 * s**2 + 1) * p[k]
        + (s**3 - 2 * s**2 + s) * h * v[k]
        + (-2 * s**3 + 3 * s**2) * p[k + 1]
        + (s**3 - s**2) * h * v[k + 1]
    )
    velocity = (
        (6 * s**2 - 6 * s) / h * p[k]
        + (3 * s**2 - 4 * s + 1) * v[k]
        + (-6 * s**2 + 6 * s) / h * p[k + 1]
        + (3 * s**2 - 2 * s) * v[k + 1]
    )
    return position, velocity


def measure_lines(products: list[dict], orbit, day: datetime) -> pd.DataFrame:
    """
    Gets the roll angles (positive left) of the near-range and far-range points of each
    grid line, measured from the geodetic nadir in the plane of the line, with the
    argument of latitude of the satellite from the element set.
    """
    rows = []
    for k, product in enumerate(products):
        times = product["time"]
        position, velocity = interpolate_orbit(product, times)
        lon, lat, _ = TO_GEODETIC.transform(*position.T)
        nadir = -up(lon, lat)
        left = normalize(np.cross(velocity, nadir))
        angles = []
        for column in (0, -1):
            los = normalize(
                ecef(product["lon"][:, column], product["lat"][:, column]) - position
            )
            angles.append(
                np.degrees(
                    np.arctan2(
                        np.sum(los * left, axis=-1), np.sum(los * nadir, axis=-1)
                    )
                )
            )
        track = orbit.to_gp_orbit().get_orbit_track(to_utc(day, times))
        rows.append(
            pd.DataFrame(
                {
                    "product": k,
                    "mode": product["mode"],
                    "argument_of_latitude": compute_argument_of_latitude(
                        track.position.m, track.velocity.m_per_s
                    ),
                    "near": angles[0],
                    "far": angles[1],
                }
            )
        )
    lines = pd.concat(rows, ignore_index=True)
    lines["roll"] = (lines.near + lines.far) / 2
    lines["field_of_view"] = np.abs(lines.near - lines.far)
    return lines


def fit_profile(lines: pd.DataFrame) -> list[tuple[float, float]]:
    """Gets the median roll angles of the lines in bins of argument of latitude."""
    centers = (
        lines.argument_of_latitude // PROFILE_BIN
    ) * PROFILE_BIN + PROFILE_BIN / 2
    medians = lines.groupby(centers).roll.median()
    return [(float(u), round(float(r), 4)) for u, r in medians.items()]


def score_edges(products, orbit, day, instrument, look: str) -> dict[str, np.ndarray]:
    """
    Gets the horizontal distances (km) between the near-range and far-range points of
    the grid lines and the rays at the edges of the view at the times of the lines, for
    a view that looks to the `"right"` (negative roll angles, with near range on the
    left) or to the `"left"` (positive roll angles, with near range on the right).
    """
    near = 1 if look == "right" else -1
    distances = {"near": [], "far": []}
    for product in products:
        times = product["time"]
        track = orbit.to_gp_orbit().get_orbit_track(to_utc(day, times))
        roll = instrument.get_roll_angle(track)
        half = instrument.cross_track_field_of_view / 2
        for edge, column, sign in (("near", 0, near), ("far", -1, -near)):
            ray = compute_projected_ray_position(track, 0, 0, roll + sign * half)
            lon, lat = product["lon"][:, column], product["lat"][:, column]
            offset = ecef(ray.longitude.degrees, ray.latitude.degrees) - ecef(lon, lat)
            vertical = up(lon, lat)
            horizontal = offset - np.sum(offset * vertical, axis=-1)[:, None] * vertical
            distances[edge].append(np.linalg.norm(horizontal, axis=-1) / 1e3)
    return {edge: np.concatenate(values) for edge, values in distances.items()}


def print_edges(edges: dict[str, np.ndarray]) -> None:
    """Prints the median and 95th percentile distances of the swath edges."""
    for edge, values in edges.items():
        print(
            f"swath edge, {edge} range: median {np.median(values):.2f} km, p95 "
            f"{np.percentile(values, 95):.2f} km",
            flush=True,
        )


def observe_passes(points, products, day) -> pd.DataFrame:
    """Finds the passes of points inside grids, one per point and acquisition."""
    lon, lat = points.geometry.x.to_numpy(), points.geometry.y.to_numpy()
    # geocentric directions, as for the grids
    unit = normalize(ecef(lon, lat))
    tree = cKDTree(unit)
    rows = []
    for product in products:
        cells = normalize(ecef(product["lon"], product["lat"]))
        center = normalize(cells.reshape(-1, 3).mean(axis=0))
        nearby = np.array(tree.query_ball_point(center, PRODUCT_RADIUS), dtype=int)
        if len(nearby) == 0:
            continue
        inside = in_lattice_cells(
            cells, np.ones(len(cells) - 1, dtype=bool), unit[nearby]
        )
        found = nearby[inside]
        if len(found) == 0:
            continue
        # the time of the nearest grid line
        _, nearest = cKDTree(cells.reshape(-1, 3)).query(unit[found])
        line = nearest // cells.shape[1]
        rows += [
            {
                "point_id": int(i),
                "datatake": product["datatake"],
                "time": pd.Timestamp(day) + pd.Timedelta(seconds=product["time"][j]),
            }
            for i, j in zip(found, line, strict=True)
        ]
    observed = pd.DataFrame(rows).drop_duplicates(["point_id", "datatake"])
    return observed.reset_index(drop=True)


def score_passes(points, satellite, products, day) -> dict:
    """Compares the passes predicted with those of lattice points inside grids."""
    observed = observe_passes(points, products, day)
    ids = np.sort(observed.point_id.unique())
    observed["point_id"] = np.searchsorted(ids, observed.point_id)
    subset = points.iloc[ids].reset_index(drop=True)
    start = to_utc(day, products[0]["time"][0])[0] - SEARCH_WINDOW
    end = to_utc(day, products[-1]["time"][-1])[0] + SEARCH_WINDOW
    print(
        f"predicting passes of {len(subset)} points; {len(observed)} observed",
        flush=True,
    )
    predicted = predict_passes(subset, satellite, 0, start, end, [])
    return match_passes(predicted, observed)


def epoch_points(products: list[dict]) -> np.ndarray:
    """
    Gets interior grid points of a sample of products, as rows of latitude, longitude,
    and zero-Doppler seconds from the start of the day.
    """
    stride = max(1, len(products) // EPOCH_PRODUCTS)
    rows = []
    for product in products[::stride]:
        n_line, n_sample = product["lon"].shape
        for i in np.linspace(1, n_line - 2, EPOCH_LINES).astype(int):
            for j in (EPOCH_SAMPLES * (n_sample - 1)).round().astype(int):
                rows.append(
                    [product["lat"][i, j], product["lon"][i, j], product["time"][i]]
                )
    return np.array(rows)


def score_epochs(points: np.ndarray, satellite, day) -> np.ndarray:
    """Gets the differences (s) of epochs from the zero-Doppler times of points."""

    def observe(row):
        latitude, longitude, seconds = row
        time = day + timedelta(seconds=float(seconds))
        found = collect_observations(
            Point(latitude=float(latitude), longitude=float(longitude)),
            satellite,
            time - SEARCH_WINDOW,
            time + SEARCH_WINDOW,
        )
        if found.empty:
            return np.nan
        differences = (found.epoch - pd.Timestamp(time)).dt.total_seconds().to_numpy()
        return differences[np.argmin(np.abs(differences))]

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        return np.array(list(pool.map(observe, points)))


def summarize_epochs(differences: np.ndarray) -> dict:
    """Summarizes and prints the differences of the epochs of the points observed."""
    found = np.isfinite(differences)
    epochs = summarize(differences[found])
    print(
        f"epochs: {int(found.sum())} of {len(differences)} points observed, median "
        f"{epochs['median']:+.3f} s, p95 {epochs['p95']:.3f} s",
        flush=True,
    )
    return epochs


def get_metrics(passes, edges, area: pd.DataFrame, epochs) -> list[dict]:
    """
    Gets the metrics of a configuration from its passes, swath edges, swath areas (a
    table of IoUs by time step, of chunks; none if it is empty), and epochs.
    """
    metrics = [
        {"name": "pass_recall", "value": float(passes["pass_recall"])},
        *[
            {
                "name": "swath_edge_p95",
                "value": float(np.percentile(values, 95)),
                "scope": f"{edge} range",
            }
            for edge, values in edges.items()
        ],
    ]
    if not area.empty:
        metrics += [
            {"name": name, "value": float(value), "scope": f"time step {time_step} s"}
            for time_step in TIME_STEPS
            for name, value in (
                ("footprint_iou_median", area[f"iou_{time_step}"].median()),
                ("footprint_iou_p5", area[f"iou_{time_step}"].quantile(0.05)),
            )
        ]
    metrics += [
        {"name": "epoch_p95", "value": float(epochs["p95"])},
        {"name": "epoch_median", "value": float(epochs["median"])},
    ]
    return metrics


def pick_points(points: np.ndarray, differences: np.ndarray) -> np.ndarray:
    """Picks observed grid points (see `epoch_points`) for a test fixture, at random."""
    observed = points[np.isfinite(differences)]
    rng = np.random.default_rng(1)
    return observed[
        rng.choice(len(observed), min(FIXTURE_PASSES, len(observed)), False)
    ]
