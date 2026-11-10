"""
Validates the ASCAT catalog entry against the Metop-C ASCAT Level 2 25 km ocean
surface wind vector product of the EUMETSAT OSI SAF (KNMI), distributed by PO.DAAC,
using the exact values of the entry through public TAT-C methods. Each file is an
orbit of rows of 42 wind vector cells (21 on each side of the nadir gap, 25 km
apart). The cells' positions are given everywhere, also where winds are not
retrieved (over land and sea ice), so they show the swaths' geometry:

* The swath edges are half a cell beyond the outermost and innermost cells of each
  side, extrapolated from the cells beside them. They are compared with the rays at
  the edges of the views at the times of the rows.
* The swath area is compared on segments of 360 s, along lines across track through
  the rows (sampled every 5 km, through the nadir gap and 50 km beyond the edges), at
  several time steps. Segments that reach latitudes beyond 70 deg are not compared,
  as swaths there are not simple in longitude and latitude.
* Passes are those of lattice points within half a cell (12.5 km) of a cell, at the
  time of its row.

The product stores times as 32-bit floating-point seconds from 1990, which resolve
128 s. Rows are 3.75 s apart, so each file's row times are recovered as the line
nearest (in the largest difference) to the stored times, beyond their rounding by
the excess printed for each file.

The left and right swaths' roll angles and fields of view come from `--fit`, which
measures the look angles of the edges on a day of data. Running on a day at least 30
days from it tests them on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_ascat.py <data directory> 2026-10-05 --fit
    python docs/validation/validate_catalog_ascat.py <data directory> 2026-11-08

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
import netCDF4
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.optimize import linprog
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
    get_with_retries,
    match_passes,
    normalize,
    predict_passes,
    read_element_set_file,
    recompute_passes,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "ascat.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "ascat.npz"

SHORT_NAME = "ASCATC-L2-25km"
NORAD_ID = 43689
FIT_DAY = datetime(2026, 10, 5, tzinfo=timezone.utc)
PLATFORM = "metop-c"
CONFIGURATION = "swaths"
# element sets of the fitting day, and public element sets from CelesTrak, with
# epochs after each day of data
ELEMENT_SET_FILE = REPOSITORY / "docs" / "validation" / "metop_2026-10-05.3le"
ELEMENT_SETS = {
    "2026-11-08": [
        "1 43689U 18087A   26313.63049069  .00000132  00000+0  80485-4 0  9990",
        "2 43689  98.7253  10.9624 0001551 125.2403 234.8919 14.21484744415495",
    ],
}
EPOCH = datetime(1990, 1, 1, tzinfo=timezone.utc)
CELLS = 42
# (edge cell, the cell beside it inward) of each side's outer and inner edges
SIDES = ({"outer": (0, 1), "inner": (20, 19)}, {"outer": (41, 40), "inner": (21, 22)})
ROW_GAP = 5.0  # s, the largest step between consecutive rows of a segment
TIME_STEPS = (10, 30, 120)
POINT_SPACING = 100e3
SEARCH_WINDOW = timedelta(minutes=10)
EDGE_STRIDE = 10  # rows
CELL_CHORD = 2 * np.sin(12.5 / 6371 / 2)  # half the spacing of cells
SEGMENT = 360  # s
LINE_STEP, LINE_MARGIN = 5e3, 50e3  # m
MAX_LATITUDE = 70
FIXTURE_PASSES = 150
EDGES = ("left inner", "left outer", "right inner", "right outer")
GEOD = Geod(ellps="WGS84")


def reconstruct_times(stored: np.ndarray) -> tuple[np.ndarray, float]:
    """
    Gets the times of rows (seconds from 1990) on the line nearest to their stored
    times, and the largest difference beyond the stored times' rounding (s).
    """
    valid = np.flatnonzero(stored < 1e29)
    middle = valid[len(valid) // 2]
    q = stored[valid] - stored[middle]
    half = float(np.spacing(np.float32(stored[valid].max()))) / 2
    rows = np.c_[np.ones(len(valid)), valid - middle]
    result = linprog(
        [0, 0, 1],
        A_ub=np.vstack(
            [np.c_[rows, -np.ones(len(valid))], np.c_[-rows, -np.ones(len(valid))]]
        ),
        b_ub=np.r_[q + half, -q + half],
        bounds=[(None, None), (None, None), (0, None)],
    )
    if result.status != 0:
        raise RuntimeError(f"no line fits the stored times: {result.message}")
    offset, period, excess = result.x
    times = stored[middle] + offset + period * (np.arange(len(stored)) - middle)
    return times, float(excess)


def read_orbit(path: pathlib.Path, day: datetime) -> dict:
    """Reads an orbit's rows: their times (s from the day), positions, and excess."""
    with netCDF4.Dataset(path) as f:
        f.set_auto_maskandscale(False)
        stored = f["time"][:, CELLS // 2].astype(float)
        lat = f["lat"][:].astype(float)
        lon = f["lon"][:].astype(float)
    times, excess = reconstruct_times(stored)
    valid = (stored < 1e29) & np.all(np.abs(lat) <= 90, axis=1)
    return {
        "time": (times - (day - EPOCH).total_seconds())[valid],
        "lat": lat[valid],
        "lon": (lon[valid] + 180) % 360 - 180,
        "excess": excess,
    }


def load_day(day: datetime, directory: pathlib.Path) -> dict:
    """Downloads and reduces the rows of a day, or loads them from a cache."""
    cache = directory / f"ascat_{day:%Y%m%d}.npz"
    if not cache.exists():
        items = get_with_retries(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params={
                "short_name": SHORT_NAME,
                "temporal": f"{day - timedelta(hours=2):%Y-%m-%dT%H:%M:%SZ},"
                f"{day + timedelta(days=1):%Y-%m-%dT%H:%M:%SZ}",
                "page_size": 100,
            },
        ).json()["items"]
        urls = sorted(
            u["URL"]
            for item in items
            for u in item["umm"]["RelatedUrls"]
            if u["URL"].startswith("https://archive.podaac")
            and u["URL"].endswith(".nc")
        )
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        folder = directory / f"ascat_{day:%Y%m%d}"
        folder.mkdir(exist_ok=True)

        def download(url: str) -> pathlib.Path:
            path = folder / url.split("/")[-1]
            if not path.exists():
                response = session.get(url, timeout=600)
                response.raise_for_status()
                path.write_bytes(response.content)
            return path

        with concurrent.futures.ThreadPoolExecutor(4) as pool:
            paths = list(pool.map(download, urls))
        orbits = [read_orbit(path, day) for path in paths]
        for path, orbit in zip(paths, orbits):
            print(f"{path.name}: times within {orbit['excess']:.3f} s of rounding")
        data = {
            k: np.concatenate([o[k] for o in orbits]) for k in ("time", "lat", "lon")
        }
        # rows of the day, without those that consecutive orbits repeat
        _, first = np.unique(np.round(data["time"], 1), return_index=True)
        keep = first[(data["time"][first] >= 0) & (data["time"][first] < 86400)]
        np.savez_compressed(
            cache,
            excess=np.array([o["excess"] for o in orbits]),
            **{k: v[keep] for k, v in data.items()},
        )
    with np.load(cache, allow_pickle=False) as cached:
        data = {key: cached[key] for key in cached.files}
    print(
        f"{len(data['time'])} rows; times within {data['excess'].max():.3f} s of "
        "rounding (largest over orbits)",
        flush=True,
    )
    return data


def to_utc(day: datetime, seconds) -> list[datetime]:
    """Converts seconds from the start of a day to times."""
    return [day + timedelta(seconds=float(s)) for s in np.atleast_1d(seconds)]


def boundary(data: dict, rows: np.ndarray, cell: int, beside: int) -> np.ndarray:
    """Gets the Earth-fixed positions (m) half a cell beyond a cell, away from another."""
    near = ecef(data["lon"][rows, cell], data["lat"][rows, cell])
    far = ecef(data["lon"][rows, beside], data["lat"][rows, beside])
    lon, lat, _ = TO_GEODETIC.transform(*(near + (near - far) / 2).T)
    return ecef(lon, lat)


def look_angles(orbit, day, seconds, points) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the look angles across and along track (deg, positive left and forward) of
    points from the propagated orbit, in the Earth-fixed frame.
    """
    track = orbit.to_gp_orbit().get_orbit_track(to_utc(day, seconds))
    position, velocity = track.frame_xyz_and_velocity(itrs)
    position, velocity = np.array(position.m).T, np.array(velocity.m_per_s).T
    lon, lat, _ = TO_GEODETIC.transform(*position.T)
    nadir = -up(lon, lat)
    left = normalize(np.cross(velocity, nadir))
    forward = normalize(np.cross(nadir, left))
    los = normalize(points - position)
    down = np.sum(los * nadir, axis=-1)
    return (
        np.degrees(np.arctan2(np.sum(los * left, axis=-1), down)),
        np.degrees(np.arctan2(np.sum(los * forward, axis=-1), down)),
    )


def edge_rows(data: dict) -> np.ndarray:
    """Gets the rows whose edges are measured."""
    return np.arange(0, len(data["time"]), EDGE_STRIDE)


def side_labels(data: dict, orbit, day) -> tuple[str, str]:
    """Gets whether each side of the rows is left or right of the ground track."""
    rows = edge_rows(data)[:50]
    labels = []
    for side in SIDES:
        across, _ = look_angles(
            orbit, day, data["time"][rows], boundary(data, rows, *side["outer"])
        )
        labels.append("left" if np.median(across) > 0 else "right")
    if sorted(labels) != ["left", "right"]:
        raise RuntimeError(f"sides are not on either side of the track: {labels}")
    return tuple(labels)


def fit(data: dict, orbit, day) -> None:
    """Measures the look angles of the edges of each swath."""
    rows = edge_rows(data)
    labels = side_labels(data, orbit, day)
    for side, label in zip(SIDES, labels):
        medians = {}
        for edge, cells in side.items():
            across, along = look_angles(
                orbit, day, data["time"][rows], boundary(data, rows, *cells)
            )
            medians[edge] = float(np.median(across))
            print(
                f"{label} {edge}: {len(rows)} rows; look angle median "
                f"{medians[edge]:.4f} deg, p5 {np.percentile(across, 5):.4f}, p95 "
                f"{np.percentile(across, 95):.4f}; along track median "
                f"{np.median(along):.4f} deg, p5 {np.percentile(along, 5):.4f}, p95 "
                f"{np.percentile(along, 95):.4f}"
            )
        print(
            f"{label}: roll {(medians['inner'] + medians['outer']) / 2:.4f} deg, field "
            f"of view {abs(medians['outer'] - medians['inner']):.4f} deg"
        )


def score_edges(data: dict, orbit, day, satellite, labels) -> dict[str, np.ndarray]:
    """
    Gets the horizontal distances (km) between the measured edges and the rays at the
    edges of the views at the times of the rows.
    """
    rows = edge_rows(data)
    track = orbit.to_gp_orbit().get_orbit_track(to_utc(day, data["time"][rows]))
    distances = {}
    for side, label in zip(SIDES, labels):
        instrument = next(
            i for i in satellite.instruments if (i.roll_angle > 0) == (label == "left")
        )
        half = instrument.cross_track_field_of_view / 2
        sign = 1 if instrument.roll_angle > 0 else -1
        for edge, offset in (("inner", -sign * half), ("outer", sign * half)):
            ray = compute_projected_ray_position(
                track, 0, 0, instrument.roll_angle + offset
            )
            point = boundary(data, rows, *side[edge])
            vertical = normalize(point)
            difference = ecef(ray.longitude.degrees, ray.latitude.degrees) - point
            horizontal = (
                difference - np.sum(difference * vertical, axis=-1)[:, None] * vertical
            )
            distances[f"{label} {edge}"] = np.linalg.norm(horizontal, axis=-1) / 1e3
    return {edge: distances[edge] for edge in EDGES}


def find_segments(data: dict) -> list[np.ndarray]:
    """Gets the rows of segments of 360 s without gaps, below 70 deg latitude."""
    breaks = np.flatnonzero(np.diff(data["time"]) > ROW_GAP) + 1
    segments = []
    for run in np.split(np.arange(len(data["time"])), breaks):
        count = int((data["time"][run[-1]] - data["time"][run[0]]) // SEGMENT)
        for k in range(count):
            rows = run[
                (data["time"][run] >= data["time"][run[0]] + k * SEGMENT)
                & (data["time"][run] < data["time"][run[0]] + (k + 1) * SEGMENT)
            ]
            if np.all(np.abs(data["lat"][rows]) <= MAX_LATITUDE):
                segments.append(rows)
    return segments


def line_points(data: dict, rows: np.ndarray):
    """
    Gets points along lines across track through rows, and whether each is within a
    swath (between the edges of a side).
    """
    lons, lats, observed = [], [], []
    for row in rows:
        ends = [boundary(data, np.array([row]), *side["outer"])[0] for side in SIDES]
        (lon0, lon1), (lat0, lat1), _ = TO_GEODETIC.transform(*np.array(ends).T)
        azimuth, _, length = GEOD.inv(lon0, lat0, lon1, lat1)
        start_lon, start_lat, _ = GEOD.fwd(lon0, lat0, azimuth + 180, LINE_MARGIN)
        total = length + 2 * LINE_MARGIN
        steps = np.arange(0, total + LINE_STEP / 2, LINE_STEP)
        lon, lat, _ = GEOD.fwd(
            np.full(len(steps), start_lon),
            np.full(len(steps), start_lat),
            np.full(len(steps), azimuth),
            steps,
        )
        inside = np.zeros(len(steps), dtype=bool)
        for side in SIDES:
            limits = []
            for cells in (side["outer"], side["inner"]):
                point = boundary(data, np.array([row]), *cells)[0]
                plon, plat, _ = TO_GEODETIC.transform(*point)
                limits.append(GEOD.inv(start_lon, start_lat, plon, plat)[2])
            inside |= (steps >= min(limits)) & (steps <= max(limits))
        lons.append(lon)
        lats.append(lat)
        observed.append(inside)
    return np.concatenate(lons), np.concatenate(lats), np.concatenate(observed)


def score_segment(data: dict, rows: np.ndarray, day, satellite_for) -> dict:
    """Gets the intersection over union of a segment's swaths, by time step."""
    lon, lat, observed = line_points(data, rows)
    first, last = data["time"][rows[0]], data["time"][rows[-1]]
    results = {}
    for time_step in TIME_STEPS:
        n = max(1, int(round((last - first) / time_step)))
        step = (last - first) / n
        times = to_utc(day, first + step * (np.arange(n) + 0.5))
        satellite = satellite_for(step)
        geometry = unary_union(
            [
                compute_ground_track(
                    satellite, times, instrument_index=i
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


def observe_passes(points, data: dict, day) -> pd.DataFrame:
    """Finds the passes of lattice points within half a cell of a cell."""
    unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
    cells = up(data["lon"], data["lat"]).reshape(-1, 3)
    tree = cKDTree(cells)
    rows = []
    for point, nearby in enumerate(tree.query_ball_point(unit, CELL_CHORD)):
        if not nearby:
            continue
        # one pass for each run of rows, at the nearest cell of the run
        found = np.array(sorted(nearby))
        row = found // CELLS
        distance = np.linalg.norm(cells[found] - unit[point], axis=-1)
        runs = np.split(np.arange(len(found)), np.flatnonzero(np.diff(row) > 2) + 1)
        for run in runs:
            nearest = found[run[np.argmin(distance[run])]] // CELLS
            rows.append(
                {
                    "point_id": point,
                    "time": pd.Timestamp(day)
                    + pd.Timedelta(seconds=float(data["time"][nearest])),
                }
            )
    return pd.DataFrame(rows)


def validate(day, tle, entry, orbit, data, write) -> None:
    """Scores the configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, time_step=time_step, platform=PLATFORM)

    satellite = satellite_for(TIME_STEPS[0])
    labels = side_labels(data, orbit, day)
    edges = score_edges(data, orbit, day, satellite, labels)
    for edge, values in edges.items():
        print(
            f"swath edge, {edge}: {len(values)} rows, median {np.median(values):.2f} "
            f"km, p95 {np.percentile(values, 95):.2f} km",
            flush=True,
        )
    lattice = generate_points_fibonacci_lattice(POINT_SPACING)
    observed = observe_passes(lattice, data, day)
    ids = np.sort(observed.point_id.unique())
    subset = lattice.iloc[ids].reset_index(drop=True)
    observed["point_id"] = np.searchsorted(ids, observed.point_id)
    start = to_utc(day, data["time"][0])[0] - SEARCH_WINDOW
    end = to_utc(day, data["time"][-1])[0] + SEARCH_WINDOW
    print(f"predicting passes of {len(subset)} points; {len(observed)} observed")
    predicted = predict_passes(subset, satellite, None, start, end, [])
    results = match_passes(predicted, observed)
    print(json.dumps(results), flush=True)
    segments = find_segments(data)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        area = pd.DataFrame(
            list(
                pool.map(
                    lambda rows: score_segment(data, rows, day, satellite_for),
                    segments,
                )
            )
        )
    for time_step in TIME_STEPS:
        print(
            f"swath IoU at {time_step} s over {len(area)} segments: median "
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
        "ascat",
        CONFIGURATION,
        {
            "digest": entry.configurations[CONFIGURATION].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_ascat.py",
            "reference": f"{SHORT_NAME} (wind vector cell positions and row times)",
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
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    data = load_day(day, arguments.directory)
    if arguments.day in ELEMENT_SETS:
        tle = ELEMENT_SETS[arguments.day]
    else:
        tle = read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    if arguments.fit:
        fit(data, orbit, day)
        return
    validate(day, tle, get_entry("ascat"), orbit, data, arguments.write)


if __name__ == "__main__":
    main()
