"""
Validates the TROPOMI catalog entry against the pixel bounds of Sentinel-5P's aerosol
index product, whose 5.5 km by 3.5 km pixels are those of the instrument's ultraviolet
and visible bands, using the exact values of the entry through public TAT-C methods.

* The swath width is compared with the distance between the ends of each line.
* The swath area is compared with the union of the lines over one sunlit arc of an
  orbit. TROPOMI's product covers the sunlit part of each orbit only.
* Passes are those of lattice points inside the lines, and their times are those of
  the lines, interpolated to the points.

TROPOMI's lines curve aft toward the swath edges, which it images about 20 s after
nadir, in a shape that does not change with latitude in the Earth-fixed velocity frame.
No single rectangular view can represent these times, so the entry is a union of views,
each rolled and pitched to a part of the line. Their angles come from `--fit`, on a day
of data, from the satellite positions recorded in the product. Running on a day at
least 30 days from it tests them on data held out from the fit. The tests need an
element set near the data, and none is used for the fit.

Needs the `validation` extra:

    python docs/validation/validate_catalog_tropomi.py <data directory> 2026-10-05 --fit
    python docs/validation/validate_catalog_tropomi.py <data directory> 2026-11-05
    python docs/validation/validate_catalog_tropomi.py <data directory> 2026-11-10

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
import re
from datetime import datetime, timedelta, timezone

import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import tatc
from tatc.catalog import get_entry
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from tatc.utils import compute_projected_ray_position
from validation_io import (
    RemoteFile,
    build_grid,
    ecef,
    get_holdout_gap,
    get_with_retries,
    in_lattice_cells,
    match_passes,
    normalize,
    predict_passes,
    recompute_passes,
    score_swath_area,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "tropomi.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "tropomi.npz"

BUCKET = "https://meeo-s5p.s3.amazonaws.com/"
PREFIX = "NRTI/L2__AER_AI/"
FIT_DAY = datetime(2026, 10, 5, tzinfo=timezone.utc)
PLATFORM = "sentinel-5p"
ELEMENT_SETS = {
    "2026-11-04": [
        "1 42969U 17064A   26311.28829078  .00000048  00000+0  44053-4 0  9996",
        "2 42969  98.8082 252.4800 0000829  99.6231 260.5039 14.19499222469868",
    ],
    "2026-11-05": [
        "1 42969U 17064A   26311.28829078  .00000048  00000+0  44053-4 0  9996",
        "2 42969  98.8082 252.4800 0000829  99.6231 260.5039 14.19499222469868",
    ],
    # a fresh held-out day for swept swaths, not used before
    "2026-11-10": [
        "1 42969U 17064A   26315.30604670  .00000068  00000+0  53520-4 0  9993",
        "2 42969  98.8091 256.4758 0000925 101.4950 258.6329 14.19500041470438",
    ],
}

KEY = "swath"
N_ROWS = 450
NODES = np.r_[np.arange(0, N_ROWS, 16), N_ROWS]
LINE_PERIOD = 0.84
ANALYSIS_TIME_STEP = 10
TIME_STEPS = (ANALYSIS_TIME_STEP, 30, 120)
VIEW_TOLERANCE = 2.0  # s, within a view, from its center to the curve of the line
VIEW_OVERLAP = 0.04  # deg, so that adjacent views leave no gaps
POINT_SPACING = 500e3
ORBIT_MARGIN = 120
AREA_RADIUS = 2 * np.sin(300 / 6371 / 2)  # chord, on the unit sphere
WIDTH_LINE_STRIDE = 20
# lines of the cells that test whether points are observed, which are then about as
# long as they are wide, so that the nearest vertices of a point are on its cell
CELL_LINE_STRIDE = 10
CELL_NEIGHBORS = 8
FIT_LINE_STRIDE = 10
FIXTURE_PASSES = 200
EARTH_ROTATION = np.array([0, 0, 7.2921150e-5])
PRODUCT = "PRODUCT/"
GEOLOCATION = PRODUCT + "SUPPORT_DATA/GEOLOCATIONS/"


def list_granules(day: datetime) -> list[str]:
    """Lists the keys of a day's granules."""
    prefix = f"{PREFIX}{day:%Y/%m/%d}/"
    keys: list[str] = []
    token = None
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            params["continuation-token"] = token
        listing = get_with_retries(BUCKET, params=params).text
        keys += re.findall(r"<Key>(.*?)</Key>", listing)
        following = re.search(r"<NextContinuationToken>(.*?)<", listing)
        if following is None:
            return sorted(key for key in keys if key.endswith(".nc"))
        token = following.group(1)


def reduce_granule(key: str, directory: pathlib.Path) -> pathlib.Path:
    """
    Reduces a granule to its line times, satellite positions, and the forward corners
    of every sixteenth pixel of its lines, or finds it in a cache.
    """
    path = directory / (pathlib.PurePosixPath(key).stem + ".npz")
    if path.exists():
        return path
    with h5py.File(RemoteFile(BUCKET + key, block_size=2**20), "r") as f:
        times = pd.to_datetime(
            [t.decode() for t in f[PRODUCT + "time_utc"][0]], utc=True, format="ISO8601"
        )

        def nodes(name: str) -> np.ndarray:
            bounds = f[GEOLOCATION + name][0]
            # bounds run around each pixel, with the first two on the side of the
            # previous line, in the order of the pixels
            starts = bounds[:, NODES[:-1], 0]
            return np.concatenate([starts, bounds[:, [N_ROWS - 1], 1]], axis=1)

        arrays = {
            "time": times.as_unit("ms").astype("int64").to_numpy(),
            "node_lat": nodes("latitude_bounds"),
            "node_lon": nodes("longitude_bounds"),
            "satellite_lat": f[GEOLOCATION + "satellite_latitude"][0],
            "satellite_lon": f[GEOLOCATION + "satellite_longitude"][0],
            "satellite_alt": f[GEOLOCATION + "satellite_altitude"][0],
        }
    # written whole, so that an interrupted download leaves no partial granule
    temporary = path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)
    return path


def load_day(day: datetime, directory: pathlib.Path) -> dict:
    """Downloads and reduces the geolocation of a day, or loads it from a cache."""
    path = directory / f"tropomi_{day:%Y%m%d}.npz"
    if not path.exists():
        keys = list_granules(day)
        granules = directory / f"tropomi_{day:%Y%m%d}"
        granules.mkdir(parents=True, exist_ok=True)
        with concurrent.futures.ThreadPoolExecutor(6) as pool:
            paths = list(pool.map(lambda key: reduce_granule(key, granules), keys))
        arrays: dict[str, list] = {}
        for granule in sorted(paths):
            with np.load(granule) as data:
                for name in data.files:
                    arrays.setdefault(name, []).append(data[name])
        print(f"{len(keys)} granules", flush=True)
        np.savez_compressed(path, **{k: np.concatenate(v) for k, v in arrays.items()})
    return dict(np.load(path))


class Reference:
    """A day of TROPOMI lines, with what is derived from them."""

    def __init__(self, data: dict[str, np.ndarray]):
        good = (
            np.all(np.isfinite(data["node_lat"]), axis=1)
            & np.all(np.isfinite(data["node_lon"]), axis=1)
            & (np.abs(data["node_lat"]).max(axis=1) <= 90)
            & (data["satellite_alt"] > 7e5)
        )
        self.n_dropped = int(np.sum(~good))
        self.milliseconds = data["time"][good]
        self.time = pd.to_datetime(self.milliseconds, unit="ms", utc=True)
        self.lat = data["node_lat"][good].astype(float)
        self.lon = data["node_lon"][good].astype(float)
        self.altitude = data["satellite_alt"][good].astype(float)
        self.satellite = ecef(
            data["satellite_lon"][good].astype(float),
            data["satellite_lat"][good].astype(float),
            self.altitude,
        )
        self.satellite_up = up(
            data["satellite_lon"][good].astype(float),
            data["satellite_lat"][good].astype(float),
        )
        self.unit = up(self.lon, self.lat)
        self.xyz = ecef(self.lon, self.lat)
        self.seconds = (self.time - self.time[0]).total_seconds().to_numpy()
        period = np.diff(self.seconds)
        self.line_period = float(np.median(period))
        self.is_next_contiguous = period < 1.5 * LINE_PERIOD
        ends = np.nonzero(~self.is_next_contiguous)[0]
        self.arcs = list(
            zip(np.r_[0, ends + 1], np.r_[ends, len(self.seconds) - 1], strict=True)
        )
        self.gaps = [
            (self.time[i], self.time[i + 1]) for i in np.nonzero(period > 60)[0]
        ]

    def in_cells(self, points: np.ndarray, first: int, last: int) -> np.ndarray:
        """
        Gets whether unit vectors lie in cells between every tenth line of a range of
        lines, with the last line of the range, and between lines with gaps.
        """
        lines = np.r_[np.arange(first, last, CELL_LINE_STRIDE), last]
        contiguous = np.array(
            [
                bool(np.all(self.is_next_contiguous[a:b]))
                for a, b in zip(lines[:-1], lines[1:], strict=True)
            ]
        )
        return in_lattice_cells(
            self.unit[lines], contiguous, points, neighbors=CELL_NEIGHBORS
        )

    def observe_passes(self, points) -> pd.DataFrame:
        """
        Finds the passes on which TROPOMI observed points inside the cells between its
        lines, with the times of the lines interpolated to the points. A line's time is
        that of the middle of the line, half a period after its first edge.
        """
        unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        point_xyz = ecef(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        rows = []
        n_nodes = len(NODES)
        for first, last in self.arcs:
            if last - first < 2:
                continue
            cells = self.unit[first : last + 1]
            found = np.nonzero(self.in_cells(unit, first, last))[0]
            if len(found) == 0:
                continue
            _, nearest = cKDTree(cells.reshape(-1, 3)).query(unit[found], workers=-1)
            lines, columns = np.divmod(nearest, n_nodes)
            lines = np.minimum(lines, last - first - 1)
            a = self.xyz[first + lines, columns]
            b = self.xyz[first + lines + 1, columns]
            fraction = np.sum((point_xyz[found] - a) * (b - a), axis=-1) / np.sum(
                (b - a) ** 2, axis=-1
            )
            seconds = (
                self.seconds[first + lines]
                + fraction * self.line_period
                - self.line_period / 2
            )
            rows.append(
                pd.DataFrame(
                    {
                        "point_id": found,
                        "time": self.time[0] + pd.to_timedelta(seconds, unit="s"),
                    }
                )
            )
        return pd.concat(rows, ignore_index=True).astype({"point_id": "int64"})


def measure_rays(reference: Reference, frame: str = "earth_fixed") -> dict:
    """
    Measures the directions to the middles of lines from the recorded satellite
    positions, as tangents of the along-track and cross-track (left) view angles.
    """
    regular = np.nonzero(
        reference.is_next_contiguous[:-1] & reference.is_next_contiguous[1:]
    )[0]
    lines = regular[::FIT_LINE_STRIDE] + 1
    position = reference.satellite[lines]
    velocity = (reference.satellite[lines + 1] - reference.satellite[lines - 1]) / (
        (reference.seconds[lines + 1] - reference.seconds[lines - 1])[:, None]
    )
    if frame == "inertial":
        velocity = velocity + np.cross(EARTH_ROTATION, position)
    zenith = reference.satellite_up[lines]
    along = normalize(velocity - np.sum(velocity * zenith, axis=-1)[:, None] * zenith)
    left = np.cross(zenith, along)
    nodes = (reference.xyz[lines] + reference.xyz[lines + 1]) / 2
    rays = nodes - position[:, None]
    down = -np.sum(rays * zenith[:, None], axis=-1)
    center = len(NODES) // 2
    advance = np.linalg.norm(
        reference.xyz[lines + 1, center] - reference.xyz[lines, center], axis=-1
    ) / (reference.seconds[lines + 1] - reference.seconds[lines])
    return {
        "latitude": reference.lat[lines, center],
        "along": np.sum(rays * along[:, None], axis=-1) / down,
        "cross": np.sum(rays * left[:, None], axis=-1) / down,
        "down": down,
        "ground_speed": float(np.median(advance)),
    }


def to_ray(along: float, cross: float) -> np.ndarray:
    """Gets the unit vector of tangents, in along-track, left, and nadir coordinates."""
    return normalize(np.array([along, cross, 1.0]))


def design_views(rays: dict) -> list[dict]:
    """
    Divides the swath into views, each as wide as it can be with the line within a
    tolerance of its center, in time. Each is centered between its ends: a roll about
    the along-track axis, then a pitch about the rolled cross-track axis.
    """
    along = np.median(rays["along"], axis=0)
    cross = np.median(rays["cross"], axis=0)
    down = np.median(rays["down"], axis=0)

    def build(i: int, j: int) -> tuple[dict, float]:
        boresight = normalize(to_ray(along[i], cross[i]) + to_ray(along[j], cross[j]))
        roll = np.arctan2(boresight[1], boresight[2])
        pitch = np.arcsin(boresight[0])
        cross_axis = np.array([0, np.cos(roll), -np.sin(roll)])
        along_axis = np.cos(pitch) * np.array([1.0, 0, 0]) - np.sin(pitch) * np.array(
            [0, np.sin(roll), np.cos(roll)]
        )
        end = to_ray(along[i], cross[i])
        half = np.arctan2(np.abs(end @ cross_axis), end @ boresight)
        worst = 0.0
        for k in range(i, j + 1):
            ray = to_ray(along[k], cross[k])
            # the distance along track, at the ray's range
            worst = max(worst, abs(ray @ along_axis) * down[k] / ray[2])
        return (
            {
                "first": int(NODES[i]),
                "last": int(NODES[j]),
                "roll": float(np.degrees(roll)),
                "pitch": float(np.degrees(pitch)),
                "field_of_view": float(np.degrees(2 * half)) + VIEW_OVERLAP,
            },
            worst / rays["ground_speed"],
        )

    views = []
    i = 0
    while i < len(NODES) - 1:
        j = i + 1
        while j < len(NODES) - 1 and build(i, j + 1)[1] <= VIEW_TOLERANCE:
            j += 1
        views.append(build(i, j)[0] | {"largest_time": build(i, j)[1]})
        i = j
    return views


def fit_geometry(reference: Reference) -> None:
    """Measures the swath's view angles and the shape of its lines."""
    print(f"line period {reference.line_period:.4f} s")
    for frame in ("earth_fixed", "inertial"):
        rays = measure_rays(reference, frame)
        along = np.degrees(np.arctan(rays["along"]))
        cross = np.degrees(np.arctan(rays["cross"]))
        latitude = rays["latitude"]
        first, last = np.median(cross[:, 0]), np.median(cross[:, -1])
        print(
            f"{frame}: {len(latitude)} lines, ends at {first:.4f} and {last:.4f} deg, "
            f"along {np.median(along[:, 0]):.4f} and {np.median(along[:, -1]):.4f} deg; "
            f"ground speed {rays['ground_speed']:.0f} m/s"
        )
        for band in ((-90, -30), (-30, 30), (30, 90)):
            chosen = (latitude >= band[0]) & (latitude < band[1])
            print(
                f"  latitudes {band}: along at the ends "
                f"{np.median(along[chosen, 0]):+.3f} and "
                f"{np.median(along[chosen, -1]):+.3f} deg, center "
                f"{np.median(along[chosen, len(NODES) // 2]):+.3f} deg"
            )
    views = design_views(measure_rays(reference))
    print(f"{len(views)} views, tolerance {VIEW_TOLERANCE} s:")
    for view in views:
        print(
            f"  pixels {view['first']:3d} to {view['last']:3d}: roll "
            f"{view['roll']:+.4f}, pitch {view['pitch']:+.4f}, field of view "
            f"{view['field_of_view']:.4f} deg, line within {view['largest_time']:.2f} s"
        )
    print(json.dumps(views))


def score_swath_width(reference: Reference, satellite) -> float:
    """Gets the median absolute difference (km) from the distance between line ends."""
    n_lines = len(reference.seconds)
    lines = np.arange(0, n_lines - 1, WIDTH_LINE_STRIDE)
    middle = (reference.xyz[lines] + reference.xyz[lines + 1]) / 2
    observed = np.linalg.norm(middle[:, 0] - middle[:, -1], axis=-1)
    track = satellite.orbit.to_gp_orbit().get_orbit_track(
        [t.floor("us").to_pydatetime() for t in reference.time[lines]]
    )
    ends = []
    # the first view is on the left of the motion, and the last on the right; the
    # views are given a small along-track extent because, before a fix in the same
    # release, a ray at 180 deg of a rectangular view without one projected to its
    # center (kept so that the script reproduces the recorded results)
    for view, angle in (
        (satellite.instruments[0], 0),
        (satellite.instruments[-1], 180),
    ):
        position = compute_projected_ray_position(
            track,
            view.cross_track_field_of_view,
            1e-3,
            view.roll_angle,
            view.pitch_angle,
            True,
            angle,
            velocity_frame=view.velocity_frame,
        )
        ends.append(ecef(position.longitude.degrees, position.latitude.degrees))
    predicted = np.linalg.norm(ends[0] - ends[1], axis=-1)
    print(
        f"swath width: observed {np.median(observed) / 1e3:.2f} km, "
        f"predicted {np.median(predicted) / 1e3:.2f} km",
        flush=True,
    )
    return float(np.median(np.abs(predicted - observed))) / 1e3


def observe_arc(reference: Reference, grid):
    """Gets the cells compared over the sunlit arc nearest noon, and those observed."""
    lon, lat, _ = grid
    unit = up(lon, lat)
    noon = reference.time[0].floor("D") + pd.Timedelta(hours=12)
    lengths = [
        (abs((reference.time[(a + b) // 2] - noon).total_seconds()), a, b)
        for a, b in reference.arcs
        if b - a > 1000
    ]
    _, first, last = min(lengths)
    cells = reference.unit[first : last + 1]
    distance, nearest = cKDTree(cells.reshape(-1, 3)).query(unit, workers=-1)
    elapsed = (
        reference.seconds[first + nearest // len(NODES)] - reference.seconds[first]
    )
    compared = (
        (distance < AREA_RADIUS)
        & (elapsed > ORBIT_MARGIN)
        & (elapsed < reference.seconds[last] - reference.seconds[first] - ORBIT_MARGIN)
    )
    observed = np.zeros(len(unit), dtype=bool)
    observed[compared] = reference.in_cells(unit[compared], first, last)
    print(
        f"arc of {reference.seconds[last] - reference.seconds[first]:.0f} s from "
        f"{reference.time[first]}: {int(compared.sum())} cells compared, "
        f"{int(observed.sum())} observed",
        flush=True,
    )
    return (
        reference.time[first].floor("us").to_pydatetime(),
        reference.time[last].floor("us").to_pydatetime(),
        compared,
        observed,
    )


def validate(day, tle, entry, orbit, reference, write) -> None:
    """Scores the configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, KEY, time_step=time_step, platform=PLATFORM)

    satellite = satellite_for(ANALYSIS_TIME_STEP)
    points = generate_points_fibonacci_lattice(POINT_SPACING)
    observed = reference.observe_passes(points)
    start = reference.time[0].floor("us").to_pydatetime() - timedelta(minutes=1)
    end = reference.time[-1].floor("us").to_pydatetime() + timedelta(minutes=1)
    predicted = predict_passes(points, satellite, None, start, end, reference.gaps)
    results = match_passes(predicted, observed)
    results["swath_width_median"] = score_swath_width(reference, satellite)
    print(json.dumps(results), flush=True)
    grid = build_grid()
    first, last, compared, covered = observe_arc(reference, grid)
    results["footprint_iou_median"] = score_swath_area(
        satellite_for, None, first, last, TIME_STEPS, grid, compared, covered
    )
    print(json.dumps(results), flush=True)
    if not write:
        return
    rng = np.random.default_rng(1)
    picks = observed.iloc[rng.choice(len(observed), FIXTURE_PASSES, replace=False)]
    passes = np.stack(
        [
            points.geometry.y.to_numpy()[picks.point_id],
            points.geometry.x.to_numpy()[picks.point_id],
            (picks.time - pd.Timestamp(day)).dt.total_seconds().to_numpy(),
        ],
        axis=-1,
    )
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=passes
    )
    names = (
        "pass_recall",
        "epoch_p95",
        "epoch_median",
        "swath_width_median",
        "footprint_iou_median",
    )
    update_record(
        RECORD_PATH,
        "tropomi",
        KEY,
        {
            "digest": entry.configurations[KEY].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_tropomi.py",
            "reference": "Sentinel-5P TROPOMI aerosol index product (near real time)",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": [
                {"name": name, "value": float(results[name])} for name in names
            ],
            "fixture_metrics": recompute_passes(satellite, None, day, passes),
        },
    )
    print(f"wrote {RECORD_PATH}")


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
    reference = Reference(load_day(day, arguments.directory))
    print(
        f"{arguments.day}: {len(reference.seconds)} lines, "
        f"{reference.n_dropped} dropped, {len(reference.arcs)} arcs, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        fit_geometry(reference)
        return
    tle = ELEMENT_SETS[arguments.day]
    validate(
        day,
        tle,
        get_entry("tropomi"),
        GeneralPerturbationsOrbit.from_tle(tle),
        reference,
        arguments.write,
    )


if __name__ == "__main__":
    main()
