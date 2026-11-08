"""
Validates the OLI catalog entry against the scenes of Landsat 9 and Landsat 8, using the
exact values of the entry through public TAT-C methods. Scenes are frames of a swath
whose ends are arbitrary, so consecutive scenes are compared together:

* The swath width is compared with the width of each scene through its center.
* The swath area is compared with the union of consecutive scenes, from the center of
  the first to the center of the last, so that the arbitrary frame ends are left out.
* Passes are those of lattice points inside the scenes.
* Epochs are compared with the imaging times of points across the swath, from the line
  times and polynomial coefficients of the scenes' angle files. Those files are only
  available about a week after acquisition, so on days without them, epochs are
  compared at the centers of scenes only.

The entry's angles and velocity frame come from `--fit`, on a day of Landsat 9 scenes.
Running on a day at least 30 days from it tests them on data held out from the fit. As
element sets for the fit's period came first, the holdout day precedes the fit day. On
Landsat 8, results are reported as a transfer to a sibling platform, under metric names
that do not count toward the entry's validation.

Needs the `validation` extra:

    python docs/validation/validate_catalog_oli.py <data directory> 2026-11-05 --fit
    python docs/validation/validate_catalog_oli.py <data directory> 2026-10-03
    python docs/validation/validate_catalog_oli.py <data directory> 2026-10-03 \
        --platform landsat-8

With `--write`, it writes the entry's validation record and, on Landsat 9, the reduced
fixture from which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
import re
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, shape
from shapely.ops import unary_union
from skyfield.framelib import itrs

import tatc
from tatc.analysis import collect_observations, compute_ground_track
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit, Point
from validation_io import (
    TO_GEODETIC,
    add_transfer_metrics,
    ecef,
    get_holdout_gap,
    get_with_retries,
    match_passes,
    normalize,
    post_with_retries,
    predict_passes,
    recompute_passes,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "oli.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "oli.npz"

SCENES_URL = "https://landsatlook.usgs.gov/stac-server/search"
ANGLES_URL = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
ANGLES_TOKEN_URL = (
    "https://planetarycomputer.microsoft.com/api/sas/v1/token/landsat-c2-l2"
)

FIT_DAY = datetime(2026, 11, 5, tzinfo=timezone.utc)
FIT_PLATFORM = "landsat-9"
PLATFORMS = {"landsat-8": "LANDSAT_8", "landsat-9": "LANDSAT_9"}
_OCTOBER = {
    "landsat-8": [
        "1 39084U 13008A   26277.62187981  .00000150  00000+0  43377-4 0  9995",
        "2 39084  98.2196 346.3945 0001312  94.4520 265.6828 14.57108036713830",
    ],
    "landsat-9": [
        "1 49260U 21088A   26277.17558559  .00000158  00000+0  45237-4 0  9997",
        "2 49260  98.2176 345.9720 0001438  89.7555 270.3809 14.57106411266879",
    ],
}
_NOVEMBER = {
    "landsat-8": [
        "1 39084U 13008A   26311.27004377  .00000333  00000+0  84030-4 0  9996",
        "2 39084  98.2158  19.5861 0001266  90.4879 269.6465 14.57106512718737",
    ],
    "landsat-9": [
        "1 49260U 21088A   26311.16701518  .00000257  00000+0  67047-4 0  9993",
        "2 49260  98.2139  19.4951 0001142  86.0967 274.0362 14.57106672271821",
    ],
}
ELEMENT_SETS = {
    "2026-10-02": _OCTOBER,
    "2026-10-03": _OCTOBER,
    "2026-11-05": _NOVEMBER,
}

KEY = "swath"
BAND = "BAND04"  # the red band, as for Sentinel-2 MSI
ANGLE_SCENE_STRIDE = 5
ALONG_FRACTIONS = (0.25, 0.5, 0.75)
CROSS_FRACTIONS = np.linspace(0.02, 0.98, 9)
MAX_SCENE_GAP = (
    30  # s, between consecutive scenes of a pass (they are about 24 s apart)
)
CHUNK_SCENES = 20
MIN_CHUNK_SCENES = 3
TIME_STEPS = (10, 30, 120)
POINT_SPACING = 100e3
SCENE_RADIUS = 2 * np.sin(160 / 6371 / 2)  # chord, on the unit sphere
SEARCH_WINDOW = timedelta(minutes=10)
FIXTURE_PASSES = 300
EARTH_ROTATION = np.array([0, 0, 7.2921150e-5])


def load_scenes(platform: str, day: datetime, directory: pathlib.Path) -> pd.DataFrame:
    """Gets the Level-1 OLI scenes acquired on a day, viewing nadir, in time order."""
    path = directory / f"scenes_{platform}_{day:%Y%m%d}.json"
    if not path.exists():
        body = {
            "collections": ["landsat-c2l1"],
            "datetime": f"{day:%Y-%m-%d}T00:00:00Z/{day:%Y-%m-%d}T23:59:59.999Z",
            "limit": 500,
            "query": {"platform": {"eq": PLATFORMS[platform]}},
            "fields": {
                "include": [
                    "id",
                    "geometry",
                    "properties.datetime",
                    "properties.view:off_nadir",
                    "properties.landsat:wrs_path",
                    "properties.landsat:wrs_row",
                ],
                "exclude": ["links", "assets"],
            },
        }
        features = []
        while True:
            page = post_with_retries(SCENES_URL, body).json()
            features += page["features"]
            following = [
                link for link in page.get("links", []) if link.get("rel") == "next"
            ]
            if not following or not page["features"]:
                break
            body = following[0]["body"]
        path.write_text(json.dumps(features), encoding="utf-8")
    features = json.loads(path.read_text(encoding="utf-8"))
    rows = [
        {
            "id": feature["id"],
            "time": feature["properties"]["datetime"],
            "off_nadir": feature["properties"].get("view:off_nadir") or 0,
            "path": feature["properties"]["landsat:wrs_path"],
            "row": feature["properties"]["landsat:wrs_row"],
            "geometry": shape(feature["geometry"]),
        }
        for feature in features
        # LT scenes are of the thermal instrument alone
        if feature["id"][:2] in ("LC", "LO")
    ]
    scenes = pd.DataFrame(rows)
    scenes["time"] = pd.to_datetime(scenes.time, utc=True, format="ISO8601")
    n_off = int((scenes.off_nadir > 0).sum())
    scenes = scenes[scenes.off_nadir == 0].sort_values("time").reset_index(drop=True)
    print(f"{len(scenes)} scenes viewing nadir, {n_off} off nadir left out", flush=True)
    return scenes


def frames_at(orbit, times, frame: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gets a satellite's Earth-fixed position, and along-track and left directions."""
    track = orbit.to_gp_orbit().get_orbit_track(
        [time.floor("us").to_pydatetime() for time in pd.DatetimeIndex(times)]
    )
    position, velocity = track.frame_xyz_and_velocity(itrs)
    position = position.m.T
    velocity = velocity.m_per_s.T
    if frame == "inertial":
        velocity = velocity + np.cross(EARTH_ROTATION, position)
    lon, lat, _ = TO_GEODETIC.transform(*position.T)
    nadir_up = up(lon, lat)
    along = normalize(
        velocity - np.sum(velocity * nadir_up, axis=-1)[:, None] * nadir_up
    )
    return position, along, np.cross(nadir_up, along)


def outline_xyz(geometry) -> np.ndarray:
    """Gets the Earth-fixed vertices of a footprint, which may be split in longitude."""
    parts = getattr(geometry, "geoms", [geometry])
    lon_lat = np.concatenate([np.array(part.exterior.coords)[:, :2] for part in parts])
    return ecef(lon_lat[:, 0], lon_lat[:, 1])


def add_corners(scenes: pd.DataFrame, orbit) -> np.ndarray:
    """
    Gets the corners of scenes as rows of forward left, forward right, aft right, and aft
    left. Vertices that split a footprint in longitude lie between its corners, which
    are its extremes along the diagonals of the along-track and left directions.
    """
    _, along, left = frames_at(orbit, scenes.time, "earth_fixed")
    corners = np.empty((len(scenes), 4, 3))
    for k, geometry in enumerate(scenes.geometry):
        xyz = outline_xyz(geometry)
        a, c = xyz @ along[k], xyz @ left[k]
        corners[k] = xyz[
            [
                np.argmax(a + c),
                np.argmax(a - c),
                np.argmin(a + c),
                np.argmin(a - c),
            ]
        ]
    return corners


def get_intervals(scenes: pd.DataFrame) -> np.ndarray:
    """Numbers the passes to which scenes belong, from the times between them."""
    gaps = scenes.time.diff().dt.total_seconds().fillna(0).to_numpy()
    return np.cumsum(gaps > MAX_SCENE_GAP)


def fit_geometry(scenes, corners, orbit) -> None:
    """Measures the swath's view angles and the orientation of its lines."""
    left_middle = (corners[:, 0] + corners[:, 3]) / 2
    right_middle = (corners[:, 1] + corners[:, 2]) / 2
    lon, lat, _ = TO_GEODETIC.transform(*((left_middle + right_middle) / 2).T)
    latitude = np.asarray(lat)
    for frame in ("earth_fixed", "inertial"):
        position, along, left = frames_at(orbit, scenes.time, frame)
        lon, lat, _ = TO_GEODETIC.transform(*position.T)
        nadir = -up(lon, lat)

        def view(point):
            line = point - position
            return np.degrees(
                np.arctan2(np.sum(line * left, axis=-1), np.sum(line * nadir, axis=-1))
            )

        to_left, to_right = view(left_middle), -view(right_middle)
        # the orientation of the ends of scenes, from their forward and aft edges
        yaws = []
        for first, second in ((1, 0), (2, 3)):
            edge = corners[:, second] - corners[:, first]
            yaws.append(
                np.degrees(
                    np.arctan2(
                        np.sum(edge * along, axis=-1), np.sum(edge * left, axis=-1)
                    )
                )
            )
        yaw = np.concatenate(yaws)
        print(
            f"{frame}: left {np.median(to_left):.4f} deg, right "
            f"{np.median(to_right):.4f} deg, field of view "
            f"{np.median(to_left + to_right):.4f} deg, roll "
            f"{np.median(to_left - to_right) / 2:.4f} deg; lines forward at the left by "
            f"{np.median(yaw):.3f} deg (p5 {np.percentile(yaw, 5):.3f}, "
            f"p95 {np.percentile(yaw, 95):.3f})"
        )
    width = np.linalg.norm(left_middle - right_middle, axis=-1) / 1e3
    print(
        f"scene width: median {np.median(width):.2f} km, p5 "
        f"{np.percentile(width, 5):.2f}, p95 {np.percentile(width, 95):.2f}; "
        f"latitudes {latitude.min():.1f} to {latitude.max():.1f}"
    )


def score_swath_width(scenes, corners, orbit, instrument) -> float:
    """Gets the median absolute difference (km) from the width through scene centers."""
    observed = np.linalg.norm(
        (corners[:, 0] + corners[:, 3]) / 2 - (corners[:, 1] + corners[:, 2]) / 2,
        axis=-1,
    )
    position, _, _ = frames_at(orbit, scenes.time, "earth_fixed")
    altitude = TO_GEODETIC.transform(*position.T)[2]
    predicted = np.vectorize(instrument.get_swath_width)(altitude)
    print(
        f"swath width: observed {np.median(observed) / 1e3:.2f} km, "
        f"predicted {np.median(predicted) / 1e3:.2f} km",
        flush=True,
    )
    return float(np.median(np.abs(predicted - observed))) / 1e3


def inside_convex(unit: np.ndarray, corners: np.ndarray) -> np.ndarray:
    """Checks whether unit vectors lie within a convex outline of ordered corners."""
    outline = normalize(corners)
    normals = np.cross(outline, np.roll(outline, -1, axis=0))
    sides = unit @ normals.T
    return np.all(sides >= 0, axis=1) | np.all(sides <= 0, axis=1)


def observe_passes(points, scenes, corners, intervals) -> pd.DataFrame:
    """Finds the passes of points inside scenes, one per point and sequence of scenes."""
    lon, lat = points.geometry.x.to_numpy(), points.geometry.y.to_numpy()
    # geocentric directions, as for the corners
    unit = normalize(ecef(lon, lat))
    tree = cKDTree(unit)
    centers = normalize(corners.mean(axis=1))
    rows = []
    for k, nearby in enumerate(tree.query_ball_point(centers, SCENE_RADIUS)):
        nearby = np.array(nearby, dtype=int)
        if len(nearby) == 0:
            continue
        found = nearby[inside_convex(unit[nearby], corners[k])]
        rows += [
            {"point_id": int(i), "interval": intervals[k], "time": scenes.time[k]}
            for i in found
        ]
    observed = pd.DataFrame(rows).drop_duplicates(["point_id", "interval"])
    return observed.reset_index(drop=True)


def score_passes(points, satellite, scenes, corners, intervals) -> dict:
    """Compares the passes predicted with those of lattice points inside scenes."""
    observed = observe_passes(points, scenes, corners, intervals)
    ids = np.sort(observed.point_id.unique())
    observed["point_id"] = np.searchsorted(ids, observed.point_id)
    subset = points.iloc[ids].reset_index(drop=True)
    start = scenes.time.iloc[0].floor("us").to_pydatetime() - SEARCH_WINDOW
    end = scenes.time.iloc[-1].floor("us").to_pydatetime() + SEARCH_WINDOW
    print(
        f"predicting passes of {len(subset)} points; {len(observed)} observed",
        flush=True,
    )
    predicted = predict_passes(subset, satellite, 0, start, end, [])
    return match_passes(predicted, observed)


def to_local(geometry, transformer):
    """Transforms a geometry to a local projection, as one valid geometry."""
    parts = [
        Polygon(
            np.column_stack(
                transformer.transform(*np.array(part.exterior.coords).T[:2])
            )
        ).buffer(0)
        for part in getattr(geometry, "geoms", [geometry])
    ]
    return unary_union(parts)


def get_chunks(intervals: np.ndarray) -> list[np.ndarray]:
    """Divides the sequences of scenes into chunks of comparable lengths."""
    chunks = []
    for interval in np.unique(intervals):
        members = np.nonzero(intervals == interval)[0]
        for k in range(0, len(members), CHUNK_SCENES):
            chunk = members[k : k + CHUNK_SCENES]
            if len(chunk) >= MIN_CHUNK_SCENES:
                chunks.append(chunk)
    return chunks


def score_chunk(chunk, scenes, corners, satellite_for) -> dict[str, float]:
    """Gets the intersection over union of a chunk, by time step."""
    middle = corners[chunk[len(chunk) // 2]].mean(axis=0)
    lon, lat, _ = TO_GEODETIC.transform(*middle)
    project = Transformer.from_crs(
        "EPSG:4978", f"+proj=aeqd +lat_0={lat} +lon_0={lon} +units=m", always_xy=True
    )
    to_plane = Transformer.from_crs(
        "EPSG:4326", f"+proj=aeqd +lat_0={lat} +lon_0={lon} +units=m", always_xy=True
    )

    def plane(xyz):
        x, y, _ = project.transform(*xyz.T)
        return np.column_stack([x, y])

    parts = []
    for position, k in enumerate(chunk):
        c = corners[k]
        center_left, center_right = (c[0] + c[3]) / 2, (c[1] + c[2]) / 2
        if position == 0:
            # the forward half of the first scene, from its center
            outline = np.array([c[0], c[1], center_right, center_left])
        elif position == len(chunk) - 1:
            outline = np.array([center_left, center_right, c[2], c[3]])
        else:
            outline = c
        parts.append(Polygon(plane(outline)).buffer(0))
    observed = unary_union(parts)
    first, last = scenes.time[chunk[0]], scenes.time[chunk[-1]]
    duration = (last - first).total_seconds()
    results = {
        "scenes": len(chunk),
        "latitude": float(lat),
    }
    for time_step in TIME_STEPS:
        n = max(1, int(round(duration / time_step)))
        step = duration / n
        times = [
            (first + pd.Timedelta(seconds=step * (i + 0.5))).floor("us").to_pydatetime()
            for i in range(n)
        ]
        track = compute_ground_track(satellite_for(step), times)
        predicted = to_local(track.geometry.iloc[0], to_plane)
        results[f"iou_{time_step}"] = (
            predicted.intersection(observed).area / predicted.union(observed).area
        )
        if time_step == TIME_STEPS[0]:
            results["observed_only"] = (
                observed.difference(predicted).area / observed.area
            )
            results["predicted_only"] = (
                predicted.difference(observed).area / predicted.area
            )
    return results


def score_area(scenes, corners, intervals, satellite_for) -> tuple[float, float]:
    """Gets the median and 5th percentile of the lowest IoU over time steps."""
    chunks = get_chunks(intervals)

    def score(chunk):
        return score_chunk(chunk, scenes, corners, satellite_for)

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        table = pd.DataFrame(list(pool.map(score, chunks)))
    print(f"swath area over {len(table)} chunks of 3 to {CHUNK_SCENES} scenes:")
    for time_step in TIME_STEPS:
        print(
            f"  IoU at {time_step} s: median {table[f'iou_{time_step}'].median():.5f}, "
            f"p5 {table[f'iou_{time_step}'].quantile(0.05):.5f}"
        )
    print(
        f"  observed only {100 * table.observed_only.median():.2f}%, "
        f"predicted only {100 * table.predicted_only.median():.2f}%",
        flush=True,
    )
    medians = [table[f"iou_{step}"].median() for step in TIME_STEPS]
    lows = [table[f"iou_{step}"].quantile(0.05) for step in TIME_STEPS]
    return float(min(medians)), float(min(lows))


def score_scene_centers(scenes, corners, satellite) -> dict[str, float]:
    """Compares epochs with the center times of scenes, at their centers."""
    lon, lat, _ = TO_GEODETIC.transform(*corners.mean(axis=1).T)

    def observe(k):
        time = scenes.time[k].floor("us").to_pydatetime()
        found = collect_observations(
            Point(latitude=float(lat[k]), longitude=float(lon[k])),
            satellite,
            time - SEARCH_WINDOW,
            time + SEARCH_WINDOW,
        )
        if found.empty:
            return np.nan
        differences = (found.epoch - scenes.time[k]).dt.total_seconds().to_numpy()
        return differences[np.argmin(np.abs(differences))]

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        differences = np.array(list(pool.map(observe, range(len(scenes)))))
    found = np.isfinite(differences)
    results = {
        "epoch_scene_center_p95": summarize(differences[found])["p95"],
        "epoch_scene_center_median": summarize(differences[found])["median"],
        "scene_centers_observed": 100 * float(found.mean()),
    }
    print(json.dumps(results), flush=True)
    return results


def parse_angle_file(text: str) -> dict:
    """Parses the values of a Landsat angle coefficient file."""
    values = {}
    for name, value in re.findall(
        r"^\s*(\w+)\s*=\s*(\([^)]*\)|\"[^\"]*\"|\S+)", text, re.MULTILINE
    ):
        if value.startswith("("):
            values[name] = np.array(
                [float(item) for item in value[1:-1].replace("\n", " ").split(",")]
            )
        elif value.startswith('"'):
            values[name] = value[1:-1]
        elif name not in ("GROUP", "END_GROUP"):
            try:
                values[name] = float(value)
            except ValueError:
                values[name] = value
    return values


def imaging_times(values: dict) -> np.ndarray | None:
    """
    Locates sample points of a scene and gets the times its band imaged them, using the
    line time and the polynomials from its map coordinates to image lines. Gets rows of
    latitude, longitude, and seconds from the start of the file's day of ephemeris, or
    `None` if the scene is not in a UTM projection.
    """
    if values["MAP_PROJECTION"] != "UTM":
        return None
    size = values[f"{BAND}_PIXEL_SIZE"]
    lines = values[f"{BAND}_L1T_IMAGE_CORNER_LINES"]
    samples = values[f"{BAND}_L1T_IMAGE_CORNER_SAMPS"]
    corner_x, corner_y = values["UL_CORNER"]
    to_geographic = Transformer.from_crs(
        f"+proj=utm +zone={int(values['UTM_ZONE'])} +datum=WGS84",
        "EPSG:4326",
        always_xy=True,
    )
    rows = []
    for along in ALONG_FRACTIONS:
        for cross in CROSS_FRACTIONS:
            # corners are upper left, upper right, lower right, and lower left
            weights = np.array(
                [
                    (1 - along) * (1 - cross),
                    (1 - along) * cross,
                    along * cross,
                    along * (1 - cross),
                ]
            )
            line, sample = float(weights @ lines), float(weights @ samples)
            times = []
            for sca in values[f"{BAND}_SCA_LIST"].astype(int):
                prefix = f"{BAND}_SCA{sca:02d}"
                mean_r = values[f"{prefix}_MEAN_L1R_LINE_SAMP"]
                mean_t = values[f"{prefix}_MEAN_L1T_LINE_SAMP"]
                dl, ds = line - mean_t[0], sample - mean_t[1]
                terms = np.array([1, dl, ds, 0, dl * ds])

                def evaluate(name):
                    numerator = terms @ values[f"{prefix}_{name}_NUM_COEF"]
                    denominator = 1 + terms[1:] @ values[f"{prefix}_{name}_DEN_COEF"]
                    return numerator / denominator

                image_line = mean_r[0] + evaluate("LINE")
                image_sample = mean_r[1] + evaluate("SAMP")
                if (
                    0 <= image_sample < values[f"{BAND}_NUM_L1R_SAMPS"]
                    and 0 <= image_line < values[f"{BAND}_NUM_L1R_LINES"]
                ):
                    times.append(
                        values["EPHEMERIS_EPOCH_SECONDS"]
                        + values[f"{BAND}_START_TIME"]
                        + image_line * values[f"{BAND}_LINE_TIME"]
                    )
            if times:
                lon, lat = to_geographic.transform(
                    corner_x + sample * size, corner_y - line * size
                )
                rows.append([lat, lon, along, cross, times[0]])
    return np.array(rows)


def load_imaging_times(platform, day, directory) -> pd.DataFrame | None:
    """Gets imaging times of points across a sample of scenes, if files are available."""
    path = directory / f"imaging_times_{platform}_{day:%Y%m%d}.csv"
    if path.exists():
        return pd.read_csv(path)
    body = {
        "collections": ["landsat-c2-l2"],
        "datetime": f"{day:%Y-%m-%d}T00:00:00Z/{day:%Y-%m-%d}T23:59:59.999Z",
        "limit": 500,
        "query": {"platform": {"eq": platform}},
    }
    items = []
    while True:
        page = post_with_retries(ANGLES_URL, body).json()
        items += page["features"]
        following = [
            link for link in page.get("links", []) if link.get("rel") == "next"
        ]
        if not following or not page["features"]:
            break
        body = following[0]["body"]
    if not items:
        return None
    items.sort(key=lambda item: item["properties"]["datetime"])
    chosen = items[::ANGLE_SCENE_STRIDE]
    print(f"{len(items)} scenes with angle files, reading {len(chosen)}", flush=True)
    token = get_with_retries(ANGLES_TOKEN_URL).json()["token"]

    def read(item):
        response = get_with_retries(f"{item['assets']['ang']['href']}?{token}")
        values = parse_angle_file(response.text)
        rows = imaging_times(values)
        if rows is None or len(rows) == 0:
            return None
        year = int(values["EPHEMERIS_EPOCH_YEAR"])
        origin = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(
            days=int(values["EPHEMERIS_EPOCH_DAY"]) - 1
        )
        frame = pd.DataFrame(
            rows, columns=["latitude", "longitude", "along", "cross", "seconds"]
        )
        frame["time"] = pd.Timestamp(origin) + pd.to_timedelta(frame.seconds, unit="s")
        frame["center_time"] = pd.Timestamp(item["properties"]["datetime"])
        frame["scene"] = item["id"]
        return frame.drop(columns="seconds")

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        frames = [frame for frame in pool.map(read, chosen) if frame is not None]
    table = pd.concat(frames, ignore_index=True)
    table.to_csv(path, index=False)
    return pd.read_csv(path)


def check_imaging_times(table: pd.DataFrame) -> None:
    """Checks the times decoded from angle files against the center times of scenes."""
    middle = table[(table.along == 0.5) & (np.isclose(table["cross"], 0.5))]
    difference = (
        pd.to_datetime(middle.time, utc=True, format="ISO8601")
        - pd.to_datetime(middle.center_time, utc=True, format="ISO8601")
    ).dt.total_seconds()
    print(
        f"decoded times at scene centers minus catalog center times: median "
        f"{difference.median():+.3f} s, largest {difference.abs().max():.3f} s, "
        f"{len(difference)} scenes",
        flush=True,
    )


def score_imaging_times(table: pd.DataFrame, satellite) -> dict[str, float]:
    """Compares epochs with the imaging times of points across the swath."""
    times = pd.to_datetime(table.time, utc=True, format="ISO8601")

    def observe(k):
        time = times[k].floor("us").to_pydatetime()
        found = collect_observations(
            Point(latitude=table.latitude[k], longitude=table.longitude[k]),
            satellite,
            time - SEARCH_WINDOW,
            time + SEARCH_WINDOW,
        )
        if found.empty:
            return np.nan
        differences = (found.epoch - times[k]).dt.total_seconds().to_numpy()
        return differences[np.argmin(np.abs(differences))]

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        differences = np.array(list(pool.map(observe, range(len(table)))))
    found = np.isfinite(differences)
    print(f"imaging times: {int(found.sum())} of {len(table)} points observed")
    for cross in np.unique(table["cross"]):
        chosen = found & (table["cross"].to_numpy() == cross)
        print(
            f"  across the scene at {cross:.2f}: median "
            f"{np.median(differences[chosen]):+.3f} s"
        )
    results = {
        "epoch_p95": summarize(differences[found])["p95"],
        "epoch_median": summarize(differences[found])["median"],
        "imaging_points_observed": 100 * float(found.mean()),
    }
    print(json.dumps(results), flush=True)
    return results


def validate(
    day, tle, entry, orbit, scenes, corners, platform, directory, write
) -> None:
    """Scores the configuration, and writes its results if requested."""

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, KEY, time_step=time_step, platform=platform)

    satellite = satellite_for(TIME_STEPS[0])
    intervals = get_intervals(scenes)
    print(f"{intervals.max() + 1} sequences of scenes", flush=True)
    results = {
        "swath_width_median": score_swath_width(
            scenes, corners, orbit, satellite.instruments[0]
        )
    }
    centers = score_scene_centers(scenes, corners, satellite)
    results.update(centers)
    table = load_imaging_times(platform, day, directory)
    if table is None:
        print("no angle files: epochs are compared at the centers of scenes only")
        results["epoch_p95"] = centers["epoch_scene_center_p95"]
        results["epoch_median"] = centers["epoch_scene_center_median"]
    else:
        check_imaging_times(table)
        results.update(score_imaging_times(table, satellite))
    points = generate_points_fibonacci_lattice(POINT_SPACING)
    passes = score_passes(points, satellite, scenes, corners, intervals)
    print(json.dumps(passes), flush=True)
    results["pass_recall"] = passes["pass_recall"]
    results["observed"] = passes["observed"]
    results["matched"] = passes["matched"]
    median, low = score_area(scenes, corners, intervals, satellite_for)
    results["footprint_iou_median"] = median
    results["footprint_iou_p5"] = low
    print(json.dumps(results), flush=True)
    if not write:
        return
    if table is None:
        raise ValueError("a validation record needs imaging times across the swath")
    names = (
        "pass_recall",
        "epoch_p95",
        "epoch_median",
        "epoch_scene_center_p95",
        "swath_width_median",
        "footprint_iou_median",
        "footprint_iou_p5",
    )
    if platform != FIT_PLATFORM:
        add_transfer_metrics(
            RECORD_PATH,
            KEY,
            f"{platform}, {day:%Y-%m-%d}",
            {name: results[name] for name in names},
        )
        return
    rng = np.random.default_rng(1)
    picks = table.iloc[rng.choice(len(table), FIXTURE_PASSES, replace=False)]
    fixture = np.stack(
        [
            picks.latitude.to_numpy(),
            picks.longitude.to_numpy(),
            (pd.to_datetime(picks.time, utc=True, format="ISO8601") - pd.Timestamp(day))
            .dt.total_seconds()
            .to_numpy(),
        ],
        axis=-1,
    )
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=fixture
    )
    update_record(
        RECORD_PATH,
        "oli",
        KEY,
        {
            "digest": entry.configurations[KEY].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_oli.py",
            "reference": "USGS Landsat Collection 2 scene footprints and center times; "
            + "angle coefficient files of Level-2 products",
            "period": f"{day:%Y-%m-%d}",
            "platform": platform,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": [
                {"name": name, "value": float(results[name])} for name in names
            ],
            "fixture_metrics": recompute_passes(satellite, 0, day, fixture),
        },
    )
    print(f"wrote {RECORD_PATH}")


def main() -> None:
    """Fits or validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--platform", choices=sorted(PLATFORMS), default=FIT_PLATFORM)
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    tle = ELEMENT_SETS[arguments.day][arguments.platform]
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    scenes = load_scenes(arguments.platform, day, arguments.directory)
    corners = add_corners(scenes, orbit)
    print(f"{arguments.platform}, {arguments.day}: {tle[0][:32]}", flush=True)
    if arguments.fit:
        fit_geometry(scenes, corners, orbit)
        return
    validate(
        day,
        tle,
        get_entry("oli"),
        orbit,
        scenes,
        corners,
        arguments.platform,
        arguments.directory,
        arguments.write,
    )


if __name__ == "__main__":
    main()
