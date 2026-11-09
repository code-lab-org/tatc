"""
Validates the C-SAR catalog entry against the Level-1 Ground Range Detected (GRD)
products of Sentinel-1C, using the exact values of the entry through public TAT-C
methods. The annotation of each product holds a geolocation grid, with the zero-Doppler
time of each of its lines, and the restituted orbit:

* The swath edges are compared with the near-range and far-range points of the grid
  lines, against the rays at the edges of the view at the times of the lines.
* The swath area is compared with the union of consecutive products of an acquisition,
  from the middle of the first to the middle of the last, so that the arbitrary ends of
  products are left out, at several time steps.
* Passes are those of lattice points inside the grids.
* Epochs are compared with the zero-Doppler times of interior grid points.

Each mode is a configuration, whose field of view and roll profile (the roll steering,
as a function of the argument of latitude) come from `--fit`, on a day of products.
Running on a day at least 30 days from it tests them on data held out from the fit.

Needs the `validation` extra:

    python docs/validation/validate_catalog_csar.py <data directory> 2026-10-02 --fit
    python docs/validation/validate_catalog_csar.py <data directory> 2026-08-12

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

With `--platform sentinel-1d`, it tests the same values on a sibling satellite, and
adds the results to the record as a transfer, which does not count toward the tier.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
from shapely.ops import unary_union

import tatc
from tatc.analysis import collect_observations, compute_ground_track
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit, Point, PointedInstrument
from tatc.utils import compute_argument_of_latitude, compute_projected_ray_position
from validation_io import (
    add_transfer_metrics,
    TO_GEODETIC,
    ecef,
    get_holdout_gap,
    get_with_retries,
    in_lattice_cells,
    match_passes,
    normalize,
    post_with_retries,
    predict_passes,
    read_element_set_file,
    recompute_passes,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "csar.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "csar.npz"
ELEMENT_SET_FILE = REPOSITORY / "docs" / "validation" / "sentinel_2026-08-12.3le"

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
TOKEN_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-1-grd"

FIT_DAY = datetime(2026, 10, 2, tzinfo=timezone.utc)
PLATFORM = FIT_PLATFORM = "sentinel-1c"
# sibling platforms, on which results are reported as a transfer, with public element
# sets from CelesTrak, with epochs after each day of data
TRANSFER_ELEMENTS = {
    ("sentinel-1d", "2026-11-05"): [
        "1 66315U 25251A   26312.30420366  .00000154  00000+0  42399-4 0  9994",
        "2 66315  98.1818 317.6134 0001246  87.7733 272.3609 14.59200353 53741",
    ],
}
NORAD_ID = 62261
# the element set of the validation notebook, for the fit day
ELEMENT_SETS = {
    "2026-10-02": [
        "1 62261U 24235A   26275.78977191 -.00000047  00000+0 -21660-6 0  9997",
        "2 62261  98.1826 281.5711 0001343  90.3841 269.7512 14.59198718 97106",
    ],
}
MODES = {"iw": "IW", "ew": "EW"}
SPEED_OF_LIGHT = 299792458

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


def get_element_set(day: datetime) -> list[str]:
    """Gets the element set of a day."""
    key = f"{day:%Y-%m-%d}"
    if PLATFORM != FIT_PLATFORM:
        return TRANSFER_ELEMENTS[(PLATFORM, key)]
    if key in ELEMENT_SETS:
        return ELEMENT_SETS[key]
    return read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)


def to_utc(day: datetime, seconds) -> list[datetime]:
    """Converts seconds from the start of a day to times."""
    return [day + timedelta(seconds=float(s)) for s in np.atleast_1d(seconds)]


def search_items(day: datetime) -> list[dict]:
    """Gets the catalog items of the IW and EW products of a day."""
    body = {
        "collections": ["sentinel-1-grd"],
        "datetime": f"{day:%Y-%m-%d}T00:00:00Z/{day:%Y-%m-%d}T23:59:59Z",
        "limit": 1000,
        "query": {"platform": {"eq": PLATFORM}},
    }
    items = []
    while True:
        page = post_with_retries(STAC_URL, body).json()
        items += page["features"]
        following = [
            link for link in page.get("links", []) if link.get("rel") == "next"
        ]
        if not following or not page["features"]:
            break
        body = following[0].get("body", body)
    return [
        item
        for item in items
        if item["properties"]["sar:instrument_mode"] in MODES.values()
    ]


def read_annotation(item: dict, token: str, day: datetime) -> dict:
    """Reads the geolocation grid and orbit of a product from its annotation."""
    base = item["assets"]["safe-manifest"]["href"].rsplit("/", 1)[0]
    mode = item["properties"]["sar:instrument_mode"]
    polarization = item["properties"]["sar:polarizations"][0].lower()
    response = get_with_retries(
        f"{base}/annotation/{mode.lower()}-{polarization}.xml?{token}"
    )
    root = ET.fromstring(response.content)

    def seconds(text):
        return (pd.Timestamp(text, tz="UTC") - pd.Timestamp(day)).total_seconds()

    points = root.findall(".//geolocationGridPoint")
    lines = sorted({int(p.find("line").text) for p in points})
    samples = sorted({int(p.find("pixel").text) for p in points})

    def grid(tag, convert=float):
        return np.array([convert(p.find(tag).text) for p in points]).reshape(
            len(lines), len(samples)
        )

    orbits = root.findall(".//orbitList/orbit")
    return {
        "id": item["id"],
        "mode": mode,
        "datatake": str(item["properties"]["s1:datatake_id"]),
        "time": grid("azimuthTime", seconds),
        "slant_range": grid("slantRangeTime")[0] * SPEED_OF_LIGHT / 2,
        "lon": grid("longitude"),
        "lat": grid("latitude"),
        "orbit_time": np.array([seconds(o.find("time").text) for o in orbits]),
        "orbit_position": np.array(
            [[float(o.find(f"position/{c}").text) for c in "xyz"] for o in orbits]
        ),
        "orbit_velocity": np.array(
            [[float(o.find(f"velocity/{c}").text) for c in "xyz"] for o in orbits]
        ),
    }


def load_products(day: datetime, directory: pathlib.Path) -> list[dict]:
    """Gets the reduced annotations of a day's IW and EW products, in time order."""
    path = directory / f"csar_{PLATFORM}_{day:%Y%m%d}.npz"
    if not path.exists():
        items = search_items(day)
        print(f"reading the annotations of {len(items)} products", flush=True)
        token = get_with_retries(TOKEN_URL).json()["token"]
        with concurrent.futures.ThreadPoolExecutor(16) as pool:
            records = list(pool.map(lambda i: read_annotation(i, token, day), items))
        np.savez_compressed(
            path,
            id=[r["id"] for r in records],
            mode=[r["mode"] for r in records],
            datatake=[r["datatake"] for r in records],
            n_line=[r["time"].shape[0] for r in records],
            n_sample=[r["time"].shape[1] for r in records],
            n_orbit=[len(r["orbit_time"]) for r in records],
            **{
                k: np.concatenate([r[k].ravel() for r in records])
                for k in ("time", "lon", "lat")
            },
            **{
                k: np.concatenate([r[k] for r in records])
                for k in ("orbit_time", "orbit_position", "orbit_velocity")
            },
        )
    with np.load(path, allow_pickle=False) as data:
        cache = {key: data[key] for key in data.files}
    size = cache["n_line"] * cache["n_sample"]

    def cut(key, counts):
        return np.split(cache[key], np.cumsum(counts)[:-1])

    products = []
    for k, (t, lon, lat, ot, op, ov) in enumerate(
        zip(
            cut("time", size),
            cut("lon", size),
            cut("lat", size),
            cut("orbit_time", cache["n_orbit"]),
            cut("orbit_position", cache["n_orbit"]),
            cut("orbit_velocity", cache["n_orbit"]),
            strict=True,
        )
    ):
        shape = (int(cache["n_line"][k]), int(cache["n_sample"][k]))
        products.append(
            {
                "id": str(cache["id"][k]),
                "mode": str(cache["mode"][k]),
                "datatake": str(cache["datatake"][k]),
                # the samples of a line differ by microseconds
                "time": t.reshape(shape).mean(axis=1),
                "lon": lon.reshape(shape),
                "lat": lat.reshape(shape),
                "orbit_time": ot,
                "orbit_position": op,
                "orbit_velocity": ov,
            }
        )
    products.sort(key=lambda p: p["time"][0])
    print(
        f"{len(products)} products: "
        + ", ".join(
            f"{sum(p['mode'] == mode for p in products)} {mode}"
            for mode in MODES.values()
        ),
        flush=True,
    )
    return products


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


def score_edges(products, orbit, day, instrument) -> dict[str, np.ndarray]:
    """
    Gets the horizontal distances (km) between the near-range and far-range points of
    the grid lines and the rays at the edges of the view at the times of the lines.
    """
    distances = {"near": [], "far": []}
    for product in products:
        times = product["time"]
        track = orbit.to_gp_orbit().get_orbit_track(to_utc(day, times))
        roll = instrument.get_roll_angle(track)
        half = instrument.cross_track_field_of_view / 2
        # right-looking views have negative roll angles, with near range on the left
        for edge, column, sign in (("near", 0, 1), ("far", -1, -1)):
            ray = compute_projected_ray_position(track, 0, 0, roll + sign * half)
            lon, lat = product["lon"][:, column], product["lat"][:, column]
            offset = ecef(ray.longitude.degrees, ray.latitude.degrees) - ecef(lon, lat)
            vertical = up(lon, lat)
            horizontal = offset - np.sum(offset * vertical, axis=-1)[:, None] * vertical
            distances[edge].append(np.linalg.norm(horizontal, axis=-1) / 1e3)
    return {edge: np.concatenate(values) for edge, values in distances.items()}


def fit(products, orbit, day) -> None:
    """Measures the field of view and roll profile of each mode."""
    lines = measure_lines(products, orbit, day)
    for key, mode in MODES.items():
        chosen = lines[lines["mode"] == mode]
        field_of_view = float(chosen.field_of_view.median())
        roll = float(chosen.roll.median())
        print(
            f"{key}: {chosen['product'].nunique()} products, {len(chosen)} lines; near "
            f"{chosen.near.median():.4f} deg, far {chosen.far.median():.4f} deg, roll "
            f"{roll:.4f} deg (p5 {chosen.roll.quantile(0.05):.4f}, p95 "
            f"{chosen.roll.quantile(0.95):.4f}), field of view {field_of_view:.4f} deg "
            f"(p5 {chosen.field_of_view.quantile(0.05):.4f}, p95 "
            f"{chosen.field_of_view.quantile(0.95):.4f})"
        )
        profile = fit_profile(chosen)
        print(f"  roll profile: {json.dumps(profile)}")
        print(f"  far look angle, largest: {-chosen.far.min():.4f} deg")
        # a check within the fit day: profiles of even products, on odd products
        even = chosen[chosen["product"] % 2 == 0]
        odd = sorted(set(chosen["product"][chosen["product"] % 2 == 1]))
        for label, parameters in (
            ("fixed roll", {"roll_angle": roll}),
            ("roll profile", {"roll_angle_profile": fit_profile(even)}),
        ):
            instrument = PointedInstrument(
                name="check",
                field_of_regard=120,
                cross_track_field_of_view=field_of_view,
                along_track_field_of_view=1,
                is_rectangular=True,
                **parameters,
            )
            edges = score_edges([products[k] for k in odd], orbit, day, instrument)
            print(
                f"  {label}, odd products: edges p95 near "
                f"{np.percentile(edges['near'], 95):.2f} km, far "
                f"{np.percentile(edges['far'], 95):.2f} km"
            )


def outline(product: dict, first: int, last: int) -> np.ndarray:
    """Gets the Earth-fixed outline of grid lines from a first to a last."""
    lon, lat = product["lon"][first : last + 1], product["lat"][first : last + 1]
    ring_lon = np.concatenate([lon[0], lon[1:, -1], lon[-1, -2::-1], lon[-2:0:-1, 0]])
    ring_lat = np.concatenate([lat[0], lat[1:, -1], lat[-1, -2::-1], lat[-2:0:-1, 0]])
    return ecef(ring_lon, ring_lat)


def get_chunks(products: list[dict]) -> list[list[int]]:
    """Divides the consecutive products of each acquisition into chunks."""
    sequences, current = [], []
    for k, product in enumerate(products):
        if current:
            previous = products[current[-1]]
            if (
                product["datatake"] != previous["datatake"]
                or product["time"][0] - previous["time"][-1] > MAX_PRODUCT_GAP
            ):
                sequences.append(current)
                current = []
        current.append(k)
    sequences.append(current)
    chunks = []
    for sequence in sequences:
        for k in range(0, len(sequence), CHUNK_PRODUCTS):
            chunk = sequence[k : k + CHUNK_PRODUCTS]
            if len(chunk) >= MIN_CHUNK_PRODUCTS:
                chunks.append(chunk)
    return chunks


def score_chunk(chunk, products, day, satellite_for) -> dict[str, float]:
    """Gets the intersection over union of a chunk, by time step."""
    first, last = products[chunk[0]], products[chunk[-1]]
    middle_first, middle_last = len(first["time"]) // 2, len(last["time"]) // 2
    center = products[chunk[len(chunk) // 2]]
    lon0, lat0 = (
        float(np.mean(center["lon"])),
        float(np.mean(center["lat"])),
    )
    project = Transformer.from_crs(
        "EPSG:4978",
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +units=m",
        always_xy=True,
    )
    to_plane = Transformer.from_crs(
        "EPSG:4326",
        f"+proj=aeqd +lat_0={lat0} +lon_0={lon0} +units=m",
        always_xy=True,
    )

    def plane(xyz):
        x, y, _ = project.transform(*xyz.T)
        return Polygon(np.column_stack([x, y])).buffer(0)

    parts = []
    for position, k in enumerate(chunk):
        product = products[k]
        n = len(product["time"])
        if position == 0:
            parts.append(plane(outline(product, middle_first, n - 1)))
        elif position == len(chunk) - 1:
            parts.append(plane(outline(product, 0, middle_last)))
        else:
            parts.append(plane(outline(product, 0, n - 1)))
    observed = unary_union(parts)
    start = first["time"][middle_first]
    duration = last["time"][middle_last] - start
    results = {"products": len(chunk), "latitude": lat0}
    for time_step in TIME_STEPS:
        n = max(1, int(round(duration / time_step)))
        step = duration / n
        times = to_utc(day, start + step * (np.arange(n) + 0.5))
        track = compute_ground_track(satellite_for(step), times)
        geometry = track.geometry.iloc[0]
        predicted = unary_union(
            [
                Polygon(
                    np.column_stack(
                        to_plane.transform(*np.array(part.exterior.coords).T[:2])
                    )
                ).buffer(0)
                for part in getattr(geometry, "geoms", [geometry])
            ]
        )
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


def score_area(products, day, satellite_for) -> pd.DataFrame:
    """Gets the intersection over union of each chunk of products, by time step."""
    chunks = get_chunks(products)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        table = pd.DataFrame(
            list(
                pool.map(
                    lambda chunk: score_chunk(chunk, products, day, satellite_for),
                    chunks,
                )
            )
        )
    print(f"swath area over {len(table)} chunks of 3 to {CHUNK_PRODUCTS} products:")
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
    return table


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


def validate_mode(key, day, tle, entry, orbit, products, write) -> None:
    """Scores the configuration of a mode, and writes its results if requested."""
    chosen = [p for p in products if p["mode"] == MODES[key]]
    print(f"--- {key}: {len(chosen)} products", flush=True)

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, key, time_step=time_step, platform=PLATFORM)

    satellite = satellite_for(TIME_STEPS[0])
    instrument = satellite.instruments[0]
    edges = score_edges(chosen, orbit, day, instrument)
    for edge, values in edges.items():
        print(
            f"swath edge, {edge} range: median {np.median(values):.2f} km, p95 "
            f"{np.percentile(values, 95):.2f} km",
            flush=True,
        )
    points = epoch_points(chosen)
    differences = score_epochs(points, satellite, day)
    found = np.isfinite(differences)
    epochs = summarize(differences[found])
    print(
        f"epochs: {int(found.sum())} of {len(points)} points observed, median "
        f"{epochs['median']:+.3f} s, p95 {epochs['p95']:.3f} s",
        flush=True,
    )
    lattice = generate_points_fibonacci_lattice(POINT_SPACING)
    passes = score_passes(lattice, satellite, chosen, day)
    print(json.dumps(passes), flush=True)
    area = score_area(chosen, day, satellite_for)
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
        *[
            {
                "name": name,
                "value": float(value),
                "scope": f"time step {time_step} s",
            }
            for time_step in TIME_STEPS
            for name, value in (
                ("footprint_iou_median", area[f"iou_{time_step}"].median()),
                ("footprint_iou_p5", area[f"iou_{time_step}"].quantile(0.05)),
            )
        ],
        {"name": "epoch_p95", "value": float(epochs["p95"])},
        {"name": "epoch_median", "value": float(epochs["median"])},
    ]
    print(json.dumps(metrics), flush=True)
    if not write:
        return
    if PLATFORM != FIT_PLATFORM:
        add_transfer_metrics(
            RECORD_PATH,
            key,
            f"{PLATFORM}, {day:%Y-%m-%d}",
            {
                metric["name"]
                + (
                    "_" + metric["scope"].replace(" ", "_").replace(".", "_")
                    if "scope" in metric
                    else ""
                ): metric["value"]
                for metric in metrics
            },
        )
        print(f"added {key} transfer metrics to {RECORD_PATH}")
        return
    rng = np.random.default_rng(1)
    fixture = points[found][
        rng.choice(int(found.sum()), min(FIXTURE_PASSES, int(found.sum())), False)
    ]
    update_fixture(
        FIXTURE_PATH,
        tle=np.array(tle),
        day=np.array(day.isoformat()),
        **{f"passes_{key}": fixture},
    )
    update_record(
        RECORD_PATH,
        "csar",
        key,
        {
            "digest": entry.configurations[key].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_csar.py",
            "reference": "Geolocation grids and restituted orbits of the annotations "
            + "of Sentinel-1 Level-1 GRD products",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": metrics,
            "fixture_metrics": recompute_passes(satellite, 0, day, fixture),
        },
    )
    print(f"wrote {key} to {RECORD_PATH}")


def main() -> None:
    """Fits or validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--mode", choices=sorted(MODES), action="append")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="read and cache the reference data, without scoring",
    )
    parser.add_argument("--platform", default=FIT_PLATFORM)
    arguments = parser.parse_args()
    global PLATFORM  # pylint: disable=global-statement
    PLATFORM = arguments.platform
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    products = load_products(day, arguments.directory)
    if arguments.download_only:
        return
    tle = get_element_set(day)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    print(f"{PLATFORM}, {arguments.day}: {tle[0][:32]}", flush=True)
    if arguments.fit:
        fit(products, orbit, day)
        return
    entry = get_entry("csar")
    for key in arguments.mode or list(MODES):
        validate_mode(key, day, tle, entry, orbit, products, arguments.write)


if __name__ == "__main__":
    main()
