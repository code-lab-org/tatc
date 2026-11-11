"""
Validates the L-SAR catalog entry against NISAR's Level 1 range-Doppler single-look
complex (RSLC) products, using the exact values of the entry through public TAT-C
methods. Only metadata are read, with HTTP range requests, as in `ValidateSARNISAR`:
the zero-Doppler times and bounding polygon of every product of a day, and the
ground layer of the geolocation grid and the orbit of a sample of them.

* The swath edges are compared with the near-range and far-range points of the grid
  lines, against the rays at the edges of the view at the times of the lines.
* The swath area is compared with the union of the bounding polygons of consecutive
  full frames of an acquisition, at several time steps. Partial frames are left out,
  as their polygons outline only the part with valid data.
* Passes are those of lattice points inside the grids.
* Epochs are compared with the zero-Doppler times of interior grid points.

The view of each group of modes (the 77 MHz mode images a narrower swath than the
others) and its roll profile come from `--fit`, on a day of products, as for C-SAR.
Running on a day at least 30 days from it tests them on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_lsar.py <data directory> 2026-10-02 --fit
    python docs/validation/validate_catalog_lsar.py <data directory> 2026-11-04

Several days pool their results, as for the 77 MHz mode, of which few products are
acquired, each with its own element set; the fixture keeps the passes of the first day,
with its element set:

    python docs/validation/validate_catalog_lsar.py <data directory> 2026-11-01 \
        2026-11-02 2026-11-03 2026-11-05 --configuration swath_77mhz

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import io
import json
import pathlib
import re
from datetime import datetime, timedelta, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd
import requests
import shapely
from pyproj import Transformer
from shapely.ops import transform as transform_geometry
from shapely.ops import unary_union

import tatc
from tatc.analysis import compute_ground_track
from tatc.catalog import get_entry
from tatc.catalog.metrics import compute_pass_scores, summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from tatc.utils import compute_projected_ray_position
from validate_catalog_csar import (
    epoch_points,
    fit_profile,
    measure_lines,
    observe_passes,
    score_epochs,
    to_utc,
)
from validation_io import (
    ecef,
    get_holdout_gap,
    match_passes,
    predict_passes,
    read_element_set_file,
    recompute_passes,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "lsar.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "lsar.npz"
SEARCH_URL = "https://api.daac.asf.alaska.edu/services/search/param"

FIT_DAY = datetime(2026, 10, 2, tzinfo=timezone.utc)
PLATFORM = "nisar"
# public element sets from CelesTrak, with epochs after each day of data (that of the
# fit day is the one of ValidateSARNISAR)
ELEMENT_SETS = {
    "2026-10-02": [
        "1 65053U 25163A   26276.92799516  .00000171  00000+0  64585-4 0  9997",
        "2 65053  98.4055 103.0181 0001222  89.8720 270.2610 14.42503663 62063",
    ],
    "2026-11-04": [
        "1 65053U 25163A   26312.30372745  .00000140  00000+0  55244-4 0  9999",
        "2 65053  98.4064 137.8903 0001231  89.5513 270.5818 14.42501786 67160",
    ],
}
# the 77 MHz mode's held-out days, pooled, whose products became available on
# 2026-11-09 and 10, each with the element set nearest the middle of the day from
# docs/validation/nisar_2026-11-01.3le (they were first tested with the element set of
# 2026-11-04 above, whose epoch is 3 to 7 days after them)
POOLED_DAYS = ("2026-11-01", "2026-11-02", "2026-11-03", "2026-11-05")
ELEMENT_SET_FILE = pathlib.Path(__file__).parent / "nisar_2026-11-01.3le"
NORAD_ID = 65053


def get_element_set(day: str) -> list[str]:
    """Gets the element set with which a day of data is tested."""
    if day in POOLED_DAYS:
        return read_element_set_file(
            ELEMENT_SET_FILE,
            NORAD_ID,
            datetime.fromisoformat(day).replace(tzinfo=timezone.utc),
        )
    return ELEMENT_SETS[day]


# groups of modes, by the bandwidth codes of the product names, whose swaths differ: the
# modes with a 20 or 40 MHz band, those with only a 5 MHz band, and the 77 MHz mode
CONFIGURATIONS = ("swath_20_40mhz", "swath_5mhz", "swath_77mhz")
VALIDATED = CONFIGURATIONS
N_RANGE = 41  # resampled slant ranges across the image
SAMPLE_PRODUCTS = 150  # products of each configuration whose grids are read
TIME_STEPS = (10, 30, 120)
CHUNK_PRODUCTS = 8
MIN_CHUNK_PRODUCTS = 3
MAX_PRODUCT_GAP = 5  # s
POINT_SPACING = 100e3
SEARCH_WINDOW = timedelta(minutes=10)
FIXTURE_PASSES = 150


class RemoteFile(io.RawIOBase):
    """A read-only remote file, read in blocks with HTTP range requests on demand."""

    def __init__(self, session, url, block_size=2**16, prefetch=4):
        super().__init__()
        self.session, self.block_size, self.blocks, self.position = (
            session,
            block_size,
            {},
            0,
        )
        response = self._get(url, 0, prefetch * block_size)
        self.url = response.url
        self.size = int(response.headers["Content-Range"].split("/")[1])
        self._store(0, response.content)

    def _get(self, url, start, end):
        for attempt in range(5):
            try:
                response = self.session.get(
                    url, headers={"Range": f"bytes={start}-{end - 1}"}, timeout=300
                )
                response.raise_for_status()
                return response
            except requests.RequestException:
                if attempt == 4:
                    raise
        raise RuntimeError("unreachable")

    def _store(self, first, content):
        for i in range(0, len(content), self.block_size):
            self.blocks[first + i // self.block_size] = content[i : i + self.block_size]

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.position

    def seek(self, offset, whence=io.SEEK_SET):
        self.position = {
            io.SEEK_SET: 0,
            io.SEEK_CUR: self.position,
            io.SEEK_END: self.size,
        }[whence] + offset
        return self.position

    def readinto(self, buffer):
        n = min(len(buffer), self.size - self.position)
        if n <= 0:
            return 0
        first = self.position // self.block_size
        last = (self.position + n - 1) // self.block_size
        missing = [b for b in range(first, last + 1) if b not in self.blocks]
        if missing:
            end = min((missing[-1] + 1) * self.block_size, self.size)
            self._store(
                missing[0],
                self._get(self.url, missing[0] * self.block_size, end).content,
            )
        data = b"".join(self.blocks[b] for b in range(first, last + 1))
        offset = self.position - first * self.block_size
        buffer[:n] = data[offset : offset + n]
        self.position += n
        return n


def configuration_of(mode: str) -> str:
    """Gets the configuration of a mode (bandwidth and polarization codes)."""
    if mode.startswith("77"):
        return "swath_77mhz"
    return "swath_5mhz" if mode.startswith("0") else "swath_20_40mhz"


def search(day: datetime) -> pd.DataFrame:
    """Lists the RSLC products that start on a day."""
    for attempt in range(5):
        try:
            response = requests.get(
                SEARCH_URL,
                params={
                    "dataset": "NISAR",
                    "processingLevel": "RSLC",
                    "start": f"{day:%Y-%m-%d}T00:00:00Z",
                    "end": f"{day:%Y-%m-%d}T23:59:59Z",
                    "output": "geojson",
                    "maxResults": "5000",
                },
                timeout=600,
            )
            response.raise_for_status()
            break
        except requests.RequestException:
            if attempt == 4:
                raise
    catalog = pd.DataFrame(
        [
            {
                "granule": f["properties"]["fileID"],
                "url": f["properties"]["url"],
                "start": pd.Timestamp(f["properties"]["startTime"]),
                "orbit": str(f["properties"].get("orbit")),
                "coverage": str(f["properties"].get("frameCoverage", "")).lower(),
                "mode": "_".join(f["properties"]["fileID"].split("_")[8:10]),
            }
            for f in response.json()["features"]
        ]
    )
    catalog = catalog[catalog.start >= pd.Timestamp(day)]
    catalog = catalog.drop_duplicates("granule").sort_values("start", ignore_index=True)
    catalog["configuration"] = catalog["mode"].map(configuration_of)
    return catalog


def seconds_since(dataset, day: datetime) -> np.ndarray:
    """Gets the times of a dataset in seconds since the start of a day."""
    epoch = re.search(r"since (\S+)", dataset.attrs["units"].decode()).group(1)
    return (
        dataset[:] + (pd.Timestamp(epoch, tz="UTC") - pd.Timestamp(day)).total_seconds()
    )


def read_bounds(session, url: str, day: datetime) -> dict:
    """Reads the bounding polygon and zero-Doppler times of a product."""
    with h5py.File(RemoteFile(session, url)) as f:
        identification = f["science/LSAR/identification"]

        def seconds(key):
            text = identification[key][()].decode()
            return (pd.Timestamp(text, tz="UTC") - pd.Timestamp(day)).total_seconds()

        return {
            "url": url,
            "polygon": identification["boundingPolygon"][()].decode(),
            "zero_doppler_start": seconds("zeroDopplerStartTime"),
            "zero_doppler_end": seconds("zeroDopplerEndTime"),
            "look_direction": identification["lookDirection"][()].decode(),
        }


def read_grid(session, url: str, day: datetime) -> dict:
    """Reads the ground layer of the geolocation grid and the orbit of a product."""
    with h5py.File(RemoteFile(session, url, block_size=2**20, prefetch=2)) as f:
        rslc = f["science/LSAR/RSLC"]
        grid = rslc["metadata/geolocationGrid"]
        layer = int(np.argmin(np.abs(grid["heightAboveEllipsoid"][:])))
        grid_range = grid["slantRange"][:]
        frequency = sorted(k for k in rslc["swaths"] if k.startswith("frequency"))[0]
        image_range = rslc[f"swaths/{frequency}/slantRange"][[0, -1]]
        ranges = np.linspace(image_range[0], image_range[1], N_RANGE)

        def resample(x):
            return np.array([np.interp(ranges, grid_range, row) for row in x])

        return {
            "time": seconds_since(grid["zeroDopplerTime"], day),
            "lon": (resample(np.unwrap(grid["coordinateX"][layer], period=360)) + 180)
            % 360
            - 180,
            "lat": resample(grid["coordinateY"][layer]),
            "orbit_time": seconds_since(rslc["metadata/orbit/time"], day),
            "orbit_position": rslc["metadata/orbit/position"][:],
            "orbit_velocity": rslc["metadata/orbit/velocity"][:],
        }


def load_day(day: datetime, directory: pathlib.Path) -> tuple[pd.DataFrame, list]:
    """Reads the bounds of a day's products and the grids of a sample, or loads them."""
    bounds_path = directory / f"lsar_bounds_{day:%Y%m%d}.csv"
    grids_path = directory / f"lsar_grids_{day:%Y%m%d}.npz"
    if not (bounds_path.exists() and grids_path.exists()):
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        catalog = search(day)
        print(f"reading {len(catalog)} products", flush=True)
        with concurrent.futures.ThreadPoolExecutor(16) as pool:
            bounds = pd.DataFrame(
                pool.map(lambda url: read_bounds(session, url, day), catalog.url)
            )
        catalog = catalog.merge(bounds, on="url")
        catalog.to_csv(bounds_path, index=False)
        sample = pd.concat(
            [
                group.iloc[
                    np.unique(
                        np.linspace(0, len(group) - 1, SAMPLE_PRODUCTS)
                        .round()
                        .astype(int)
                    )
                ]
                for _, group in catalog.groupby("configuration")
            ]
        ).sort_values("start")
        print(f"reading the grids of {len(sample)} products", flush=True)
        with concurrent.futures.ThreadPoolExecutor(16) as pool:
            records = list(
                pool.map(lambda url: read_grid(session, url, day), sample.url)
            )
        np.savez_compressed(
            grids_path,
            granule=sample.granule.to_numpy().astype(str),
            n_line=[len(r["time"]) for r in records],
            n_orbit=[len(r["orbit_time"]) for r in records],
            **{
                k: np.concatenate([r[k] for r in records])
                for k in ("time", "lon", "lat")
            },
            **{
                k: np.concatenate([r[k] for r in records])
                for k in ("orbit_time", "orbit_position", "orbit_velocity")
            },
        )
    catalog = pd.read_csv(bounds_path, dtype={"orbit": str, "mode": str})
    catalog["configuration"] = catalog["mode"].map(configuration_of)
    with np.load(grids_path, allow_pickle=False) as data:
        cache = {key: data[key] for key in data.files}

    def cut(key, counts):
        return np.split(cache[key], np.cumsum(counts)[:-1])

    modes = catalog.set_index("granule")
    products = []
    for k, (t, lon, lat, ot, op, ov) in enumerate(
        zip(
            cut("time", cache["n_line"]),
            cut("lon", cache["n_line"]),
            cut("lat", cache["n_line"]),
            cut("orbit_time", cache["n_orbit"]),
            cut("orbit_position", cache["n_orbit"]),
            cut("orbit_velocity", cache["n_orbit"]),
            strict=True,
        )
    ):
        granule = str(cache["granule"][k])
        products.append(
            {
                "id": granule,
                "mode": modes.loc[granule, "mode"],
                "configuration": modes.loc[granule, "configuration"],
                "datatake": modes.loc[granule, "orbit"],
                "time": t,
                "lon": lon,
                "lat": lat,
                "orbit_time": ot,
                "orbit_position": op,
                "orbit_velocity": ov,
            }
        )
    products.sort(key=lambda p: p["time"][0])
    print(
        f"{len(catalog)} products ("
        + ", ".join(f"{n} {k}" for k, n in catalog.configuration.value_counts().items())
        + f"); grids of {len(products)}",
        flush=True,
    )
    return catalog, products


def fit(products, orbit, day) -> None:
    """Measures the field of view and roll profile of each configuration."""
    lines = measure_lines(products, orbit, day)
    lines["configuration"] = lines["product"].map(
        {k: p["configuration"] for k, p in enumerate(products)}
    )
    for key in CONFIGURATIONS:
        chosen = lines[lines.configuration == key]
        print(
            f"{key}: {chosen['product'].nunique()} products, {len(chosen)} lines; "
            f"near {chosen.near.median():.4f} deg, far {chosen.far.median():.4f} deg, "
            f"roll {chosen.roll.median():.4f} deg (p5 {chosen.roll.quantile(0.05):.4f}, "
            f"p95 {chosen.roll.quantile(0.95):.4f}), field of view "
            f"{chosen.field_of_view.median():.4f} deg (p5 "
            f"{chosen.field_of_view.quantile(0.05):.4f}, p95 "
            f"{chosen.field_of_view.quantile(0.95):.4f})"
        )
        print(f"  roll profile: {json.dumps(fit_profile(chosen))}")
        print(f"  far look angle, largest: {chosen.far.max():.4f} deg")


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
        # left-looking views have positive roll angles, with near range on the right
        for edge, column, sign in (("near", 0, -1), ("far", -1, 1)):
            ray = compute_projected_ray_position(track, 0, 0, roll + sign * half)
            lon, lat = product["lon"][:, column], product["lat"][:, column]
            offset = ecef(ray.longitude.degrees, ray.latitude.degrees) - ecef(lon, lat)
            vertical = up(lon, lat)
            horizontal = offset - np.sum(offset * vertical, axis=-1)[:, None] * vertical
            distances[edge].append(np.linalg.norm(horizontal, axis=-1) / 1e3)
    return {edge: np.concatenate(values) for edge, values in distances.items()}


def get_chunks(frames: pd.DataFrame) -> list[pd.DataFrame]:
    """Divides the consecutive full frames of each acquisition into chunks."""
    frames = frames.sort_values("zero_doppler_start")
    chunks, current = [], []
    for _, frame in frames.iterrows():
        if current and (
            frame.orbit != current[-1].orbit
            or frame["mode"] != current[-1]["mode"]
            or frame.zero_doppler_start - current[-1].zero_doppler_end > MAX_PRODUCT_GAP
        ):
            chunks.append(current)
            current = []
        current.append(frame)
    chunks.append(current)
    result = []
    for sequence in chunks:
        for k in range(0, len(sequence), CHUNK_PRODUCTS):
            chunk = sequence[k : k + CHUNK_PRODUCTS]
            if len(chunk) >= MIN_CHUNK_PRODUCTS:
                result.append(pd.DataFrame(chunk))
    return result


def unwrap(geometry):
    """Shifts the longitudes of a polygon across the antimeridian to be continuous."""
    lon = np.array(geometry.exterior.coords)[:, 0]
    if lon.max() - lon.min() > 180:
        return shapely.ops.transform(
            lambda x, y, z=None: (
                np.where(np.asarray(x) < 0, np.asarray(x) + 360, x),
                y,
            ),
            geometry,
        )
    return geometry


def score_chunk(chunk, day, satellite_for) -> dict[str, float]:
    """Gets the intersection over union of a chunk of frames, by time step."""
    polygons = [
        unwrap(shapely.force_2d(shapely.from_wkt(text))) for text in chunk.polygon
    ]
    observed = unary_union(polygons)
    center = observed.centroid
    project = Transformer.from_crs(
        "EPSG:4326",
        f"+proj=aeqd +lat_0={center.y} +lon_0={center.x} +units=m",
        always_xy=True,
    ).transform
    observed = transform_geometry(project, observed).buffer(0)
    start = chunk.zero_doppler_start.min()
    duration = chunk.zero_doppler_end.max() - start
    results = {"frames": len(chunk), "latitude": center.y}
    for time_step in TIME_STEPS:
        n = max(1, int(round(duration / time_step)))
        step = duration / n
        times = to_utc(day, start + step * (np.arange(n) + 0.5))
        geometry = compute_ground_track(satellite_for(step), times).geometry.iloc[0]
        predicted = unary_union(
            [
                transform_geometry(project, unwrap(part)).buffer(0)
                for part in getattr(geometry, "geoms", [geometry])
            ]
        )
        results[f"iou_{time_step}"] = (
            predicted.intersection(observed).area / predicted.union(observed).area
        )
    return results


def measure(key, day, entry, orbit, catalog, products) -> dict | None:
    """Measures a configuration on a day of data."""
    chosen = [p for p in products if p["configuration"] == key]
    frames = catalog[(catalog.configuration == key) & (catalog.coverage == "full")]
    print(
        f"--- {key}, {day:%Y-%m-%d}: {len(chosen)} grids, {len(frames)} full frames",
        flush=True,
    )
    if not chosen:
        return None

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, key, time_step=time_step, platform=PLATFORM)

    satellite = satellite_for(TIME_STEPS[0])
    edges = score_edges(chosen, orbit, day, satellite.instruments[0])
    points = epoch_points(chosen)
    differences = score_epochs(points, satellite, day)
    lattice = generate_points_fibonacci_lattice(POINT_SPACING)
    observed = observe_passes(lattice, chosen, day)
    ids = np.sort(observed.point_id.unique())
    observed["point_id"] = np.searchsorted(ids, observed.point_id)
    subset = lattice.iloc[ids].reset_index(drop=True)
    start = to_utc(day, chosen[0]["time"][0])[0] - SEARCH_WINDOW
    end = to_utc(day, chosen[-1]["time"][-1])[0] + SEARCH_WINDOW
    passes = match_passes(
        predict_passes(subset, satellite, 0, start, end, []), observed
    )
    print(json.dumps(passes), flush=True)
    chunks = get_chunks(frames)
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        area = pd.DataFrame(
            list(pool.map(lambda c: score_chunk(c, day, satellite_for), chunks))
        )
    return {
        "satellite": satellite,
        "edges": edges,
        "points": points,
        "differences": differences,
        "passes": passes,
        "area": area,
    }


def validate(key, days, tles, entry, catalogs, write) -> None:
    """Scores a configuration on one or more days, and writes its results if asked."""
    results = {}
    for day in days:
        catalog, products = catalogs[day]
        orbit = GeneralPerturbationsOrbit.from_tle(tles[day])
        result = measure(key, day, entry, orbit, catalog, products)
        if result is not None:
            results[day] = result
    if not results:
        return
    edges = {
        edge: np.concatenate([r["edges"][edge] for r in results.values()])
        for edge in ("near", "far")
    }
    for edge, values in edges.items():
        print(
            f"swath edge, {edge} range: median {np.median(values):.2f} km, p95 "
            f"{np.percentile(values, 95):.2f} km",
            flush=True,
        )
    differences = np.concatenate([r["differences"] for r in results.values()])
    found = np.isfinite(differences)
    epochs = summarize(differences[found])
    print(
        f"epochs: {int(found.sum())} of {len(differences)} points observed, median "
        f"{epochs['median']:+.3f} s, p95 {epochs['p95']:.3f} s",
        flush=True,
    )
    counts = {
        name: sum(r["passes"][name] for r in results.values())
        for name in ("matched", "observed", "predicted")
    }
    passes = {**compute_pass_scores(**counts), **counts}
    print(json.dumps(passes), flush=True)
    area = pd.concat([r["area"] for r in results.values()], ignore_index=True)
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
    if area.empty:
        print("swath IoU: no chunks of full frames", flush=True)
    else:
        for time_step in TIME_STEPS:
            print(
                f"swath IoU at {time_step} s over {len(area)} chunks: median "
                f"{area[f'iou_{time_step}'].median():.5f}, p5 "
                f"{area[f'iou_{time_step}'].quantile(0.05):.5f}",
                flush=True,
            )
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
    print(json.dumps(metrics), flush=True)
    if write:
        write_results(key, results, tles, entry, metrics)


def write_results(key, results, tles, entry, metrics) -> None:
    """Writes a configuration's record and fixture."""
    days = list(results)
    # the passes of the first day, with its element set: shared by the
    # configurations of a single day, and kept apart for pooled days, whose element
    # sets differ
    fixture_day = days[0]
    points = results[fixture_day]["points"][
        np.isfinite(results[fixture_day]["differences"])
    ]
    rng = np.random.default_rng(1)
    fixture = points[rng.choice(len(points), min(FIXTURE_PASSES, len(points)), False)]
    suffix = f"_{key}" if len(days) > 1 else ""
    update_fixture(
        FIXTURE_PATH,
        **{
            f"tle{suffix}": np.array(tles[fixture_day]),
            f"day{suffix}": np.array(fixture_day.isoformat()),
            f"passes_{key}": fixture,
        },
    )
    gaps = [get_holdout_gap(day, FIT_DAY) for day in days]
    update_record(
        RECORD_PATH,
        "lsar",
        key,
        {
            "digest": entry.configurations[key].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_lsar.py",
            "reference": "Geolocation grids, orbits, and bounding polygons of NISAR "
            "Level 1 RSLC products",
            "period": ", ".join(f"{day:%Y-%m-%d}" for day in days),
            "platform": PLATFORM,
            "holdout_gap_days": None if None in gaps else min(gaps),
            "metrics": metrics,
            "fixture_metrics": recompute_passes(
                results[days[0]]["satellite"], 0, fixture_day, fixture
            ),
        },
    )
    print(f"wrote {key} to {RECORD_PATH}")


def main() -> None:
    """Fits or validates the entry on one or more days of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("days", nargs="+")
    parser.add_argument("--configuration", choices=VALIDATED, action="append")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--download-only", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    days = [
        datetime.fromisoformat(day).replace(tzinfo=timezone.utc)
        for day in arguments.days
    ]
    catalogs = {day: load_day(day, arguments.directory) for day in days}
    if arguments.download_only:
        return
    tles = {day: get_element_set(name) for day, name in zip(days, arguments.days)}
    for day, tle in tles.items():
        print(f"{day:%Y-%m-%d}: element set {tle[0][18:32]}", flush=True)
    if arguments.fit:
        if len(days) > 1:
            raise ValueError("fit on one day")
        fit(
            catalogs[days[0]][1],
            GeneralPerturbationsOrbit.from_tle(tles[days[0]]),
            days[0],
        )
        return
    entry = get_entry("lsar")
    # the 77 MHz mode has its own held-out days
    for key in arguments.configuration or ["swath_20_40mhz", "swath_5mhz"]:
        validate(key, days, tles, entry, catalogs, arguments.write)


if __name__ == "__main__":
    main()
