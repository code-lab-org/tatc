"""
Helpers shared by the scripts that validate catalog entries against NOAA-20 data:
reading the public archive, geometry, and scoring through public TAT-C methods.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import io
import json
import pathlib
import re
import time as clock
from datetime import datetime, timedelta

import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
import requests
from pyproj import Transformer
from shapely import contains_xy, prepare

from tatc.analysis import (
    collect_observations,
    compute_ground_track,
    compute_sample_sizes,
)
from tatc.catalog import STANDARD_VERSION
from tatc.catalog.metrics import compute_mask_iou, compute_pass_scores, summarize
from tatc.schemas import CircularOrbit, Point
from tatc.utils import VelocityFrame, compute_projected_ray_position

BUCKET = "https://noaa-nesdis-n20-pds.s3.amazonaws.com/"
IET_EPOCH = pd.Timestamp("1958-01-01", tz="UTC")
LEAP_SECONDS = 37
PUBLISHED_EPOCH = pd.Timestamp("2026-10-01", tz="UTC").to_pydatetime()
SIZE_TOLERANCE = 0.03
GRID_STEP = 0.1

TO_ECEF = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
TO_GEODETIC = Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)

_SESSION = requests.Session()


def get_with_retries(url: str, **kwargs) -> requests.Response:
    """Sends a GET request, retrying transient failures with a backoff."""
    for attempt in range(6):
        try:
            response = _SESSION.get(url, timeout=180, **kwargs)
            response.raise_for_status()
            return response
        except requests.RequestException:
            if attempt == 5:
                raise
            clock.sleep(2**attempt)
    raise RuntimeError("unreachable")


def list_keys(prefix: str) -> list[str]:
    """Lists the keys of all objects in the archive with a prefix."""
    keys, token = [], None
    while True:
        params = {"list-type": "2", "prefix": prefix}
        if token:
            params["continuation-token"] = token
        listing = get_with_retries(BUCKET, params=params).text
        keys += re.findall(r"<Key>(.*?)</Key>", listing)
        match = re.search(
            r"<NextContinuationToken>(.*?)</NextContinuationToken>", listing
        )
        if match is None:
            return keys
        token = match.group(1)


def fetch(key: str) -> bytes:
    """Downloads an object of the archive."""
    return get_with_retries(BUCKET + key).content


def load_element_set(day: datetime, directory: pathlib.Path) -> list[str]:
    """Gets the operational element set with an epoch at the start of a day."""
    path = directory / f"n20_tle_{day:%Y%m%d}.txt"
    if not path.exists():
        key = next(
            k
            for k in list_keys(f"TLE-AUX/{day:%Y/%m/%d}/")
            if f"Z_{day:%Y%m%d}000000Z_" in k
        )
        with h5py.File(io.BytesIO(fetch(key))) as f:
            text = bytes(f["All_Data/TLE-AUX_All/Dataset_Array_0"][:]).decode()
        path.write_text(text.strip() + "\n", encoding="utf-8")
    return path.read_text(encoding="utf-8").splitlines()


def iet_to_utc(iet) -> pd.DatetimeIndex:
    """Converts IDPS Epoch Times (microseconds since 1958, with leap seconds) to UTC."""
    return IET_EPOCH + pd.to_timedelta(
        np.asarray(iet) - LEAP_SECONDS * 1_000_000, unit="us"
    )


def ecef(lon, lat, height=0.0) -> np.ndarray:
    """Gets the Earth-fixed position (m) of geodetic coordinates."""
    return np.stack(
        TO_ECEF.transform(lon, lat, np.broadcast_to(height, np.shape(lon))), axis=-1
    )


def up(lon, lat) -> np.ndarray:
    """Gets the geodetic up unit vector."""
    lon, lat = np.radians(lon), np.radians(lat)
    return np.stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1
    )


def normalize(x: np.ndarray) -> np.ndarray:
    """Gets unit vectors along the last axis."""
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def get_holdout_gap(day: datetime, fit_day: datetime) -> int | None:
    """Gets the days from the fit to a day, to which the standard applies its gap."""
    return None if day == fit_day else abs((day - fit_day).days)


def to_observations(satellite, instruments, epochs, lon, lat) -> gpd.GeoDataFrame:
    """Builds observations of points by named instruments of a satellite."""
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": instruments,
            "epoch": pd.DatetimeIndex(epochs),
        },
        geometry=gpd.points_from_xy(np.ravel(lon), np.ravel(lat)),
        crs="EPSG:4326",
    )


def compute_published_errors(entry, key: str, footprints) -> np.ndarray:
    """
    Gets the relative errors of sample sizes against published footprints, given as
    rows of instrument name, altitude (km), inclination (deg), scan angle (deg), and
    the sizes (km) along and across track. Gives two errors for each row.
    """
    errors = []
    for name, altitude, inclination, angle, along, cross in footprints:
        orbit = CircularOrbit(
            altitude=altitude * 1e3, inclination=inclination, epoch=PUBLISHED_EPOCH
        )
        satellite = entry.to_satellite(orbit, key)
        time = PUBLISHED_EPOCH + timedelta(minutes=10)
        target = compute_projected_ray_position(
            orbit.to_gp_orbit().get_orbit_track(time),
            0,
            0,
            angle,
            0,
            velocity_frame=VelocityFrame.INERTIAL,
        )
        sizes = compute_sample_sizes(
            to_observations(
                satellite,
                [name],
                [time],
                [target.longitude.degrees],
                [target.latitude.degrees],
            ),
            satellite,
        )
        errors += [
            sizes.sample_along.iloc[0] / 1e3 / along - 1,
            sizes.sample_cross.iloc[0] / 1e3 / cross - 1,
        ]
        print(f"{name}, {angle} deg: {errors[-2]:+.2%} along, {errors[-1]:+.2%} across")
    return np.array(errors)


def score_published(errors: np.ndarray) -> tuple[float, float]:
    """Gets the percent of dimensions within tolerance and the largest error (%)."""
    return (
        100 * float(np.mean(np.abs(errors) <= SIZE_TOLERANCE)),
        100 * float(np.abs(errors).max()),
    )


def compute_spacing_ratios(satellite, names, epochs, rows) -> tuple[float, float]:
    """
    Gets the median ratios of predicted to measured spacing along and across track,
    for rows of latitude, longitude, along-track (km), and cross-track spacing (km).
    """
    sizes = compute_sample_sizes(
        to_observations(satellite, names, epochs, rows[:, 1], rows[:, 0]), satellite
    )
    return (
        float(np.nanmedian(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 2])),
        float(np.nanmedian(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
    )


def score_passes(points, satellite, index, observed, start, end, gaps) -> dict:
    """
    Compares predicted and observed passes of points, given observations with
    `point_id` and `time`, and the periods without data as pairs of times.
    """

    def predict(item):
        point_id, geometry = item
        return collect_observations(
            Point(latitude=geometry.y, longitude=geometry.x),
            satellite,
            start,
            end,
            index,
        ).assign(point_id=point_id)

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        predicted = pd.concat(list(pool.map(predict, enumerate(points.geometry))))
    predicted = predicted.reset_index(drop=True)
    predicted["point_id"] = predicted.point_id.astype("int64")
    in_gap = np.zeros(len(predicted), dtype=bool)
    for first, last in gaps:
        in_gap |= ((predicted.epoch >= first) & (predicted.epoch <= last)).to_numpy()
    predicted = predicted[~in_gap]
    merged = pd.merge_asof(
        observed.assign(time=observed.time.dt.as_unit("ns")).sort_values("time"),
        predicted[["point_id", "epoch"]]
        .assign(epoch=predicted.epoch.dt.as_unit("ns"))
        .sort_values("epoch"),
        left_on="time",
        right_on="epoch",
        by="point_id",
        direction="nearest",
        tolerance=pd.Timedelta(minutes=10),
    )
    hit = merged.epoch.notna()
    differences = (merged.epoch - merged.time).dt.total_seconds()[hit]
    return {
        **compute_pass_scores(int(hit.sum()), len(observed), len(predicted)),
        "epoch_p95": summarize(differences)["p95"],
        "epoch_median": summarize(differences)["median"],
        "observed": len(observed),
        "predicted": len(predicted),
        "matched": int(hit.sum()),
    }


def recompute_passes(satellite, index, day, passes) -> dict[str, float]:
    """
    Recomputes observed passes, given as rows of latitude, longitude, and seconds from
    the start of a day. Gets the percent found and the epochs' 95th percentile difference.
    """
    differences = []
    for latitude, longitude, seconds in passes:
        time = day + timedelta(seconds=float(seconds))
        found = collect_observations(
            Point(latitude=float(latitude), longitude=float(longitude)),
            satellite,
            time - timedelta(minutes=10),
            time + timedelta(minutes=10),
            index,
        )
        differences.append(
            min(
                ((e - pd.Timestamp(time)).total_seconds() for e in found.epoch), key=abs
            )
            if not found.empty
            else np.nan
        )
    differences = np.array(differences, dtype=float)
    return {
        "passes_found": 100 * float(np.isfinite(differences).mean()),
        "epoch_p95": float(summarize(differences)["p95"]),
    }


def build_grid() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gets the longitudes, latitudes, and areas of the cells of a global grid."""
    lon, lat = np.meshgrid(
        np.arange(-180 + GRID_STEP / 2, 180, GRID_STEP),
        np.arange(-90 + GRID_STEP / 2, 90, GRID_STEP),
    )
    weights = np.cos(np.radians(lat)) * (GRID_STEP * 111.32) ** 2
    return lon.ravel(), lat.ravel(), weights.ravel()


def score_swath_area(
    satellite_for, index, first, last, time_steps, grid, compared, observed
) -> float:
    """Gets the lowest IoU of a swath with the observed area, over time steps."""
    lon, lat, weights = grid
    scores = []
    for time_step in time_steps:
        count = int((last - first).total_seconds() // time_step)
        times = [first + timedelta(seconds=time_step * (i + 0.5)) for i in range(count)]
        geometry = compute_ground_track(
            satellite_for(time_step), times, instrument_index=index
        ).geometry.iloc[0]
        prepare(geometry)
        predicted = np.zeros(len(lon), dtype=bool)
        predicted[compared] = contains_xy(geometry, lon[compared], lat[compared])
        scores.append(compute_mask_iou(observed, predicted, weights))
        print(f"  swath IoU at {time_step:.4g} s: {scores[-1]:.5f}", flush=True)
    return min(scores)


def update_fixture(path: pathlib.Path, **arrays) -> None:
    """Writes arrays to a fixture, keeping the others."""
    existing = {}
    if path.exists():
        with np.load(path, allow_pickle=False) as data:
            existing = {key: data[key] for key in data.files}
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **{**existing, **arrays})


def update_record(path: pathlib.Path, entry_id: str, key: str, results: dict) -> None:
    """Writes the results of a configuration to a record, keeping the others."""
    record = {"entry_id": entry_id, "configurations": {}}
    if path.exists():
        record = json.loads(path.read_text(encoding="utf-8"))
    record["standard_version"] = STANDARD_VERSION
    record["configurations"][key] = results
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
