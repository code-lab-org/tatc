"""
Validates the AHI catalog entry against Himawari-9's full-disk cloud products of
NOAA's Enterprise algorithms, using the exact values of the entry through public
TAT-C methods, as `validate_catalog_abi.py` does for ABI. The cloud top height is
retrieved only within a local zenith angle that the product does not declare, so it
is measured on one scan (run with `--fit`), with the nominal position of the fixed
grid, and tested on a scan at least 30 days from it, with a public element set.

Reads anonymous public data (about 1.1 GB a scan), and needs the `validation` extra:

    python docs/validation/validate_catalog_ahi.py <data directory> --fit
    python docs/validation/validate_catalog_ahi.py <data directory> --write

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta, timezone

import numpy as np
import requests
import xarray as xr
from skyfield.framelib import itrs

import tatc
from tatc.analysis import collect_observations
from tatc.catalog import STANDARD_VERSION, get_entry
from tatc.catalog.metrics import compute_product_limit_error
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit, Point
from tatc.utils import compute_field_of_regard
from validation_io import add_earth_orientation, get_holdout_gap

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "ahi.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "ahi.npz"

BUCKET = "https://noaa-himawari9.s3.amazonaws.com/"
PREFIX = "AHI-L2-FLDK-Clouds"
FIT_TIME = datetime(2026, 10, 1, 12, tzinfo=timezone.utc)
TIME = datetime(2026, 11, 8, 12, tzinfo=timezone.utc)
# the public element set from CelesTrak, with an epoch within a day of TIME
TLE = (
    "1 41836U 16064A   26313.19704546 -.00000266  00000+0  00000+0 0  9993",
    "2 41836   0.0374  36.5534 0000882 218.8169   4.7319  1.00274517 36669",
)
# the nominal position of the fixed grid: longitude (deg) and height (m)
GRID_LONGITUDE, GRID_HEIGHT = 140.7, 35785863.0
WGS84_A, WGS84_F = 6378137.0, 1 / 298.257223563
SAMPLES = 2000
FIXTURE_SAMPLES = 60
CHUNK_ROWS = 250


def download(product: str, time: datetime, directory: pathlib.Path) -> pathlib.Path:
    """Downloads the full-disk file of a cloud product of the scan starting at a time."""
    prefix = f"{PREFIX}/{time:%Y/%m/%d/%H%M}/"
    root = ET.fromstring(
        requests.get(
            BUCKET, params={"list-type": "2", "prefix": prefix}, timeout=120
        ).content
    )
    namespace = "{http://s3.amazonaws.com/doc/2006-03-01/}"
    key = next(
        element.text
        for element in root.iter(f"{namespace}Key")
        if f"/AHI-{product}_" in element.text
    )
    path = directory / key.split("/")[-1]
    if not path.exists():
        with requests.get(BUCKET + key, timeout=600, stream=True) as response:
            response.raise_for_status()
            with open(path, "wb") as file:
                for block in response.iter_content(2**22):
                    file.write(block)
    return path


def satellite_position(orbit: GeneralPerturbationsOrbit | None) -> np.ndarray:
    """
    Gets the Earth-fixed position (m) of the satellite at TIME from an orbit, or of
    the fixed grid's nominal position without one.
    """
    if orbit is None:
        radius = WGS84_A + GRID_HEIGHT
        longitude = np.radians(GRID_LONGITUDE)
        return np.array([radius * np.cos(longitude), radius * np.sin(longitude), 0.0])
    return np.array(
        orbit.elements[0]
        .to_skyfield()
        .at(timescale.from_datetime(TIME))
        .frame_xyz(itrs)
        .m
    )


def compute_zenith(position: np.ndarray, latitude, longitude) -> np.ndarray:
    """Gets the local zenith angle (deg) of a satellite at geodetic positions."""
    lat, lon = np.radians(latitude), np.radians(longitude)
    e2 = WGS84_F * (2 - WGS84_F)
    n = WGS84_A / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    point = np.stack(
        [
            n * np.cos(lat) * np.cos(lon),
            n * np.cos(lat) * np.sin(lon),
            n * (1 - e2) * np.sin(lat),
        ],
        axis=-1,
    )
    up = np.stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1
    )
    sight = position - point
    return 90 - np.degrees(
        np.arcsin(np.sum(up * sight, axis=-1) / np.linalg.norm(sight, axis=-1))
    )


def read_pixels(path: pathlib.Path, variable: str, position: np.ndarray, step: int):
    """
    Gets the latitude, longitude, local zenith angle, and whether a variable is
    retrieved, of the pixels of the Earth's disk, every `step` rows and columns.
    """
    lats, lons, zeniths, retrieved = [], [], [], []
    with xr.open_dataset(path) as dataset:
        rows = dataset.sizes["Rows"]
        for start in range(0, rows, CHUNK_ROWS * step):
            part = dataset.isel(
                Rows=slice(start, min(start + CHUNK_ROWS * step, rows), step),
                Columns=slice(None, None, step),
            )
            lat = part.Latitude.values.ravel()
            lon = part.Longitude.values.ravel()
            value = part[variable].values.ravel()
            disk = np.isfinite(lat) & (np.abs(lat) <= 90) & np.isfinite(lon)
            lats.append(lat[disk])
            lons.append(lon[disk])
            zeniths.append(compute_zenith(position, lat[disk], lon[disk]))
            retrieved.append(np.isfinite(value[disk]))
    return tuple(np.concatenate(x) for x in (lats, lons, zeniths, retrieved))


def predict(satellite, latitude, longitude) -> tuple[np.ndarray, np.ndarray]:
    """Gets whether each pixel is observed in 10 minutes, and the satellite's elevation."""

    def evaluate(index: int) -> tuple[bool, float]:
        observations = collect_observations(
            Point(latitude=float(latitude[index]), longitude=float(longitude[index])),
            satellite,
            TIME,
            TIME + timedelta(minutes=10),
        )
        if observations.empty:
            return False, np.nan
        return True, float(observations.sat_alt.iloc[0])

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        results = list(pool.map(evaluate, range(len(latitude))))
    return np.array([r[0] for r in results]), np.array([r[1] for r in results])


def fit(directory: pathlib.Path) -> None:
    """
    Prints the greatest local zenith angle at which the cloud top height is retrieved
    on the fitting scan, at full resolution, with the fixed grid's nominal position,
    and the field of regard it gives.
    """
    path = download("CHGT", FIT_TIME, directory)
    _, _, zenith, retrieved = read_pixels(
        path, "CldTopHght", satellite_position(None), 1
    )
    bound = float(zenith[retrieved].max())
    print(
        f"{FIT_TIME:%Y-%m-%dT%H:%MZ}: cloud top height retrieved in "
        f"{retrieved.mean():.1%} of the disk's pixels, up to a local zenith angle of "
        f"{bound:.4f} deg (99.9th percentile "
        f"{np.percentile(zenith[retrieved], 99.9):.4f}); field of regard "
        f"{compute_field_of_regard(GRID_HEIGHT, 90 - bound):.6f} deg"
    )


def validate(directory: pathlib.Path, write: bool) -> None:
    """Validates the entry on the held-out scan, and writes its record and fixture."""
    rng = np.random.default_rng(1)
    entry = get_entry("ahi")
    orbit = GeneralPerturbationsOrbit.from_tle(list(TLE))
    position = satellite_position(orbit)
    common = {
        "evidence": "rerun",
        "tatc_version": tatc.__version__,
        "source": "docs/validation/validate_catalog_ahi.py",
        "period": f"{TIME:%Y-%m-%dT%H:%MZ}",
        "platform": "himawari-9",
    }
    fixture = {"tle": np.array(TLE), "time": np.array(TIME.isoformat())}
    records = {}

    # earth disk: elevation angles, and which imaged pixels are observed
    lat, lon, zenith, imaged = read_pixels(
        download("CMSK", TIME, directory), "CloudMask", position, 8
    )
    picks = rng.choice(np.flatnonzero(imaged), SAMPLES, replace=False)
    satellite = entry.to_satellite(orbit, "full_disk", platform="himawari-9")
    observed, elevation = predict(satellite, lat[picks], lon[picks])
    error = np.abs(elevation - (90 - zenith[picks]))
    records["full_disk"] = {
        **common,
        "digest": entry.configurations["full_disk"].digest,
        "reference": "AHI-CMSK cloud mask, every 8th pixel",
        "metrics": [
            {"name": "imaged_pixels_observed", "value": 100 * float(observed.mean())},
            {"name": "elevation_error_max", "value": float(np.nanmax(error))},
        ],
    }
    keep = rng.choice(SAMPLES, FIXTURE_SAMPLES, replace=False)
    fixture["full_disk"] = np.stack(
        [lat[picks][keep], lon[picks][keep], zenith[picks][keep]], axis=-1
    )
    records["full_disk"]["fixture_metrics"] = {
        "imaged_pixels_observed": 100 * float(observed[keep].mean()),
        "elevation_error_max": float(np.nanmax(error[keep])),
    }
    print(json.dumps(records["full_disk"]["metrics"]), flush=True)

    # cloud top height: the boundary of predicted observations near the greatest
    # local zenith angle at which the product retrieves it on this scan
    key = "cloud_top_height"
    lat, lon, zenith, retrieved = read_pixels(
        download("CHGT", TIME, directory), "CldTopHght", position, 1
    )
    bound = float(zenith[retrieved].max())
    picks = rng.choice(
        np.flatnonzero(np.abs(zenith - bound) <= 5), SAMPLES, replace=False
    )
    satellite = entry.to_satellite(orbit, key, platform="himawari-9")
    observed, _ = predict(satellite, lat[picks], lon[picks])
    picked_zenith = zenith[picks]
    records[key] = {
        **common,
        "digest": entry.configurations[key].digest,
        "reference": "AHI-CHGT cloud top height, the greatest local zenith angle "
        "at which it is retrieved",
        "holdout_gap_days": get_holdout_gap(TIME, FIT_TIME),
        "metrics": [
            {
                "name": "product_limit",
                "value": compute_product_limit_error(picked_zenith, observed, bound),
            },
            {"name": "retrieval_bound", "value": bound},
        ],
    }
    order = np.argsort(np.abs(picked_zenith - bound))[: FIXTURE_SAMPLES * 4]
    keep = rng.choice(order, FIXTURE_SAMPLES, replace=False)
    fixture[key] = np.stack(
        [lat[picks][keep], lon[picks][keep], picked_zenith[keep]], axis=-1
    )
    fixture[f"{key}_bound"] = np.array(bound)
    records[key]["fixture_metrics"] = {
        "product_limit": compute_product_limit_error(
            picked_zenith[keep], observed[keep], bound
        )
    }
    print(json.dumps(records[key]["metrics"]), flush=True)
    if not write:
        return
    FIXTURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(FIXTURE_PATH, **add_earth_orientation(fixture))
    record = {
        "entry_id": "ahi",
        "standard_version": STANDARD_VERSION,
        "configurations": records,
    }
    RECORD_PATH.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


def main() -> None:
    """Fits the bound of the cloud top height, or validates the entry."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    if arguments.fit:
        fit(arguments.directory)
    else:
        validate(arguments.directory, arguments.write)


if __name__ == "__main__":
    main()
