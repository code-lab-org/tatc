"""
Validates the ABI catalog entry against GOES-19 products, using the exact values
of the entry through public TAT-C methods, as `ValidateGeostationaryGOES` does
with values of its own. Writes the entry's validation record and the reduced
fixture from which unit tests recompute the results.

Reads anonymous public data, and needs the `validation` extra:

    python docs/validation/validate_catalog_abi.py <empty data directory>

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import pathlib
import sys
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
from validation_io import add_earth_orientation

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "abi.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "abi.npz"

BUCKET = "https://noaa-goes19.s3.amazonaws.com/"
TIME = datetime(2026, 10, 3, 12, tzinfo=timezone.utc)
TLE = (
    "1 60133U 24119A   26277.37537957 -.00000263  00000+0  00000+0 0  9990",
    "2 60133   0.0344  24.5665 0000927 249.5134 158.9494  1.00270254  8054",
)
SAMPLES = 2000
FIXTURE_SAMPLES = 60
PRODUCTS = {
    "sst_good_quality": ("ABI-L2-SSTF", "SST", 4),
    "lst_good_quality": ("ABI-L2-LSTF", "LST", 1),
}


def download(product: str, directory: pathlib.Path) -> pathlib.Path:
    """Downloads the first full-disk file of a product starting in the hour of TIME."""
    prefix = f"{product}/{TIME:%Y}/{TIME.timetuple().tm_yday:03d}/{TIME:%H}/"
    root = ET.fromstring(
        requests.get(
            BUCKET, params={"list-type": "2", "prefix": prefix}, timeout=120
        ).content
    )
    namespace = "{http://s3.amazonaws.com/doc/2006-03-01/}"
    key = root.find(f"{namespace}Contents/{namespace}Key").text
    path = directory / key.split("/")[-1]
    if not path.exists():
        response = requests.get(BUCKET + key, timeout=600)
        response.raise_for_status()
        path.write_bytes(response.content)
    return path


def fixed_grid(dataset: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """Gets the geodetic latitude and longitude (deg) of ABI fixed grid pixels."""
    projection = dataset.goes_imager_projection.attrs
    r_eq, r_pol = projection["semi_major_axis"], projection["semi_minor_axis"]
    h = projection["perspective_point_height"] + r_eq
    lon_0 = np.radians(projection["longitude_of_projection_origin"])
    x, y = np.meshgrid(dataset.x.values, dataset.y.values)
    a = np.sin(x) ** 2 + np.cos(x) ** 2 * (
        np.cos(y) ** 2 + (r_eq / r_pol) ** 2 * np.sin(y) ** 2
    )
    b = -2 * h * np.cos(x) * np.cos(y)
    c = h**2 - r_eq**2
    with np.errstate(invalid="ignore"):
        r_s = (-b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)
        s_x = r_s * np.cos(x) * np.cos(y)
        s_y = -r_s * np.sin(x)
        s_z = r_s * np.cos(x) * np.sin(y)
        latitude = np.degrees(
            np.arctan((r_eq / r_pol) ** 2 * s_z / np.sqrt((h - s_x) ** 2 + s_y**2))
        )
        longitude = np.degrees(lon_0 - np.arctan(s_y / (h - s_x)))
    return latitude, longitude


def compute_zenith(
    orbit: GeneralPerturbationsOrbit, latitude: np.ndarray, longitude: np.ndarray
) -> np.ndarray:
    """Gets the local zenith angle (deg) of the satellite at geodetic positions."""
    satellite = np.array(
        orbit.elements[0]
        .to_skyfield()
        .at(timescale.from_datetime(TIME))
        .frame_xyz(itrs)
        .m
    )
    lat, lon = np.radians(latitude), np.radians(longitude)
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    position = np.stack(
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
    line_of_sight = satellite - position
    return 90 - np.degrees(
        np.arcsin(
            np.sum(up * line_of_sight, axis=-1) / np.linalg.norm(line_of_sight, axis=-1)
        )
    )


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


def main(directory: pathlib.Path) -> None:
    """Validates the entry, then writes its record and fixture."""
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    entry = get_entry("abi")
    orbit = GeneralPerturbationsOrbit.from_tle(list(TLE))
    common = {
        "evidence": "rerun",
        "tatc_version": tatc.__version__,
        "source": "docs/validation/validate_catalog_abi.py",
        "period": f"{TIME:%Y-%m-%dT%H:%MZ}",
        "platform": "goes-19",
    }
    fixture = {"tle": np.array(TLE), "time": np.array(TIME.isoformat())}
    records = {}

    # earth disk: elevation angles, and which imaged pixels are observed
    with xr.open_dataset(download("ABI-L2-ACMF", directory)) as dataset:
        reduced = dataset.isel(x=slice(0, None, 8), y=slice(0, None, 8))
        latitude, longitude = fixed_grid(reduced)
        on_grid = np.isfinite(latitude)
        imaged = np.isfinite(reduced.ACM.values)
    picks = rng.choice(np.flatnonzero(on_grid & imaged), SAMPLES, replace=False)
    lat, lon = latitude.flat[picks], longitude.flat[picks]
    zenith = compute_zenith(orbit, lat, lon)
    satellite = entry.to_satellite(orbit, "full_disk", platform="goes-19")
    observed, elevation = predict(satellite, lat, lon)
    records["full_disk"] = {
        **common,
        "digest": entry.configurations["full_disk"].digest,
        "reference": "ABI-L2-ACMF clear sky mask, every 8th pixel",
        "metrics": [
            {"name": "imaged_pixels_observed", "value": 100 * float(observed.mean())},
            {
                "name": "elevation_error_max",
                "value": float(np.nanmax(np.abs(elevation - (90 - zenith)))),
            },
        ],
    }
    keep = rng.choice(SAMPLES, FIXTURE_SAMPLES, replace=False)
    fixture["full_disk"] = np.stack([lat[keep], lon[keep], zenith[keep]], axis=-1)
    records["full_disk"]["fixture_metrics"] = {
        "imaged_pixels_observed": 100 * float(observed[keep].mean()),
        "elevation_error_max": float(
            np.nanmax(np.abs(elevation[keep] - (90 - zenith[keep])))
        ),
    }

    # product coverage: the boundary of predicted observations near declared bounds
    for key, (product, variable, step) in PRODUCTS.items():
        with xr.open_dataset(download(product, directory)) as dataset:
            reduced = dataset.isel(x=slice(0, None, step), y=slice(0, None, step))
            latitude, longitude = fixed_grid(reduced)
            valid = np.isfinite(latitude)
            retrieved = np.isfinite(reduced[variable].values) & valid
            good = retrieved & (reduced.DQF.values == 0)
            bound = float(reduced.quantitative_local_zenith_angle_bounds[1])
        expected = compute_field_of_regard(35786e3, 90 - bound)
        configured = entry.configurations[key].instruments[0]
        if abs(configured.parameters["field_of_regard"] - expected) > 1e-9:
            raise ValueError(
                f"{key}: the entry does not match the bound of {bound} deg"
            )
        zenith = np.full(latitude.shape, np.nan)
        zenith[valid] = compute_zenith(orbit, latitude[valid], longitude[valid])
        picks = rng.choice(
            np.flatnonzero(retrieved & (np.abs(zenith - bound) <= 5)),
            SAMPLES,
            replace=False,
        )
        lat, lon = latitude.flat[picks], longitude.flat[picks]
        picked_zenith, picked_good = zenith.flat[picks], good.flat[picks]
        satellite = entry.to_satellite(orbit, key, platform="goes-19")
        observed, _ = predict(satellite, lat, lon)
        records[key] = {
            **common,
            "digest": entry.configurations[key].digest,
            "reference": f"{product} {variable} and quality flags",
            "metrics": [
                {
                    "name": "product_limit",
                    "value": compute_product_limit_error(
                        picked_zenith, observed, bound
                    ),
                },
                {
                    "name": "good_pixels_observed",
                    "value": 100 * float(observed[picked_good].mean()),
                },
            ],
        }
        # the pixels closest to the bound on each side define its error
        order = np.argsort(np.abs(picked_zenith - bound))[: FIXTURE_SAMPLES * 4]
        keep = rng.choice(order, FIXTURE_SAMPLES, replace=False)
        fixture[key] = np.stack([lat[keep], lon[keep], picked_zenith[keep]], axis=-1)
        fixture[f"{key}_bound"] = np.array(bound)
        records[key]["fixture_metrics"] = {
            "product_limit": compute_product_limit_error(
                picked_zenith[keep], observed[keep], bound
            )
        }

    FIXTURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(FIXTURE_PATH, **add_earth_orientation(fixture))
    RECORD_PATH.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "entry_id": "abi",
        "standard_version": STANDARD_VERSION,
        "configurations": records,
    }
    RECORD_PATH.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: validate_catalog_abi.py <data directory>")
    main(pathlib.Path(sys.argv[1]))
