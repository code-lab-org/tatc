"""
Validates the ABI catalog entry against GOES-19 full-disk products, using the exact
values of the entry through public TAT-C methods: whether the full disk's imaged pixels
are observed, with the satellite's elevation angle, and the local zenith angles to
which the sea and land surface temperature products are of good quality. The values
are derived from the products' declarations, not fitted.

Reads anonymous public data, and needs the `validation` extra:

    python validation/catalog/validate_catalog_abi.py <data directory> [TIME]

It compares the products of the full-disk scan starting in the hour of TIME (UTC, as
2026-10-03T12:00Z), by default the time of the entry's record. `--configuration`
tests only some configurations, on the same samples. With `--write`, it writes the
entry's validation record and the reduced fixture from which unit tests recompute
the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import pathlib
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta

import numpy as np
import xarray as xr
from skyfield.framelib import itrs

from tatc.analysis import collect_observations
from tatc.catalog.metrics import compute_product_limit_error
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit, Point
from tatc.utils import compute_field_of_regard
from validation_io import get_with_retries, parse_run, write_results

BUCKET = "https://noaa-goes19.s3.amazonaws.com/"
SAMPLES = 2000
FIXTURE_SAMPLES = 60
PRODUCTS = {
    "sst_good_quality": ("ABI-L2-SSTF", "SST", 4),
    "lst_good_quality": ("ABI-L2-LSTF", "LST", 1),
}
WGS84_A, WGS84_F = 6378137.0, 1 / 298.257223563


def download(product: str, time: datetime, directory: pathlib.Path) -> pathlib.Path:
    """Downloads the first full-disk file of a product starting in the hour of a time."""
    prefix = f"{product}/{time:%Y}/{time.timetuple().tm_yday:03d}/{time:%H}/"
    root = ET.fromstring(
        get_with_retries(BUCKET, params={"list-type": "2", "prefix": prefix}).content
    )
    namespace = "{http://s3.amazonaws.com/doc/2006-03-01/}"
    key = root.find(f"{namespace}Contents/{namespace}Key").text
    path = directory / key.split("/")[-1]
    if not path.exists():
        path.write_bytes(get_with_retries(BUCKET + key).content)
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


def satellite_position(orbit: GeneralPerturbationsOrbit, time: datetime) -> np.ndarray:
    """Gets the Earth-fixed position (m) of the satellite at a time."""
    return np.array(
        orbit.elements[0]
        .to_skyfield()
        .at(timescale.from_datetime(time))
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


def predict(satellite, time, latitude, longitude) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets whether each pixel is observed in 10 minutes from a time, and the
    satellite's elevation.
    """

    def evaluate(index: int) -> tuple[bool, float]:
        observations = collect_observations(
            Point(latitude=float(latitude[index]), longitude=float(longitude[index])),
            satellite,
            time,
            time + timedelta(minutes=10),
        )
        if observations.empty:
            return False, np.nan
        return True, float(observations.sat_alt.iloc[0])

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        results = list(pool.map(evaluate, range(len(latitude))))
    return np.array([r[0] for r in results]), np.array([r[1] for r in results])


def pick_near_bound(rng, zenith: np.ndarray, bound: float) -> np.ndarray:
    """
    Picks the samples of a fixture among those closest to a bound of the local
    zenith angle on each side, which define its error.
    """
    order = np.argsort(np.abs(zenith - bound))[: FIXTURE_SAMPLES * 4]
    return rng.choice(order, FIXTURE_SAMPLES, replace=False)


def main() -> None:
    """Validates the entry on a full-disk scan, and writes its results if requested."""
    run = parse_run("abi", __doc__)
    if run.fit:
        sys.exit("ABI's values are derived from its products' declarations, not fitted")
    time = run.day
    tle = run.element_set(time.replace(hour=0, minute=0, second=0, microsecond=0))
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    position = satellite_position(orbit, time)
    # the samples of every configuration are drawn in turn, so that a configuration
    # is tested on the same samples whichever others are tested
    rng = np.random.default_rng(1)
    fixture = {"tle": np.array(tle), "time": np.array(time.isoformat())}
    results = {}

    # earth disk: elevation angles, and which imaged pixels are observed
    key = "full_disk"
    with xr.open_dataset(download("ABI-L2-ACMF", time, run.directory)) as dataset:
        reduced = dataset.isel(x=slice(0, None, 8), y=slice(0, None, 8))
        latitude, longitude = fixed_grid(reduced)
        on_grid = np.isfinite(latitude)
        imaged = np.isfinite(reduced.ACM.values)
    picks = rng.choice(np.flatnonzero(on_grid & imaged), SAMPLES, replace=False)
    lat, lon = latitude.flat[picks], longitude.flat[picks]
    zenith = compute_zenith(position, lat, lon)
    keep = rng.choice(SAMPLES, FIXTURE_SAMPLES, replace=False)
    if key in run.configurations:
        satellite = run.entry.to_satellite(
            orbit, configuration=key, platform=run.platform
        )
        observed, elevation = predict(satellite, time, lat, lon)
        error = np.abs(elevation - (90 - zenith))
        fixture[key] = np.stack([lat[keep], lon[keep], zenith[keep]], axis=-1)
        results[key] = {
            "reference": "ABI-L2-ACMF clear sky mask, every 8th pixel",
            "metrics": [
                {"name": "imaged_pixels_observed", "value": 100 * float(observed.mean())},
                {"name": "elevation_error_max", "value": float(np.nanmax(error))},
            ],
            "fixture_metrics": {
                "imaged_pixels_observed": 100 * float(observed[keep].mean()),
                "elevation_error_max": float(np.nanmax(error[keep])),
            },
        }
        print(json.dumps(results[key]["metrics"]), flush=True)

    # product coverage: the boundary of predicted observations near declared bounds
    for key, (product, variable, step) in PRODUCTS.items():
        with xr.open_dataset(download(product, time, run.directory)) as dataset:
            reduced = dataset.isel(x=slice(0, None, step), y=slice(0, None, step))
            latitude, longitude = fixed_grid(reduced)
            valid = np.isfinite(latitude)
            retrieved = np.isfinite(reduced[variable].values) & valid
            good = retrieved & (reduced.DQF.values == 0)
            bound = float(reduced.quantitative_local_zenith_angle_bounds[1])
        expected = compute_field_of_regard(35786e3, 90 - bound)
        configured = run.entry.configurations[key].instruments[0]
        if abs(configured.parameters["field_of_regard"] - expected) > 1e-9:
            raise ValueError(f"{key}: the entry does not match the bound of {bound} deg")
        zenith = np.full(latitude.shape, np.nan)
        zenith[valid] = compute_zenith(position, latitude[valid], longitude[valid])
        picks = rng.choice(
            np.flatnonzero(retrieved & (np.abs(zenith - bound) <= 5)),
            SAMPLES,
            replace=False,
        )
        lat, lon = latitude.flat[picks], longitude.flat[picks]
        picked_zenith, picked_good = zenith.flat[picks], good.flat[picks]
        keep = pick_near_bound(rng, picked_zenith, bound)
        if key not in run.configurations:
            continue
        satellite = run.entry.to_satellite(
            orbit, configuration=key, platform=run.platform
        )
        observed, _ = predict(satellite, time, lat, lon)
        fixture[key] = np.stack([lat[keep], lon[keep], picked_zenith[keep]], axis=-1)
        fixture[f"{key}_bound"] = np.array(bound)
        results[key] = {
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
            "fixture_metrics": {
                "product_limit": compute_product_limit_error(
                    picked_zenith[keep], observed[keep], bound
                )
            },
        }
        print(json.dumps(results[key]["metrics"]), flush=True)

    if not run.write:
        return
    # the fixture, then the records
    for index, (key, result) in enumerate(results.items()):
        write_results(
            run,
            key,
            **result,
            fixture=fixture if index == 0 else None,
            period=f"{time:%Y-%m-%dT%H:%MZ}",
        )


if __name__ == "__main__":
    main()
