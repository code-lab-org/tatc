"""
Validates the AHI catalog entry against Himawari-9's full-disk cloud products of
NOAA's Enterprise algorithms, using the exact values of the entry through public
TAT-C methods, as `validate_catalog_abi.py` does for ABI. The cloud top height is
retrieved only within a local zenith angle that the product does not declare, so it
is measured on one scan (run with `--fit`), with the nominal position of the fixed
grid, and tested on a scan at least 30 days from it, with a public element set.

Reads anonymous public data (about 1.1 GB a scan), and needs the `validation` extra:

    python validation/catalog/validate_catalog_ahi.py <data directory> [TIME] --fit
    python validation/catalog/validate_catalog_ahi.py <data directory> [TIME]

It reads the scan starting at TIME (UTC, as 2026-11-08T12:00Z), by default the time
of the entry's record, or with `--fit` that of the fit. `--configuration` tests only
some configurations, on the same samples. With `--write`, it writes the entry's
validation record and the reduced fixture from which unit tests recompute the
results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pathlib
import xml.etree.ElementTree as ET
from datetime import datetime

import numpy as np
import xarray as xr

from tatc.catalog.metrics import compute_product_limit_error
from tatc.schemas import GeneralPerturbationsOrbit
from tatc.utils import compute_field_of_regard
from validate_catalog_abi import (
    FIXTURE_SAMPLES,
    SAMPLES,
    WGS84_A,
    compute_zenith,
    pick_near_bound,
    predict,
    satellite_position,
)
from validation_io import get_with_retries, parse_run, write_results

BUCKET = "https://noaa-himawari9.s3.amazonaws.com/"
PREFIX = "AHI-L2-FLDK-Clouds"
# the nominal position of the fixed grid: longitude (deg) and height (m)
GRID_LONGITUDE, GRID_HEIGHT = 140.7, 35785863.0
CHUNK_ROWS = 250


def download(product: str, time: datetime, directory: pathlib.Path) -> pathlib.Path:
    """Downloads the full-disk file of a cloud product of the scan starting at a time."""
    prefix = f"{PREFIX}/{time:%Y/%m/%d/%H%M}/"
    root = ET.fromstring(
        get_with_retries(BUCKET, params={"list-type": "2", "prefix": prefix}).content
    )
    namespace = "{http://s3.amazonaws.com/doc/2006-03-01/}"
    key = next(
        element.text
        for element in root.iter(f"{namespace}Key")
        if f"/AHI-{product}_" in element.text
    )
    path = directory / key.split("/")[-1]
    if not path.exists():
        with get_with_retries(BUCKET + key, stream=True) as response:
            with open(path, "wb") as file:
                for block in response.iter_content(2**22):
                    file.write(block)
    return path


def grid_position() -> np.ndarray:
    """Gets the Earth-fixed position (m) of the fixed grid's nominal position."""
    radius = WGS84_A + GRID_HEIGHT
    longitude = np.radians(GRID_LONGITUDE)
    return np.array([radius * np.cos(longitude), radius * np.sin(longitude), 0.0])


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


def fit(time: datetime, directory: pathlib.Path) -> None:
    """
    Prints the greatest local zenith angle at which the cloud top height is retrieved
    on the fitting scan, at full resolution, with the fixed grid's nominal position,
    and the field of regard it gives.
    """
    path = download("CHGT", time, directory)
    _, _, zenith, retrieved = read_pixels(path, "CldTopHght", grid_position(), 1)
    bound = float(zenith[retrieved].max())
    print(
        f"{time:%Y-%m-%dT%H:%MZ}: cloud top height retrieved in "
        f"{retrieved.mean():.1%} of the disk's pixels, up to a local zenith angle of "
        f"{bound:.4f} deg (99.9th percentile "
        f"{np.percentile(zenith[retrieved], 99.9):.4f}); field of regard "
        f"{compute_field_of_regard(GRID_HEIGHT, 90 - bound):.6f} deg"
    )


def validate(run) -> None:
    """Validates the entry on a held-out scan, and writes its results if requested."""
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
    lat, lon, zenith, imaged = read_pixels(
        download("CMSK", time, run.directory), "CloudMask", position, 8
    )
    picks = rng.choice(np.flatnonzero(imaged), SAMPLES, replace=False)
    keep = rng.choice(SAMPLES, FIXTURE_SAMPLES, replace=False)
    if key in run.configurations:
        satellite = run.entry.to_satellite(
            orbit, configuration=key, platform=run.platform
        )
        observed, elevation = predict(satellite, time, lat[picks], lon[picks])
        error = np.abs(elevation - (90 - zenith[picks]))
        fixture[key] = np.stack(
            [lat[picks][keep], lon[picks][keep], zenith[picks][keep]], axis=-1
        )
        results[key] = {
            "reference": "AHI-CMSK cloud mask, every 8th pixel",
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

    # cloud top height: the boundary of predicted observations near the greatest
    # local zenith angle at which the product retrieves it on this scan
    key = "cloud_top_height"
    lat, lon, zenith, retrieved = read_pixels(
        download("CHGT", time, run.directory), "CldTopHght", position, 1
    )
    bound = float(zenith[retrieved].max())
    picks = rng.choice(
        np.flatnonzero(np.abs(zenith - bound) <= 5), SAMPLES, replace=False
    )
    picked_zenith = zenith[picks]
    keep = pick_near_bound(rng, picked_zenith, bound)
    if key in run.configurations:
        satellite = run.entry.to_satellite(
            orbit, configuration=key, platform=run.platform
        )
        observed, _ = predict(satellite, time, lat[picks], lon[picks])
        fixture[key] = np.stack(
            [lat[picks][keep], lon[picks][keep], picked_zenith[keep]], axis=-1
        )
        fixture[f"{key}_bound"] = np.array(bound)
        results[key] = {
            "reference": "AHI-CHGT cloud top height, the greatest local zenith angle "
            "at which it is retrieved",
            "metrics": [
                {
                    "name": "product_limit",
                    "value": compute_product_limit_error(
                        picked_zenith, observed, bound
                    ),
                },
                {"name": "retrieval_bound", "value": bound},
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


def main() -> None:
    """Fits the bound of the cloud top height, or validates the entry."""
    run = parse_run("ahi", __doc__)
    if run.fit:
        fit(run.day, run.directory)
    else:
        validate(run)


if __name__ == "__main__":
    main()
