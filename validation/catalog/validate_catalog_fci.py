"""
Validates the FCI catalog entry against Meteosat-12's (MTG-I1's) full-disk products,
using the exact values of the entry through public TAT-C methods, as for ABI (see
`validate_catalog_abi.py`, whose geometry this script shares): whether the full disk's
imaged pixels (of the cloud mask, EO:EUM:DAT:0678) are observed, with the satellite's
elevation angle; and the local zenith angles to which the land surface temperature
(LSA SAF, EO:EUM:DAT:1088) is retrieved and the sea surface temperature (OSI SAF Level
3C, EO:EUM:DAT:0694) is of acceptable quality or better (quality level 4 or 5).

The products do not declare these bounds, which are sharp in their data: `--fit`
measures them on a scan, from the satellite's nominal position (and, for the sea
surface temperature, the zenith angles the product gives), and running on a scan at
least 30 days from the fit tests them on data held out from it.

Needs the `validation` extra and a login to the EUMETSAT Data Store
(`eumdac auth --interactive`):

    python validation/catalog/validate_catalog_fci.py <data directory> 2026-10-16T12:00Z --fit
    python validation/catalog/validate_catalog_fci.py <data directory> 2026-11-17T12:00Z

It compares the products of the scan starting at TIME (UTC), by default the time of
the entry's record. The element set is the one nearest the time in the files of
`elements/` (or in a file given with `--elements`). `--configuration` tests only some
configurations, on the same samples. With `--write`, it writes the entry's validation
record and the reduced fixture from which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pathlib
from datetime import datetime, timedelta

import netCDF4
import numpy as np
from pyproj import Proj

from tatc.catalog.metrics import compute_product_limit_error
from tatc.schemas import GeneralPerturbationsOrbit
from tatc.utils import compute_field_of_regard
from validate_catalog_abi import (
    FIXTURE_SAMPLES,
    SAMPLES,
    compute_zenith,
    pick_near_bound,
    predict,
    satellite_position,
)
from validation_io import download_eumetsat, parse_run, write_results

# the products of each configuration: collection, the text of the names of a scan's
# products (from its start), and the file of each product
CLOUD_MASK = ("EO:EUM:DAT:0678", "OPE_{:%Y%m%d%H%M%S}_", ".nc")
LST = ("EO:EUM:DAT:1088", "LPMG_{:%Y%m%d%H%M%S}", ".nc")
SST = ("EO:EUM:DAT:0694", "{:%Y%m%d%H%M%S}-OSISAF", ".nc")
# the least quality level of the sea surface temperature of acceptable quality
SST_QUALITY = 4
# the bounds of local zenith angle (deg) of the configurations, measured by `--fit`
BOUNDS = {"lst_retrieved": 70.0, "sst_acceptable_quality": 65.0}
# the nominal position of the satellite, over the equator at 0 deg longitude
NOMINAL_HEIGHT = 35786400.0


def download(product: tuple, time: datetime, directory: pathlib.Path) -> pathlib.Path:
    """Downloads the file of a product of the scan starting at a time."""
    collection, match, entry = product
    (path,) = download_eumetsat(
        collection,
        time,
        time + timedelta(minutes=5),
        directory,
        match=match.format(time),
        entry=entry,
    )
    return path


def fixed_grid(dataset, step: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the geodetic latitude and longitude (deg) of every `step`-th pixel of the
    MTG fixed grid (NaN off the Earth's disk), whose angles sweep about the y axis.
    """
    name = next(
        v
        for v in dataset.variables
        if getattr(dataset[v], "grid_mapping_name", "") == "geostationary"
    )
    projection = dataset[name]
    geos = Proj(
        proj="geos",
        h=float(projection.perspective_point_height),
        lon_0=float(projection.longitude_of_projection_origin),
        sweep=str(projection.sweep_angle_axis),
        a=float(projection.semi_major_axis),
        b=float(projection.semi_minor_axis),
    )
    height = float(projection.perspective_point_height)
    x = dataset["x"][::step].astype(float) * height
    y = dataset["y"][::step].astype(float) * height
    longitude, latitude = geos(*np.meshgrid(x, y), inverse=True)
    off_disk = (np.abs(longitude) > 360) | (np.abs(latitude) > 90)
    return np.where(off_disk, np.nan, latitude), np.where(off_disk, np.nan, longitude)


def read_lst(path: pathlib.Path, step: int = 4):
    """Gets the positions of a land surface temperature scan, and which are retrieved."""
    with netCDF4.Dataset(path) as dataset:
        latitude, longitude = fixed_grid(dataset, step)
        values = np.ma.squeeze(dataset["LST"][..., ::step, ::step])
    return latitude, longitude, ~np.ma.getmaskarray(values) & np.isfinite(latitude)


def read_sst(path: pathlib.Path, step: int = 2):
    """
    Gets the positions of a sea surface temperature product's cells, which are
    retrieved, which are of acceptable quality or better, and their zenith angles.
    """
    with netCDF4.Dataset(path) as dataset:
        latitude, longitude = np.meshgrid(
            dataset["lat"][::step].astype(float),
            dataset["lon"][::step].astype(float),
            indexing="ij",
        )
        values = np.ma.squeeze(dataset["sea_surface_temperature"][..., ::step, ::step])
        quality = np.ma.squeeze(dataset["quality_level"][..., ::step, ::step])
        zenith = np.ma.squeeze(dataset["satellite_zenith_angle"][..., ::step, ::step])
    retrieved = ~np.ma.getmaskarray(values)
    good = retrieved & (np.ma.filled(quality, 0) >= SST_QUALITY)
    return latitude, longitude, retrieved, good, np.ma.filled(zenith, np.nan)


def nominal_position() -> np.ndarray:
    """Gets the Earth-fixed position (m) of the satellite's nominal location."""
    return np.array([6378137.0 + NOMINAL_HEIGHT, 0.0, 0.0])


def fit_bounds(run, time: datetime) -> None:
    """Prints the greatest local zenith angles of the products' retrievals."""
    latitude, longitude, retrieved = read_lst(download(LST, time, run.directory), 2)
    lst = compute_zenith(nominal_position(), latitude[retrieved], longitude[retrieved])
    *_, retrieved, good, zenith = read_sst(download(SST, time, run.directory), 1)
    print(
        json.dumps(
            {
                "lst_zenith_max": float(lst.max()),
                "lst_zenith_p999": float(np.percentile(lst, 99.9)),
                "sst_acceptable_zenith_max": float(np.nanmax(zenith[good])),
                "sst_retrieved_zenith_max": float(np.nanmax(zenith[retrieved])),
            }
        )
    )


def main() -> None:
    """Validates the entry on a full-disk scan, or fits its bounds."""
    run = parse_run("fci", __doc__)
    time = run.day
    if run.fit:
        fit_bounds(run, time)
        return
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
    with netCDF4.Dataset(download(CLOUD_MASK, time, run.directory)) as dataset:
        latitude, longitude = fixed_grid(dataset, 8)
        imaged = ~np.ma.getmaskarray(dataset["cloud_state"][::8, ::8])
    picks = rng.choice(
        np.flatnonzero(np.isfinite(latitude) & imaged), SAMPLES, replace=False
    )
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
            "reference": "FCI Level 2 cloud mask (EO:EUM:DAT:0678), every 8th pixel",
            "metrics": [
                {
                    "name": "imaged_pixels_observed",
                    "value": 100 * float(observed.mean()),
                },
                {"name": "elevation_error_max", "value": float(np.nanmax(error))},
            ],
            "fixture_metrics": {
                "imaged_pixels_observed": 100 * float(observed[keep].mean()),
                "elevation_error_max": float(np.nanmax(error[keep])),
            },
        }
        print(json.dumps(results[key]["metrics"]), flush=True)

    # product limits: the boundary of predicted observations near the bounds, among
    # pixels on the disk within 5 deg of them, and the share of the product's pixels
    # (retrieved, or of acceptable quality) that are observed
    latitude, longitude, retrieved = read_lst(download(LST, time, run.directory))
    on_disk = np.isfinite(latitude)
    products = {
        "lst_retrieved": (latitude, longitude, on_disk, retrieved, LST),
    }
    latitude, longitude, _, good, _ = read_sst(download(SST, time, run.directory))
    products["sst_acceptable_quality"] = (
        latitude,
        longitude,
        np.ones(latitude.shape, dtype=bool),
        good,
        SST,
    )
    references = {
        "lst_retrieved": "LSA SAF land surface temperature (EO:EUM:DAT:1088)",
        "sst_acceptable_quality": "OSI SAF Level 3C sea surface temperature "
        "(EO:EUM:DAT:0694), quality levels 4 and 5",
    }
    for key, (latitude, longitude, candidate, good, _) in products.items():
        bound = BOUNDS[key]
        configured = run.entry.configurations[key].instruments[0]
        expected = compute_field_of_regard(35786e3, 90 - bound)
        if abs(configured.parameters["field_of_regard"] - expected) > 1e-9:
            raise ValueError(
                f"{key}: the entry does not match the bound of {bound} deg"
            )
        zenith = np.full(latitude.shape, np.nan)
        zenith[candidate] = compute_zenith(
            position, latitude[candidate], longitude[candidate]
        )
        picks = rng.choice(
            np.flatnonzero(candidate & (np.abs(zenith - bound) <= 5)),
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
            "reference": references[key],
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
