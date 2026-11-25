"""
Validates the CERES catalog entry against the NOAA-20 CERES FM6 ES-8 product (Level 2
instantaneous footprints), using the exact values of the entry through public TAT-C
methods. In its cross-track mode, CERES sweeps the Earth twice in each 6.6 s record,
in opposite directions; the first sweep of each record is taken as a scan, with each
footprint at its own time (10 ms apart). The product places footprints at the top of the atmosphere
(30 km above the ellipsoid), so each is projected to the ellipsoid along the line of
sight from the spacecraft's position, which the product records. The scan's angles
are fitted to one day (run with `--fit`), which needs no element set: running on a
day at least 30 days from it tests the fitted values on data held out from the fit,
with NOAA's operational element set of the day.

Needs the `validation` extra, `pyhdf`, and an Earthdata token in `EARTHDATA_TOKEN`.
Each day's file (about 280 MB, HDF4) is downloaded whole:

    python validation/catalog/validate_catalog_ceres.py <data directory> 2026-08-25 --fit
    python validation/catalog/validate_catalog_ceres.py <data directory> 2026-09-30

The product is released about six weeks after the data are taken. A day is reduced
to the first sweep of each record.

With `--platform suomi-npp`, it validates the configuration of Suomi NPP's CERES
(FM5) with the CER_ES8_NPP-FM5 Edition2 product, which ends on 2026-07-31, before
the satellite was retired. NOAA's archive has Suomi NPP's element sets only from
2026-08-19, so earlier days read them from a file of the history of its element sets:

    python validation/catalog/validate_catalog_ceres.py <data directory> 2026-07-25 --platform suomi-npp --fit
    python validation/catalog/validate_catalog_ceres.py <data directory> 2026-06-20 --platform suomi-npp

It also validates the platform's sample size configuration (`beams_noaa_20` or
`beams_suomi_npp`), following `ValidateSampleSizes`: the projected pixels are
compared with the footprints' positions, and their spacing with that measured. The
only published footprints, of the optical field of view at nadir in the CERES ATBD
(for TRMM and EOS), are compared without the gate.

Without a day, it takes that of the platform's record (or, with `--fit`, of the
fit).

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
from dataclasses import dataclass
from datetime import datetime

import earthaccess
import numpy as np
import pyhdf.VS  # pylint: disable=unused-import
from pyhdf.HDF import HDF
from pyhdf.SD import SD, SDC

from tatc.schemas import GeneralPerturbationsOrbit
from scanner import Reference, Setup, fit, fit_yaw, validate_beams, validate_swath
from validation_io import (
    FIXTURES,
    TO_GEODETIC,
    get_with_retries,
    login_earthdata,
    parse_run,
)


@dataclass(frozen=True)
class Product:
    """
    The ES-8 product of a platform's CERES, the configuration it validates, its
    fixture, and the prefix of its cached files.
    """

    short_name: str
    version: str
    configuration: str
    beams: str
    fixture_path: pathlib.Path
    prefix: str


PRODUCTS = {
    "noaa-20": Product(
        "CER_ES8_NOAA20-FM6",
        "Edition1",
        "swath_noaa_20",
        "beams_noaa_20",
        FIXTURES / "ceres.npz",
        "",
    ),
    "suomi-npp": Product(
        "CER_ES8_NPP-FM5",
        "Edition2",
        "swath_suomi_npp",
        "beams_suomi_npp",
        FIXTURES / "ceres_suomi_npp.npz",
        "npp_",
    ),
}
RECORD_SAMPLES, SAMPLE_PERIOD, RECORD_PERIOD = 660, 0.01, 6.6
# the samples of the first Earth sweep of each record; the second sweeps in the
# opposite direction, so that its footprints and the first's, taken as consecutive
# scans, would cross and leave the cells between them twisted
SWEEPS = (np.arange(72, 263),)
SAMPLE_STRIDE = 1
TOA_HEIGHT = 30e3  # m
WGS84_A, WGS84_F = 6378137.0, 1 / 298.257223563
JULIAN_UNIX_EPOCH = 2440587.5
# published footprints of the optical field of view at nadir: instrument, altitude
# (km), inclination (deg), scan angle (deg), and sizes along and across track (km),
# from Table 4.4-1 of the CERES ATBD Subsystem 4.4 (for TRMM and EOS, whose sizes
# along the scan are across track in the cross-track mode); none are published at
# other angles or for this orbit, so they are compared without the gate
PUBLISHED_FOOTPRINTS = (
    ("CERES", 350, 35, 0, 16, 8),
    ("CERES", 705, 98.2, 0, 32, 16),
)


def download(day: datetime, directory: pathlib.Path, product: Product) -> pathlib.Path:
    """Downloads the product's file of a day, if it is not already."""
    items = get_with_retries(
        "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
        params={
            "short_name": product.short_name,
            "version": product.version,
            "temporal": f"{day:%Y-%m-%dT00:00:00Z},{day:%Y-%m-%dT23:59:59Z}",
            "page_size": 10,
        },
    ).json()["items"]
    url = next(
        u["URL"]
        for item in items
        for u in item["umm"]["RelatedUrls"]
        if u["URL"].startswith("https://data.asdc")
        and u["URL"].endswith(f"{day:%Y%m%d}")
    )
    path = directory / url.split("/")[-1]
    if not path.exists():
        login_earthdata()
        session = earthaccess.get_requests_https_session()
        with session.get(url, stream=True, timeout=600) as response:
            response.raise_for_status()
            with open(path, "wb") as file:
                for block in response.iter_content(2**22):
                    file.write(block)
    return path


def project_to_ellipsoid(satellite: np.ndarray, colatitude, longitude) -> np.ndarray:
    """
    Gets the Earth-fixed positions (m) at which the lines of sight from the satellite
    through points at the top of the atmosphere (geocentric colatitude and longitude,
    in degrees) meet the ellipsoid.
    """
    b = WGS84_A * (1 - WGS84_F)
    latitude, lon = np.radians(90 - colatitude), np.radians(longitude)
    unit = np.stack(
        [
            np.cos(latitude) * np.cos(lon),
            np.cos(latitude) * np.sin(lon),
            np.sin(latitude),
        ],
        axis=-1,
    )
    radius = WGS84_A * b / np.hypot(b * np.cos(latitude), WGS84_A * np.sin(latitude))
    top = unit * (radius + TOA_HEIGHT)[..., np.newaxis]
    direction = top - satellite
    scale = np.array([1 / WGS84_A, 1 / WGS84_A, 1 / b])
    p, d = satellite * scale, direction * scale
    qa = np.sum(d * d, axis=-1)
    qb = 2 * np.sum(p * d, axis=-1)
    qc = np.sum(p * p, axis=-1) - 1
    root = (-qb - np.sqrt(np.maximum(qb**2 - 4 * qa * qc, 0))) / (2 * qa)
    return satellite + root[..., np.newaxis] * direction


def load_day(day: datetime, directory: pathlib.Path, product: Product) -> dict:
    """Reads and reduces the sweeps of a day, or loads them from a cache."""
    cache = directory / f"{product.prefix}ceres_es8_{day:%Y%m%d}.npz"
    if not cache.exists():
        path = download(day, directory, product)
        vs = HDF(str(path)).vstart()

        def table(name: str) -> np.ndarray:
            vdata = vs.attach(name)
            values = np.array(vdata[:], dtype=float).ravel()
            vdata.detach()
            return values

        julian = table("Time of observation")
        start, end = (
            np.stack(
                [
                    table(f"{c} component of satellite position at record {edge}")
                    for c in "XYZ"
                ],
                axis=-1,
            )
            for edge in ("start", "end")
        )
        sd = SD(str(path), SDC.READ)
        colatitude = sd.select("Colatitude of CERES FOV at TOA")[:]
        longitude = sd.select("Longitude of CERES FOV at TOA")[:]
        seconds = (julian - JULIAN_UNIX_EPOCH) * 86400 - day.timestamp()
        parts = []
        for sweep in SWEEPS:
            columns = sweep[
                np.unique(np.r_[0 : len(sweep) : SAMPLE_STRIDE, len(sweep) - 1])
            ]
            fraction = columns * SAMPLE_PERIOD / RECORD_PERIOD
            satellite = (
                start[:, np.newaxis]
                + fraction[:, np.newaxis] * (end - start)[:, np.newaxis]
            )
            surface = project_to_ellipsoid(
                satellite, colatitude[:, columns], longitude[:, columns]
            )
            lon, lat, _ = TO_GEODETIC.transform(*np.moveaxis(surface, -1, 0))
            middle = columns[len(columns) // 2] * SAMPLE_PERIOD / RECORD_PERIOD
            sc_lon, sc_lat, sc_alt = TO_GEODETIC.transform(
                *(start + middle * (end - start)).T
            )
            valid = np.all(colatitude[:, columns] <= 180, axis=1)
            parts.append(
                {
                    "time": seconds + columns[len(columns) // 2] * SAMPLE_PERIOD,
                    "beam_time": seconds[:, np.newaxis] + columns * SAMPLE_PERIOD,
                    "lat": np.where(valid[:, np.newaxis], lat, np.nan),
                    "lon": np.where(valid[:, np.newaxis], lon, np.nan),
                    "sc_lat": np.asarray(sc_lat),
                    "sc_lon": np.asarray(sc_lon),
                    "sc_alt": np.asarray(sc_alt),
                }
            )
        data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        order = np.argsort(data["time"])
        np.savez_compressed(cache, **{k: v[order] for k, v in data.items()})
    return dict(np.load(cache))


def make_setup(platform: str) -> Setup:
    """Gets what a run on a platform validates, and where it writes its results."""
    product = PRODUCTS[platform]
    return Setup(
        swath=product.configuration,
        beams=product.beams,
        size_instrument="CERES",
        published_footprints=PUBLISHED_FOOTPRINTS,
        fixture_path=product.fixture_path,
        reference=f"{product.short_name} {product.version} (footprint positions, "
        "projected to the ellipsoid, and times)",
        published_reference="Nadir footprints of the CERES ATBD Subsystem 4.4 "
        f"(not gated); footprint positions of {product.short_name} "
        f"{product.version}, projected to the ellipsoid",
        time_steps=(RECORD_PERIOD, 30, 120),
        pass_radius=60,
        point_spacing=200e3,
        gate_published=False,
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run("ceres", __doc__)
    data = load_day(run.day, run.directory, PRODUCTS[run.platform])
    stride = np.unique(np.r_[0 : len(SWEEPS[0]) : SAMPLE_STRIDE, len(SWEEPS[0]) - 1])
    data.update(pixel=stride, pixels=len(SWEEPS[0]))
    reference = Reference(run.day, data, RECORD_PERIOD)
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} sweeps from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        offsets = reference.beam_seconds - reference.seconds[:, np.newaxis]
        print(f"footprint times from the sweep's middle: {offsets[0, [0, -1]]}")
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    # a transfer tests the configuration of the default platform
    entry = run.entry
    default = entry.configurations[entry.default_configuration].platforms[0]
    setup = make_setup(default if run.transfer else run.platform)
    if setup.beams in run.configurations:
        validate_beams(setup, run, reference, orbit, tle)
    if setup.swath in run.configurations:
        validate_swath(setup, run, reference, orbit, tle)


if __name__ == "__main__":
    main()
