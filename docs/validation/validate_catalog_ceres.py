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

    python docs/validation/validate_catalog_ceres.py <data directory> 2026-08-25 --fit
    python docs/validation/validate_catalog_ceres.py <data directory> 2026-09-30

The product is released about six weeks after the data are taken. A day is reduced
to the first sweep of each record.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import pathlib
from datetime import datetime, timezone

import earthaccess
import numpy as np
import pyhdf.VS  # pylint: disable=unused-import
from pyhdf.HDF import HDF
from pyhdf.SD import SD, SDC

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import Reference, Setup, fit, fit_yaw, validate_swath
from validation_io import (
    TO_GEODETIC,
    get_with_retries,
    load_element_set,
    use_archive,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 8, 25, tzinfo=timezone.utc)
SHORT_NAME, VERSION = "CER_ES8_NOAA20-FM6", "Edition1"
RECORD_SAMPLES, SAMPLE_PERIOD, RECORD_PERIOD = 660, 0.01, 6.6
# the samples of the first Earth sweep of each record; the second sweeps in the
# opposite direction, so that its footprints and the first's, taken as consecutive
# scans, would cross and leave the cells between them twisted
SWEEPS = (np.arange(72, 263),)
SAMPLE_STRIDE = 1
TOA_HEIGHT = 30e3  # m
WGS84_A, WGS84_F = 6378137.0, 1 / 298.257223563
JULIAN_UNIX_EPOCH = 2440587.5


def download(day: datetime, directory: pathlib.Path) -> pathlib.Path:
    """Downloads the product's file of a day, if it is not already."""
    items = get_with_retries(
        "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
        params={
            "short_name": SHORT_NAME,
            "version": VERSION,
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
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
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


def load_day(day: datetime, directory: pathlib.Path) -> dict:
    """Reads and reduces the sweeps of a day, or loads them from a cache."""
    cache = directory / f"ceres_es8_{day:%Y%m%d}.npz"
    if not cache.exists():
        path = download(day, directory)
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


SETUP = Setup(
    entry_id="ceres",
    platform="noaa-20",
    swath="swath",
    beams="",
    size_instrument="CERES",
    published_footprints=(),
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "ceres.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "ceres.npz",
    source="docs/validation/validate_catalog_ceres.py",
    reference="CER_ES8_NOAA20-FM6 Edition1 (footprint positions, projected to the "
    "ellipsoid, and times)",
    published_reference="",
    fit_day=FIT_DAY,
    time_steps=(RECORD_PERIOD, 30, 120),
    pass_radius=60,
    point_spacing=200e3,
)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    data = load_day(day, arguments.directory)
    stride = np.unique(np.r_[0 : len(SWEEPS[0]) : SAMPLE_STRIDE, len(SWEEPS[0]) - 1])
    data.update(pixel=stride, pixels=len(SWEEPS[0]))
    reference = Reference(day, data, RECORD_PERIOD)
    print(
        f"{arguments.day}: {len(reference.time)} sweeps from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        offsets = reference.beam_seconds - reference.seconds[:, np.newaxis]
        print(f"footprint times from the sweep's middle: {offsets[0, [0, -1]]}")
        return
    use_archive("noaa-20")
    tle = load_element_set(day, arguments.directory)
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate_swath(
        SETUP, day, tle, get_entry("ceres"), orbit, reference, arguments.write
    )
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
