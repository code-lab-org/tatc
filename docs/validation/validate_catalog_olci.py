"""
Validates the OLCI catalog entry against the full-resolution Level 1 OLCI products of
Sentinel-3A (`S3A_OL_1_EFR`, distributed by NASA LAADS DAAC), using the exact values
of the entry through public TAT-C methods. OLCI images a fan of five cameras, tilted
away from the sun, as a pushbroom. The entry's angles are fitted to the tie points of
the products on one day (run with `--fit`): the products record no spacecraft
position, which is found for each image line as the point nearest the lines of sight
of its tie points (from their viewing zenith and azimuth angles), so the fit needs no
element set. Running on a day at least 30 days from it tests the fitted angles on data
held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the tie
points and times are read from each product's archive, with range requests (about
1.3 MB a product):

    python docs/validation/validate_catalog_olci.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_olci.py <data directory> 2026-11-04

A day is reduced to one image line in 23 (about 1 s), at the 77 tie points of each
line (every 64th of 4865 columns). The products cover the sunlit part of each orbit,
and only some of its products are distributed, whose gaps are left out.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import io
import pathlib
import zipfile
from datetime import datetime, timezone

import h5py
import numpy as np
import pandas as pd

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import (
    Product,
    Reference,
    Setup,
    fit,
    fit_yaw,
    load_day,
    validate_swath,
)
from validation_io import RemoteFile, TO_GEODETIC, ecef, normalize, up

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-04": (
        "1 41335U 16011A   26312.33341724  .00000144  00000+0  77522-4 0  9991",
        "2 41335  98.6306  17.6324 0001215  89.3094 270.8226 14.26734893558614",
    ),
}
EPOCH_2000 = pd.Timestamp("2000-01-01", tz="UTC")
GPS_UTC_OFFSET = 18  # s
LINE_STRIDE = 23
LINE_PERIOD = 0.044  # s, approximately, of an image line


def nearest_point(ground: np.ndarray, sight: np.ndarray) -> np.ndarray:
    """Gets the point nearest lines from ground points along unit directions."""
    projection = np.eye(3) - sight[:, :, None] * sight[:, None, :]
    return np.linalg.solve(
        projection.sum(axis=0), (projection @ ground[:, :, None]).sum(axis=0)
    )[:, 0]


def read_olci(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the tie points, times, and viewing angles of a product from its archive with
    range requests, reduces them to one line in 23, and finds the spacecraft's position
    for each line from the lines of sight of its tie points.
    """
    archive = zipfile.ZipFile(
        io.BufferedReader(RemoteFile(url, session, block_size=2**18), 2**18)
    )

    def read(name: str, fields: tuple[str, ...]) -> dict:
        member = next(n for n in archive.namelist() if n.endswith(name))
        with h5py.File(io.BytesIO(archive.read(member))) as f:
            values = {}
            for field in fields:
                raw = f[field][:]
                fill = f[field].attrs.get("_FillValue")
                scale = f[field].attrs.get("scale_factor", [1.0])[0]
                value = raw.astype(float) * scale
                if fill is not None:
                    value[raw == fill[0]] = np.nan
                values[field] = value
            return values

    geo = read("tie_geo_coordinates.nc", ("latitude", "longitude"))
    angles = read("tie_geometries.nc", ("OZA", "OAA"))
    stamps = read("time_coordinates.nc", ("time_stamp",))["time_stamp"]
    keep = np.arange(0, len(stamps), LINE_STRIDE)
    lat, lon = geo["latitude"][keep], geo["longitude"][keep]
    zenith, azimuth = np.radians(angles["OZA"][keep]), np.radians(angles["OAA"][keep])
    # the direction from each tie point to the spacecraft, from east, north, and up
    lam, phi = np.radians(lon), np.radians(lat)
    east = np.stack([-np.sin(lam), np.cos(lam), np.zeros_like(lam)], axis=-1)
    north = np.stack(
        [-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)], axis=-1
    )
    sight = (
        (np.sin(zenith) * np.sin(azimuth))[..., None] * east
        + (np.sin(zenith) * np.cos(azimuth))[..., None] * north
        + np.cos(zenith)[..., None] * up(lon, lat)
    )
    ground = ecef(np.nan_to_num(lon), np.nan_to_num(lat))
    position = np.full((len(keep), 3), np.nan)
    for row in range(len(keep)):
        valid = (
            np.isfinite(lat[row]) & np.isfinite(zenith[row]) & np.isfinite(azimuth[row])
        )
        if valid.sum() >= 10:
            position[row] = nearest_point(ground[row, valid], sight[row, valid])
    sc_lon, sc_lat, sc_alt = TO_GEODETIC.transform(*position.T)
    # the times are in GPS time, 18 s ahead of UTC since 2017
    seconds = (
        EPOCH_2000 + pd.to_timedelta(stamps[keep], unit="us") - pd.Timestamp(day)
    ).total_seconds() - GPS_UTC_OFFSET
    return {
        "time": np.asarray(seconds, dtype=float),
        "lat": np.nan_to_num(lat, nan=-999.0),
        "lon": np.nan_to_num(lon, nan=-999.0),
        "sc_lat": np.nan_to_num(np.asarray(sc_lat), nan=-999.0),
        "sc_lon": np.nan_to_num(np.asarray(sc_lon), nan=-999.0),
        "sc_alt": np.nan_to_num(np.asarray(sc_alt), nan=0.0),
        "quality": (~np.all(np.isfinite(position), axis=1)).astype(np.uint8),
    }


PRODUCT = Product(
    short_name="S3A_OL_1_EFR",
    version="1",
    read=read_olci,
    host="https://data.laadsdaac",
)
SETUP = Setup(
    entry_id="olci",
    platform="sentinel-3a",
    swath="swath",
    beams="",
    size_instrument="OLCI",
    published_footprints=(),
    record_path=REPOSITORY / "src" / "tatc" / "resources" / "validation" / "olci.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "olci.npz",
    source="docs/validation/validate_catalog_olci.py",
    reference="S3A_OL_1_EFR (tie point positions and times)",
    published_reference="",
    fit_day=FIT_DAY,
    time_steps=(10, 30, 120),
    pass_radius=40,
    point_spacing=200e3,
)


def load_reference(day: datetime, directory: pathlib.Path) -> Reference:
    """Loads a reduced day."""
    return Reference(day, load_day(PRODUCT, day, directory), LINE_STRIDE * LINE_PERIOD)


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
    reference = load_reference(day, arguments.directory)
    print(
        f"{arguments.day}: {len(reference.time)} lines from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute, median altitude "
        f"{np.median(reference.altitude) / 1e3:.1f} km, line period "
        f"{np.median(np.diff(reference.seconds)) / LINE_STRIDE:.5f} s",
        flush=True,
    )
    if arguments.fit:
        for frame in ("earth_fixed", "inertial"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        # the along-track angle of each tie point, which a single view holds constant
        lon, lat, _ = TO_GEODETIC.transform(*reference.position.T)
        nadir = -up(lon, lat)
        motion = reference.velocity("earth_fixed")
        forward = normalize(motion - np.sum(motion * nadir, axis=-1)[:, None] * nadir)
        sight = normalize(
            ecef(reference.lon, reference.lat) - reference.position[:, np.newaxis]
        )
        along = np.degrees(
            np.arctan2(
                np.sum(sight * forward[:, np.newaxis], axis=-1),
                np.sum(sight * nadir[:, np.newaxis], axis=-1),
            )
        )
        print(
            "along-track angle (deg) by tie point, every 4th: "
            + " ".join(f"{a:+.3f}" for a in np.median(along, axis=0)[::4])
        )
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("olci")
    validate_swath(SETUP, day, tle, entry, orbit, reference, arguments.write)
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
