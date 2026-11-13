"""
Validates the SPEXone catalog entry against the PACE SPEXone Level 1C product, using
the exact values of the entry through public TAT-C methods. SPEXone is a multi-angle
spectropolarimeter that views each place at five angles along track (-50, -20, 0,
20, and 50 deg); the entry models its nadir view, from the positions of the
product's bins (about 5 km) that it observes (6 to 24 of the grid's 29) and the
times at which it observed them.
The product gives no spacecraft position, which is read from the navigation of the
Level 1B product (Earth-fixed, every 4 s). The angles are fitted to one day (run
with `--fit`), which needs no element set: running on a day at least 30 days from it
tests the fitted values on data held out from the fit.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields are read, with range requests:

    python docs/validation/validate_catalog_spexone.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_spexone.py <data directory> 2026-11-07

The product covers the sunlit part of each orbit, whose dark part is a gap.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import pathlib
from datetime import datetime, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd

from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import (
    Product,
    Reference,
    Setup,
    find_granules,
    fit,
    fit_yaw,
    load_day,
    validate_swath,
)
from validation_io import TO_GEODETIC, RemoteFile

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
PACE_ELEMENTS = (
    "1 58928U 24025A   26313.27483190  .00000501  00000+0  10091-3 0  9993",
    "2 58928  98.1040 242.4857 0002666 129.8006 230.3433 14.65644467147180",
)
ELEMENTS = {"2026-11-07": PACE_ELEMENTS, "2026-11-08": PACE_ELEMENTS}
NADIR_VIEW = 2
# the bins of the product's grid that the nadir view observes
COLUMNS = np.arange(6, 25)


def read_spexone(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads a Level 1C granule with range requests: each row's nadir-view time and each
    bin's (seconds from the start of the day), and the bins' positions.
    """
    with h5py.File(RemoteFile(url, session)) as f:
        start = pd.Timestamp(f.attrs["time_coverage_start"].decode())
        seconds = f["bin_attributes/nadir_view_time"][:].astype(float)
        offsets = f["bin_attributes/view_time_offsets"][:, :, NADIR_VIEW]
        offsets = offsets[:, COLUMNS].astype(float)
        lat = f["geolocation_data/latitude"][:][:, COLUMNS]
        lon = f["geolocation_data/longitude"][:][:, COLUMNS]
    # the times are seconds of the day on which the granule starts
    seconds = seconds + (start.normalize() - pd.Timestamp(day)).total_seconds()
    keep = (seconds > -1e4) & np.all(np.abs(offsets) < 200, axis=1)
    return {
        "time": seconds[keep],
        "beam_time": seconds[keep, np.newaxis] + offsets[keep],
        "lat": lat[keep],
        "lon": lon[keep],
    }


def read_navigation(url: str, session, day: datetime) -> tuple[np.ndarray, np.ndarray]:
    """Reads a Level 1B granule's spacecraft positions (m, Earth-fixed) and times."""
    with h5py.File(RemoteFile(url, session)) as f:
        start = pd.Timestamp(f.attrs["time_coverage_start"].decode()[:10], tz="UTC")
        seconds = f["navigation_data/orb_time"][:].astype(float)
        position = f["navigation_data/orb_pos"][:].astype(float)
    seconds = seconds + (start - pd.Timestamp(day)).total_seconds()
    keep = (np.abs(position) < 7.2e6).all(axis=1) & (seconds > -1e4)
    return seconds[keep], position[keep]


PRODUCT = Product(
    short_name="PACE_SPEXONE_L1C_SCI",
    version="4",
    read=read_spexone,
    host="https://obdaac-tea",
)
NAVIGATION = Product(
    short_name="PACE_SPEXONE_L1B_SCI", version="4", host="https://obdaac-tea"
)
SETUP = Setup(
    entry_id="spexone",
    platform="pace",
    swath="swath",
    beams="",
    size_instrument="SPEXone",
    published_footprints=(),
    record_path=REPOSITORY
    / "src"
    / "tatc"
    / "resources"
    / "validation"
    / "spexone.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "spexone.npz",
    source="docs/validation/validate_catalog_spexone.py",
    reference="PACE_SPEXONE_L1C_SCI version 4 (positions of bins and times of the "
    "nadir view)",
    published_reference="",
    fit_day=FIT_DAY,
    time_steps=(10, 30, 120),
    pass_radius=10,
    point_spacing=100e3,
)


def load_navigation(day: datetime, directory: pathlib.Path):
    """Reads the spacecraft positions of a day, or loads them from a cache."""
    cache = directory / f"spexone_navigation_{day:%Y%m%d}.npz"
    if not cache.exists():
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        with concurrent.futures.ThreadPoolExecutor(4) as pool:
            parts = list(
                pool.map(
                    lambda url: read_navigation(url, session, day),
                    find_granules(NAVIGATION, day),
                )
            )
        seconds = np.concatenate([p[0] for p in parts])
        position = np.concatenate([p[1] for p in parts])
        _, first = np.unique(seconds, return_index=True)
        np.savez_compressed(cache, seconds=seconds[first], position=position[first])
    data = np.load(cache)
    return data["seconds"], data["position"]


def load_reference(day: datetime, directory: pathlib.Path) -> Reference:
    """Loads a day of bins, with the spacecraft's position at each row's time."""
    data = load_day(PRODUCT, day, directory)
    seconds, position = load_navigation(day, directory)
    inside = (data["time"] >= seconds[0]) & (data["time"] <= seconds[-1])
    data = {k: v[inside] for k, v in data.items() if k != "granules"}
    xyz = np.stack(
        [np.interp(data["time"], seconds, position[:, i]) for i in range(3)], axis=-1
    )
    sc_lon, sc_lat, sc_alt = TO_GEODETIC.transform(*xyz.T)
    data.update(sc_lat=np.asarray(sc_lat), sc_lon=np.asarray(sc_lon))
    data["sc_alt"] = np.asarray(sc_alt)
    return Reference(day, data, float(np.median(np.diff(data["time"]))))


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
        f"{arguments.day}: {len(reference.time)} rows from {reference.time[0]}, "
        f"{reference.beams} bins, {len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit:
        for frame in ("earth_fixed", "inertial"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        offsets = reference.beam_seconds - reference.seconds[:, np.newaxis]
        print(f"bin times from the row's: {np.round(np.median(offsets, axis=0), 3)}")
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate_swath(
        SETUP, day, tle, get_entry("spexone"), orbit, reference, arguments.write
    )
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
