"""
Validates the ECOSTRESS catalog entry against the ECOSTRESS Level 1B geolocation
product, using the exact values of the entry through public TAT-C methods. ECOSTRESS
is a thermal radiometer on the International Space Station whose mirror sweeps the
Earth across track every 1.181 s, imaging 128 lines at once; each sweep is taken as a
scan, from its middle line. The station's position is read from the attitude and
ephemeris product (inertial, rotated to Earth-fixed). The scan's angles are fitted
to one day (run with `--fit`), which needs no element set: running on a day at
least 30 days from it tests the fitted values on data held out from the fit, with a
public element set.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. The
geolocation files are large (about 700 MB each), so a day is reduced to blocks of
eight consecutive scenes in every 32 (spans of about seven minutes, over which the
swath is compared) and every fourth sweep, whose middle lines are read with range
requests (about 5 GB a day), and every 50th sample (and the last). The fit of
2026-10-04 read every fourth scene instead.

    python docs/validation/validate_catalog_ecostress.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_ecostress.py <data directory> 2026-10-04 --fit-view
    python docs/validation/validate_catalog_ecostress.py <data directory> 2026-11-08

ECOSTRESS images scenes over land and coastal waters, whose gaps are not compared.

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
from skyfield.framelib import itrs

from tatc.catalog import get_entry
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import (
    Product,
    Reference,
    Setup,
    find_granules,
    fit,
    fit_view,
    fit_yaw,
    validate_swath,
)
from validation_io import TO_GEODETIC, RemoteFile

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-07": (
        "1 25544U 98067A   26313.20747964  .00008201  00000+0  15534-3 0  9995",
        "2 25544  51.6321 302.4538 0007297   5.2479 354.8585 15.49270212594231",
    ),
    "2026-11-08": (
        "1 25544U 98067A   26314.17503267  .00009346  00000+0  17575-3 0  9996",
        "2 25544  51.6322 297.6644 0007292   8.9462 351.1656 15.49289712594382",
    ),
}
GEOLOCATION = Product(
    short_name="ECO_L1B_GEO", version="002", host="https://data.lpdaac"
)
EPHEMERIS = Product(short_name="ECO_L1B_ATT", version="002", host="https://data.lpdaac")
LINES_PER_SWEEP, SWEEP_PERIOD = 128, 1.181
# blocks of consecutive scenes, so that the swath is compared over spans of minutes
SCENE_BLOCK, BLOCK_STRIDE, SWEEP_STRIDE, SAMPLE_STRIDE = 8, 4, 4, 50
SAMPLES = 5400
COLUMNS = np.r_[0:SAMPLES:SAMPLE_STRIDE, SAMPLES - 1]
J2000 = 2451545.0


def read_scene(url: str, session, day: datetime) -> dict:
    """
    Reads the middle line of every fourth sweep of a scene with range requests: its
    time (seconds from the start of the day) and the positions of its samples.
    """
    with h5py.File(RemoteFile(url, session, block_size=2**20)) as f:
        g = f["Geolocation"]
        start = g["line_start_time_j2000"][:]
        first = np.flatnonzero(np.r_[True, np.diff(start) > 0.5])
        rows = (first + LINES_PER_SWEEP // 2)[::SWEEP_STRIDE]
        rows = rows[rows < len(start)]
        lat = g["latitude"][list(rows), :][:, COLUMNS]
        lon = g["longitude"][list(rows), :][:, COLUMNS]
    utc = timescale.tt_jd(J2000 + start[rows] / 86400).utc_datetime()
    seconds = np.array([(t - day).total_seconds() for t in utc])
    return {"time": seconds, "lat": lat, "lon": lon}


def read_ephemeris(url: str, session, day: datetime) -> tuple[np.ndarray, np.ndarray]:
    """Reads an orbit's station positions (m, Earth-fixed) and their times."""
    with h5py.File(RemoteFile(url, session)) as f:
        seconds = f["Ephemeris/time_j2000"][:]
        position = f["Ephemeris/eci_position"][:]
    t = timescale.tt_jd(J2000 + seconds / 86400)
    fixed = np.einsum("ijk,kj->ki", itrs.rotation_at(t), position)
    utc = t.utc_datetime()
    return np.array([(u - day).total_seconds() for u in utc]), fixed


def load_day(day: datetime, directory: pathlib.Path) -> dict:
    """Reads and reduces the scenes of a day, or loads them from a cache."""
    cache = directory / f"ecostress_{day:%Y%m%d}.npz"
    if not cache.exists():
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        scenes = [
            url
            for i, url in enumerate(find_granules(GEOLOCATION, day))
            if (i // SCENE_BLOCK) % BLOCK_STRIDE == 0
        ]
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            parts = list(pool.map(lambda url: read_scene(url, session, day), scenes))
            orbits = list(
                pool.map(
                    lambda url: read_ephemeris(url, session, day),
                    find_granules(EPHEMERIS, day),
                )
            )
        data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        seconds = np.concatenate([o[0] for o in orbits])
        position = np.concatenate([o[1] for o in orbits])
        order = np.argsort(seconds)
        seconds, position = seconds[order], position[order]
        keep = (data["time"] >= max(seconds[0], 0)) & (
            data["time"] <= min(seconds[-1], 86400)
        )
        data = {k: v[keep] for k, v in data.items()}
        xyz = np.stack(
            [np.interp(data["time"], seconds, position[:, i]) for i in range(3)],
            axis=-1,
        )
        sc_lon, sc_lat, sc_alt = TO_GEODETIC.transform(*xyz.T)
        order = np.argsort(data["time"])
        np.savez_compressed(
            cache,
            **{k: v[order] for k, v in data.items()},
            sc_lat=np.asarray(sc_lat)[order],
            sc_lon=np.asarray(sc_lon)[order],
            sc_alt=np.asarray(sc_alt)[order],
            scenes=len(scenes),
        )
    data = dict(np.load(cache))
    data.pop("scenes")
    return data


SETUP = Setup(
    entry_id="ecostress",
    platform="iss",
    swath="swath",
    beams="",
    size_instrument="ECOSTRESS",
    published_footprints=(),
    record_path=REPOSITORY
    / "src"
    / "tatc"
    / "resources"
    / "validation"
    / "ecostress.json",
    fixture_path=REPOSITORY / "tests" / "catalog" / "fixtures" / "ecostress.npz",
    source="docs/validation/validate_catalog_ecostress.py",
    reference="ECO_L1B_GEO version 002 (sample positions and times of every fourth "
    "sweep, in blocks of eight consecutive scenes of every 32)",
    published_reference="",
    fit_day=FIT_DAY,
    time_steps=(10, 30, 120),
    pass_radius=20,
    point_spacing=100e3,
)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument(
        "--fit-view",
        action="store_true",
        help="center the view in TAT-C's frame, with the entry's yaw and pitch",
    )
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    data = load_day(day, arguments.directory)
    data.update(pixel=COLUMNS, pixels=SAMPLES, line_stride=SWEEP_STRIDE)
    reference = Reference(day, data, SWEEP_STRIDE * SWEEP_PERIOD)
    print(
        f"{arguments.day}: {len(reference.time)} sweeps from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if arguments.fit_view:
        attitude = get_entry("ecostress").configurations["swath"].attitude
        fit_view(reference, "inertial", attitude["yaw_angle"], attitude["pitch_angle"])
        return
    if arguments.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate_swath(
        SETUP, day, tle, get_entry("ecostress"), orbit, reference, arguments.write
    )
    if arguments.write:
        print(f"wrote {SETUP.record_path} and {SETUP.fixture_path}")


if __name__ == "__main__":
    main()
