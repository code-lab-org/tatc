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

    python validation/catalog/validate_catalog_ecostress.py <data directory> 2026-10-04 --fit
    python validation/catalog/validate_catalog_ecostress.py <data directory> 2026-10-04 --fit --fit-view
    python validation/catalog/validate_catalog_ecostress.py <data directory> 2026-11-08

It also validates the entry's sample size configuration (`beams`), following
`ValidateSampleSizes`: the projected pixels of the middle line of a sweep are
compared with the samples' positions, and their spacing with that measured. The only
published pixel, at nadir in the L1 Calibration ATBD, is compared without the gate.

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
With `--fit-view`, it centers the view in TAT-C's frame, with the entry's yaw and
pitch, instead of fitting the scan's angles. ECOSTRESS images scenes over land and
coastal waters, whose gaps are not compared.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import pathlib
from datetime import datetime

import earthaccess
import h5py
import numpy as np
from skyfield.framelib import itrs

from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit
from scanner import (
    Product,
    Reference,
    Setup,
    find_granules,
    fit,
    fit_view,
    fit_yaw,
    validate_beams,
    validate_swath,
)
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    FIXTURES,
    login_earthdata,
    parse_run,
)

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
# the published pixel at nadir: instrument, altitude (km), inclination (deg), scan
# angle (deg), sizes along and across track (km), and angle along track (deg), from
# Table 1 of the L1 Calibration ATBD, compared without the gate, as no other is
# published; at the center of the scan line, which the body's pitch places behind
# nadir
PUBLISHED_FOOTPRINTS = (("ECOSTRESS", 400, 51.6, 0, 0.0753, 0.06851, -7.2427),)


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
        login_earthdata()
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
    swath="swath",
    beams="beams",
    size_instrument="ECOSTRESS",
    published_footprints=PUBLISHED_FOOTPRINTS,
    fixture_path=FIXTURES / "ecostress.npz",
    reference="ECO_L1B_GEO version 002 (sample positions and times of every fourth "
    "sweep, in blocks of eight consecutive scenes of every 32)",
    published_reference="Nadir pixel of the ECOSTRESS L1 Calibration ATBD (not "
    "gated); sample positions of ECO_L1B_GEO version 002 (the middle line of every "
    "fourth sweep)",
    time_steps=(10, 30, 120),
    pass_radius=20,
    point_spacing=100e3,
    gate_published=False,
)


def main() -> None:
    """Validates the entry on a day of data, or fits its angles."""
    run = parse_run(
        "ecostress",
        __doc__,
        add_arguments=lambda parser: parser.add_argument(
            "--fit-view",
            action="store_true",
            help="center the view in TAT-C's frame, with the entry's yaw and pitch",
        ),
    )
    data = load_day(run.day, run.directory)
    # the middle line of each kept sweep, of the lines a sweep images at once
    data.update(
        pixel=COLUMNS,
        pixels=SAMPLES,
        line_stride=SWEEP_STRIDE * LINES_PER_SWEEP,
        line=LINES_PER_SWEEP // 2,
    )
    reference = Reference(run.day, data, SWEEP_STRIDE * SWEEP_PERIOD)
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} sweeps from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute",
        flush=True,
    )
    if run.fit_view:
        attitude = run.entry.configurations[SETUP.swath].attitude
        fit_view(reference, "inertial", attitude["yaw_angle"], attitude["pitch_angle"])
        return
    if run.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        return
    tle = run.element_set()
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    if SETUP.beams in run.configurations:
        validate_beams(SETUP, run, reference, orbit, tle)
    if SETUP.swath in run.configurations:
        validate_swath(SETUP, run, reference, orbit, tle)


if __name__ == "__main__":
    main()
