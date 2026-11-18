"""
Validates the GEDI catalog entry against the GEDI Level 1B geolocated waveform product
(GEDI01_B, version 003), using the exact values of the entry through public TAT-C
methods, as `validate_catalog_atlas.py` does for ATLAS. GEDI is a lidar on the
International Space Station whose four lasers are dithered across track into eight
beam tracks about 600 m apart. GEDI points its beams off nadir by plan; most of the
time it points within about a degree of nadir, where the entry models them, and it
slews several degrees on some passes, which are left out: a shot is compared where
the center of the beams (their median off-nadir angle across track, in bins of 10 s
of each file) lies within 1 deg of nadir.

The entry models each beam as a narrow view. Their centers and widths are fitted to
the cross-track positions of the beams' shots relative to TAT-C's ground track on one
day (run with `--fit`): each view is centered on the median offset of its beam, and
is as wide as twice the 99.5th percentile of the offsets' distance from it, and the
body is pitched to place them at the shots' median offset along track. The widths
include the pointing's spread and the error of the element set (the one nearest the
middle of the day), so they hold for element sets of a similar age. Running on a day
at least 30 days from the fit tests the fitted values on data held out from the fit.

Needs the `validation` extra, an Earthdata token in `EARTHDATA_TOKEN`, and the ISS's
element sets of the days in `elements/iss_2026-04-01.3le` (or in a file given with
`--elements`). Every 200th shot (about 12 km apart) is read, with range requests:

    python validation/catalog/validate_catalog_gedi.py <data directory> 2026-05-01 --fit
    python validation/catalog/validate_catalog_gedi.py <data directory> 2026-06-08

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import pathlib
from datetime import datetime, timedelta, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd

from tatc.schemas import GeneralPerturbationsOrbit
from validate_catalog_atlas import add_offsets, score_views, track_bearings, write_views
from validation_io import (
    RemoteFile,
    ecef,
    get_with_retries,
    login_earthdata,
    parse_run,
    up,
)

CONFIGURATION = "beams"
SHORT_NAME, VERSION = "GEDI01_B", "003"
# the beams, as the product names them, from left to right of the ground track (as
# measured on the fitting day)
BEAMS = (
    "BEAM0000",
    "BEAM0001",
    "BEAM0010",
    "BEAM0011",
    "BEAM0101",
    "BEAM0110",
    "BEAM1000",
    "BEAM1011",
)
GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
GPS_UTC = 18  # s
SHOT_STRIDE = 200
NADIR_LIMIT = 1.0  # deg, of the beams' center
CENTER_BIN = 10  # s
SAMPLE_PER_BEAM = 250
FIT_PERCENTILE = 99.5


def read_granule(url: str, session) -> list[pd.DataFrame]:
    """
    Reads every 200th shot of each beam of a granule: its position (at the last bin
    of the waveform), time (GPS seconds), and off-nadir angle across track.
    """
    frames = []
    with h5py.File(RemoteFile(url, session, block_size=2**22)) as f:
        for beam in BEAMS:
            if f"{beam}/geolocation" not in f:
                continue
            g = f[f"{beam}/geolocation"]
            epoch = float(f[f"{beam}/ancillary/master_time_epoch"][0])
            data = {
                "latitude": g["latitude_lastbin"][::SHOT_STRIDE],
                "longitude": g["longitude_lastbin"][::SHOT_STRIDE],
                "gps_seconds": epoch + g["delta_time"][::SHOT_STRIDE],
                "inst_lat": g["latitude_instrument"][::SHOT_STRIDE],
                "inst_lon": g["longitude_instrument"][::SHOT_STRIDE],
                "inst_alt": g["altitude_instrument"][::SHOT_STRIDE],
            }
            frame = pd.DataFrame(data)
            frame = frame[
                (frame.latitude.abs() <= 90)
                & (frame.inst_lat.abs() <= 90)
                & np.isfinite(frame.inst_alt)
            ]
            if len(frame) < 3:
                continue
            frames.append(
                frame.assign(
                    file=url.rsplit("/", 1)[-1],
                    beam=beam,
                    off_nadir=off_nadir(frame),
                )
            )
    return frames


def off_nadir(frame: pd.DataFrame) -> np.ndarray:
    """Gets the shots' off-nadir angles across track (deg, positive to the left)."""
    instrument = ecef(
        frame.inst_lon.to_numpy(), frame.inst_lat.to_numpy(), frame.inst_alt.to_numpy()
    )
    shot = ecef(frame.longitude.to_numpy(), frame.latitude.to_numpy())
    nadir = -up(frame.inst_lon.to_numpy(), frame.inst_lat.to_numpy())
    sight = shot - instrument
    motion = np.gradient(instrument, axis=0)
    left = np.cross(motion, nadir)
    left /= np.linalg.norm(left, axis=1)[:, None]
    return np.degrees(
        np.arctan2(np.sum(sight * left, axis=1), np.sum(sight * nadir, axis=1))
    )


def load_day(day: datetime, directory: pathlib.Path) -> pd.DataFrame:
    """
    Reads every 200th shot of each beam from the day's files, with range requests, or
    loads them from a cache, and keeps those of the day where the beams' center lies
    within 1 deg of nadir.
    """
    cache = directory / f"gedi01b_{day:%Y%m%d}.csv"
    if not cache.exists():
        items = get_with_retries(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params={
                "short_name": SHORT_NAME,
                "version": VERSION,
                "temporal": f"{day - timedelta(hours=2):%Y-%m-%dT%H:%M:%SZ},"
                f"{day + timedelta(days=1):%Y-%m-%dT%H:%M:%SZ}",
                "page_size": 200,
            },
        ).json()["items"]
        urls = [
            next(
                u["URL"]
                for u in item["umm"]["RelatedUrls"]
                if u["URL"].startswith("https://data.lpdaac")
                and u["URL"].endswith(".h5")
            )
            for item in items
        ]
        login_earthdata()
        session = earthaccess.get_requests_https_session()
        with concurrent.futures.ThreadPoolExecutor(8) as pool:
            parts = [
                frame
                for frames in pool.map(lambda u: read_granule(u, session), urls)
                for frame in frames
            ]
        pd.concat(parts, ignore_index=True).to_csv(cache, index=False)
    shots = pd.read_csv(cache)
    shots["time"] = pd.to_datetime(GPS_EPOCH) + pd.to_timedelta(
        shots.gps_seconds - GPS_UTC, unit="s"
    )
    shots = shots[
        (shots.time >= pd.Timestamp(day))
        & (shots.time < pd.Timestamp(day) + pd.Timedelta(days=1))
    ]
    shots = shots.assign(bin=(shots.gps_seconds // CENTER_BIN).astype("int64"))
    center = shots.groupby(["file", "bin"]).off_nadir.transform("median")
    nominal = center.abs() <= NADIR_LIMIT
    print(
        f"{len(shots)} shots; {100 * nominal.mean():.1f}% where the beams' center is "
        f"within {NADIR_LIMIT} deg of nadir",
        flush=True,
    )
    # the views are analyzed as alternatives, indexed like the beams
    shots = shots[nominal].assign(pair=lambda s: s.beam)
    return shots.reset_index(drop=True)


def angle(distance_m: np.ndarray, altitude: float) -> np.ndarray:
    """Off-nadir angle (deg) of a signed ground distance (m) from nadir."""
    radius = 6371008.8
    theta = np.asarray(distance_m) / radius
    return np.degrees(
        np.arctan2(radius * np.sin(theta), radius + altitude - radius * np.cos(theta))
    )


def along_offsets(shots: pd.DataFrame, orbit) -> np.ndarray:
    """Gets the shots' offsets (km, positive ahead) along the ground track."""
    bearing, length = track_bearings(shots, orbit)
    return length * np.cos(np.radians(bearing)) / 1e3


def fit(shots: pd.DataFrame, orbit) -> None:
    """
    Prints the fitted roll angle and field of view of each beam's view, and the
    body's pitch that places the views at the shots' median offset along track.
    """
    altitude = float(np.median(shots.inst_alt))
    print(f"altitude: median {altitude / 1e3:.1f} km")
    along = along_offsets(shots, orbit)
    pitch = np.degrees(np.arctan(np.median(along) * 1e3 / altitude))
    print(
        f"along track: median {np.median(along):+.3f} km (p5 "
        f"{np.percentile(along, 5):+.3f}, p95 {np.percentile(along, 95):+.3f}); "
        f"pitch {pitch:+.4f} deg"
    )
    for beam in BEAMS:
        offsets = shots.cross_km[shots.beam == beam].to_numpy()
        center = np.median(offsets)
        half = np.percentile(np.abs(offsets - center), FIT_PERCENTILE)
        near = angle(1e3 * (center - half), altitude)
        far = angle(1e3 * (center + half), altitude)
        print(
            f"{beam}: {len(offsets)} shots, center {center:.3f} km, half-width "
            f"{half:.3f} km; roll {(near + far) / 2:.4f} deg, cross-track field of "
            f"view {far - near:.4f} deg"
        )


def score_spacing(entry, shots: pd.DataFrame) -> float:
    """
    Gets the median error (km) of the spacing across track between adjacent beams,
    measured from the beams' median offsets in bins of 10 s of each file, against the
    spacing of the views' centers.
    """
    altitude = float(np.median(shots.inst_alt))
    rolls = np.array(
        [
            spec.parameters["roll_angle"]
            for spec in entry.configurations[CONFIGURATION].instruments
        ]
    )
    # the ground distance of each roll angle, by inverting `angle` numerically
    grid = np.linspace(-20e3, 20e3, 40001)
    centers = np.interp(rolls, angle(grid, altitude), grid) / 1e3
    binned = shots.groupby(["file", "bin", "beam"]).cross_km.median().unstack("beam")
    errors = np.concatenate(
        [
            np.abs(
                (binned[BEAMS[i]] - binned[BEAMS[i + 1]]).dropna().to_numpy()
                - (centers[i] - centers[i + 1])
            )
            for i in range(len(BEAMS) - 1)
        ]
    )
    print(f"  spacing: {len(errors)} bins, median error {np.median(errors):.4f} km")
    return float(np.median(errors))


def main() -> None:
    """Validates the entry on a day of data, or fits its views."""
    run = parse_run("gedi", __doc__)
    tle = run.element_set()
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    shots = add_offsets(load_day(run.day, run.directory), orbit)
    print(
        f"{run.day:%Y-%m-%d}: {tle[0][:32]}; {len(shots)} shots in "
        f"{shots.file.nunique()} files",
        flush=True,
    )
    if run.fit:
        fit(shots, orbit)
        return
    scored = score_views(run, CONFIGURATION, orbit, shots, BEAMS, SAMPLE_PER_BEAM)
    spacing = score_spacing(run.entry, shots)
    if run.write:
        write_views(
            run,
            CONFIGURATION,
            orbit,
            tle,
            scored,
            spacing,
            BEAMS,
            f"{SHORT_NAME} version {VERSION} (every 200th shot, where the beams' "
            "center lies within 1 deg of nadir)",
        )


if __name__ == "__main__":
    main()
