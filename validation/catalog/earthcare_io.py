"""
Reads the Level 1 products of EarthCARE from ESA's MAAP, and validates the catalog
entries of its profilers (ATLID and CPR), whose profiles lie along a narrow track.

Products are found in the MAAP catalogue (STAC, no authentication) and read with range
requests, with an access token exchanged for the MAAP offline token in `ESA_MAAP_KEY`.
Only the needed fields of each frame (an eighth of an orbit) are read, and a day is
cached in a directory.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import os
import pathlib
import threading
import time as clock
from collections.abc import Callable
from datetime import datetime, timedelta

import h5py
import numpy as np
import pandas as pd
import requests

from tatc.analysis import collect_observations
from tatc.catalog.metrics import summarize
from tatc.schemas import Point
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    ecef,
    get_with_retries,
    normalize,
    recompute_passes,
    up,
    update_fixture,
    write_record,
)

CATALOGUE = "https://catalog.maap.eo.esa.int/catalogue/collections"
TOKEN_URL = "https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token"
# the public client of MAAP's offline tokens, as in ESA's examples of data access
CLIENT = {
    "client_id": "offline-token",
    "client_secret": "p1eL7uonXs6MDxtGbgKdPVRAmnGxHpVE",
}
VALIDATED = "EarthCAREL1Validated_MAAP"
EPOCH = pd.Timestamp("2000-01-01", tz="UTC")
# the products' times count elapsed seconds since 2000-01-01 00:00 UTC, including the
# leap seconds since (five, from 2017), so they run 5 s ahead of UTC: the spacecraft
# positions they record match an element set's orbit 5.0 s before their times, and
# a frame's first time follows the UTC sensing start of its header by 5.1 s
LEAP_SECONDS_SINCE_2000 = 5
EARTH_ROTATION = 7.2921150e-5  # rad/s
TIME_STEPS = (10, 30, 120)


class MaapSession:
    """A session with an access token, renewed before it expires."""

    def __init__(self):
        if "ESA_MAAP_KEY" not in os.environ:
            raise RuntimeError("MAAP access needs an offline token in ESA_MAAP_KEY")
        self.session = requests.Session()
        self.expires = 0.0
        self.lock = threading.Lock()

    def get(self) -> requests.Session:
        """Gets the session, with a current access token."""
        with self.lock:
            if clock.time() > self.expires - 600:
                response = requests.post(
                    TOKEN_URL,
                    data={
                        **CLIENT,
                        "grant_type": "refresh_token",
                        "refresh_token": os.environ["ESA_MAAP_KEY"],
                        "scope": "offline_access openid",
                    },
                    timeout=60,
                )
                response.raise_for_status()
                token = response.json()
                self.session.headers["Authorization"] = (
                    f"Bearer {token['access_token']}"
                )
                self.expires = clock.time() + float(token.get("expires_in", 3600))
            return self.session


def find_frames(product_type: str, day: datetime, collection: str = VALIDATED):
    """
    Gets the identifiers and URLs of the frames of a product from an hour before to
    after a day, keeping the latest processing of each frame.
    """
    url = f"{CATALOGUE}/{collection}/items"
    params = {
        "productType": product_type,
        "datetime": f"{day - timedelta(hours=1):%Y-%m-%dT%H:%M:%SZ}/"
        f"{day + timedelta(days=1, hours=1):%Y-%m-%dT%H:%M:%SZ}",
        "limit": 500,
    }
    frames = {}
    while url:
        page = get_with_retries(url, params=params).json()
        for feature in page.get("features", []):
            properties = feature["properties"]
            key = (properties["sat:absolute_orbit"], properties["grid:code"])
            created = properties.get("created", "")
            if key not in frames or created > frames[key][2]:
                frames[key] = (
                    feature["id"],
                    feature["assets"]["enclosure_h5"]["href"],
                    created,
                )
        url = next(
            (link["href"] for link in page.get("links", []) if link["rel"] == "next"),
            None,
        )
        params = None
    return sorted((name, href) for name, href, _ in frames.values())


def to_seconds(seconds_since_2000, day: datetime) -> np.ndarray:
    """
    Converts the products' times (elapsed seconds since 2000, with leap seconds) to
    UTC seconds from the start of a day.
    """
    return (
        np.asarray(seconds_since_2000, dtype=float)
        - LEAP_SECONDS_SINCE_2000
        - (pd.Timestamp(day) - EPOCH).total_seconds()
    )


def load_day(
    product_type: str,
    day: datetime,
    directory: pathlib.Path,
    read: Callable[[h5py.File, datetime], dict],
    name: str | None = None,
    collection: str = VALIDATED,
    threads: int = 4,
    orbit_stride: int = 1,
) -> dict:
    """
    Reads fields of the frames of a day, with range requests, or loads them from a
    cache. `read` gets a frame's fields, including `time` (seconds from the start of
    the day) along its first axis. Keeps each time once, in order, within the day.
    With `orbit_stride`, reads only the frames of orbits whose number it divides.
    """
    cache = directory / f"{(name or product_type).lower()}_{day:%Y%m%d}.npz"
    if cache.exists():
        return dict(np.load(cache))
    frames = [
        frame
        for frame in find_frames(product_type, day, collection)
        if int(frame[0].rsplit("_", 1)[1][:-1]) % orbit_stride == 0
    ]
    maap = MaapSession()

    def read_frame(frame):
        _, href = frame
        for attempt in range(3):
            try:
                with h5py.File(RemoteFile(href, maap.get(), block_size=2**20)) as f:
                    return read(f, day)
            except (OSError, requests.RequestException):
                if attempt == 2:
                    raise
                clock.sleep(5)
        return None

    with concurrent.futures.ThreadPoolExecutor(threads) as pool:
        parts = list(pool.map(read_frame, frames))
    reduced = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    _, first = np.unique(np.round(reduced["time"], 4), return_index=True)
    keep = first[(reduced["time"][first] >= 0) & (reduced["time"][first] < 86400)]
    np.savez_compressed(
        cache, frames=len(frames), **{k: v[keep] for k, v in reduced.items()}
    )
    print(f"cached {len(frames)} frames of {product_type} in {cache.name}", flush=True)
    return dict(np.load(cache))


class Profiles:
    """
    A day of profiles of a single-beam profiler: each profile's time (seconds from
    the start of the day), its position on the ground, and the spacecraft's position.
    """

    def __init__(self, day: datetime, data: dict):
        valid = (
            np.isfinite(data["time"])
            & (np.abs(data["lat"]) <= 90)
            & (np.abs(data["lon"]) <= 360)
            & np.all(np.isfinite(data["sc_position"]), axis=-1)
            & (np.linalg.norm(data["sc_position"], axis=-1) > 6e6)
            & (np.linalg.norm(data["sc_position"], axis=-1) < 8e6)
        )
        order = np.argsort(data["time"][valid])
        self.day = day
        self.seconds = data["time"][valid][order].astype(float)
        self.lat = data["lat"][valid][order].astype(float)
        self.lon = (data["lon"][valid][order].astype(float) + 180) % 360 - 180
        self.position = data["sc_position"][valid][order].astype(float)
        self.time = pd.DatetimeIndex(
            pd.Timestamp(day) + pd.to_timedelta(self.seconds, unit="s")
        )
        period = np.diff(self.seconds)
        self.gaps = [
            (self.time[i].to_pydatetime(), self.time[i + 1].to_pydatetime())
            for i in np.nonzero(period > 60)[0]
        ]

    def velocity(self, frame: str) -> np.ndarray:
        """Gets the spacecraft's velocity (m/s), Earth-fixed or inertial."""
        velocity = np.gradient(self.position, self.seconds, axis=0)
        if frame == "inertial":
            velocity = velocity + np.cross(
                [0, 0, EARTH_ROTATION], self.position, axis=-1
            )
        return velocity

    def pointing(self, frame: str) -> tuple[np.ndarray, np.ndarray]:
        """
        Gets each profile's angles (deg) from the geodetic nadir along track (positive
        forward of the motion in the frame) and across track (positive to the left).
        """
        lon, lat, _ = TO_GEODETIC.transform(*self.position.T)
        zenith = up(lon, lat)
        motion = self.velocity(frame)
        forward = normalize(motion - np.sum(motion * zenith, axis=-1)[:, None] * zenith)
        left = np.cross(zenith, forward)
        ray = ecef(self.lon, self.lat) - self.position
        down = -np.sum(ray * zenith, axis=-1)
        along = np.degrees(np.arctan2(np.sum(ray * forward, axis=-1), down))
        cross = np.degrees(np.arctan2(np.sum(ray * left, axis=-1), down))
        return along, cross


def fit_pointing(profiles: Profiles) -> None:
    """
    Prints the profiler's fitted pointing: the median angles along and across track
    from the geodetic nadir, with their spread and their medians by latitude and
    direction, in the inertial and Earth-fixed frames, and the period of profiles.
    """
    ascending = np.gradient(profiles.lat) > 0
    for frame in ("inertial", "earth_fixed"):
        along, cross = profiles.pointing(frame)
        print(
            f"{frame}: along track {np.median(along):+.4f} deg (p5 "
            f"{np.percentile(along, 5):+.4f}, p95 {np.percentile(along, 95):+.4f}); "
            f"across track {np.median(cross):+.4f} deg (p5 "
            f"{np.percentile(cross, 5):+.4f}, p95 {np.percentile(cross, 95):+.4f})"
        )
        for name, chosen in (("ascending", ascending), ("descending", ~ascending)):
            bands = []
            for low in range(-90, 90, 30):
                band = chosen & (profiles.lat >= low) & (profiles.lat < low + 30)
                if np.any(band):
                    bands.append(
                        f"{low:+d}..{low + 30:+d}: {np.median(along[band]):+.3f}/"
                        f"{np.median(cross[band]):+.3f}"
                    )
            print(f"  {name} along/across by latitude: " + ", ".join(bands))
    altitude = TO_GEODETIC.transform(*profiles.position.T)[2]
    print(
        f"profile period: median {np.median(np.diff(profiles.seconds)):.4f} s; "
        f"altitude median {np.median(altitude) / 1e3:.2f} km "
        f"(p5 {np.percentile(altitude, 5) / 1e3:.2f}, "
        f"p95 {np.percentile(altitude, 95) / 1e3:.2f})"
    )


def track_offsets(satellite, sample: pd.DataFrame) -> np.ndarray:
    """
    Gets the distances (km) across track of profiles from the predicted track of the
    view's center (the ground track, for a view at nadir), near each profile's time.
    """
    instrument = satellite.instruments[0]
    orbit = satellite.orbit.to_gp_orbit()
    offsets = np.arange(-30, 30.5, 0.5)
    result = []
    for row in sample.itertuples():
        times = [
            row.time.floor("us").to_pydatetime() + timedelta(seconds=float(s))
            for s in offsets
        ]
        center = instrument.compute_footprint_center(orbit.get_orbit_track(times))
        xyz = ecef(center.longitude.degrees, center.latitude.degrees)
        target = ecef(row.longitude, row.latitude)
        k = int(
            np.clip(np.argmin(np.linalg.norm(xyz - target, axis=-1)), 1, len(xyz) - 2)
        )
        direction = normalize(xyz[k + 1] - xyz[k - 1])
        relative = target - xyz[k]
        result.append(
            np.linalg.norm(relative - np.dot(relative, direction) * direction) / 1e3
        )
    return np.array(result)


def find_epochs(satellite, sample: pd.DataFrame) -> np.ndarray:
    """Gets the difference (s) of each profile's nearest predicted epoch, or NaN."""

    def find(row):
        time = row.time.floor("us").to_pydatetime()
        found = collect_observations(
            Point(latitude=row.latitude, longitude=row.longitude),
            satellite,
            time - timedelta(minutes=2),
            time + timedelta(minutes=2),
        )
        if found.empty:
            return np.nan
        return min(
            ((e - pd.Timestamp(time)).total_seconds() for e in found.epoch), key=abs
        )

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        return np.array(list(pool.map(find, sample.itertuples())), dtype=float)


def validate_profiler(
    entry,
    configuration: str,
    platform: str,
    profiles: Profiles,
    orbit,
    tle,
    day: datetime,
    record_path: pathlib.Path,
    fixture_path: pathlib.Path,
    source: str,
    reference: str,
    write: bool,
    sample_size: int = 1000,
    fixture_passes: int = 60,
) -> dict:
    """
    Scores a profiler's configuration on a sample of a day's profiles: the share
    within the view at their times (segment recall), at time steps up to 120 s, the
    epochs' differences, and the distance across track from the predicted track of
    the view's center. Writes the record and the reduced fixture if requested.
    """
    rng = np.random.default_rng(1)
    picks = np.sort(rng.choice(len(profiles.seconds), sample_size, replace=False))
    sample = pd.DataFrame(
        {
            "time": profiles.time[picks],
            "latitude": profiles.lat[picks],
            "longitude": profiles.lon[picks],
        }
    )
    satellite = entry.to_satellite(
        orbit, configuration=configuration, time_step=TIME_STEPS[0], platform=platform
    )
    offsets = track_offsets(satellite, sample)
    results = {"track_offset_median": float(np.median(offsets))}
    print(
        f"track offset: median {np.median(offsets):.3f} km, 95th percentile "
        f"{np.percentile(offsets, 95):.3f} km",
        flush=True,
    )
    by_step = {}
    for time_step in TIME_STEPS:
        satellite = entry.to_satellite(
            orbit, configuration=configuration, time_step=time_step, platform=platform
        )
        differences = find_epochs(satellite, sample)
        found = np.isfinite(differences)
        by_step[time_step] = {
            "segment_recall": 100 * float(found.mean()),
            "epoch_p95": float(summarize(differences[found])["p95"]),
            "epoch_median": float(summarize(differences[found])["median"]),
        }
        print(
            f"  time step {time_step} s: {json.dumps(by_step[time_step])}", flush=True
        )
    # the worst over the time steps
    results["segment_recall"] = min(r["segment_recall"] for r in by_step.values())
    results["epoch_p95"] = max(r["epoch_p95"] for r in by_step.values())
    results["epoch_median"] = by_step[TIME_STEPS[0]]["epoch_median"]
    print(json.dumps(results), flush=True)
    if not write:
        return results
    chosen = sample.iloc[rng.choice(len(sample), fixture_passes, replace=False)]
    passes = np.stack(
        [
            chosen.latitude.to_numpy(),
            chosen.longitude.to_numpy(),
            (chosen.time - pd.Timestamp(day)).dt.total_seconds().to_numpy(),
        ],
        axis=-1,
    )
    update_fixture(
        fixture_path, tle=np.array(tle), day=np.array(day.isoformat()), passes=passes
    )
    satellite = entry.to_satellite(
        orbit, configuration=configuration, time_step=TIME_STEPS[0], platform=platform
    )
    write_record(
        record_path,
        entry,
        configuration,
        source=source,
        reference=reference,
        period=f"{day:%Y-%m-%d}",
        platforms=[platform],
        metrics=[
            {"name": name, "value": float(results[name])}
            for name in (
                "segment_recall",
                "track_offset_median",
                "epoch_p95",
                "epoch_median",
            )
        ],
        fixture_metrics=recompute_passes(satellite, 0, day, passes),
    )
    print(f"wrote {record_path} and {fixture_path}")
    return results


def check_time_scale(
    product_type: str, day: datetime, time_field: str, collection: str = VALIDATED
) -> None:
    """
    Prints the difference (s) between the first time of a frame of a day and the
    sensing start in its header, which reveals a time scale other than UTC (such as
    GPS, 18 s ahead).
    """
    frames = find_frames(product_type, day, collection)
    _, href = frames[len(frames) // 2]
    with h5py.File(RemoteFile(href, MaapSession().get(), block_size=2**20)) as f:
        first = EPOCH + pd.Timedelta(
            seconds=float(np.nanmin(f[time_field][:])) - LEAP_SECONDS_SINCE_2000
        )
        header = f["HeaderData/VariableProductHeader/MainProductHeader"]
        start = header["sensingStartTime"][()]
    # the header gives the time scale before the time, as in "UTC=2026-09-27T11:49:32"
    start = (start.decode() if isinstance(start, bytes) else str(start)).split("=")
    print(f"sensing start's time scale: {start[0] if len(start) > 1 else 'unstated'}")
    start = pd.Timestamp(start[-1])
    if start.tzinfo is None:
        start = start.tz_localize("UTC")
    print(
        f"first time {first.isoformat()} vs sensing start {start.isoformat()}: "
        f"{(first - start).total_seconds():+.3f} s"
    )


__all__ = [
    "MaapSession",
    "Profiles",
    "check_time_scale",
    "find_epochs",
    "find_frames",
    "fit_pointing",
    "load_day",
    "to_seconds",
    "track_offsets",
    "validate_profiler",
]
