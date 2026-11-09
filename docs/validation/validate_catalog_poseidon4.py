"""
Validates the Poseidon-4 catalog entry against the Sentinel-6 Michael Freilich Level 2
low-resolution near-real-time product (reduced), using the exact values of the entry
through public TAT-C methods. Poseidon-4 is a nadir radar altimeter, whose
measurements (1 Hz) lie along the ground track. A sample of the measurements of a day
is scored: whether the nadir view observes each within two minutes of its time, the
epoch's difference, and its distance across track from the ground track of a public
element set. The entry has no fitted values, so any day tests it.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
measurements' times and positions are read, with range requests:

    python docs/validation/validate_catalog_poseidon4.py <data directory> 2026-11-07

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
from datetime import datetime, timedelta, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd
from pyproj import Geod
from skyfield.api import wgs84

import tatc
from tatc.analysis import collect_observations
from tatc.catalog.metrics import summarize
from tatc.catalog import get_entry
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit, Point
from validation_io import (
    RemoteFile,
    get_with_retries,
    recompute_passes,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = (
    REPOSITORY / "src" / "tatc" / "resources" / "validation" / "poseidon4.json"
)
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "poseidon4.npz"
SHORT_NAME = "JASON_CS_S6A_L2_ALT_LR_RED_OST_NRT_F"
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-07": (
        "1 46984U 20086A   26313.07249362 -.00000031  00000+0  13085-3 0  9999",
        "2 46984  66.0401 194.1205 0007726 269.8082  90.2048 12.80930425278889",
    ),
}
EPOCH = pd.Timestamp("2000-01-01", tz="UTC")
TIME_STEPS = (10, 30, 120)
SAMPLE = 1000
FIXTURE_PASSES = 60
GEOD = Geod(ellps="WGS84")


def read_granule(url: str, session) -> pd.DataFrame:
    """Reads the 1 Hz measurements' times and positions of a granule."""
    with h5py.File(RemoteFile(url, session)) as f:
        g = f["data_01"]
        seconds = g["time"][:].astype(float)
        latitude = g["latitude"][:] * g["latitude"].attrs.get("scale_factor", 1.0)
        longitude = g["longitude"][:] * g["longitude"].attrs.get("scale_factor", 1.0)
    return pd.DataFrame(
        {
            "time": EPOCH + pd.to_timedelta(seconds, unit="s"),
            "latitude": np.asarray(latitude, dtype=float),
            "longitude": (np.asarray(longitude, dtype=float) + 180) % 360 - 180,
        }
    )


def load_day(day: datetime, directory: pathlib.Path) -> pd.DataFrame:
    """Reads the measurements of a day, or loads them from a cache."""
    cache = directory / f"poseidon4_{day:%Y%m%d}.npz"
    start = pd.Timestamp(day)
    if not cache.exists():
        items = get_with_retries(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params={
                "short_name": SHORT_NAME,
                "temporal": f"{day - timedelta(hours=2):%Y-%m-%dT%H:%M:%SZ},"
                f"{day + timedelta(days=1, hours=2):%Y-%m-%dT%H:%M:%SZ}",
                "page_size": 200,
            },
        ).json()["items"]
        urls = [
            next(
                u["URL"]
                for u in item["umm"]["RelatedUrls"]
                if u["URL"].startswith("https://archive.podaac")
                and u["URL"].endswith(".nc")
            )
            for item in items
        ]
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        with concurrent.futures.ThreadPoolExecutor(4) as pool:
            parts = list(pool.map(lambda url: read_granule(url, session), urls))
        data = pd.concat(parts).drop_duplicates("time").sort_values("time")
        data = data[(data.time >= start) & (data.time < start + pd.Timedelta(days=1))]
        data = data[data.latitude.abs() <= 90]
        np.savez_compressed(
            cache,
            seconds=(data.time - start).dt.total_seconds().to_numpy(),
            latitude=data.latitude.to_numpy(),
            longitude=data.longitude.to_numpy(),
        )
    cached = np.load(cache)
    return pd.DataFrame(
        {
            "time": start + pd.to_timedelta(cached["seconds"], unit="s"),
            "latitude": cached["latitude"],
            "longitude": cached["longitude"],
        }
    )


def track_offsets(orbit, sample: pd.DataFrame) -> np.ndarray:
    """Gets the distances (km) across track of measurements from the ground track."""
    satellite = orbit.elements[0].to_skyfield()
    times = timescale.from_datetimes(list(sample.time.dt.floor("us")))
    now = wgs84.subpoint_of(satellite.at(times))
    later = wgs84.subpoint_of(satellite.at(timescale.tt_jd(times.tt + 1 / 86400)))
    heading = GEOD.inv(
        now.longitude.degrees,
        now.latitude.degrees,
        later.longitude.degrees,
        later.latitude.degrees,
    )[0]
    azimuth, _, length = GEOD.inv(
        now.longitude.degrees,
        now.latitude.degrees,
        sample.longitude.to_numpy(),
        sample.latitude.to_numpy(),
    )
    return length * np.sin(np.radians(azimuth - heading)) / 1e3


def find_epochs(satellite, sample: pd.DataFrame) -> np.ndarray:
    """Gets the difference (s) of each measurement's nearest predicted epoch, or NaN."""

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


def main() -> None:
    """Validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    data = load_day(day, arguments.directory)
    print(f"{arguments.day}: {len(data)} measurements", flush=True)
    tle = list(ELEMENTS[arguments.day])
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("poseidon4")
    rng = np.random.default_rng(1)
    sample = data.iloc[rng.choice(len(data), SAMPLE, replace=False)].reset_index(
        drop=True
    )
    offsets = np.abs(track_offsets(orbit, sample))
    results = {"track_offset_median": float(np.median(offsets))}
    print(
        f"track offset: median {np.median(offsets):.3f} km, 95th percentile "
        f"{np.percentile(offsets, 95):.3f} km",
        flush=True,
    )
    by_step = {}
    for time_step in TIME_STEPS:
        satellite = entry.to_satellite(
            orbit, time_step=time_step, platform="sentinel-6a"
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
    if not arguments.write:
        return
    picks = sample.iloc[rng.choice(len(sample), FIXTURE_PASSES, replace=False)]
    passes = np.stack(
        [
            picks.latitude.to_numpy(),
            picks.longitude.to_numpy(),
            (picks.time - pd.Timestamp(day)).dt.total_seconds().to_numpy(),
        ],
        axis=-1,
    )
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=passes
    )
    satellite = entry.to_satellite(
        orbit, time_step=TIME_STEPS[0], platform="sentinel-6a"
    )
    update_record(
        RECORD_PATH,
        "poseidon4",
        "nadir",
        {
            "digest": entry.configurations["nadir"].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_poseidon4.py",
            "reference": f"{SHORT_NAME} (1 Hz measurement positions and times)",
            "period": f"{day:%Y-%m-%d}",
            "platform": "sentinel-6a",
            "metrics": [
                {"name": name, "value": float(results[name])}
                for name in (
                    "segment_recall",
                    "track_offset_median",
                    "epoch_p95",
                    "epoch_median",
                )
            ],
            "fixture_metrics": recompute_passes(satellite, 0, day, passes),
        },
    )
    print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()
