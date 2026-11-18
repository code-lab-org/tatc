"""
Validates the Poseidon-4 catalog entry against the Sentinel-6 Michael Freilich Level 2
low-resolution near-real-time product (reduced), using the exact values of the entry
through public TAT-C methods. Poseidon-4 is a nadir radar altimeter, whose
measurements (1 Hz) lie along the ground track. A sample of the measurements of a day
is scored: whether the nadir view observes each within two minutes of its time, the
epoch's difference, and its distance across track from the ground track of a public
element set. The entry has no fitted values, so any day tests it. With `--platform
sentinel-6b`, it scores the same configuration on Sentinel-6B, as a transfer.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
measurements' times and positions are read, with range requests:

    python validation/catalog/validate_catalog_poseidon4.py <data directory> 2026-11-07
    python validation/catalog/validate_catalog_poseidon4.py <data directory> 2026-11-08 \
        --platform sentinel-6b

Without a day, it takes that of the entry's record (or of the transfer's). With
`--write`, it writes the entry's validation record and the reduced fixture from which
unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import concurrent.futures
import json
import pathlib
from datetime import datetime, timedelta

import earthaccess
import h5py
import numpy as np
import pandas as pd
from pyproj import Geod
from skyfield.api import wgs84

from tatc.analysis import collect_observations
from tatc.catalog.metrics import summarize
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit, Point
from validation_io import (
    RemoteFile,
    get_with_retries,
    login_earthdata,
    parse_run,
    recompute_passes,
    write_results,
)

CONFIGURATION = "nadir"
SHORT_NAMES = {
    "sentinel-6a": "JASON_CS_S6A_L2_ALT_LR_RED_OST_NRT_F",
    "sentinel-6b": "JASON_CS_S6B_L2_ALT_LR_RED_OST_NRT_G",
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


def load_day(day: datetime, directory: pathlib.Path, platform: str) -> pd.DataFrame:
    """Reads the measurements of a day, or loads them from a cache."""
    suffix = "" if platform == "sentinel-6a" else f"_{platform}"
    cache = directory / f"poseidon4_{day:%Y%m%d}{suffix}.npz"
    start = pd.Timestamp(day)
    if not cache.exists():
        items = get_with_retries(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params={
                "short_name": SHORT_NAMES[platform],
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
        login_earthdata()
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
    run = parse_run("poseidon4", __doc__)
    if run.fit:
        raise SystemExit("the entry has no fitted values")
    day, platform, entry = run.day, run.platform, run.entry
    data = load_day(day, run.directory, platform)
    print(f"{day:%Y-%m-%d}: {len(data)} measurements", flush=True)
    tle = run.element_set()
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
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
        satellite = entry.to_satellite(orbit, time_step=time_step, platform=platform)
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
    if not run.write:
        return
    reference = f"{SHORT_NAMES[platform]} (1 Hz measurement positions and times)"
    if run.transfer:
        write_results(
            run,
            CONFIGURATION,
            reference=reference,
            metrics={
                name: results[name]
                for name in ("segment_recall", "track_offset_median", "epoch_p95")
            },
        )
        return
    # the fixture's passes are drawn after the sample, from the same generator
    picks = sample.iloc[rng.choice(len(sample), FIXTURE_PASSES, replace=False)]
    passes = np.stack(
        [
            picks.latitude.to_numpy(),
            picks.longitude.to_numpy(),
            (picks.time - pd.Timestamp(day)).dt.total_seconds().to_numpy(),
        ],
        axis=-1,
    )
    satellite = entry.to_satellite(orbit, time_step=TIME_STEPS[0], platform=platform)
    write_results(
        run,
        CONFIGURATION,
        reference=reference,
        metrics={
            name: results[name]
            for name in (
                "segment_recall",
                "track_offset_median",
                "epoch_p95",
                "epoch_median",
            )
        },
        fixture={
            "tle": np.array(tle),
            "day": np.array(day.isoformat()),
            "passes": passes,
        },
        fixture_metrics=recompute_passes(satellite, 0, day, passes),
    )


if __name__ == "__main__":
    main()
