"""
Validates the EMIT catalog entry against the EMIT Level 1B observation geometry
product, using the exact values of the entry through public TAT-C methods. EMIT is
an imaging spectrometer on the International Space Station, a pushbroom of 1,242
pixels across track, which images scenes (about 12 s each) by plan. The station's
position is read from the attitude and ephemeris product (inertial, rotated to
Earth-fixed). The view's angles are fitted to one day (run with `--fit`), which
needs no element set: running on a day at least 30 days from it tests the fitted
values on data held out from the fit, with a public element set.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
geometry of every 16th line (about 0.15 s apart) and every 25th pixel (and the last)
is read, with range requests, and the time of every 8th of those lines, interpolated
between them:

    python docs/validation/validate_catalog_emit.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_emit.py <data directory> 2026-11-07

The swath, about 75 km wide, is imaged in runs of scenes lasting about a minute (at
most about five), too short and sparse for a global lattice of points or grid, so a
sample of pixels is scored instead: whether the view observes each within two minutes
of its time, and the epoch's difference. The swath width is compared, line by line,
with the distance between the outermost pixels, and the swath with the strip between
them over segments of the runs (120 s each, as long as the longest time step).

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

import tatc
from tatc.analysis import compute_ground_track
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.schemas import GeneralPerturbationsOrbit
from gesdisc_scanner import Product, Reference, find_granules, fit, fit_yaw
from validate_catalog_ecostress import read_ephemeris
from validate_catalog_oco2 import area, find_epochs, observed_strip, score_width
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    get_holdout_gap,
    recompute_passes,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "emit.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "emit.npz"
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
PLATFORM = "iss"
CONFIGURATION = "swath"
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-11-07": (
        "1 25544U 98067A   26313.20747964  .00008201  00000+0  15534-3 0  9995",
        "2 25544  51.6321 302.4538 0007297   5.2479 354.8585 15.49270212594231",
    ),
}
SCENES = Product(short_name="EMITL1BRAD", version="002", host="https://data.lpdaac")
EPHEMERIS = Product(short_name="EMITL1BATT", version="002", host="https://data.lpdaac")
PIXELS, LINE_STRIDE, PIXEL_STRIDE = 1242, 16, 25
COLUMNS = np.r_[0:PIXELS:PIXEL_STRIDE, PIXELS - 1]
TIME_BAND = 9  # UTC time (decimal hours)
TIME_STRIDE = 8  # of the lines read
TIME_STEPS = (10, 30, 120)
SAMPLE = 1000
FIXTURE_PASSES = 60
# segments of the runs of scenes (s, a multiple of the time steps), and the lines of
# a segment's strip (every 4th line read, about 0.6 s apart)
SEGMENT, SEGMENT_LINE_STRIDE = 120, 4


def read_scene(url: str, session, day: datetime) -> dict:
    """
    Reads every 16th line of a scene's observation geometry with range requests: its
    time (seconds from the start of the day), and the positions of its pixels.
    """
    folder, name = url.rsplit("/", 1)
    url = f"{folder}/{name.replace('EMIT_L1B_RAD_', 'EMIT_L1B_OBS_')}"
    with h5py.File(RemoteFile(url, session, block_size=2**20)) as f:
        lines = f["location/lat"].shape[0]
        rows = np.arange(LINE_STRIDE // 2, lines, LINE_STRIDE)
        lat = f["location/lat"][list(rows), :][:, COLUMNS]
        lon = f["location/lon"][list(rows), :][:, COLUMNS]
        # the times of every 8th line read, interpolated between them (lines are
        # imaged at a constant rate), as each is a separate request
        timed = np.unique(np.r_[rows[::TIME_STRIDE], rows[-1]])
        sampled = f["obs"][list(timed), PIXELS // 2, TIME_BAND].astype(float)
        valid = sampled >= 0
        hours = (
            np.interp(rows, timed[valid], sampled[valid])
            if valid.sum() >= 2
            else np.full(len(rows), -1.0)
        )
    start = datetime(day.year, day.month, day.day, tzinfo=timezone.utc)
    # the hours of the scene's day, which may precede the analysis day
    stamp = url.split("_")[-1]
    scene_day = datetime.fromisoformat(stamp[:8]).replace(tzinfo=timezone.utc)
    # hours after midnight, in a scene that started before it
    hours = np.where(hours < int(stamp[9:11]) - 1, hours + 24, hours)
    seconds = hours * 3600 + (scene_day - start).total_seconds()
    seconds[hours < 0] = np.nan
    return {"time": seconds, "lat": lat, "lon": lon}


def load_day(day: datetime, directory: pathlib.Path) -> dict:
    """Reads and reduces the scenes of a day, or loads them from a cache."""
    cache = directory / f"emit_{day:%Y%m%d}.npz"
    if not cache.exists():
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        scenes = [
            url
            for url in find_granules(SCENES, day)
            if url.split("/")[-1].startswith("EMIT_L1B_RAD_")
        ]
        folder = directory / f"emit_{day:%Y%m%d}"
        folder.mkdir(exist_ok=True)
        done = []

        def read_cached(url: str) -> dict:
            # each scene is kept as it is read, so that an interrupted run resumes
            path = folder / f"{url.rsplit('/', 1)[-1][:-3]}.npz"
            if not path.exists():
                np.savez(path, **read_scene(url, session, day))
            done.append(url)
            if len(done) % 25 == 0:
                print(f"  {len(done)} of {len(scenes)} scenes read", flush=True)
            with np.load(path) as scene:
                return {key: scene[key] for key in scene.files}

        with concurrent.futures.ThreadPoolExecutor(4) as pool:
            parts = list(pool.map(read_cached, scenes))
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
        keep = (
            np.isfinite(data["time"])
            & (data["time"] >= max(seconds[0], 0))
            & (data["time"] <= min(seconds[-1], 86400))
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
    print(f"{int(data.pop('scenes'))} scenes", flush=True)
    return data


def find_segments(reference: Reference) -> list[np.ndarray]:
    """
    Gets the lines of segments of the runs of scenes, each of SEGMENT seconds without
    a gap, that do not cross the antimeridian.
    """
    seconds = reference.seconds
    breaks = np.nonzero(np.diff(seconds) > 1)[0] + 1
    segments = []
    for run in np.split(np.arange(len(seconds)), breaks):
        start = seconds[run[0]]
        while start + SEGMENT <= seconds[run[-1]]:
            chosen = run[(seconds[run] >= start) & (seconds[run] <= start + SEGMENT)]
            if np.ptp(reference.lon[chosen]) < 180:
                segments.append(chosen)
            start += SEGMENT
    return segments


def score_strips(satellite_for, reference: Reference, segments) -> tuple[float, dict]:
    """
    Gets the lowest median, over time steps, of the intersection over union of the
    predicted swath with the strip between the outermost pixels, over segments.
    """
    medians = {}
    for time_step in TIME_STEPS:
        satellite = satellite_for(time_step)
        scores = []
        for lines in segments:
            first = reference.time[lines[0]].to_pydatetime()
            times = [
                first + timedelta(seconds=time_step * (i + 0.5))
                for i in range(int(SEGMENT // time_step))
            ]
            predicted = compute_ground_track(satellite, times).geometry.iloc[0]
            observed = observed_strip(reference, lines[::SEGMENT_LINE_STRIDE])
            union = area(predicted.union(observed))
            scores.append(area(predicted.intersection(observed)) / union)
        medians[time_step] = float(np.median(scores))
        print(
            f"  strip IoU at {time_step} s over {len(scores)} segments: median "
            f"{medians[time_step]:.4f}, lowest {min(scores):.4f}",
            flush=True,
        )
    return min(medians.values()), medians


def validate(day, tle, orbit, reference: Reference, write: bool) -> None:
    """
    Scores the configuration on a sample of pixels and segments of the runs of scenes,
    and writes its results if requested.
    """
    entry = get_entry("emit")

    def satellite_for(time_step):
        return entry.to_satellite(orbit, time_step=time_step, platform=PLATFORM)

    rng = np.random.default_rng(1)
    lines = rng.choice(len(reference.time), SAMPLE, replace=False)
    pixels = rng.integers(0, reference.beams, SAMPLE)
    latitude = reference.lat[lines, pixels]
    longitude = reference.lon[lines, pixels]
    seconds = reference.seconds[lines]
    times = [day + timedelta(seconds=float(s)) for s in seconds]
    differences = find_epochs(satellite_for(TIME_STEPS[0]), latitude, longitude, times)
    found = np.isfinite(differences)
    results = {
        "pass_recall": 100 * float(found.mean()),
        "epoch_p95": float(summarize(differences[found])["p95"]),
        "epoch_median": float(summarize(differences[found])["median"]),
        "swath_width_median": score_width(
            satellite_for(TIME_STEPS[0]), reference, np.sort(lines[:200])
        ),
    }
    print(json.dumps(results), flush=True)
    results["footprint_iou_median"], _ = score_strips(
        satellite_for, reference, find_segments(reference)
    )
    print(json.dumps(results), flush=True)
    if not write:
        return
    picks = rng.choice(SAMPLE, FIXTURE_PASSES, replace=False)
    passes = np.stack([latitude[picks], longitude[picks], seconds[picks]], axis=-1)
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=passes
    )
    update_record(
        RECORD_PATH,
        "emit",
        CONFIGURATION,
        {
            "digest": entry.configurations[CONFIGURATION].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_emit.py",
            "reference": "EMITL1BRAD version 002, observation geometry (pixel "
            "positions and line times of every 16th line)",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": [
                {"name": name, "value": float(results[name])}
                for name in (
                    "pass_recall",
                    "epoch_p95",
                    "epoch_median",
                    "swath_width_median",
                    "footprint_iou_median",
                )
            ],
            "fixture_metrics": recompute_passes(
                satellite_for(TIME_STEPS[0]), 0, day, passes
            ),
        },
    )
    print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


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
    data.update(pixel=COLUMNS, pixels=PIXELS, line_stride=LINE_STRIDE)
    period = float(np.median(np.diff(data["time"])))
    reference = Reference(day, data, period)
    print(
        f"{arguments.day}: {len(reference.time)} lines from {reference.time[0]}, "
        f"{len(reference.gaps)} gaps over a minute; line spacing {period:.4f} s",
        flush=True,
    )
    if arguments.fit:
        for frame in ("inertial", "earth_fixed"):
            fit(reference, frame)
            fit_yaw(reference, frame)
        return
    tle = list(ELEMENTS[arguments.day])
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    validate(day, tle, orbit, reference, arguments.write)


if __name__ == "__main__":
    main()
