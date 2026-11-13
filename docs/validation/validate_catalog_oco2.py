"""
Validates the OCO-2 catalog entry against the OCO-2 Level 1B science product, using
the exact values of the entry through public TAT-C methods. OCO-2 observes either at
nadir or at the sun glint, alternating by orbit; each frame (about 0.333 s) has eight
soundings across its narrow swath, whose positions and times the product gives with
the spacecraft's position. The nadir view's angles, and the glint pointing's offset
and limits, are fitted to one day (run with `--fit`), which needs no element set:
running on a day at least 30 days from it tests the fitted values on data held out
from the fit, with a public element set.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
geometry is read, with range requests:

    python docs/validation/validate_catalog_oco2.py <data directory> 2026-10-04 --fit
    python docs/validation/validate_catalog_oco2.py <data directory> 2026-11-07

The swath, about 10 km wide, is too narrow for a global lattice of points or grid, so
a sample of soundings is scored instead: whether the view observes each within two
minutes of its time, and the epoch's difference. The swath width is compared, frame
by frame, with the distance between the first and eighth soundings, and the swath
with the strip between them over segments of the track (360 s each).

With `--write`, it writes the entry's validation record and the reduced fixtures from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
from datetime import datetime, timedelta, timezone

import h5py
import numpy as np
import pandas as pd
from pyproj import Geod, Transformer
from shapely.geometry import Polygon
from shapely.ops import unary_union
from skyfield.positionlib import ICRF
from skyfield.framelib import itrs
from skyfield.units import Distance

import tatc
from tatc.analysis import collect_observations, compute_ground_track
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit, Point
from tatc.utils import compute_projected_ray_position, compute_specular_point
from gesdisc_scanner import Product, Reference, load_day, tai93_to_seconds
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    ecef,
    get_holdout_gap,
    recompute_passes,
    transform_geometry,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
FIT_DAY = datetime(2026, 10, 4, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
OCO2_ELEMENTS = (
    "1 40059U 14035A   26313.25381781  .00000489  00000+0  11850-3 0  9999",
    "2 40059  98.2167 252.4180 0000990  88.3344 271.7968 14.57127601657173",
)
ELEMENTS = {"2026-11-07": OCO2_ELEMENTS}
MODES = {b"ND": 0, b"GL": 1, b"TG": 2}
CONFIGURATIONS = {0: "nadir", 1: "glint"}
FRAME_PERIOD = 0.333


def read_oco2(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads a granule's soundings with range requests: each frame's time and each
    sounding's (seconds from the start of the day), their positions, the spacecraft's
    position, and the frame's mode.
    """
    with h5py.File(RemoteFile(url, session, block_size=2**18)) as f:
        g = f["SoundingGeometry"]
        beam_time = tai93_to_seconds(g["sounding_time_tai93"][:], day)
        mode = np.array(
            [MODES.get(m, -1) for m in g["sounding_operation_mode"][:, 0]], dtype=int
        )
        frame = f["FrameGeometry"]
        return {
            "time": np.nanmedian(np.where(beam_time > -1e8, beam_time, np.nan), axis=1),
            "beam_time": beam_time,
            "lat": g["sounding_latitude"][:],
            "lon": g["sounding_longitude"][:],
            "sc_lat": frame["spacecraft_lat"][:],
            "sc_lon": frame["spacecraft_lon"][:],
            "sc_alt": frame["spacecraft_alt"][:],
            "mode": mode,
        }


PRODUCT = Product(short_name="OCO2_L1B_Science", version="11.2", read=read_oco2)
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "oco2.json"
FIXTURES = REPOSITORY / "tests" / "catalog" / "fixtures"
TIME_STEPS = (10, 30, 120)
SAMPLE = 1000
FIXTURE_PASSES = 60
# segments of the track (s, a multiple of the time steps) whose swaths are compared
SEGMENT, SEGMENTS, SEGMENT_FRAME_STRIDE = 360, 8, 3
GEOD = Geod(ellps="WGS84")
TO_EQUAL_AREA = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)


def find_epochs(satellite, latitude, longitude, times) -> np.ndarray:
    """
    Gets the difference (s) of the nearest predicted epoch of each sounding, within
    two minutes of its time, or NaN.
    """

    def find(index: int) -> float:
        time = times[index]
        found = collect_observations(
            Point(latitude=float(latitude[index]), longitude=float(longitude[index])),
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
        return np.array(list(pool.map(find, range(len(times)))), dtype=float)


def score_width(satellite, reference: Reference, frames: np.ndarray) -> float:
    """
    Gets the median absolute difference (km) of the distance between the projected
    ends of the slit's view and that between the first and eighth soundings.
    """
    instrument = satellite.get_instrument(0)
    track = satellite.orbit.to_gp_orbit().get_orbit_track(
        list(reference.time[frames].to_pydatetime())
    )
    differences = []
    for k, frame in enumerate(frames):
        ends = [
            compute_projected_ray_position(
                track[k],
                instrument.cross_track_field_of_view,
                1e-6,
                is_rectangular=True,
                angle=angle,
                attitude=instrument.platform_attitude,
            )
            for angle in (0, 180)
        ]
        predicted = GEOD.inv(
            float(ends[0].longitude.degrees),
            float(ends[0].latitude.degrees),
            float(ends[1].longitude.degrees),
            float(ends[1].latitude.degrees),
        )[2]
        observed = GEOD.inv(
            reference.lon[frame, 0],
            reference.lat[frame, 0],
            reference.lon[frame, -1],
            reference.lat[frame, -1],
        )[2]
        differences.append(abs(predicted - observed) / 1e3)
    return float(np.median(differences))


def find_segments(reference: Reference) -> list[np.ndarray]:
    """
    Gets the frames of segments of the track, each of SEGMENT seconds without a gap,
    spread over the day, that do not cross the antimeridian.
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
    picks = np.linspace(0, len(segments) - 1, min(SEGMENTS, len(segments)))
    return [segments[int(round(p))] for p in picks]


def observed_strip(reference: Reference, frames: np.ndarray):
    """Gets the strip between the first and eighth soundings of frames."""
    frames = frames[::SEGMENT_FRAME_STRIDE]
    quads = [
        Polygon(
            [
                (reference.lon[a, 0], reference.lat[a, 0]),
                (reference.lon[a, -1], reference.lat[a, -1]),
                (reference.lon[b, -1], reference.lat[b, -1]),
                (reference.lon[b, 0], reference.lat[b, 0]),
            ]
        ).buffer(0)
        for a, b in zip(frames[:-1], frames[1:])
    ]
    return unary_union(quads)


def area(geometry) -> float:
    """Gets the area (m2) of a geometry in longitude and latitude."""
    return transform_geometry(TO_EQUAL_AREA.transform, geometry).area


def score_strips(satellite_for, reference: Reference, segments) -> tuple[float, dict]:
    """
    Gets the lowest median, over time steps, of the intersection over union of the
    predicted swath with the strip of soundings, over segments of the track.
    """
    medians = {}
    for time_step in TIME_STEPS:
        satellite = satellite_for(time_step)
        scores = []
        for frames in segments:
            first = reference.time[frames[0]].to_pydatetime()
            times = [
                first + timedelta(seconds=time_step * (i + 0.5))
                for i in range(int(SEGMENT // time_step))
            ]
            predicted = compute_ground_track(
                satellite, times, time_step=timedelta(seconds=time_step)
            ).geometry.iloc[0]
            observed = observed_strip(reference, frames)
            union = area(predicted.union(observed))
            scores.append(area(predicted.intersection(observed)) / union)
        medians[time_step] = float(np.median(scores))
        print(f"  strip IoU at {time_step} s: {medians[time_step]:.4f}", flush=True)
    return min(medians.values()), medians


def validate_mode(day, tle, orbit, mode: int, reference: Reference, write: bool):
    """
    Scores a mode's configuration on a sample of soundings and segments of the track,
    and writes its results if requested.
    """
    key = CONFIGURATIONS[mode]
    entry = get_entry("oco2")

    def satellite_for(time_step):
        return entry.to_satellite(orbit, key, time_step=time_step, platform="oco-2")

    rng = np.random.default_rng(1)
    frames = rng.choice(len(reference.time), SAMPLE, replace=False)
    soundings = rng.integers(0, reference.beams, SAMPLE)
    latitude = reference.lat[frames, soundings]
    longitude = reference.lon[frames, soundings]
    times = [
        datetime.fromtimestamp(
            day.timestamp() + float(reference.beam_seconds[f, s]), tz=timezone.utc
        )
        for f, s in zip(frames, soundings)
    ]
    differences = find_epochs(satellite_for(TIME_STEPS[0]), latitude, longitude, times)
    found = np.isfinite(differences)
    results = {
        "pass_recall": 100 * float(found.mean()),
        "epoch_p95": float(summarize(differences[found])["p95"]),
        "epoch_median": float(summarize(differences[found])["median"]),
        "swath_width_median": score_width(
            satellite_for(TIME_STEPS[0]), reference, np.sort(frames[:200])
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
    passes = np.stack(
        [
            latitude[picks],
            longitude[picks],
            reference.beam_seconds[frames[picks], soundings[picks]],
        ],
        axis=-1,
    )
    fixture = FIXTURES / f"oco2_{key}.npz"
    update_fixture(
        fixture, tle=np.array(tle), day=np.array(day.isoformat()), passes=passes
    )
    update_record(
        RECORD_PATH,
        "oco2",
        key,
        {
            "digest": entry.configurations[key].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_oco2.py",
            "reference": f"OCO2_L1B_Science version 11.2 (sounding positions and "
            f"times, {key} orbits)",
            "period": f"{day:%Y-%m-%d}",
            "platform": "oco-2",
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
    print(f"wrote {RECORD_PATH} and {fixture}")


def load_reference(day: datetime, directory: pathlib.Path, mode: int) -> Reference:
    """Loads the soundings of a day in a mode."""
    data = load_day(PRODUCT, day, directory)
    chosen = data["mode"] == mode
    data = {k: v[chosen] for k, v in data.items() if k != "granules"}
    return Reference(day, data, FRAME_PERIOD)


def fit_glint(reference: Reference) -> None:
    """
    Prints the pointing of the glint orbits: the angle from nadir of the soundings'
    center, the angle between it and the specular point seen from the spacecraft,
    and the Sun's elevation at the specular point.
    """
    t = timescale.from_datetimes(list(reference.time.to_pydatetime()))
    # the spacecraft's position, from its Earth-fixed coordinates
    rotation = itrs.rotation_at(t)
    position = np.einsum("jik,kj->ik", rotation, reference.position) / 1.0
    track = ICRF(Distance(m=position).au, t=t, center=399)
    specular = compute_specular_point(track)
    spot = ecef(specular.longitude.degrees, specular.latitude.degrees)
    center = np.mean(ecef(reference.lon, reference.lat), axis=1)
    lon, lat, _ = TO_GEODETIC.transform(*reference.position.T)
    nadir = -up(lon, np.asarray(lat))

    def angle(a, b):
        cosine = np.sum(a * b, axis=-1) / (
            np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1)
        )
        return np.degrees(np.arccos(np.clip(cosine, -1, 1)))

    to_center = center - reference.position
    to_spot = spot - reference.position
    off_nadir = angle(to_center, nadir)
    offset = angle(to_center, to_spot)
    spot_off_nadir = angle(to_spot, nadir)
    # the Sun's elevation at the specular point equals the satellite's there
    elevation = 90 - angle(
        reference.position - spot,
        up(specular.longitude.degrees, specular.latitude.degrees),
    )
    finite = np.isfinite(offset)
    print(
        f"glint frames: {len(offset)}, with a specular point: {finite.mean():.1%}; "
        f"soundings' center from nadir: max {np.nanmax(off_nadir):.3f} deg, "
        f"99th percentile {np.nanpercentile(off_nadir, 99):.3f} deg"
    )
    print(
        f"angle from the specular point: median {np.nanmedian(offset):.4f} deg, 95th "
        f"percentile {np.nanpercentile(offset, 95):.4f} deg"
    )
    # the offset toward nadir, by the specular point's angle from nadir
    pointed = finite & (offset < 15)
    toward = spot_off_nadir - off_nadir
    bins = (spot_off_nadir[pointed] // 2) * 2 + 1
    profile = [
        (float(b), round(float(np.median(toward[pointed][bins == b])), 3))
        for b in np.unique(bins)
        if np.sum(bins == b) >= 20
    ]
    print(f"offset_angle_profile: {profile}")
    print(
        f"pointed frames: {pointed.mean():.1%}; the specular point's angle from nadir "
        f"up to {np.nanmax(spot_off_nadir[pointed]):.3f} deg, the Sun's elevation at it "
        f"down to {np.nanmin(elevation[pointed]):.3f} deg; residual of the profile: "
        f"median {np.median(np.abs(toward[pointed] - np.interp(spot_off_nadir[pointed], *zip(*profile)))):.3f} deg"
    )


def argument_of_latitude(reference: Reference) -> np.ndarray:
    """Gets the approximate argument of latitude (deg) of each frame."""
    _, lat, _ = TO_GEODETIC.transform(*reference.position.T)
    lat = np.asarray(lat)
    ascending = np.gradient(lat, reference.seconds) > 0
    u = np.degrees(
        np.arcsin(np.clip(np.sin(np.radians(lat)) / np.sin(np.radians(98.2)), -1, 1))
    )
    return np.where(ascending, u % 360, 180 - u)


def fit_nadir(reference: Reference) -> None:
    """
    Prints the slit of the nadir orbits: the angle between the first and eighth
    soundings seen from the spacecraft, and the yaw angle profile that turns the
    slit from the cross-track direction (the negative of the slit's angle from it,
    toward the front, in TAT-C's view angles), by argument of latitude.
    """
    lon, lat, _ = TO_GEODETIC.transform(*reference.position.T)
    zenith = up(lon, np.asarray(lat))
    motion = reference.velocity("earth_fixed")
    along_axis = motion - np.sum(motion * zenith, axis=-1)[:, None] * zenith
    along_axis /= np.linalg.norm(along_axis, axis=-1)[:, None]
    left_axis = np.cross(zenith, along_axis)
    rays = ecef(reference.lon, reference.lat) - reference.position[:, np.newaxis]
    down = -np.sum(rays * zenith[:, np.newaxis], axis=-1)
    along = np.sum(rays * along_axis[:, np.newaxis], axis=-1) / down
    cross = np.sum(rays * left_axis[:, np.newaxis], axis=-1) / down
    first, last = rays[:, 0], rays[:, -1]
    extent = np.degrees(
        np.arccos(
            np.sum(first * last, axis=-1)
            / (np.linalg.norm(first, axis=-1) * np.linalg.norm(last, axis=-1))
        )
    )
    slit = np.degrees(
        np.arctan2(along[:, -1] - along[:, 0], cross[:, -1] - cross[:, 0])
    )
    middle = np.degrees(np.arctan(np.mean(along, axis=1)))
    print(
        f"nadir frames: {len(slit)}; first to eighth sounding: median "
        f"{np.median(extent):.4f} deg; center along track: median {np.median(middle):+.4f} deg"
    )
    u = argument_of_latitude(reference)
    bins = (u // 4) * 4 + 2
    profile = [
        (float(b), round(float(-np.median(slit[bins == b])), 2))
        for b in np.unique(bins)
        if np.sum(bins == b) >= 50
    ]
    print(f"yaw_angle_profile: {profile}")


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
    for mode, key in CONFIGURATIONS.items():
        reference = load_reference(day, arguments.directory, mode)
        print(
            f"{arguments.day} ({key}): {len(reference.time)} frames from "
            f"{reference.time[0]}, {len(reference.gaps)} gaps over a minute",
            flush=True,
        )
        if arguments.fit:
            if mode == 0:
                fit_nadir(reference)
            else:
                fit_glint(reference)
            continue
        tle = list(ELEMENTS[arguments.day])
        orbit = GeneralPerturbationsOrbit.from_tle(tle)
        validate_mode(day, tle, orbit, mode, reference, arguments.write)


if __name__ == "__main__":
    main()
