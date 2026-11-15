"""
Validates the TMI catalog entry against the GPM-format Level 1C TMI product of TRMM,
using the exact values of the entry through public TAT-C methods. The entry's scan is
fitted to the positions of the samples of each group (S1, 10 GHz; S2, 19 to 37 GHz;
S3, 85.5 GHz) on one day (run with `--fit`), from the spacecraft positions the product
records, so the fit needs no element set; running on a day at least 30 days from it
tests the fitted scan on data held out from the fit. Both days are in TRMM's forward
orientation (the product's `SCorientation` of 0) throughout. Element sets are read
from a file of the history of TRMM's element sets. The tests follow those of SSMIS:

* A point is observed by a feed on a pass if it lies in a cell of the lattice of its
  group's samples, at the time it crosses the cone, from the recorded spacecraft positions.
* The swath width is compared with the distance between the first and last samples.
* The swath area over a span of the day (from noon) is compared with the union of
  the cells.

It also validates the entry's sample size configuration, as for AMSR2: sizes at the
center of the scan are compared with the footprints published in the TC PRIMED
documentation, which the gate tests, and the spacing of samples along the arc and
between scans with that measured.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`. Only the
needed fields of the product's files are read, with range requests:

    python validation/catalog/validate_catalog_tmi.py <data directory> 2013-06-01 --fit
    python validation/catalog/validate_catalog_tmi.py <data directory> 2013-06-01
    python validation/catalog/validate_catalog_tmi.py <data directory> 2013-07-10

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import json
import pathlib
from datetime import datetime, timedelta, timezone

import h5py
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.spatial import cKDTree

import tatc
from tatc.analysis import compute_sample_sizes
from tatc.catalog import get_entry
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import CircularOrbit, GeneralPerturbationsOrbit
from gesdisc_scanner import Product, Reference, load_day, observe_orbit
from validation_io import (
    TO_GEODETIC,
    RemoteFile,
    PUBLISHED_EPOCH,
    build_grid,
    compute_spacing_ratios,
    ecef,
    get_holdout_gap,
    normalize,
    read_element_set_file,
    recompute_passes,
    score_passes,
    score_published,
    score_swath_area,
    to_observations,
    up,
    update_fixture,
    update_record,
    RECORDS,
    FIXTURES,
    source_of,
)
from validate_catalog_trmm_pr import ELEMENT_SET_FILE, FIT_DAY, NORAD_ID

RECORD_PATH = RECORDS / "tmi.json"
FIXTURE_PATH = FIXTURES / "tmi.npz"
PLATFORM = "trmm"
CONFIGURATION = "feeds"
SAMPLE_CONFIGURATION = "sample_size"
GROUPS = ("S1", "S2", "S3")


def read_tmi(_product, url: str, session, day: datetime, _directory) -> dict:
    """
    Reads the times of a granule's scans, the spacecraft's positions, and the
    positions of the samples of each group, with range requests.
    """
    with h5py.File(RemoteFile(url, session)) as f:
        g = f["S1"]
        date = pd.to_datetime(
            pd.DataFrame(
                {
                    "year": g["ScanTime/Year"][:],
                    "month": g["ScanTime/Month"][:],
                    "day": g["ScanTime/DayOfMonth"][:],
                }
            ),
            utc=True,
        )
        part = {
            "time": (date - pd.Timestamp(day)).dt.total_seconds().to_numpy()
            + g["ScanTime/SecondOfDay"][:].astype(float),
            "sc_lat": g["SCstatus/SClatitude"][:],
            "sc_lon": g["SCstatus/SClongitude"][:],
            "sc_alt": g["SCstatus/SCaltitude"][:],
            "orientation": g["SCstatus/SCorientation"][:],
        }
        for group in GROUPS:
            part[f"lat_{group}"] = f[group]["Latitude"][:]
            part[f"lon_{group}"] = f[group]["Longitude"][:]
    return part


PRODUCT = Product(short_name="GPM_1CTRMMTMI", version="08", read=read_tmi)
# the scan period is measured on the day of the fit
SCAN_PERIOD = 1.899
TIME_STEPS = (SCAN_PERIOD, 10, 30)
POINT_SPACING = 300e3
PASS_RADIUS = 60  # km
FIXTURE_PASSES = 80
# published footprints: instrument, altitude (km), inclination (deg), scan angle
# (deg), and sizes along and across track (km), from the TC PRIMED documentation,
# after TRMM's orbit boost of August 2001, at the center of the scan
PUBLISHED_FOOTPRINTS = (
    ("TMI 10.65 GHz", 403, 35, 0, 68.1, 41.3),
    ("TMI 19.35 GHz", 403, 35, 0, 34.6, 20.9),
    ("TMI 21.3 GHz", 403, 35, 0, 31.1, 18.9),
    ("TMI 37.0 GHz", 403, 35, 0, 18.3, 11.1),
    ("TMI 85.5 GHz", 403, 35, 0, 7.8, 4.8),
)
# an instrument of each group, whose spacing is measured
SPACING_INSTRUMENTS = {
    "S1": "TMI 10.65 GHz",
    "S2": "TMI 19.35 GHz",
    "S3": "TMI 85.5 GHz",
}
SPACING_SCAN_STRIDE = 200
FIXTURE_SPACING_SAMPLES = 100
GEOD = Geod(ellps="WGS84")


def load_references(day: datetime, directory: pathlib.Path) -> dict[str, Reference]:
    """Loads a day, as a reference for each group of samples."""
    data = load_day(PRODUCT, day, directory)
    orientation = data["orientation"]
    print(
        "orientation of the scans: "
        + ", ".join(
            f"{value}: {count}"
            for value, count in zip(*np.unique(orientation, return_counts=True))
        ),
        flush=True,
    )
    references = {}
    for group in GROUPS:
        part = {
            "time": data["time"],
            "sc_lat": data["sc_lat"],
            "sc_lon": data["sc_lon"],
            "sc_alt": data["sc_alt"],
            "lat": data[f"lat_{group}"],
            "lon": data[f"lon_{group}"],
        }
        references[group] = Reference(day, part, SCAN_PERIOD)
    return references


def fit(reference: Reference, group: str) -> None:
    """
    Prints the fitted scan of a group: the median angle of the samples from nadir
    (the cone angle), and the azimuths of the first and last samples from the
    velocity, in the inertial and Earth-fixed frames, with their center by latitude,
    and the spacing of samples along the scan.
    """
    lon, lat, _ = TO_GEODETIC.transform(*reference.position.T)
    nadir = -up(lon, lat)
    sight = normalize(
        ecef(reference.lon, reference.lat) - reference.position[:, np.newaxis]
    )
    cone = np.degrees(np.arccos(np.sum(sight * nadir[:, np.newaxis], axis=-1)))
    print(
        f"{group} ({reference.beams} samples): cone angle median "
        f"{np.median(cone):.4f} deg, p5 {np.percentile(cone, 5):.4f}, p95 "
        f"{np.percentile(cone, 95):.4f}"
    )
    for frame in ("inertial", "earth_fixed"):
        motion = reference.velocity(frame)
        forward = normalize(motion - np.sum(motion * nadir, axis=-1)[:, None] * nadir)
        # to the left of the motion
        left = normalize(np.cross(forward, nadir))
        azimuth = np.degrees(
            np.arctan2(
                np.sum(sight * left[:, None], axis=-1),
                np.sum(sight * forward[:, None], axis=-1),
            )
        )
        first, last = np.median(azimuth[:, 0]), np.median(azimuth[:, -1])
        center = (azimuth[:, 0] + azimuth[:, -1]) / 2
        bands = ", ".join(
            f"{low:+d}..{low + 15:+d}: "
            f"{np.median(center[(np.asarray(lat) >= low) & (np.asarray(lat) < low + 15)]):+.3f}"
            for low in range(-45, 45, 15)
            if np.any((np.asarray(lat) >= low) & (np.asarray(lat) < low + 15))
        )
        step = np.median(np.abs(np.diff(azimuth, axis=1)))
        print(
            f"  {frame}: first sample {first:+.4f} deg, last {last:+.4f} deg; center "
            f"{np.median(center):+.4f} (p5 {np.percentile(center, 5):+.4f}, p95 "
            f"{np.percentile(center, 95):+.4f}), half width "
            f"{np.median(np.abs(azimuth[:, -1] - azimuth[:, 0]) / 2):.4f} deg; "
            f"azimuth step {step:.5f} deg; center by latitude: {bands}"
        )
    print(
        f"  scan period: median {np.median(np.diff(reference.seconds)):.6f} s; "
        f"altitude median {np.median(reference.altitude) / 1e3:.2f} km"
    )


def cone_crossing(reference: Reference, target, k: int, cone_angle: float) -> float:
    """Gets the time (s) a target crosses the cone near scan k, by bisection."""

    def angle(t: float) -> float:
        seconds = reference.seconds
        i = int(np.clip(np.searchsorted(seconds, t) - 1, 0, len(seconds) - 2))
        f = (t - seconds[i]) / (seconds[i + 1] - seconds[i])
        position = reference.position[i] + f * (
            reference.position[i + 1] - reference.position[i]
        )
        lon, lat, _ = TO_GEODETIC.transform(*position)
        return (
            np.degrees(np.arccos(np.dot(normalize(target - position), -up(lon, lat))))
            - cone_angle
        )

    lower, upper = reference.seconds[k] - 8, reference.seconds[k] + 8
    if angle(lower) * angle(upper) > 0:
        return np.nan
    for _ in range(30):
        middle = (lower + upper) / 2
        lower, upper = (middle, upper) if angle(middle) > 0 else (lower, middle)
    return (lower + upper) / 2


def observe_passes(reference: Reference, points, cone_angle: float) -> pd.DataFrame:
    """
    Finds the passes on which the samples observed points (in a cell of the lattice of
    samples), and their times (when they cross the cone).
    """
    point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
    flat = reference.unit.reshape(-1, 3)
    scan_of = np.repeat(np.arange(len(reference.seconds)), reference.beams)
    rows = []
    for p, hits in enumerate(
        cKDTree(flat).query_ball_point(
            point_unit, r=2 * np.sin(PASS_RADIUS / 6371 / 2), workers=-1
        )
    ):
        if not hits:
            continue
        hits = np.array(hits)
        hits = hits[np.argsort(reference.seconds[scan_of[hits]])]
        times = reference.seconds[scan_of[hits]]
        for group in np.split(hits, np.nonzero(np.diff(times) > 600)[0] + 1):
            k = int(
                scan_of[
                    group[
                        np.argmin(np.linalg.norm(flat[group] - point_unit[p], axis=-1))
                    ]
                ]
            )
            first, last = max(k - 2, 0), min(k + 2, len(reference.seconds) - 1)
            if not reference.in_cells(point_unit[[p]], first, last)[0]:
                continue
            crossing = cone_crossing(
                reference,
                ecef(points.geometry.x.iloc[p], points.geometry.y.iloc[p]),
                k,
                cone_angle,
            )
            if np.isfinite(crossing):
                rows.append({"point_id": p, "seconds": float(crossing)})
    observed = pd.DataFrame(rows).astype({"point_id": "int64"})
    observed["time"] = pd.Timestamp(reference.day) + pd.to_timedelta(
        observed.seconds, unit="s"
    )
    return observed


def score_swath_width(reference: Reference, instrument) -> float:
    """
    Gets the median absolute difference (km) of the swath width from the distance
    along the surface between the first and last samples of scans.
    """
    scans = np.arange(0, len(reference.seconds), 20)
    distance = GEOD.inv(
        reference.lon[scans, 0],
        reference.lat[scans, 0],
        reference.lon[scans, -1],
        reference.lat[scans, -1],
    )[2]
    width = np.vectorize(instrument.get_swath_width)(reference.altitude[scans])
    return float(np.median(np.abs(width - distance))) / 1e3


def measure_spacing(reference: Reference) -> np.ndarray:
    """
    Measures the spacing of a group's samples as `ValidateSampleSizes` does, from
    half the distance between the samples on either side. Returns rows of latitude,
    longitude, spacing between scans (km), and spacing along the arc (km), with the
    seconds of each row.
    """
    scans = reference.contiguous_scans()[::SPACING_SCAN_STRIDE]
    position = ecef(reference.lon, reference.lat)
    arc = np.linalg.norm(position[scans, 2:] - position[scans, :-2], axis=-1) / 2e3
    between = (
        np.linalg.norm(position[scans + 1, 1:-1] - position[scans - 1, 1:-1], axis=-1)
        / 2e3
    )
    rows = np.stack(
        [
            reference.lat[scans, 1:-1].ravel(),
            reference.lon[scans, 1:-1].ravel(),
            between.ravel(),
            arc.ravel(),
        ],
        axis=-1,
    )
    seconds = np.repeat(reference.seconds[scans], reference.beams - 2)
    return rows, seconds


def compute_published_errors(entry) -> np.ndarray:
    """
    Gets the relative errors of sample sizes, at the center of the scan, against the
    published footprints. Gives two errors for each footprint.
    """
    errors = []
    for name, altitude, inclination, _, along, cross in PUBLISHED_FOOTPRINTS:
        orbit = CircularOrbit(
            altitude=altitude * 1e3, inclination=inclination, epoch=PUBLISHED_EPOCH
        )
        satellite = entry.to_satellite(
            orbit, SAMPLE_CONFIGURATION, time_step=SCAN_PERIOD, platform=PLATFORM
        )
        names = [instrument.name for instrument in satellite.instruments]
        time = PUBLISHED_EPOCH + timedelta(minutes=10)
        target = satellite.get_instrument(names.index(name)).compute_footprint_center(
            orbit.to_gp_orbit().get_orbit_track(time)
        )
        sizes = compute_sample_sizes(
            to_observations(
                satellite,
                [name],
                [time],
                [target.longitude.degrees],
                [target.latitude.degrees],
            ),
            satellite,
        )
        errors += [
            sizes.sample_along.iloc[0] / 1e3 / along - 1,
            sizes.sample_cross.iloc[0] / 1e3 / cross - 1,
        ]
        print(f"{name}: {errors[-2]:+.2%} along, {errors[-1]:+.2%} across")
    return np.array(errors)


def validate_sample_sizes(day, tle, entry, orbit, references, write: bool) -> None:
    """Scores the sample size configuration, and writes its results if requested."""
    within, largest = score_published(compute_published_errors(entry))
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    metrics = [
        {"name": "sample_size_within_3pct", "value": within},
        {"name": "sample_size_largest_error", "value": largest},
    ]
    satellite = entry.to_satellite(
        orbit, SAMPLE_CONFIGURATION, time_step=SCAN_PERIOD, platform=PLATFORM
    )
    rng = np.random.default_rng(2)
    fixture_metrics = {"sample_size_within_3pct": within}
    reduced = {}
    for group, name in SPACING_INSTRUMENTS.items():
        rows, seconds = measure_spacing(references[group])
        times = pd.DatetimeIndex(pd.Timestamp(day) + pd.to_timedelta(seconds, "s"))
        # the spacing between scans is along track, and along the arc across track
        between, arc = compute_spacing_ratios(
            satellite, [name] * len(rows), times, rows
        )
        print(
            f"{group} spacing, predicted / measured: {arc:.4f} arc, {between:.4f} scans"
        )
        metrics += [
            {"name": "spacing_arc_ratio_median", "value": arc, "scope": group},
            {"name": "spacing_scan_ratio_median", "value": between, "scope": group},
        ]
        picks = rng.choice(len(rows), FIXTURE_SPACING_SAMPLES, replace=False)
        between, arc = compute_spacing_ratios(
            satellite, [name] * len(picks), times[picks], rows[picks]
        )
        fixture_metrics[f"spacing_arc_ratio_median_{group}"] = arc
        fixture_metrics[f"spacing_scan_ratio_median_{group}"] = between
        reduced[f"spacing_{group}"] = rows[picks]
        reduced[f"spacing_seconds_{group}"] = seconds[picks]
    if not write:
        return
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), **reduced
    )
    update_record(
        RECORD_PATH,
        "tmi",
        SAMPLE_CONFIGURATION,
        {
            "digest": entry.configurations[SAMPLE_CONFIGURATION].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": source_of(__file__),
            "reference": "Footprints of the TC PRIMED documentation; spacing of "
            "GPM_1CTRMMTMI version 08",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": metrics,
            "fixture_metrics": fixture_metrics,
        },
    )


def validate_swath(day, tle, entry, orbit, references, write: bool) -> None:
    """
    Scores each feed of the swath configuration on the samples of its group, and
    writes the results if requested, with the fixture of the first feed's passes.
    """

    def satellite_for(time_step: float):
        return entry.to_satellite(
            orbit, CONFIGURATION, time_step=time_step, platform=PLATFORM
        )

    satellite = satellite_for(SCAN_PERIOD)
    points = generate_points_fibonacci_lattice(POINT_SPACING)
    grid = build_grid()
    metrics, first_observed = [], None
    for index, group in enumerate(GROUPS):
        reference = references[group]
        cone_angle = (
            entry.configurations[CONFIGURATION]
            .instruments[index]
            .parameters["cone_angle"]
        )
        observed = observe_passes(reference, points, cone_angle)
        start = reference.time[0].to_pydatetime()
        end = reference.time[-1].to_pydatetime() + timedelta(seconds=SCAN_PERIOD)
        results = score_passes(
            points, satellite, index, observed, start, end, reference.gaps
        )
        results["swath_width_median"] = score_swath_width(
            reference, satellite.instruments[index]
        )
        first, last, compared, covered = observe_orbit(
            reference, orbit, grid, 2 * PASS_RADIUS
        )
        results["footprint_iou_median"] = score_swath_area(
            satellite_for, index, first, last, TIME_STEPS, grid, compared, covered
        )
        print(group, json.dumps(results), flush=True)
        metrics += [
            {"name": name, "value": float(results[name]), "scope": group}
            for name in (
                "pass_recall",
                "pass_precision",
                "epoch_p95",
                "swath_width_median",
                "footprint_iou_median",
            )
        ]
        if first_observed is None:
            first_observed = observed
    if not write:
        return
    rng = np.random.default_rng(1)
    picks = first_observed.iloc[
        rng.choice(len(first_observed), FIXTURE_PASSES, replace=False)
    ]
    passes = np.stack(
        [
            points.geometry.y.to_numpy()[picks.point_id],
            points.geometry.x.to_numpy()[picks.point_id],
            picks.seconds.to_numpy(),
        ],
        axis=-1,
    )
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=passes
    )
    update_record(
        RECORD_PATH,
        "tmi",
        CONFIGURATION,
        {
            "digest": entry.configurations[CONFIGURATION].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": source_of(__file__),
            "reference": "GPM_1CTRMMTMI version 08 (positions of the S1, S2, and S3 "
            "samples)",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": metrics,
            "fixture_metrics": recompute_passes(satellite, 0, day, passes),
        },
    )


def main() -> None:
    """Validates the entry on a day of data, or fits its scan."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--fit", action="store_true")
    parser.add_argument("--sample-sizes-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    references = load_references(day, arguments.directory)
    reference = references["S1"]
    print(
        f"{arguments.day}: {len(reference.time)} scans, {len(reference.gaps)} gaps "
        f"over a minute, median altitude {np.median(reference.altitude) / 1e3:.1f} km",
        flush=True,
    )
    if arguments.fit:
        for group in GROUPS:
            fit(references[group], group)
        return
    tle = read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)
    print(f"element set: {tle[0][:32]}", flush=True)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("tmi")
    validate_sample_sizes(day, tle, entry, orbit, references, arguments.write)
    if not arguments.sample_sizes_only:
        validate_swath(day, tle, entry, orbit, references, arguments.write)
    if arguments.write:
        print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()
