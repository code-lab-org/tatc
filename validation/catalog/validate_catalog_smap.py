"""
Validates the SMAP radiometer catalog entry against SMAP's Level 1B time-ordered
brightness temperatures (`SPL1BTB`), using the exact values of the entry through
public TAT-C methods. The radiometer scans a full circle, so it observes each point
twice on a pass: once looking forward, and once looking aft. Each look is tested as a
pass of a sector scanner is for AMSR2 and GMI:

* The forward and aft halves of each scan are resampled at fixed azimuths from the
  velocity, so that consecutive scans form a lattice of footprint centers. A point is
  observed by a look if it lies in a cell of the lattice, at the time it crosses the
  cone, from the reported spacecraft positions.
* The swath width is compared with the distance between the footprints at azimuths of
  90 deg to either side.
* The swath area over an orbit is compared with the union of the cells of both halves.

The cone angle is derived from the incidence angles of a day of data, so running on
a day at least 30 days from it tests it on data held out from it.

Needs the `validation` extra and a NASA Earthdata login:

    python validation/catalog/validate_catalog_smap.py <data directory> 2026-10-24 --fit
    python validation/catalog/validate_catalog_smap.py <data directory> 2026-10-24
    python validation/catalog/validate_catalog_smap.py <data directory> 2026-09-20

Without a day, it takes that of the entry's record (or, with `--fit`, of the fit).
Element sets for days not registered in `elements/element_sets.yml` are taken from
the file of SMAP's element sets in `elements/` (or one given with `--elements`).
`--download-only` downloads and reduces the day's data, without testing.

It also validates the entry's sample size configuration (`footprints`): sizes at
the scan's center are compared with published footprints, and the spacing of the
samples nearest it with that measured. To run only these tests, use
`--configuration footprints`.

With `--fit`, it derives the cone angle instead, as `ValidateGridSMAP` did, from the
gridded brightness temperatures of the day (`SPL1CTB`, the half orbits from two hours
before the day to two hours after it): the angle from nadir of a boresight at the
median incidence angle of the forward looks, at the mean altitude of the day's element
set, on a spherical Earth. It prints it, without writing.

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pathlib
import re
from datetime import datetime, timedelta

import earthaccess
import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from tatc.catalog.metrics import compute_pass_scores, summarize
from tatc.constants import timescale
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit
from conical import (
    GEOD,
    compute_arc_spacing_ratios,
    compute_center_errors,
    cone_crossing,
    orbit_period,
    pick_cone_passes,
    score_swath_steps,
)
from validation_io import (
    TO_GEODETIC,
    build_grid,
    ecef,
    in_lattice_cells,
    login_earthdata,
    normalize,
    parse_run,
    predict_passes,
    recompute_passes,
    score_published,
    up,
    write_results,
)

CONFIGURATION = "looks"
SCAN_PERIOD = 4.11
AZIMUTHS = np.arange(-87, 87.1, 1.5)  # deg from the velocity, in each half
TIME_STEPS = (SCAN_PERIOD, 10, 30)
POINT_SPACING = 300e3
MATCH_TOLERANCE = pd.Timedelta(seconds=60)
ORBIT_MARGIN = 120  # s
FIXTURE_PASSES = 120
SAMPLE_CONFIGURATION = "footprints"
# published footprints: name, altitude (km), and sizes along and across track (km) at
# the scan's center (forward), where the plane of incidence, which holds the major
# axis, is along track: the instantaneous field of view of Piepmeier et al. (2017),
# and the 3 dB field of view of the JPL SMAP specifications
PUBLISHED_FOOTPRINTS = (
    ("SMAP radiometer", 685, 47.0, 36.0),
    ("SMAP radiometer", 685, 49.0, 38.0),
)
PUBLISHED_INCLINATION = 98.1
SPACING_SCAN_STRIDE = 25
FIXTURE_SPACING_SAMPLES = 150


# the mean radius (m) of the spherical Earth on which the cone angle is derived
FIT_EARTH_RADIUS = 6371.0088e3


def to_seconds(j2000: np.ndarray, day: datetime) -> np.ndarray:
    """Converts seconds since J2000 (noon TT on 1 January 2000) to seconds of a day."""
    result = np.full(j2000.shape, np.nan)
    valid = j2000 > 0
    utc = pd.to_datetime(
        timescale.tt_jd(2451545.0 + j2000[valid] / 86400).utc_datetime(), utc=True
    )
    result[valid] = (utc - pd.Timestamp(day)).total_seconds()
    return result


def load_day(day: datetime, directory: pathlib.Path) -> dict[str, np.ndarray]:
    """Downloads and reduces the half orbits that overlap a day, or loads them."""
    path = directory / f"smap_l1b_{day:%Y%m%d}.npz"
    if not path.exists():
        login_earthdata()
        granules = earthaccess.search_data(
            short_name="SPL1BTB",
            temporal=(
                f"{day - timedelta(hours=1):%Y-%m-%dT%H:%M:%S}",
                f"{day + timedelta(hours=25):%Y-%m-%dT%H:%M:%S}",
            ),
        )
        # the latest version of each half orbit
        latest = {}
        for granule in sorted(granules, key=lambda g: g["umm"]["GranuleUR"]):
            name = granule["umm"]["GranuleUR"]
            latest[re.search(r"_TB_(\d+_[AD])_", name).group(1)] = granule
        files = earthaccess.download(
            list(latest.values()), str(directory / f"raw_{day:%Y%m%d}")
        )
        parts = []
        for file in sorted(files):
            with h5py.File(file) as f:
                tb, sc = f["Brightness_Temperature"], f["Spacecraft_Data"]
                parts.append(
                    {
                        "time": to_seconds(tb["tb_time_seconds"][:], day),
                        "lat": tb["tb_lat"][:],
                        "lon": tb["tb_lon"][:],
                        "sc_time": to_seconds(sc["antenna_scan_time"][:], day),
                        "sc_lat": sc["sc_nadir_lat"][:],
                        "sc_lon": sc["sc_nadir_lon"][:],
                        "sc_alt": sc["sc_geodetic_alt_ellipsoid"][:],
                    }
                )
            pathlib.Path(file).unlink()
        columns = max(p["lat"].shape[1] for p in parts)
        for part in parts:
            pad = columns - part["lat"].shape[1]
            for key in ("time", "lat", "lon"):
                part[key] = np.pad(
                    part[key], ((0, 0), (0, pad)), constant_values=-9999.0
                )
        data = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        _, first = np.unique(data["sc_time"], return_index=True)
        keep = first[(data["sc_time"][first] >= 0) & (data["sc_time"][first] < 86400)]
        np.savez_compressed(path, **{k: v[keep] for k, v in data.items()})
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


class Reference:
    """A day of SMAP scans, resampled into forward and aft lattices."""

    def __init__(self, day: datetime, data: dict[str, np.ndarray]):
        self.day = day
        # scans in time order, without those repeated (within 1 s) by overlapping files
        order = np.argsort(data["sc_time"])
        keep = order[np.r_[True, np.diff(data["sc_time"][order]) > 1]]
        data = {key: value[keep] for key, value in data.items()}
        self.time = data["sc_time"]
        self.position = ecef(
            data["sc_lon"].astype(float),
            data["sc_lat"].astype(float),
            data["sc_alt"].astype(float),
        )
        velocity = np.gradient(self.position, self.time, axis=0)
        nadir = -up(data["sc_lon"].astype(float), data["sc_lat"].astype(float))
        forward = normalize(
            velocity - np.sum(velocity * nadir, axis=-1, keepdims=True) * nadir
        )
        left = np.cross(forward, nadir)
        valid = (np.abs(data["lat"]) <= 90) & (data["time"] > -1e4)
        foot = ecef(
            np.where(valid, data["lon"], 0).astype(float),
            np.where(valid, data["lat"], 0).astype(float),
        )
        horizontal = foot - self.position[:, None]
        azimuth = np.degrees(
            np.arctan2(
                np.sum(horizontal * left[:, None], axis=-1),
                np.sum(horizontal * forward[:, None], axis=-1),
            )
        )
        # the footprints as sampled, for their spacing
        self.foot, self.valid, self.azimuth = foot, valid, azimuth
        self.foot_lat, self.foot_lon = data["lat"], data["lon"]
        self.foot_time = data["time"]
        # Earth-fixed positions of the resampled footprints of each half, and their
        # geocentric directions, the lattice of cells
        self.points, self.lattice = {}, {}
        for look, offset in (("fore", 0), ("aft", 180)):
            relative = (azimuth - offset + 180) % 360 - 180
            cells = np.full((len(self.time), len(AZIMUTHS), 3), np.nan)
            for k in range(len(self.time)):
                inside = valid[k] & (np.abs(relative[k]) < 90)
                if inside.sum() < 10:
                    continue
                order = np.argsort(relative[k][inside])
                angles = relative[k][inside][order]
                if angles[0] > AZIMUTHS[0] or angles[-1] < AZIMUTHS[-1]:
                    continue
                if np.max(np.diff(angles)) > 4:
                    continue
                position = foot[k][inside][order]
                cells[k] = np.stack(
                    [np.interp(AZIMUTHS, angles, position[:, c]) for c in range(3)],
                    axis=-1,
                )
            self.points[look] = cells
            self.lattice[look] = normalize(cells)
        complete = np.all(
            [np.all(np.isfinite(c[..., 0]), axis=1) for c in self.lattice.values()],
            axis=0,
        )
        self.contiguous = np.r_[np.diff(self.time) < SCAN_PERIOD + 1, False] & complete
        self.contiguous[:-1] &= complete[1:]
        self.gaps = [
            (self.time[i], self.time[i + 1])
            for i in np.nonzero(np.diff(self.time) > 60)[0]
        ]
        if self.time[0] > 60:
            self.gaps.insert(0, (0, self.time[0]))
        if self.time[-1] < 86400 - 60:
            self.gaps.append((self.time[-1], 86400))

    def in_gap(self, times: pd.Series, margin: float = 180) -> np.ndarray:
        """Gets whether times are within a margin of gaps of the data."""
        seconds = (times - pd.Timestamp(self.day)).dt.total_seconds().to_numpy()
        result = np.zeros(len(seconds), dtype=bool)
        for a, b in self.gaps:
            result |= (seconds > a - margin) & (seconds < b + margin)
        return result

    def in_cell(self, cells, point, k: int, j: int) -> bool:
        """Gets whether a unit vector lies in the cell of a lattice at (k, j)."""
        if not (0 <= k < len(self.time) - 1 and 0 <= j < len(AZIMUTHS) - 1):
            return False
        if not self.contiguous[k]:
            return False
        corners = [cells[k, j], cells[k, j + 1], cells[k + 1, j + 1], cells[k + 1, j]]
        signs = [
            np.dot(np.cross(a, b), point)
            for a, b in zip(corners, corners[1:] + corners[:1])
        ]
        return bool(np.all(np.array(signs) >= 0) or np.all(np.array(signs) <= 0))

    def observe_looks(self, points, cone_angle: float) -> pd.DataFrame:
        """
        Finds the looks of each half at points, and their times: for each pass near a
        point (footprints within 40 km, separated from others by 10 minutes), whether
        the point lies in a cell next to one of its eight nearest footprints.
        """
        lon, lat = points.geometry.x.to_numpy(), points.geometry.y.to_numpy()
        unit = normalize(ecef(lon, lat))
        rows = []
        columns = len(AZIMUTHS)
        for look, cells in self.lattice.items():
            filled = np.where(np.isfinite(cells), cells, 0)
            flat = filled.reshape(-1, 3)
            neighbors = cKDTree(flat).query_ball_point(
                unit, r=2 * np.sin(40 / 6371 / 2), workers=-1
            )
            for p, hits in enumerate(neighbors):
                if not hits:
                    continue
                hits = np.array(hits)
                hits = hits[np.argsort(self.time[hits // columns])]
                times = self.time[hits // columns]
                for group in np.split(hits, np.nonzero(np.diff(times) > 600)[0] + 1):
                    nearest = group[
                        np.argsort(np.linalg.norm(flat[group] - unit[p], axis=-1))[:8]
                    ]
                    if not any(
                        self.in_cell(
                            filled, unit[p], n // columns + dk, n % columns + dj
                        )
                        for n in nearest
                        for dk in (-1, 0)
                        for dj in (-1, 0)
                    ):
                        continue
                    k = int(nearest[0] // columns)
                    crossing = cone_crossing(
                        self.time,
                        self.position,
                        ecef(lon[p], lat[p]),
                        k,
                        cone_angle,
                        8,
                        signed=True,
                    )
                    if np.isfinite(crossing):
                        rows.append(
                            {"point_id": int(p), "look": look, "seconds": crossing}
                        )
        observed = pd.DataFrame(rows)
        observed["time"] = pd.Timestamp(self.day) + pd.to_timedelta(
            observed.seconds, unit="s"
        )
        # a point in the overlap of consecutive passes counts once per look and pass
        observed = observed.sort_values("time")
        repeated = observed.groupby(["point_id", "look"]).time.diff() < pd.Timedelta(
            seconds=20
        )
        return observed[~repeated].reset_index(drop=True)


def score_looks(reference, points, satellite, observed) -> dict[str, float]:
    """Compares predicted and observed looks of points."""
    start, end = (
        reference.day + timedelta(seconds=float(s))
        for s in (reference.time[0], reference.time[-1])
    )
    predicted = predict_passes(points, satellite, 0, start, end, [])
    predicted = predicted[~reference.in_gap(predicted.epoch)]
    observed = observed[~reference.in_gap(observed.time)]
    merged = pd.merge_asof(
        observed.assign(time=observed.time.dt.as_unit("ns")).sort_values("time"),
        predicted[["point_id", "epoch"]]
        .assign(epoch=predicted.epoch.dt.as_unit("ns"))
        .sort_values("epoch"),
        left_on="time",
        right_on="epoch",
        by="point_id",
        direction="nearest",
        tolerance=MATCH_TOLERANCE,
    )
    hit = merged.epoch.notna()
    differences = (merged.epoch - merged.time).dt.total_seconds()[hit]
    results = {
        **compute_pass_scores(int(hit.sum()), len(observed), len(predicted)),
        "epoch_p95": summarize(differences)["p95"],
        "epoch_median": summarize(differences)["median"],
        "observed": len(observed),
        "predicted": len(predicted),
        "matched": int(hit.sum()),
    }
    for look in ("fore", "aft"):
        chosen = merged.look == look
        results[f"recall_{look}"] = 100 * float(hit[chosen].mean())
    return results


def score_swath_width(reference, instrument) -> float:
    """Gets the median absolute swath width difference (km) across the circle."""
    width, chords = [], []
    for k in np.nonzero(reference.contiguous)[0][::20]:
        # the footprints nearest 90 deg to the left and right, at the ends of the
        # forward half
        ends = reference.points["fore"][k, [-1, 0]]
        lon, lat, _ = TO_GEODETIC.transform(*ends.T)
        chords.append(GEOD.inv(lon[0], lat[0], lon[1], lat[1])[2])
        _, _, height = TO_GEODETIC.transform(*reference.position[k])
        width.append(instrument.get_swath_width(height))
    return float(np.median(np.abs(np.array(width) - np.array(chords)))) / 1e3


def score_swath_area(reference, satellite_for, orbit) -> float:
    """Gets the lowest IoU of the swath over an orbit without gaps, by time step."""
    period = orbit_period(orbit)
    start = 12 * 3600
    while True:
        scans = np.nonzero(
            (reference.time >= start) & (reference.time < start + period)
        )[0]
        # an orbit without gaps of a minute or more (single missing scans are cells
        # that the lattice leaves out)
        if len(scans) > 2 and np.all(np.diff(reference.time[scans]) < 60):
            break
        start += 600
        if start + period > reference.time[-1]:
            raise ValueError("no orbit without a gap after 12:00 UTC")
    lon, lat, weights = build_grid()
    grid = normalize(ecef(lon, lat))
    tree = cKDTree(
        np.concatenate([c[scans].reshape(-1, 3) for c in reference.lattice.values()])
    )
    distance, nearest = tree.query(grid, workers=-1)
    count = len(scans) * len(AZIMUTHS)
    nearest_time = reference.time[scans][(nearest % count) // len(AZIMUTHS)]
    first, last = reference.time[scans[0]], reference.time[scans[-1]]
    compared = (
        (distance < 2 * np.sin(300 / 6371 / 2))
        & (nearest_time > first + ORBIT_MARGIN)
        & (nearest_time < last - ORBIT_MARGIN)
    )
    observed = np.zeros(len(grid), dtype=bool)
    for cells in reference.lattice.values():
        part = cells[scans]
        observed[compared] |= in_lattice_cells(
            part, reference.contiguous[scans][:-1], grid[compared], neighbors=12
        )
    return score_swath_steps(
        satellite_for,
        0,
        reference.day,
        first,
        last,
        TIME_STEPS,
        (lon, lat, weights),
        compared,
        observed,
    )


def measure_spacing(reference: Reference) -> np.ndarray:
    """
    Measures the spacing of the samples nearest the scan's center (forward), as
    `ValidateSampleSizes` does: along the arc, from half the distance between the
    samples before and after, and between scans, from half the distance between the
    resampled footprints at the center of the scans before and after. Returns rows of
    latitude, longitude, seconds, spacing along the arc (km), and between scans (km).
    """
    center = int(np.argmin(np.abs(AZIMUTHS)))
    fore = reference.points["fore"]
    rows = []
    usable = np.nonzero(reference.contiguous[:-1] & reference.contiguous[1:])[0] + 1
    for k in usable[::SPACING_SCAN_STRIDE]:
        inside = np.nonzero(reference.valid[k] & (np.abs(reference.azimuth[k]) < 20))[0]
        inside = inside[np.argsort(reference.foot_time[k][inside])]
        if len(inside) < 3:
            continue
        i = 1 + int(np.argmin(np.abs(reference.azimuth[k][inside[1:-1]])))
        before, sample, after = inside[i - 1], inside[i], inside[i + 1]
        foot = reference.foot[k]
        rows.append(
            [
                reference.foot_lat[k][sample],
                reference.foot_lon[k][sample],
                reference.foot_time[k][sample],
                np.linalg.norm(foot[after] - foot[before]) / 2e3,
                np.linalg.norm(fore[k + 1, center] - fore[k - 1, center]) / 2e3,
            ]
        )
    return np.array(rows, dtype=float)


def validate_sample_sizes(run, orbit, tle, reference: Reference) -> None:
    """Scores the sample size configuration, and writes its results if the run writes."""
    entry = run.entry
    within, largest = score_published(
        compute_center_errors(
            entry,
            SAMPLE_CONFIGURATION,
            PUBLISHED_FOOTPRINTS,
            PUBLISHED_INCLINATION,
            time_step=SCAN_PERIOD,
            platform=run.platform,
        )
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    satellite = entry.to_satellite(
        orbit,
        configuration=SAMPLE_CONFIGURATION,
        time_step=SCAN_PERIOD,
        platform=run.platform,
    )
    name = satellite.instruments[0].name
    rows = measure_spacing(reference)
    arc, between = compute_arc_spacing_ratios(satellite, name, run.day, rows)
    print(
        f"spacing, predicted / measured: {arc:.4f} arc, {between:.4f} scans "
        f"({np.median(rows[:, 3]):.2f} and {np.median(rows[:, 4]):.2f} km measured)",
        flush=True,
    )
    if not run.write:
        return
    rng = np.random.default_rng(2)
    picks = rows[rng.choice(len(rows), FIXTURE_SPACING_SAMPLES, replace=False)]
    fixture_arc, fixture_between = compute_arc_spacing_ratios(
        satellite, name, run.day, picks
    )
    write_results(
        run,
        SAMPLE_CONFIGURATION,
        reference="Footprints of Piepmeier et al. (2017) and the JPL SMAP "
        "specifications; spacing of SPL1BTB",
        metrics=[
            {"name": "sample_size_within_3pct", "value": within},
            {"name": "sample_size_largest_error", "value": largest},
            {"name": "spacing_arc_ratio_median", "value": arc},
            {"name": "spacing_scan_ratio_median", "value": between},
        ],
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "footprints": np.array(
                [row[1:] for row in PUBLISHED_FOOTPRINTS], dtype=float
            ),
            "spacing": picks,
        },
        fixture_metrics={
            "sample_size_within_3pct": within,
            "spacing_arc_ratio_median": fixture_arc,
            "spacing_scan_ratio_median": fixture_between,
        },
    )


def validate_looks(run, orbit, tle, reference: Reference) -> None:
    """Scores the looks configuration, and writes its results if the run writes."""
    entry = run.entry

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, time_step=time_step, platform=run.platform)

    satellite = satellite_for(SCAN_PERIOD)
    cone_angle = (
        entry.configurations[CONFIGURATION].instruments[0].parameters["cone_angle"]
    )
    points = generate_points_fibonacci_lattice(POINT_SPACING)
    observed = reference.observe_looks(points, cone_angle)
    results = score_looks(reference, points, satellite, observed)
    print(json.dumps(results), flush=True)
    results["swath_width_median"] = score_swath_width(
        reference, satellite.instruments[0]
    )
    results["footprint_iou_median"] = score_swath_area(reference, satellite_for, orbit)
    print(json.dumps(results), flush=True)
    if not run.write:
        return
    passes = pick_cone_passes(
        points, observed[~reference.in_gap(observed.time)], FIXTURE_PASSES
    )
    write_results(
        run,
        CONFIGURATION,
        reference="SPL1BTB (footprint positions and times)",
        metrics=[
            {"name": name, "value": float(results[name])}
            for name in (
                "pass_recall",
                "pass_precision",
                "epoch_p95",
                "epoch_median",
                "swath_width_median",
                "footprint_iou_median",
            )
        ],
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "passes": passes,
        },
        fixture_metrics=recompute_passes(satellite, 0, run.day, passes),
    )


def load_fit_incidence(day: datetime, directory: pathlib.Path) -> np.ndarray:
    """
    Downloads the incidence angles (deg) of the forward looks of the gridded brightness
    temperatures of the half orbits around a day, or loads them from a cache.
    """
    path = directory / f"smap_l1c_incidence_{day:%Y%m%d}.npy"
    if not path.exists():
        login_earthdata()
        granules = earthaccess.search_data(
            short_name="SPL1CTB",
            temporal=(
                f"{day - timedelta(hours=2):%Y-%m-%dT%H:%M:%S}",
                f"{day + timedelta(hours=26):%Y-%m-%dT%H:%M:%S}",
            ),
        )
        # the latest version of each half orbit
        latest = {}
        for granule in sorted(granules, key=lambda g: g["umm"]["GranuleUR"]):
            name = granule["umm"]["GranuleUR"]
            latest[re.search(r"_TB_(\d+_[AD])_", name).group(1)] = granule
        files = earthaccess.download(
            list(latest.values()), str(directory / f"raw_l1c_{day:%Y%m%d}")
        )
        parts = []
        for file in sorted(files):
            with h5py.File(file) as f:
                parts.append(f["Global_Projection/cell_boresight_incidence_fore"][:])
            pathlib.Path(file).unlink()
        np.save(path, np.concatenate(parts))
    return np.load(path)


def fit(day: datetime, directory: pathlib.Path, tle: list[str]) -> None:
    """
    Prints the cone angle of a boresight at the median incidence angle of a day's
    forward looks, at the mean altitude of its element set, on a spherical Earth.
    """
    incidence = float(np.nanmedian(load_fit_incidence(day, directory)))
    altitude = GeneralPerturbationsOrbit.from_tle(tle).get_mean_altitude()
    cone = np.degrees(
        np.arcsin(
            FIT_EARTH_RADIUS
            * np.sin(np.radians(incidence))
            / (FIT_EARTH_RADIUS + altitude)
        )
    )
    print(
        f"median incidence angle {incidence:.2f} deg, mean altitude "
        f"{altitude / 1e3:.1f} km"
    )
    print(f"cone_angle: {cone:.3f}")


def main() -> None:
    """Validates the entry on a day of data, or fits its cone angle."""
    run = parse_run(
        "smap_radiometer",
        __doc__,
        add_arguments=lambda parser: parser.add_argument(
            "--download-only",
            action="store_true",
            help="download and reduce the day's data, without testing",
        ),
    )
    if run.fit:
        fit(run.day, run.directory, run.element_set())
        return
    data = load_day(run.day, run.directory)
    if run.download_only:
        return
    tle = run.element_set()
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    reference = Reference(run.day, data)
    print(
        f"{run.day:%Y-%m-%d}: {len(reference.time)} scans, gaps "
        + ", ".join(f"{a / 3600:.2f}-{b / 3600:.2f} h" for a, b in reference.gaps),
        flush=True,
    )
    if SAMPLE_CONFIGURATION in run.configurations:
        validate_sample_sizes(run, orbit, tle, reference)
    if CONFIGURATION in run.configurations:
        validate_looks(run, orbit, tle, reference)


if __name__ == "__main__":
    main()
