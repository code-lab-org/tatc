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

The cone angle was derived in `ValidateGridSMAP` from the incidence angles of a day of
data, so running on a day at least 30 days from it tests it on data held out from it.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_smap.py <data directory> 2026-10-24
    python docs/validation/validate_catalog_smap.py <data directory> 2026-09-20

With `--write`, it writes the entry's validation record and the reduced fixture from
which unit tests recompute the results.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import json
import pathlib
import re
from datetime import datetime, timedelta, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.spatial import cKDTree
from shapely import contains_xy, prepare

import tatc
from tatc.analysis import collect_observations, compute_ground_track
from tatc.catalog import get_entry
from tatc.catalog.metrics import compute_mask_iou, compute_pass_scores, summarize
from tatc.constants import timescale
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit, Point
from validation_io import (
    TO_GEODETIC,
    build_grid,
    ecef,
    get_holdout_gap,
    in_lattice_cells,
    normalize,
    read_element_set_file,
    recompute_passes,
    up,
    update_fixture,
    update_record,
)

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = (
    REPOSITORY / "src" / "tatc" / "resources" / "validation" / "smap_radiometer.json"
)
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "smap_radiometer.npz"
ELEMENT_SET_FILE = REPOSITORY / "docs" / "validation" / "smap_2026-09-20.3le"

FIT_DAY = datetime(2026, 10, 24, tzinfo=timezone.utc)
PLATFORM = "smap"
NORAD_ID = 40376
CONFIGURATION = "looks"
# the element set of ValidateGridSMAP, for the day of its cone angle
ELEMENT_SETS = {
    "2026-10-24": [
        "1 40376U 15003A   26298.67804942  .00000247  00000+0  56844-4 0  9999",
        "2 40376  98.1310 303.9506 0001331  85.8091 274.3263 14.63397245626701",
    ],
}
SCAN_PERIOD = 4.11
AZIMUTHS = np.arange(-87, 87.1, 1.5)  # deg from the velocity, in each half
TIME_STEPS = (SCAN_PERIOD, 10, 30)
POINT_SPACING = 300e3
MATCH_TOLERANCE = pd.Timedelta(seconds=60)
ORBIT_MARGIN = 120  # s
FIXTURE_PASSES = 120
GEOD = Geod(ellps="WGS84")


def get_element_set(day: datetime) -> list[str]:
    """Gets the element set of a day."""
    key = f"{day:%Y-%m-%d}"
    if key in ELEMENT_SETS:
        return ELEMENT_SETS[key]
    return read_element_set_file(ELEMENT_SET_FILE, NORAD_ID, day)


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
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
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

    def to_utc(self, seconds) -> list[datetime]:
        """Converts seconds from the start of the day to times."""
        return [self.day + timedelta(seconds=float(s)) for s in np.atleast_1d(seconds)]

    def in_gap(self, times: pd.Series, margin: float = 180) -> np.ndarray:
        """Gets whether times are within a margin of gaps of the data."""
        seconds = (times - pd.Timestamp(self.day)).dt.total_seconds().to_numpy()
        result = np.zeros(len(seconds), dtype=bool)
        for a, b in self.gaps:
            result |= (seconds > a - margin) & (seconds < b + margin)
        return result

    def cone_crossing(self, target, k: int, cone_angle: float) -> float:
        """Gets the time (s) a target crosses the cone near scan k, by bisection."""

        def angle(t: float) -> float:
            i = int(np.clip(np.searchsorted(self.time, t) - 1, 0, len(self.time) - 2))
            f = (t - self.time[i]) / (self.time[i + 1] - self.time[i])
            position = self.position[i] + f * (self.position[i + 1] - self.position[i])
            lon, lat, _ = TO_GEODETIC.transform(*position)
            return (
                np.degrees(
                    np.arccos(np.dot(normalize(target - position), -up(lon, lat)))
                )
                - cone_angle
            )

        lower, upper = self.time[k] - 8, self.time[k] + 8
        if angle(lower) * angle(upper) > 0:
            return np.nan
        for _ in range(30):
            middle = (lower + upper) / 2
            if (angle(middle) > 0) == (angle(lower) > 0):
                lower = middle
            else:
                upper = middle
        return (lower + upper) / 2

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
                    crossing = self.cone_crossing(ecef(lon[p], lat[p]), k, cone_angle)
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
    start, end = reference.to_utc([reference.time[0], reference.time[-1]])

    def predict(item):
        p, geometry = item
        return collect_observations(
            Point(latitude=geometry.y, longitude=geometry.x), satellite, start, end
        ).assign(point_id=p)

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        predicted = pd.concat(list(pool.map(predict, enumerate(points.geometry))))
    predicted = predicted.reset_index(drop=True)
    predicted["point_id"] = predicted.point_id.astype("int64")
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
    period = 86400 / (orbit.elements[0].mean_motion * 86400 / 360)
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
    scores = []
    for time_step in TIME_STEPS:
        times = reference.to_utc(
            first + time_step * (np.arange(int((last - first) // time_step)) + 0.5)
        )
        geometry = compute_ground_track(satellite_for(time_step), times).geometry.iloc[
            0
        ]
        prepare(geometry)
        predicted = np.zeros(len(grid), dtype=bool)
        predicted[compared] = contains_xy(geometry, lon[compared], lat[compared])
        scores.append(compute_mask_iou(observed, predicted, weights))
        print(f"swath IoU at {time_step:g} s: {scores[-1]:.4f}", flush=True)
    return min(scores)


def main() -> None:
    """Validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--download-only", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)
    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    data = load_day(day, arguments.directory)
    if arguments.download_only:
        return
    tle = get_element_set(day)
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("smap_radiometer")
    reference = Reference(day, data)
    print(
        f"{arguments.day}: {len(reference.time)} scans, gaps "
        + ", ".join(f"{a / 3600:.2f}-{b / 3600:.2f} h" for a, b in reference.gaps),
        flush=True,
    )

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, time_step=time_step, platform=PLATFORM)

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
    if not arguments.write:
        return
    rng = np.random.default_rng(1)
    kept = observed[~reference.in_gap(observed.time)]
    picks = kept.iloc[rng.choice(len(kept), FIXTURE_PASSES, replace=False)]
    fixture = np.column_stack(
        [
            points.geometry.y.to_numpy()[picks.point_id],
            points.geometry.x.to_numpy()[picks.point_id],
            picks.seconds.to_numpy(),
        ]
    )
    update_fixture(
        FIXTURE_PATH, tle=np.array(tle), day=np.array(day.isoformat()), passes=fixture
    )
    metrics = [
        {"name": name, "value": float(results[name])}
        for name in (
            "pass_recall",
            "pass_precision",
            "epoch_p95",
            "epoch_median",
            "swath_width_median",
            "footprint_iou_median",
        )
    ]
    update_record(
        RECORD_PATH,
        "smap_radiometer",
        CONFIGURATION,
        {
            "digest": entry.configurations[CONFIGURATION].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_smap.py",
            "reference": "SPL1BTB (footprint positions and times)",
            "period": f"{day:%Y-%m-%d}",
            "platform": PLATFORM,
            "holdout_gap_days": get_holdout_gap(day, FIT_DAY),
            "metrics": metrics,
            "fixture_metrics": recompute_passes(satellite, 0, day, fixture),
        },
    )
    print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()
