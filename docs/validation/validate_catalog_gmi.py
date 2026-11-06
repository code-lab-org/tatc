"""
Validates the GMI catalog entry against the GPM Level 1C GMI product, using the
exact values of the entry through public TAT-C methods. The tests follow
`ValidateConicalGMI`, whose day of data the entry's angles were fitted to:
running on that day checks this script against the notebook, while running on a
day at least 30 days from it tests the fitted angles on data held out from the
fit. The day must be one on which GPM flew forward, as it does for about 40 days
at a time, because the configuration represents forward flight. Element sets for
other days can be added to `ELEMENTS`.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_gmi.py <data directory> 2026-10-03

With `--write`, it writes the entry's validation record and the reduced fixture
from which unit tests recompute the results. A record of the day of the fit has
no holdout, which the entry's validation reports as a gap.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import argparse
import concurrent.futures
import io
import json
import pathlib
from datetime import datetime, timedelta, timezone

import earthaccess
import h5py
import numpy as np
import pandas as pd
import requests
from pyproj import Geod, Transformer
from scipy.spatial import cKDTree
from shapely import contains_xy, prepare

import tatc
from tatc.analysis import collect_observations, compute_ground_track
from tatc.catalog import STANDARD_VERSION, get_entry
from tatc.catalog.metrics import compute_mask_iou, compute_pass_scores, summarize
from tatc.constants import EARTH_ROTATION_RATE
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import GeneralPerturbationsOrbit, Point

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "gmi.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "gmi.npz"

FIT_DAY = datetime(2026, 10, 3, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs just after each day of data
ELEMENTS = {
    "2026-10-03": (
        "1 39574U 14009C   26277.54227627  .00004528  00000+0  10004-3 0  9997",
        "2 39574  64.9701 301.6724 0010793 268.6572  91.3332 15.45217995714704",
    ),
    "2026-11-05": (
        "1 39574U 14009C   26310.23667006  .00009405  00000+0  20531-3 0  9995",
        "2 39574  64.9705 192.1442 0011664 254.8479 105.1371 15.44674931719752",
    ),
}
FEEDS = ("S1", "S2")
SCAN_PERIOD = 1.875
TIME_STEPS = (SCAN_PERIOD, 10, 30)
GRID_STEP, ORBIT_MARGIN = 0.1, 60
FIXTURE_PASSES = 40

GEOD = Geod(ellps="WGS84")
TO_ECEF = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
TO_GEODETIC = Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)


def load_day(day: datetime, directory: pathlib.Path) -> dict[str, np.ndarray]:
    """Downloads and reduces all orbits of a day, or loads them from a cache."""
    cache = directory / f"gmi_l1c_{day:%Y%m%d}.npz"
    if not cache.exists():
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        granules = requests.get(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params={
                "short_name": "GPM_1CGPMGMI",
                "version": "08",
                "temporal": f"{day:%Y-%m-%d}T00:00:00Z,{day:%Y-%m-%d}T23:59:59Z",
                "page_size": 100,
            },
            timeout=120,
        ).json()["items"]
        parts = []
        for granule in granules:
            url = next(
                u["URL"]
                for u in granule["umm"]["RelatedUrls"]
                if u["URL"].startswith("https://data.gesdisc")
            )
            response = session.get(url, timeout=900)
            response.raise_for_status()
            with h5py.File(io.BytesIO(response.content)) as f:
                s1 = f["S1"]
                seconds = pd.to_datetime(
                    pd.DataFrame(
                        {
                            "year": s1["ScanTime/Year"][:],
                            "month": s1["ScanTime/Month"][:],
                            "day": s1["ScanTime/DayOfMonth"][:],
                            "hour": s1["ScanTime/Hour"][:],
                            "minute": s1["ScanTime/Minute"][:],
                            "second": s1["ScanTime/Second"][:],
                            "ms": s1["ScanTime/MilliSecond"][:],
                        }
                    ),
                    utc=True,
                ) - pd.Timestamp(day)
                parts.append(
                    {
                        "time": seconds.dt.total_seconds().to_numpy(),
                        "sc_lat": s1["SCstatus/SClatitude"][:],
                        "sc_lon": s1["SCstatus/SClongitude"][:],
                        "sc_alt": s1["SCstatus/SCaltitude"][:],
                        "orientation": s1["SCstatus/SCorientation"][:],
                        **{f"{k}_lat": f[f"{k}/Latitude"][:] for k in FEEDS},
                        **{f"{k}_lon": f[f"{k}/Longitude"][:] for k in FEEDS},
                    }
                )
        reduced = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        # keep each scan once (orbit files overlap), in time order, within the day
        _, first = np.unique(reduced["time"], return_index=True)
        keep = first[(reduced["time"][first] >= 0) & (reduced["time"][first] < 86400)]
        np.savez_compressed(cache, **{k: v[keep] for k, v in reduced.items()})
    return dict(np.load(cache))


def ecef(lon, lat, height=0.0) -> np.ndarray:
    """Gets the Earth-fixed position (m) of geodetic coordinates."""
    return np.stack(
        TO_ECEF.transform(lon, lat, np.broadcast_to(height, np.shape(lon))), axis=-1
    )


def up(lon, lat) -> np.ndarray:
    """Gets the geodetic up unit vector."""
    lon, lat = np.radians(lon), np.radians(lat)
    return np.stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1
    )


def normalize(x: np.ndarray) -> np.ndarray:
    """Gets unit vectors along the last axis."""
    return x / np.linalg.norm(x, axis=-1, keepdims=True)


def inertial_velocity(position: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    """Gets the inertial velocity, expressed in Earth-fixed coordinates."""
    return velocity + EARTH_ROTATION_RATE * np.stack(
        [-position[..., 1], position[..., 0], np.zeros(position.shape[:-1])], axis=-1
    )


class Reference:
    """A day of GMI scans, with what is derived from them for each feed."""

    def __init__(self, day: datetime, data: dict[str, np.ndarray]):
        self.day = day
        orientations = np.unique(data["orientation"])
        if not np.array_equal(orientations, [0]):
            raise ValueError(
                f"GPM did not fly forward all of {day:%Y-%m-%d}: {orientations}; "
                "the configuration represents forward flight"
            )
        valid = np.all(
            [np.all(np.abs(data[f"{k}_lat"]) <= 90, axis=1) for k in FEEDS], axis=0
        )
        self.time = data["time"][valid]
        self.altitude = data["sc_alt"][valid].astype(float) * 1e3
        self.position = ecef(
            data["sc_lon"][valid].astype(float),
            data["sc_lat"][valid].astype(float),
            self.altitude,
        )
        self.velocity = np.gradient(self.position, self.time, axis=0)
        self.contiguous = (
            np.r_[np.diff(self.time) < 3, False] & np.r_[False, np.diff(self.time) < 3]
        )
        self.foot = {
            k: (
                data[f"{k}_lon"][valid].astype(float),
                data[f"{k}_lat"][valid].astype(float),
            )
            for k in FEEDS
        }
        self.samples = self.foot["S1"][0].shape[1]
        self.gaps = [
            (self.time[i], self.time[i + 1])
            for i in np.nonzero(np.diff(self.time) > 60)[0]
        ]

    def to_utc(self, seconds) -> list[datetime]:
        """Converts seconds from the start of the day to times."""
        return [self.day + timedelta(seconds=float(s)) for s in np.atleast_1d(seconds)]

    def in_gap(self, times: pd.Series) -> np.ndarray:
        """Gets whether times are in gaps of the data."""
        start = pd.Timestamp(self.day)
        result = np.zeros(len(times), dtype=bool)
        for a, b in self.gaps:
            result |= (times >= start + pd.Timedelta(seconds=a)) & (
                times <= start + pd.Timedelta(seconds=b)
            )
        return result

    def in_cell(self, unit: np.ndarray, points, k, j) -> np.ndarray:
        """Gets whether points lie in the cells of footprint centers at (k, j)."""
        valid = (k >= 0) & (k < len(self.time) - 1) & (j >= 0) & (j < self.samples - 1)
        k = np.clip(k, 0, len(self.time) - 2)
        j = np.clip(j, 0, self.samples - 2)
        valid &= np.diff(self.time)[k] < 3
        corners = [unit[k, j], unit[k, j + 1], unit[k + 1, j + 1], unit[k + 1, j]]
        signs = np.stack(
            [
                np.sum(np.cross(p, q) * points, axis=-1)
                for p, q in zip(corners, corners[1:] + corners[:1])
            ]
        )
        return valid & (np.all(signs >= 0, axis=0) | np.all(signs <= 0, axis=0))

    def cone_crossing(self, target, k: int, cone_angle: float) -> float:
        """Gets the time (s) a target crosses a cone near scan k, by bisection."""

        def angle(t: float) -> float:
            i = int(np.clip(np.searchsorted(self.time, t) - 1, 0, len(self.time) - 2))
            f = (t - self.time[i]) / (self.time[i + 1] - self.time[i])
            position = self.position[i] + f * (self.position[i + 1] - self.position[i])
            lon, lat, _ = TO_GEODETIC.transform(*position)
            nadir = -up(lon, lat)
            return (
                np.degrees(np.arccos(np.dot(normalize(target - position), nadir)))
                - cone_angle
            )

        lower, upper = self.time[k] - 6, self.time[k] + 6
        if angle(lower) * angle(upper) > 0:
            return np.nan
        for _ in range(30):
            middle = (lower + upper) / 2
            lower, upper = (middle, upper) if angle(middle) > 0 else (lower, middle)
        return (lower + upper) / 2

    def observe_passes(self, points, feed: str, cone_angle: float) -> pd.DataFrame:
        """Finds the passes on which a feed observed points, and their times."""
        unit = up(*self.foot[feed])
        flat = unit.reshape(-1, 3)
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        tree = cKDTree(flat)
        rows = []
        for p, hits in enumerate(
            tree.query_ball_point(point_unit, r=2 * np.sin(40 / 6371 / 2), workers=-1)
        ):
            if not hits:
                continue
            hits = np.array(hits)
            scans = hits // self.samples
            order = np.argsort(self.time[scans])
            hits, scans = hits[order], scans[order]
            for group in np.split(
                np.arange(len(hits)), np.nonzero(np.diff(self.time[scans]) > 600)[0] + 1
            ):
                nearest = hits[group][
                    np.argsort(
                        np.linalg.norm(flat[hits[group]] - point_unit[p], axis=-1)
                    )[:8]
                ]
                if not any(
                    self.in_cell(
                        unit,
                        point_unit[[p]],
                        np.array([n // self.samples + dk]),
                        np.array([n % self.samples + dj]),
                    )[0]
                    for n in nearest
                    for dk in (-1, 0)
                    for dj in (-1, 0)
                ):
                    continue
                crossing = self.cone_crossing(
                    ecef(points.geometry.x.iloc[p], points.geometry.y.iloc[p]),
                    int(nearest[0] // self.samples),
                    cone_angle,
                )
                if np.isfinite(crossing):
                    rows.append({"point_id": p, "seconds": float(crossing)})
        observed = pd.DataFrame(rows).astype({"point_id": "int64"})
        observed["time"] = pd.Timestamp(self.day) + pd.to_timedelta(
            observed.seconds, unit="s"
        )
        return observed


def score_swath_width(reference: Reference, instrument, feed: str) -> float:
    """Gets the median absolute swath width difference (km) from the scanned arcs."""
    i = np.nonzero(reference.contiguous)[0][::20]
    lon, lat = reference.foot[feed][0][i], reference.foot[feed][1][i]
    chord = GEOD.inv(lon[:, 0], lat[:, 0], lon[:, -1], lat[:, -1])[2]
    width = np.vectorize(instrument.get_swath_width)(reference.altitude[i])
    return float(np.median(np.abs(width - chord))) / 1e3


def score_passes(reference, points, satellite, index, observed) -> dict[str, float]:
    """Compares a feed's predicted and observed passes of points."""
    start, end = reference.to_utc([reference.time[0], reference.time[-1]])

    def predict(item):
        p, geometry = item
        return collect_observations(
            Point(latitude=geometry.y, longitude=geometry.x),
            satellite,
            start,
            end,
            index,
        ).assign(point_id=p)

    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        predicted = pd.concat(list(pool.map(predict, enumerate(points.geometry))))
    predicted = predicted.reset_index(drop=True)
    predicted["point_id"] = predicted.point_id.astype("int64")
    predicted = predicted[~reference.in_gap(predicted.epoch)]
    merged = pd.merge_asof(
        observed.sort_values("time"),
        predicted[["point_id", "epoch"]].sort_values("epoch"),
        left_on="time",
        right_on="epoch",
        by="point_id",
        direction="nearest",
        tolerance=pd.Timedelta(minutes=10),
    )
    hit = merged.epoch.notna()
    differences = (merged.epoch - merged.time).dt.total_seconds()[hit]
    scores = compute_pass_scores(int(hit.sum()), len(observed), len(predicted))
    return {
        **scores,
        "epoch_p95": summarize(differences)["p95"],
        "epoch_median": summarize(differences)["median"],
        "observed": len(observed),
        "predicted": len(predicted),
        "matched": int(hit.sum()),
    }


def score_swath_area(reference, satellite_for, orbit, index, feed) -> float:
    """Gets the lowest IoU of a feed's swath over one orbit, at each time step."""
    first_scan = 12 * 3600
    period = 86400 / (orbit.elements[0].mean_motion * 86400 / 360)
    scans = np.nonzero(
        (reference.time >= first_scan) & (reference.time < first_scan + period)
    )[0]
    lon, lat = np.meshgrid(
        np.arange(-180 + GRID_STEP / 2, 180, GRID_STEP),
        np.arange(-90 + GRID_STEP / 2, 90, GRID_STEP),
    )
    weights = (np.cos(np.radians(lat)) * (GRID_STEP * 111.32) ** 2).ravel()
    grid = up(lon.ravel(), lat.ravel())
    unit = up(*reference.foot[feed])
    tree = cKDTree(unit[scans].reshape(-1, 3))
    distance, nearest = tree.query(grid, workers=-1)
    nearest_time = (
        reference.time[scans][nearest // reference.samples] - reference.time[scans[0]]
    )
    duration = reference.time[scans[-1]] - reference.time[scans[0]]
    compared = (
        (distance < 2 * np.sin(300 / 6371 / 2))
        & (nearest_time > ORBIT_MARGIN)
        & (nearest_time < duration - ORBIT_MARGIN)
    )
    observed = np.zeros(len(grid), dtype=bool)
    _, nearest_eight = tree.query(grid[compared], k=8, workers=-1)
    for column in range(8):
        k, j = np.divmod(nearest_eight[:, column], reference.samples)
        k = scans[k]
        inside = np.zeros(len(k), dtype=bool)
        for dk in (-1, 0):
            for dj in (-1, 0):
                inside |= reference.in_cell(unit, grid[compared], k + dk, j + dj)
        observed[np.flatnonzero(compared)] |= inside
    scores = []
    first, last = reference.time[scans[0]], reference.time[scans[-1]]
    for time_step in TIME_STEPS:
        times = reference.to_utc(
            first + time_step * (np.arange(int((last - first) // time_step)) + 0.5)
        )
        geometry = compute_ground_track(
            satellite_for(time_step), times, instrument_index=index
        ).geometry.iloc[0]
        prepare(geometry)
        predicted = np.zeros(len(grid), dtype=bool)
        predicted[compared] = contains_xy(
            geometry, lon.ravel()[compared], lat.ravel()[compared]
        )
        scores.append(compute_mask_iou(observed, predicted, weights))
    return min(scores)


def write_fixture_and_record(day, tle, entry, orbit, results, fixtures) -> None:
    """Writes the reduced fixture and the validation record."""
    rng = np.random.default_rng(1)
    rows, feeds = [], []
    for index, feed in enumerate(FEEDS):
        table = fixtures[feed]
        picks = rng.choice(len(table), FIXTURE_PASSES, replace=False)
        rows.append(table[picks])
        feeds.append(np.full(FIXTURE_PASSES, index))
    rows, feeds = np.concatenate(rows), np.concatenate(feeds)
    FIXTURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        FIXTURE_PATH,
        tle=np.array(tle),
        day=np.array(day.isoformat()),
        passes=rows,
        feeds=feeds,
    )
    satellite = entry.to_satellite(orbit, time_step=SCAN_PERIOD, platform="gpm")
    differences = []
    for (lat, lon, seconds), index in zip(rows, feeds):
        time = day + timedelta(seconds=float(seconds))
        found = collect_observations(
            Point(latitude=float(lat), longitude=float(lon)),
            satellite,
            time - timedelta(minutes=10),
            time + timedelta(minutes=10),
            int(index),
        )
        differences.append(
            min((e - pd.Timestamp(time)).total_seconds() for e in found.epoch)
            if not found.empty
            else np.nan
        )
    differences = np.array(differences, dtype=float)
    metrics = []
    for name in (
        "pass_recall",
        "pass_precision",
        "epoch_p95",
        "swath_width_median",
        "footprint_iou_median",
    ):
        for feed in FEEDS:
            metrics.append(
                {"name": name, "value": float(results[feed][name]), "scope": feed}
            )
    record = {
        "entry_id": "gmi",
        "standard_version": STANDARD_VERSION,
        "configurations": {
            "feeds": {
                "digest": entry.configurations["feeds"].digest,
                "evidence": "rerun",
                "tatc_version": tatc.__version__,
                "source": "docs/validation/validate_catalog_gmi.py",
                "reference": "GPM_1CGPMGMI version 08",
                "period": f"{day:%Y-%m-%d}",
                "platform": "gpm",
                # the standard applies its minimum gap to this value
                "holdout_gap_days": (
                    None if day == FIT_DAY else abs((day - FIT_DAY).days)
                ),
                "metrics": metrics,
                "fixture_metrics": {
                    "passes_found": 100 * float(np.isfinite(differences).mean()),
                    "epoch_p95": summarize(differences)["p95"],
                },
            }
        },
    }
    RECORD_PATH.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    """Validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day", choices=sorted(ELEMENTS))
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    tle = list(ELEMENTS[arguments.day])
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("gmi")
    specifications = entry.configurations["feeds"].instruments
    reference = Reference(day, load_day(day, arguments.directory))
    print(
        f"{arguments.day}: {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute"
    )

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, time_step=time_step, platform="gpm")

    # a lattice with a typical spacing of 500 km within the latitudes GMI observes
    points = generate_points_fibonacci_lattice(500e3)
    points = points[np.abs(points.geometry.y) < 72].reset_index(drop=True)

    results, fixtures = {}, {}
    for index, feed in enumerate(FEEDS):
        cone_angle = specifications[index].parameters["cone_angle"]
        satellite = satellite_for(SCAN_PERIOD)
        observed = reference.observe_passes(points, feed, cone_angle)
        results[feed] = score_passes(reference, points, satellite, index, observed)
        results[feed]["swath_width_median"] = score_swath_width(
            reference, satellite.instruments[index], feed
        )
        results[feed]["footprint_iou_median"] = score_swath_area(
            reference, satellite_for, orbit, index, feed
        )
        fixtures[feed] = np.stack(
            [
                points.geometry.y.to_numpy()[observed.point_id],
                points.geometry.x.to_numpy()[observed.point_id],
                observed.seconds.to_numpy(),
            ],
            axis=-1,
        )
        print(feed, json.dumps(results[feed]))

    if arguments.write:
        write_fixture_and_record(day, tle, entry, orbit, results, fixtures)
        print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()
