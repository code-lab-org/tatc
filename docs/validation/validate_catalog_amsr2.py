"""
Validates the AMSR2 catalog entry against the GPM Level 1C AMSR2 product, using
the exact values of the entry through public TAT-C methods. The tests follow
`ValidateConicalAMSR2`, whose day of data the entry's angles were fitted to:
running on that day checks this script against the notebook, while running on a
day at least 30 days from it tests the fitted angles on data held out from the
fit. Element sets for other days can be given in a file of the history of
GCOM-W's element sets.

Needs the `validation` extra and an Earthdata token in `EARTHDATA_TOKEN`:

    python docs/validation/validate_catalog_amsr2.py <data directory> 2026-09-26
    python docs/validation/validate_catalog_amsr2.py <data directory> 2026-11-01

It also validates the entry's sample size configuration, following
`ValidateSampleSizes`. The gated metric compares sizes with footprints published
by JAXA, which are independent of the fit of the angles. The spacing of samples on
the day is reported, as the standard does not gate it. To run only these tests,
on data already downloaded, use `--sample-sizes-only`.

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
import geopandas as gpd
import h5py
import numpy as np
import pandas as pd
import requests
from pyproj import Geod, Transformer
from scipy.spatial import cKDTree
from shapely import contains_xy, prepare

import tatc
from tatc.analysis import (
    collect_observations,
    compute_ground_track,
    compute_sample_sizes,
)
from tatc.catalog import STANDARD_VERSION, get_entry
from tatc.catalog.metrics import compute_mask_iou, compute_pass_scores, summarize
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import CircularOrbit, GeneralPerturbationsOrbit, Point

REPOSITORY = pathlib.Path(__file__).resolve().parents[2]
RECORD_PATH = REPOSITORY / "src" / "tatc" / "resources" / "validation" / "amsr2.json"
FIXTURE_PATH = REPOSITORY / "tests" / "catalog" / "fixtures" / "amsr2.npz"

FIT_DAY = datetime(2026, 9, 26, tzinfo=timezone.utc)
# public element sets from CelesTrak, with epochs after each day of data
ELEMENTS = {
    "2026-09-26": (
        "1 38337U 12025A   26277.29678146  .00000240  00000+0  63281-4 0  9990",
        "2 38337  98.2100 216.2278 0001644  77.4488  41.8708 14.57111943764929",
    ),
    "2026-11-01": (
        "1 38337U 12025A   26312.09154684  .00000538  00000+0  12959-3 0  9999",
        "2 38337  98.2164 250.5392 0001184  78.0816 292.9047 14.57095606769992",
    ),
}
# the group whose samples the scan was fitted to, and that of the 89 GHz channels
GROUPS = ("S1", "S5")
SCAN_PERIOD = 1.5
TIME_STEPS = (SCAN_PERIOD, 10, 30)
GRID_STEP, ORBIT_MARGIN = 0.1, 60
FIXTURE_PASSES = 80

SAMPLE_CONFIGURATION = "sample_size"
# published footprints: instrument, altitude (km), along and across track (km), from
# JAXA's AMSR2 channel set (ground resolution at the scan center)
PUBLISHED_FOOTPRINTS = (
    ("AMSR2 6.925 and 7.3 GHz", 700, 62, 35),
    ("AMSR2 10.65 GHz", 700, 42, 24),
    ("AMSR2 18.7 GHz", 700, 22, 14),
    ("AMSR2 23.8 GHz", 700, 26, 15),
    ("AMSR2 36.5 GHz", 700, 12, 7),
    ("AMSR2 89.0 GHz", 700, 5, 3),
)
PUBLISHED_INCLINATION = 98.2
PUBLISHED_EPOCH = datetime(2026, 10, 1, tzinfo=timezone.utc)
SIZE_TOLERANCE = 0.03
# an instrument of each group, whose spacing is measured
SPACING_INSTRUMENTS = {"S1": "AMSR2 10.65 GHz", "S5": "AMSR2 89.0 GHz"}
SPACING_SCAN_STRIDE = 200
FIXTURE_SPACING_SAMPLES = 100

GEOD = Geod(ellps="WGS84")
TO_ECEF = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
TO_GEODETIC = Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)


def load_day(day: datetime, directory: pathlib.Path) -> dict[str, np.ndarray]:
    """Downloads and reduces all orbits of a day, or loads them from a cache."""
    cache = directory / f"amsr2_l1c_{day:%Y%m%d}.npz"
    if not cache.exists():
        if earthaccess.login(strategy="environment") is None:
            raise RuntimeError("Earthdata login failed: is EARTHDATA_TOKEN set?")
        session = earthaccess.get_requests_https_session()
        granules = requests.get(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params={
                "short_name": "GPM_1CGCOMW1AMSR2",
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
                        **{f"{k}_lat": f[f"{k}/Latitude"][:] for k in GROUPS},
                        **{f"{k}_lon": f[f"{k}/Longitude"][:] for k in GROUPS},
                    }
                )
        reduced = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        # keep each scan once (orbit files overlap), in time order, within the day
        _, first = np.unique(reduced["time"], return_index=True)
        keep = first[(reduced["time"][first] >= 0) & (reduced["time"][first] < 86400)]
        np.savez_compressed(
            cache, granules=len(granules), **{k: v[keep] for k, v in reduced.items()}
        )
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


class Reference:
    """A day of AMSR2 scans, with what is derived from them."""

    def __init__(self, day: datetime, data: dict[str, np.ndarray]):
        self.day = day
        valid = np.all(
            [np.all(np.abs(data[f"{k}_lat"]) <= 90, axis=1) for k in GROUPS], axis=0
        )
        self.time = data["time"][valid]
        self.altitude = data["sc_alt"][valid].astype(float) * 1e3
        self.position = ecef(
            data["sc_lon"][valid].astype(float),
            data["sc_lat"][valid].astype(float),
            self.altitude,
        )
        self.contiguous = (
            np.r_[np.diff(self.time) < 2, False] & np.r_[False, np.diff(self.time) < 2]
        )
        self.foot = {
            k: (
                data[f"{k}_lon"][valid].astype(float),
                data[f"{k}_lat"][valid].astype(float),
            )
            for k in GROUPS
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
        valid &= np.diff(self.time)[k] < 2
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

    def observe_passes(self, points, cone_angle: float) -> pd.DataFrame:
        """Finds the passes on which the S1 samples observed points, and their times."""
        unit = up(*self.foot["S1"])
        flat = unit.reshape(-1, 3)
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        tree = cKDTree(flat)
        rows = []
        for p, hits in enumerate(
            tree.query_ball_point(point_unit, r=2 * np.sin(60 / 6371 / 2), workers=-1)
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


def score_swath_width(reference: Reference, instrument) -> float:
    """Gets the median absolute swath width difference (km) from the scanned arcs."""
    i = np.nonzero(reference.contiguous)[0][::20]
    lon, lat = reference.foot["S1"][0][i], reference.foot["S1"][1][i]
    chord = GEOD.inv(lon[:, 0], lat[:, 0], lon[:, -1], lat[:, -1])[2]
    width = np.vectorize(instrument.get_swath_width)(reference.altitude[i])
    return float(np.median(np.abs(width - chord))) / 1e3


def score_passes(reference, points, satellite, observed) -> dict[str, float]:
    """Compares predicted and observed passes of points."""
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


def get_orbit_scans(reference: Reference, period: float) -> np.ndarray:
    """
    Gets the scans of the first orbit from 12:00 UTC without a gap in the data, as
    `ValidateConicalAMSR2` used an orbit after the gap of its day.
    """
    time = reference.time
    start = 12 * 3600
    while start + period < time[-1]:
        scans = np.nonzero((time >= start) & (time < start + period))[0]
        if len(scans) < 2:
            start += 600
            continue
        breaks = np.nonzero(np.diff(time[scans]) >= 2)[0]
        if len(breaks):
            start = time[scans[breaks[-1] + 1]]
        elif time[scans[0]] - start >= 2 or time[scans[-1]] < start + period - 2:
            start += 600
        else:
            return scans
    raise ValueError("no orbit without a gap after 12:00 UTC")


def score_swath_area(reference, satellite_for, orbit) -> float:
    """Gets the lowest IoU of the swath over one orbit, at each time step."""
    period = 86400 / (orbit.elements[0].mean_motion * 86400 / 360)
    scans = get_orbit_scans(reference, period)
    lon, lat = np.meshgrid(
        np.arange(-180 + GRID_STEP / 2, 180, GRID_STEP),
        np.arange(-90 + GRID_STEP / 2, 90, GRID_STEP),
    )
    weights = (np.cos(np.radians(lat)) * (GRID_STEP * 111.32) ** 2).ravel()
    grid = up(lon.ravel(), lat.ravel())
    unit = up(*reference.foot["S1"])
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
    print(f"swath orbit from {first / 3600:.2f} to {last / 3600:.2f} h UTC")
    for time_step in TIME_STEPS:
        times = reference.to_utc(
            first + time_step * (np.arange(int((last - first) // time_step)) + 0.5)
        )
        geometry = compute_ground_track(satellite_for(time_step), times).geometry.iloc[
            0
        ]
        prepare(geometry)
        predicted = np.zeros(len(grid), dtype=bool)
        predicted[compared] = contains_xy(
            geometry, lon.ravel()[compared], lat.ravel()[compared]
        )
        score = compute_mask_iou(observed, predicted, weights)
        print(f"swath IoU at {time_step:g} s: {score:.4f}")
        scores.append(score)
    return min(scores)


def score_published_footprints(entry) -> tuple[float, float]:
    """
    Compares sample sizes at the center of the scan with published footprints.
    Returns the percent of dimensions within tolerance and the largest error (%).
    """
    errors = []
    for altitude in sorted({row[1] for row in PUBLISHED_FOOTPRINTS}):
        orbit = CircularOrbit(
            altitude=altitude * 1e3,
            inclination=PUBLISHED_INCLINATION,
            epoch=PUBLISHED_EPOCH,
        )
        satellite = entry.to_satellite(
            orbit, SAMPLE_CONFIGURATION, time_step=SCAN_PERIOD, platform="gcom-w"
        )
        names = [instrument.name for instrument in satellite.instruments]
        time = PUBLISHED_EPOCH + timedelta(minutes=10)
        track = orbit.to_gp_orbit().get_orbit_track(time)
        for name, _, along, cross in (
            row for row in PUBLISHED_FOOTPRINTS if row[1] == altitude
        ):
            target = satellite.get_instrument(
                names.index(name)
            ).compute_footprint_center(track)
            sizes = compute_sample_sizes(
                gpd.GeoDataFrame(
                    {
                        "satellite": [satellite.name],
                        "instrument": [name],
                        "epoch": pd.DatetimeIndex([time]),
                    },
                    geometry=gpd.points_from_xy(
                        [target.longitude.degrees], [target.latitude.degrees]
                    ),
                    crs="EPSG:4326",
                ),
                satellite,
            )
            pair = (
                sizes.sample_along.iloc[0] / 1e3 / along - 1,
                sizes.sample_cross.iloc[0] / 1e3 / cross - 1,
            )
            print(f"{name}, {altitude} km: {pair[0]:+.2%} along, {pair[1]:+.2%} across")
            errors.extend(abs(error) for error in pair)
    errors = np.array(errors)
    return 100 * float(np.mean(errors <= SIZE_TOLERANCE)), 100 * float(errors.max())


def measure_spacing(reference, group: str) -> np.ndarray:
    """
    Measures the spacing of a group's samples as `ValidateSampleSizes` does, from
    half the distance between the samples on either side. Returns rows of latitude,
    longitude, seconds, spacing along the arc (km), and spacing between scans (km).
    """
    period = np.diff(reference.time)
    regular = np.abs(period - SCAN_PERIOD) < 0.05
    usable = np.nonzero(regular[:-1] & regular[1:])[0] + 1
    scans = usable[::SPACING_SCAN_STRIDE]
    lon, lat = reference.foot[group]
    position = ecef(lon, lat)
    arc = np.linalg.norm(position[scans, 2:] - position[scans, :-2], axis=-1) / 2e3
    between = (
        np.linalg.norm(position[scans + 1, 1:-1] - position[scans - 1, 1:-1], axis=-1)
        / 2e3
    )
    return np.stack(
        [
            lat[scans, 1:-1].ravel(),
            lon[scans, 1:-1].ravel(),
            np.repeat(reference.time[scans], lat.shape[1] - 2),
            arc.ravel(),
            between.ravel(),
        ],
        axis=-1,
    )


def compute_spacing_ratios(day, satellite, name: str, rows) -> tuple[float, float]:
    """Gets the median ratios of predicted to measured spacing along arcs and scans."""
    sizes = compute_sample_sizes(
        gpd.GeoDataFrame(
            {
                "satellite": satellite.name,
                "instrument": name,
                "epoch": pd.DatetimeIndex(
                    [day + timedelta(seconds=float(s)) for s in rows[:, 2]]
                ),
            },
            geometry=gpd.points_from_xy(rows[:, 1], rows[:, 0]),
            crs="EPSG:4326",
        ),
        satellite,
    )
    return (
        float(np.median(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
        float(np.median(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 4])),
    )


def update_fixture(**arrays) -> None:
    """Writes arrays to the fixture, keeping the others."""
    existing = {}
    if FIXTURE_PATH.exists():
        with np.load(FIXTURE_PATH, allow_pickle=False) as data:
            existing = {key: data[key] for key in data.files}
    FIXTURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(FIXTURE_PATH, **{**existing, **arrays})


def update_record(key: str, configuration: dict) -> None:
    """Writes the results of a configuration to the record, keeping the others."""
    record = {"entry_id": "amsr2", "configurations": {}}
    if RECORD_PATH.exists():
        record = json.loads(RECORD_PATH.read_text(encoding="utf-8"))
    record["standard_version"] = STANDARD_VERSION
    record["configurations"][key] = configuration
    RECORD_PATH.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")


def get_holdout_gap(day: datetime) -> int | None:
    """Gets the days from the fit to a day, to which the standard applies its gap."""
    return None if day == FIT_DAY else abs((day - FIT_DAY).days)


def validate_sample_sizes(day, tle, entry, orbit, reference, write: bool) -> None:
    """Scores the sample size configuration, and writes its results if requested."""
    within, largest = score_published_footprints(entry)
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    metrics = [
        {"name": "sample_size_within_3pct", "value": within},
        {"name": "sample_size_largest_error", "value": largest},
    ]
    satellite = entry.to_satellite(
        orbit, SAMPLE_CONFIGURATION, time_step=SCAN_PERIOD, platform="gcom-w"
    )
    names = [instrument.name for instrument in satellite.instruments]
    rng = np.random.default_rng(2)
    reduced, indices, fixture_metrics = [], [], {}
    for group, name in SPACING_INSTRUMENTS.items():
        rows = measure_spacing(reference, group)
        arc, between = compute_spacing_ratios(day, satellite, name, rows)
        print(
            f"{group} spacing, predicted / measured: {arc:.4f} arc, {between:.4f} scans"
        )
        metrics += [
            {"name": "spacing_arc_ratio_median", "value": arc, "scope": group},
            {"name": "spacing_scan_ratio_median", "value": between, "scope": group},
        ]
        picks = rows[rng.choice(len(rows), FIXTURE_SPACING_SAMPLES, replace=False)]
        arc, between = compute_spacing_ratios(day, satellite, name, picks)
        fixture_metrics[f"spacing_arc_ratio_median_{group}"] = arc
        fixture_metrics[f"spacing_scan_ratio_median_{group}"] = between
        reduced.append(picks)
        indices.append(np.full(FIXTURE_SPACING_SAMPLES, names.index(name)))
    if not write:
        return
    update_fixture(
        tle=np.array(tle),
        day=np.array(day.isoformat()),
        footprints=np.array(
            [
                [names.index(name), altitude, along, cross]
                for name, altitude, along, cross in PUBLISHED_FOOTPRINTS
            ]
        ),
        spacing=np.concatenate(reduced),
        spacing_instruments=np.concatenate(indices),
    )
    update_record(
        SAMPLE_CONFIGURATION,
        {
            "digest": entry.configurations[SAMPLE_CONFIGURATION].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_amsr2.py",
            "reference": "Footprints of JAXA's AMSR2 channel set; spacing of "
            "GPM_1CGCOMW1AMSR2 version 08",
            "period": f"{day:%Y-%m-%d}",
            "platform": "gcom-w",
            # the angles were not fitted to the published footprints, and the
            # spacing is measured on this day
            "holdout_gap_days": get_holdout_gap(day),
            "metrics": metrics,
            "fixture_metrics": {
                "sample_size_within_3pct": within,
                **fixture_metrics,
            },
        },
    )


def write_fixture_and_record(day, tle, entry, orbit, results, table) -> None:
    """Writes the reduced fixture and the validation record."""
    rng = np.random.default_rng(1)
    rows = table[rng.choice(len(table), FIXTURE_PASSES, replace=False)]
    update_fixture(tle=np.array(tle), day=np.array(day.isoformat()), passes=rows)
    satellite = entry.to_satellite(orbit, time_step=SCAN_PERIOD, platform="gcom-w")
    differences = []
    for lat, lon, seconds in rows:
        time = day + timedelta(seconds=float(seconds))
        found = collect_observations(
            Point(latitude=float(lat), longitude=float(lon)),
            satellite,
            time - timedelta(minutes=10),
            time + timedelta(minutes=10),
        )
        differences.append(
            min(
                ((e - pd.Timestamp(time)).total_seconds() for e in found.epoch),
                key=abs,
            )
            if not found.empty
            else np.nan
        )
    differences = np.array(differences, dtype=float)
    metrics = [
        {"name": name, "value": float(results[name])}
        for name in (
            "pass_recall",
            "pass_precision",
            "epoch_p95",
            "swath_width_median",
            "footprint_iou_median",
        )
    ]
    update_record(
        "swath",
        {
            "digest": entry.configurations["swath"].digest,
            "evidence": "rerun",
            "tatc_version": tatc.__version__,
            "source": "docs/validation/validate_catalog_amsr2.py",
            "reference": "GPM_1CGCOMW1AMSR2 version 08",
            "period": f"{day:%Y-%m-%d}",
            "platform": "gcom-w",
            "holdout_gap_days": get_holdout_gap(day),
            "metrics": metrics,
            "fixture_metrics": {
                "passes_found": 100 * float(np.isfinite(differences).mean()),
                "epoch_p95": summarize(differences)["p95"],
            },
        },
    )


def read_elements(path: pathlib.Path, day: datetime) -> list[str]:
    """Reads the unique element sets with epochs within 1.5 days of the middle of a day."""
    lines = [line.rstrip() for line in path.read_text(encoding="utf-8").splitlines()]
    pairs = {}
    for first, second in zip(lines, lines[1:]):
        if first.startswith("1 ") and second.startswith("2 "):
            year, day_of_year = int(first[18:20]), float(first[20:32])
            epoch = datetime(2000 + year, 1, 1, tzinfo=timezone.utc) + timedelta(
                days=day_of_year - 1
            )
            offset = abs((epoch - day - timedelta(hours=12)).total_seconds())
            if offset <= 1.5 * 86400:
                pairs[epoch] = [first, second]
    if not pairs:
        raise ValueError(f"{path} has no element sets near {day:%Y-%m-%d}")
    return [line for epoch in sorted(pairs) for line in pairs[epoch]]


def main() -> None:
    """Validates the entry on a day of data."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("day")
    parser.add_argument("--elements", type=pathlib.Path)
    parser.add_argument("--sample-sizes-only", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    arguments.directory.mkdir(parents=True, exist_ok=True)

    day = datetime.fromisoformat(arguments.day).replace(tzinfo=timezone.utc)
    if arguments.elements is not None:
        tle = read_elements(arguments.elements, day)
    elif arguments.day in ELEMENTS:
        tle = list(ELEMENTS[arguments.day])
    else:
        parser.error(f"--elements is needed for {arguments.day}")
    print(f"{len(tle) // 2} element sets")
    orbit = GeneralPerturbationsOrbit.from_tle(tle)
    entry = get_entry("amsr2")
    reference = Reference(day, load_day(day, arguments.directory))
    print(
        f"{arguments.day}: {len(reference.time)} scans, "
        f"{len(reference.gaps)} gaps over a minute"
    )
    validate_sample_sizes(day, tle, entry, orbit, reference, arguments.write)
    if arguments.sample_sizes_only:
        return

    def satellite_for(time_step: float):
        return entry.to_satellite(orbit, time_step=time_step, platform="gcom-w")

    # a lattice with a typical spacing of 500 km, as in ValidateConicalAMSR2
    points = generate_points_fibonacci_lattice(500e3)
    cone_angle = entry.configurations["swath"].instruments[0].parameters["cone_angle"]
    satellite = satellite_for(SCAN_PERIOD)
    observed = reference.observe_passes(points, cone_angle)
    results = score_passes(reference, points, satellite, observed)
    results["swath_width_median"] = score_swath_width(
        reference, satellite.instruments[0]
    )
    results["footprint_iou_median"] = score_swath_area(reference, satellite_for, orbit)
    table = np.stack(
        [
            points.geometry.y.to_numpy()[observed.point_id],
            points.geometry.x.to_numpy()[observed.point_id],
            observed.seconds.to_numpy(),
        ],
        axis=-1,
    )
    print(json.dumps(results))

    if arguments.write:
        write_fixture_and_record(day, tle, entry, orbit, results, table)
        print(f"wrote {RECORD_PATH} and {FIXTURE_PATH}")


if __name__ == "__main__":
    main()
