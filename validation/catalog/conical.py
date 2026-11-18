"""
Helpers shared by the scripts that validate conical scanners (AMSR2, GMI, SMAP's
radiometer, SSMIS, and TMI). A point is observed on a pass if it lies in a cell of
the lattice of a scan's samples, at the time it crosses the cone, from the
spacecraft's recorded positions; the swath width is compared with the distance
between the first and last samples of scans; and sample sizes at the center of the
scan are compared with published footprints.

The GPM Level 1C products of AMSR2 and GMI are read whole (`load_day`) into a
`Reference`, with the scans' times to the millisecond, and scored with
`score_swath` and `validate_sample_sizes`. SSMIS and TMI are read as cross-track
scanners are (`scanner.Reference`), and scored with `validate_cone_swath`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import io
import json
import pathlib
from dataclasses import dataclass
from datetime import datetime, timedelta

import earthaccess
import h5py
import numpy as np
import pandas as pd
from pyproj import Geod
from scipy.spatial import cKDTree
from shapely import contains_xy, prepare

from tatc.analysis import compute_ground_track, compute_sample_sizes
from tatc.catalog.metrics import compute_mask_iou
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import CircularOrbit
import scanner
from validation_io import (
    PUBLISHED_EPOCH,
    TO_GEODETIC,
    build_grid,
    compute_spacing_ratios,
    ecef,
    get_with_retries,
    login_earthdata,
    normalize,
    recompute_passes,
    score_passes,
    score_published,
    score_swath_area,
    to_observations,
    up,
    write_results,
)

GEOD = Geod(ellps="WGS84")
GATED = (
    "pass_recall",
    "pass_precision",
    "epoch_p95",
    "swath_width_median",
    "footprint_iou_median",
)


def cone_crossing(
    seconds, position, target, k: int, cone_angle: float, window: float, signed=False
) -> float:
    """
    Gets the time (s) at which a target (Earth-fixed, m) crosses a cone about the
    geodetic nadir, by bisection within a window (s) on either side of scan k, with
    the spacecraft's position interpolated between scans. The bisection takes the
    angle from the cone to decrease through the crossing, or, if `signed`, follows
    its sign at the lower bound.
    """

    def angle(t: float) -> float:
        i = int(np.clip(np.searchsorted(seconds, t) - 1, 0, len(seconds) - 2))
        f = (t - seconds[i]) / (seconds[i + 1] - seconds[i])
        at = position[i] + f * (position[i + 1] - position[i])
        lon, lat, _ = TO_GEODETIC.transform(*at)
        return (
            np.degrees(np.arccos(np.dot(normalize(target - at), -up(lon, lat))))
            - cone_angle
        )

    lower, upper = seconds[k] - window, seconds[k] + window
    if angle(lower) * angle(upper) > 0:
        return np.nan
    for _ in range(30):
        middle = (lower + upper) / 2
        if (angle(middle) > 0) == (angle(lower) > 0 if signed else True):
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2


def to_observed(day: datetime, rows: list[dict]) -> pd.DataFrame:
    """Gets observed passes, from rows of `point_id` and `seconds` of a day."""
    observed = pd.DataFrame(rows).astype({"point_id": "int64"})
    observed["time"] = pd.Timestamp(day) + pd.to_timedelta(observed.seconds, unit="s")
    return observed


def pick_cone_passes(points, observed, count: int, seed=1) -> np.ndarray:
    """
    Picks observed passes for a test fixture, at random (from a seed or a random
    generator): rows of the point's latitude and longitude, and the seconds of the
    pass from the start of the day.
    """
    rng = np.random.default_rng(seed)
    table = np.stack(
        [
            points.geometry.y.to_numpy()[observed.point_id],
            points.geometry.x.to_numpy()[observed.point_id],
            observed.seconds.to_numpy(),
        ],
        axis=-1,
    )
    return table[rng.choice(len(table), count, replace=False)]


def score_chord_width(instrument, lon, lat, altitude, scans) -> float:
    """
    Gets the median absolute difference (km) of the swath width from the distance
    along the surface between the first and last samples of scans, at the
    spacecraft's altitudes (m).
    """
    chord = GEOD.inv(lon[scans, 0], lat[scans, 0], lon[scans, -1], lat[scans, -1])[2]
    width = np.vectorize(instrument.get_swath_width)(altitude[scans])
    return float(np.median(np.abs(width - chord))) / 1e3


def orbit_period(orbit) -> float:
    """Gets the period (s) of an orbit's first element set."""
    return 86400 / (orbit.elements[0].mean_motion * 86400 / 360)


def score_swath_steps(
    satellite_for, index: int, day, first, last, time_steps, grid, compared, observed
) -> float:
    """
    Gets the lowest IoU, over time steps, of an instrument's swath from `first` to
    `last` (s from the start of a day) with the observed area of a grid (see
    `validation_io.score_swath_area`, which counts time from a time).
    """
    lon, lat, weights = grid
    scores = []
    for time_step in time_steps:
        seconds = first + time_step * (
            np.arange(int((last - first) // time_step)) + 0.5
        )
        times = [day + timedelta(seconds=float(s)) for s in seconds]
        geometry = compute_ground_track(
            satellite_for(time_step), times, instrument_index=index
        ).geometry.iloc[0]
        prepare(geometry)
        predicted = np.zeros(len(lon), dtype=bool)
        predicted[compared] = contains_xy(geometry, lon[compared], lat[compared])
        scores.append(compute_mask_iou(observed, predicted, weights))
        print(f"swath IoU at {time_step:g} s: {scores[-1]:.4f}", flush=True)
    return min(scores)


def compute_center_errors(
    entry, key: str, footprints, inclination: float, *, time_step: float, platform
) -> np.ndarray:
    """
    Gets the relative errors of sample sizes at the center of the scan (the
    footprint center) against published footprints, given as rows of instrument
    name, altitude (km), and the sizes (km) along and across track, on circular
    orbits of an inclination (deg). Gives two errors for each row.
    """
    errors = []
    for name, altitude, along, cross in footprints:
        orbit = CircularOrbit(
            altitude=altitude * 1e3, inclination=inclination, epoch=PUBLISHED_EPOCH
        )
        satellite = entry.to_satellite(
            orbit, configuration=key, time_step=time_step, platform=platform
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
        print(
            f"{name}, {altitude} km: {errors[-2]:+.2%} along, {errors[-1]:+.2%} across"
        )
    return np.array(errors)


def measure_arc_spacing(lon, lat, seconds, scans) -> np.ndarray:
    """
    Measures the spacing of samples as `ValidateSampleSizes` does, from half the
    distance between the samples on either side, on scans whose neighbors are one
    scan period before and after. Returns rows of latitude, longitude, seconds of
    the scan, spacing along the arc (km), and spacing between scans (km), for the
    samples between the first and last.
    """
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
            np.repeat(seconds[scans], lat.shape[1] - 2),
            arc.ravel(),
            between.ravel(),
        ],
        axis=-1,
    )


def regular_scans(seconds, scan_period: float) -> np.ndarray:
    """Gets the scans whose neighbors are one scan period before and after."""
    regular = np.abs(np.diff(seconds) - scan_period) < 0.05
    return np.nonzero(regular[:-1] & regular[1:])[0] + 1


def compute_arc_spacing_ratios(satellite, name: str, day, rows) -> tuple[float, float]:
    """
    Gets the median ratios of predicted to measured spacing along arcs and between
    scans, for rows of `measure_arc_spacing` of an instrument (the spacing along the
    arc is across track, and between scans along track).
    """
    epochs = [day + timedelta(seconds=float(s)) for s in rows[:, 2]]
    between, arc = compute_spacing_ratios(
        satellite, [name] * len(rows), epochs, rows[:, [0, 1, 4, 3]]
    )
    return arc, between


# the GPM Level 1C products of AMSR2 and GMI


@dataclass(frozen=True)
class Product:
    """
    A GPM Level 1C product of a conical scanner, read whole: its groups of
    samples (the first is that whose sample count the lattices use), its scan
    period (s), the longest time (s) from a scan to the next that are contiguous,
    and the radius (km) within which samples are near a point. GMI's product also
    gives the spacecraft's orientation.
    """

    short_name: str
    cache_prefix: str
    groups: tuple[str, ...]
    scan_period: float
    contiguous_gap: float
    pass_radius: float
    orientation: bool = False
    version: str = "08"

    @property
    def time_steps(self) -> tuple[float, ...]:
        """The time steps (s) at which the swath area is compared."""
        return (self.scan_period, 10, 30)


ORBIT_MARGIN = 60  # s


def load_day(product: Product, day: datetime, directory: pathlib.Path) -> dict:
    """Downloads and reduces all orbits of a day, or loads them from a cache."""
    cache = directory / f"{product.cache_prefix}_{day:%Y%m%d}.npz"
    if not cache.exists():
        login_earthdata()
        session = earthaccess.get_requests_https_session()
        granules = get_with_retries(
            "https://cmr.earthdata.nasa.gov/search/granules.umm_json",
            params={
                "short_name": product.short_name,
                "version": product.version,
                "temporal": f"{day:%Y-%m-%d}T00:00:00Z,{day:%Y-%m-%d}T23:59:59Z",
                "page_size": 100,
            },
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
                part = {
                    "time": seconds.dt.total_seconds().to_numpy(),
                    "sc_lat": s1["SCstatus/SClatitude"][:],
                    "sc_lon": s1["SCstatus/SClongitude"][:],
                    "sc_alt": s1["SCstatus/SCaltitude"][:],
                }
                if product.orientation:
                    part["orientation"] = s1["SCstatus/SCorientation"][:]
                part.update({f"{k}_lat": f[f"{k}/Latitude"][:] for k in product.groups})
                part.update(
                    {f"{k}_lon": f[f"{k}/Longitude"][:] for k in product.groups}
                )
                parts.append(part)
        reduced = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
        # keep each scan once (orbit files overlap), in time order, within the day
        _, first = np.unique(reduced["time"], return_index=True)
        keep = first[(reduced["time"][first] >= 0) & (reduced["time"][first] < 86400)]
        counts = {} if product.orientation else {"granules": len(granules)}
        np.savez_compressed(cache, **counts, **{k: v[keep] for k, v in reduced.items()})
    return dict(np.load(cache))


class Reference:
    """
    A day of scans of a GPM Level 1C product, with what is derived from them: the
    times of scans (s from the start of the day), the spacecraft's positions and
    velocities (Earth-fixed), and the samples' positions of each group.
    """

    def __init__(self, product: Product, day: datetime, data: dict[str, np.ndarray]):
        self.day = day
        self.scan_period = product.scan_period
        valid = np.all(
            [np.all(np.abs(data[f"{k}_lat"]) <= 90, axis=1) for k in product.groups],
            axis=0,
        )
        self.time = data["time"][valid]
        self.altitude = data["sc_alt"][valid].astype(float) * 1e3
        self.position = ecef(
            data["sc_lon"][valid].astype(float),
            data["sc_lat"][valid].astype(float),
            self.altitude,
        )
        self.velocity = np.gradient(self.position, self.time, axis=0)
        self.contiguous_gap = product.contiguous_gap
        self.contiguous = (
            np.r_[np.diff(self.time) < product.contiguous_gap, False]
            & np.r_[False, np.diff(self.time) < product.contiguous_gap]
        )
        self.foot = {
            k: (
                data[f"{k}_lon"][valid].astype(float),
                data[f"{k}_lat"][valid].astype(float),
            )
            for k in product.groups
        }
        self.samples = self.foot[product.groups[0]][0].shape[1]
        self.gaps = [
            (self.time[i], self.time[i + 1])
            for i in np.nonzero(np.diff(self.time) > 60)[0]
        ]

    @property
    def span(self) -> tuple[datetime, datetime]:
        """The times of the first and last scans."""
        return tuple(
            self.day + timedelta(seconds=float(s))
            for s in (self.time[0], self.time[-1])
        )

    @property
    def gap_times(self) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
        """The periods without data over a minute, as pairs of times."""
        start = pd.Timestamp(self.day)
        return [
            (start + pd.Timedelta(seconds=a), start + pd.Timedelta(seconds=b))
            for a, b in self.gaps
        ]

    def in_cell(self, unit: np.ndarray, points, k, j) -> np.ndarray:
        """Gets whether points lie in the cells of footprint centers at (k, j)."""
        valid = (k >= 0) & (k < len(self.time) - 1) & (j >= 0) & (j < self.samples - 1)
        k = np.clip(k, 0, len(self.time) - 2)
        j = np.clip(j, 0, self.samples - 2)
        valid &= np.diff(self.time)[k] < self.contiguous_gap
        corners = [unit[k, j], unit[k, j + 1], unit[k + 1, j + 1], unit[k + 1, j]]
        signs = np.stack(
            [
                np.sum(np.cross(p, q) * points, axis=-1)
                for p, q in zip(corners, corners[1:] + corners[:1])
            ]
        )
        return valid & (np.all(signs >= 0, axis=0) | np.all(signs <= 0, axis=0))

    def observe_passes(
        self, points, group: str, cone_angle: float, radius: float
    ) -> pd.DataFrame:
        """
        Finds the passes on which a group's samples observed points (in a cell next
        to one of the eight samples nearest a point, among those within a radius, in
        km), and their times (when they cross the cone).
        """
        unit = up(*self.foot[group])
        flat = unit.reshape(-1, 3)
        point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
        tree = cKDTree(flat)
        rows = []
        for p, hits in enumerate(
            tree.query_ball_point(
                point_unit, r=2 * np.sin(radius / 6371 / 2), workers=-1
            )
        ):
            if not hits:
                continue
            hits = np.array(hits)
            scans = hits // self.samples
            order = np.argsort(self.time[scans])
            hits, scans = hits[order], scans[order]
            for group_hits in np.split(
                np.arange(len(hits)), np.nonzero(np.diff(self.time[scans]) > 600)[0] + 1
            ):
                nearest = hits[group_hits][
                    np.argsort(
                        np.linalg.norm(flat[hits[group_hits]] - point_unit[p], axis=-1)
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
                crossing = cone_crossing(
                    self.time,
                    self.position,
                    ecef(points.geometry.x.iloc[p], points.geometry.y.iloc[p]),
                    int(nearest[0] // self.samples),
                    cone_angle,
                    6,
                )
                if np.isfinite(crossing):
                    rows.append({"point_id": p, "seconds": float(crossing)})
        return to_observed(self.day, rows)

    def observe_orbit(self, group: str, scans: np.ndarray):
        """
        Gets the cells of a global grid compared over a range of scans, within 300
        km of a group's samples and more than a margin from the ends of the range,
        and which of them the samples observed (in a cell next to one of the eight
        nearest samples).

        Returns:
            tuple: The grid (longitudes, latitudes, and areas of its cells), and
            whether each cell is compared and was observed.
        """
        lon, lat, weights = build_grid()
        grid = up(lon, lat)
        unit = up(*self.foot[group])
        tree = cKDTree(unit[scans].reshape(-1, 3))
        distance, nearest = tree.query(grid, workers=-1)
        nearest_time = self.time[scans][nearest // self.samples] - self.time[scans[0]]
        duration = self.time[scans[-1]] - self.time[scans[0]]
        compared = (
            (distance < 2 * np.sin(300 / 6371 / 2))
            & (nearest_time > ORBIT_MARGIN)
            & (nearest_time < duration - ORBIT_MARGIN)
        )
        observed = np.zeros(len(grid), dtype=bool)
        _, nearest_eight = tree.query(grid[compared], k=8, workers=-1)
        for column in range(8):
            k, j = np.divmod(nearest_eight[:, column], self.samples)
            k = scans[k]
            inside = np.zeros(len(k), dtype=bool)
            for dk in (-1, 0):
                for dj in (-1, 0):
                    inside |= self.in_cell(unit, grid[compared], k + dk, j + dj)
            observed[np.flatnonzero(compared)] |= inside
        return (lon, lat, weights), compared, observed

    def measure_spacing(self, group: str, stride: int) -> np.ndarray:
        """Measures the spacing of a group's samples (see `measure_arc_spacing`)."""
        scans = regular_scans(self.time, self.scan_period)[::stride]
        return measure_arc_spacing(*self.foot[group], self.time, scans)


def score_swath(
    product: Product,
    reference: Reference,
    satellite_for,
    points,
    group: str,
    index: int,
    cone_angle: float,
    scans: np.ndarray,
) -> tuple[dict, pd.DataFrame]:
    """
    Scores an instrument's swath on a group's samples: its passes and their epochs,
    its swath width (on every 20th scan with contiguous neighbors), and its swath
    area over a range of scans.

    Returns:
        tuple: The results, and the observed passes.
    """
    satellite = satellite_for(product.scan_period)
    observed = reference.observe_passes(points, group, cone_angle, product.pass_radius)
    start, end = reference.span
    results = score_passes(
        points, satellite, index, observed, start, end, reference.gap_times
    )
    results["swath_width_median"] = score_chord_width(
        satellite.instruments[index],
        *reference.foot[group],
        reference.altitude,
        np.nonzero(reference.contiguous)[0][::20],
    )
    grid, compared, covered = reference.observe_orbit(group, scans)
    results["footprint_iou_median"] = score_swath_steps(
        satellite_for,
        index,
        reference.day,
        reference.time[scans[0]],
        reference.time[scans[-1]],
        product.time_steps,
        grid,
        compared,
        covered,
    )
    return results, observed


@dataclass(frozen=True)
class SampleSizes:
    """
    The sample size configuration of a conical scanner: its published footprints
    (rows of instrument name, altitude in km, and the sizes in km along and across
    track at the center of the scan, on orbits of an inclination in deg), an
    instrument of each group whose spacing is measured, and the reference data.
    """

    key: str
    footprints: tuple
    inclination: float
    spacing_instruments: dict[str, str]
    reference: str
    spacing_stride: int = 200
    fixture_spacing: int = 100


def validate_sample_sizes(
    sizes: SampleSizes, product: Product, run, orbit, tle, reference: Reference
) -> None:
    """
    Scores the sample size configuration of a run (see `validation_io.parse_run`),
    and writes its results if the run writes.
    """
    entry = run.entry
    within, largest = score_published(
        compute_center_errors(
            entry,
            sizes.key,
            sizes.footprints,
            sizes.inclination,
            time_step=product.scan_period,
            platform=run.platform,
        )
    )
    print(
        f"published dimensions: {within:.1f}% within 3%, largest error {largest:.2f}%"
    )
    metrics = [
        {"name": "sample_size_within_3pct", "value": within},
        {"name": "sample_size_largest_error", "value": largest},
    ]
    satellite = entry.to_satellite(
        orbit,
        configuration=sizes.key,
        time_step=product.scan_period,
        platform=run.platform,
    )
    names = [instrument.name for instrument in satellite.instruments]
    rng = np.random.default_rng(2)
    reduced, indices, fixture_metrics = [], [], {}
    for group, name in sizes.spacing_instruments.items():
        rows = reference.measure_spacing(group, sizes.spacing_stride)
        arc, between = compute_arc_spacing_ratios(satellite, name, run.day, rows)
        print(
            f"{group} spacing, predicted / measured: {arc:.4f} arc, {between:.4f} scans"
        )
        metrics += [
            {"name": "spacing_arc_ratio_median", "value": arc, "scope": group},
            {"name": "spacing_scan_ratio_median", "value": between, "scope": group},
        ]
        picks = rows[rng.choice(len(rows), sizes.fixture_spacing, replace=False)]
        arc, between = compute_arc_spacing_ratios(satellite, name, run.day, picks)
        fixture_metrics[f"spacing_arc_ratio_median_{group}"] = arc
        fixture_metrics[f"spacing_scan_ratio_median_{group}"] = between
        reduced.append(picks)
        indices.append(np.full(sizes.fixture_spacing, names.index(name)))
    if not run.write:
        return
    write_results(
        run,
        sizes.key,
        reference=sizes.reference,
        metrics=metrics,
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "footprints": np.array(
                [
                    [names.index(name), altitude, along, cross]
                    for name, altitude, along, cross in sizes.footprints
                ]
            ),
            "spacing": np.concatenate(reduced),
            "spacing_instruments": np.concatenate(indices),
        },
        fixture_metrics={"sample_size_within_3pct": within, **fixture_metrics},
    )


# the products of SSMIS and TMI, read as those of cross-track scanners


def observe_cone_passes(reference, points, cone_angle: float, radius: float):
    """
    Finds the passes on which the samples of a `scanner.Reference` observed points
    (in a cell of the lattice of samples within two scans of the sample nearest a
    point, among those within a radius, in km), and their times (when they cross the
    cone).
    """
    point_unit = up(points.geometry.x.to_numpy(), points.geometry.y.to_numpy())
    flat = reference.unit.reshape(-1, 3)
    scan_of = np.repeat(np.arange(len(reference.seconds)), reference.beams)
    rows = []
    for p, hits in enumerate(
        cKDTree(flat).query_ball_point(
            point_unit, r=2 * np.sin(radius / 6371 / 2), workers=-1
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
                reference.seconds,
                reference.position,
                ecef(points.geometry.x.iloc[p], points.geometry.y.iloc[p]),
                k,
                cone_angle,
                8,
            )
            if np.isfinite(crossing):
                rows.append({"point_id": p, "seconds": float(crossing)})
    return to_observed(reference.day, rows)


def validate_cone_swath(
    run,
    references: dict,
    orbit,
    tle,
    *,
    key: str,
    reference: str,
    scan_period: float,
    pass_radius: float,
    point_spacing: float,
    fixture_passes: int,
) -> None:
    """
    Scores the instruments of the swath configuration of a run (see
    `validation_io.parse_run`), each on the samples of its group of a
    `scanner.Reference` (in order), and writes the results if the run writes, with
    the fixture of the first instrument's passes. The metrics of several
    instruments are scoped by group.
    """
    entry = run.entry

    def satellite_for(time_step: float):
        return entry.to_satellite(
            orbit, configuration=key, time_step=time_step, platform=run.platform
        )

    satellite = satellite_for(scan_period)
    points = generate_points_fibonacci_lattice(point_spacing)
    grid = build_grid()
    metrics, first_observed = [], None
    for index, (group, lattice) in enumerate(references.items()):
        cone_angle = (
            entry.configurations[key].instruments[index].parameters["cone_angle"]
        )
        observed = observe_cone_passes(lattice, points, cone_angle, pass_radius)
        start = lattice.time[0].to_pydatetime()
        end = lattice.time[-1].to_pydatetime() + timedelta(seconds=scan_period)
        results = score_passes(
            points, satellite, index, observed, start, end, lattice.gaps
        )
        results["swath_width_median"] = score_chord_width(
            satellite.instruments[index],
            lattice.lon,
            lattice.lat,
            lattice.altitude,
            np.arange(0, len(lattice.seconds), 20),
        )
        first, last, compared, covered = scanner.observe_orbit(
            lattice, orbit, grid, 2 * pass_radius
        )
        results["footprint_iou_median"] = score_swath_area(
            satellite_for,
            index,
            first,
            last,
            (scan_period, 10, 30),
            grid,
            compared,
            covered,
        )
        print(group, json.dumps(results), flush=True)
        scope = {"scope": group} if len(references) > 1 else {}
        metrics += [
            {"name": name, "value": float(results[name]), **scope} for name in GATED
        ]
        if first_observed is None:
            first_observed = observed
    if not run.write:
        return
    passes = pick_cone_passes(points, first_observed, fixture_passes)
    write_results(
        run,
        key,
        reference=reference,
        metrics=metrics,
        fixture={
            "tle": np.array(tle),
            "day": np.array(run.day.isoformat()),
            "passes": passes,
        },
        fixture_metrics=recompute_passes(satellite, 0, run.day, passes),
    )


def fit_cone_scan(reference, bands: range, group: str | None = None) -> None:
    """
    Prints the fitted scan of a `scanner.Reference` (of a group, if given): the
    median angle of the samples from nadir (the cone angle), and the azimuths of
    the first and last samples from the velocity, in the inertial and Earth-fixed
    frames, with their center by latitude (in bands of latitude starting at each
    of `bands`, as wide as its step); for a group, also the azimuth step between
    samples and the median altitude.
    """
    lon, lat, _ = TO_GEODETIC.transform(*reference.position.T)
    lat = np.asarray(lat)
    nadir = -up(lon, lat)
    sight = normalize(
        ecef(reference.lon, reference.lat) - reference.position[:, np.newaxis]
    )
    cone = np.degrees(np.arccos(np.sum(sight * nadir[:, np.newaxis], axis=-1)))
    statistics = (
        f"{np.median(cone):.4f} deg, p5 {np.percentile(cone, 5):.4f}, p95 "
        f"{np.percentile(cone, 95):.4f}"
    )
    if group is None:
        print(f"cone angle: median {statistics}")
    else:
        print(f"{group} ({reference.beams} samples): cone angle median {statistics}")
    indent = "" if group is None else "  "
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
        by_latitude = ", ".join(
            f"{low:+d}..{low + bands.step:+d}: "
            f"{np.median(center[(lat >= low) & (lat < low + bands.step)]):+.3f}"
            for low in bands
            if np.any((lat >= low) & (lat < low + bands.step))
        )
        step = (
            ""
            if group is None
            else f"azimuth step {np.median(np.abs(np.diff(azimuth, axis=1))):.5f} deg; "
        )
        print(
            f"{indent}{frame}: first sample {first:+.4f} deg, last {last:+.4f} deg; "
            f"center {np.median(center):+.4f} (p5 {np.percentile(center, 5):+.4f}, "
            f"p95 {np.percentile(center, 95):+.4f}), half width "
            f"{np.median(np.abs(azimuth[:, -1] - azimuth[:, 0]) / 2):.4f} deg; "
            f"{step}center by latitude: {by_latitude}"
        )
    period = f"scan period: median {np.median(np.diff(reference.seconds)):.6f} s"
    if group is None:
        print(period)
    else:
        print(
            f"  {period}; altitude median {np.median(reference.altitude) / 1e3:.2f} km"
        )
