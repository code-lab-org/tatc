"""
Methods to perform coverage analysis of points (see `coverage_metrics` for the
aggregation, reduction, and gridding of observations).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import geometry as geo
from skyfield.api import wgs84
from skyfield.framelib import itrs

from ..constants import EARTH_POLAR_RADIUS, timescale
from ..schemas import (
    ConicalInstrument,
    GeneralPerturbationsOrbit,
    Instrument,
    Point,
    PointedInstrument,
    Satellite,
)
from ..utils.geometry import _get_point_coordinates, _hash_geometries
from ..utils.observation import (
    compute_max_access_time,
    compute_min_elevation_angle,
)
from ..utils.orbital import compute_apoapsis_radius
from ..utils.projection import (
    ViewGeometry,
    _compute_line_of_sight,
    _compute_nadir_axes,
    _compute_view_frame,
)
from ..utils.computation import TimeRequest, _run
from ..utils.time import _index_orbit_track, _to_time, _to_time_from_offsets
from .check import (
    _check_satellite,
    _check_satellites,
    _check_time_window,
    _combine_results,
    _get_instrument_indices,
    _is_single,
    _warn_low_perigees,
)
from .sampling import (
    _assemble_periods,
    _build_observation_frame,
    _find_crossings,
    _get_empty_coverage_frame,
)


def _get_elevation_angles(
    point: Point | geo.Point, satellite: Satellite, times: list[datetime]
) -> TimeRequest:
    """
    Get a satellite's elevation angles (degrees) from a point at times, as a
    computation (see `tatc.utils.computation.TimeRequest`).

    Args:
        point (Point | shapely.geometry.Point): The point.
        satellite (Satellite): The satellite.
        times (list[datetime.datetime]): The times.

    Returns:
        TimeRequest: the computation of the elevation angles (numpy.ndarray).
    """
    longitude, latitude, elevation = _get_point_coordinates(point)
    topos = wgs84.latlon(latitude, longitude, elevation)
    t_track, t_topos = _to_time(times), timescale.from_datetimes(times)
    yield [t_track, t_topos]
    orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track_at_time(t_track)
    return np.atleast_1d((orbit_track - topos.at(t_topos)).altaz()[0].degrees)


def _get_visible_interval_series(
    point: Point | geo.Point,
    satellite: Satellite,
    min_elevation_angle: float,
    max_altitude: float,
    start: datetime,
    end: datetime,
) -> TimeRequest:
    """
    Get the series of visible intervals based on altitude angle constraints,
    as a computation (see `tatc.utils.computation.TimeRequest`).

    Args:
        point (Point | shapely.geometry.Point): Point to observe: a TAT-C
                point or a shapely point (longitude, latitude, and optional
                elevation in meters).
        satellite (Satellite): Satellite doing the observation.
        min_elevation_angle (float): Minimum elevation angle (degrees) for valid observation.
        max_altitude (float): A conservative upper-bound satellite altitude
                (meters, e.g. the orbit's apogee altitude), used only to
                compute a generously large `max_access_time` bound for
                matching rise events to their corresponding set events.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.

    Returns:
        TimeRequest: the computation of the series of observation intervals.
    """
    # compute the maximum access time to filter bad data
    max_access_time = timedelta(
        seconds=compute_max_access_time(max_altitude, min_elevation_angle)
    )
    # find the set of observation events
    orbit = satellite.orbit.to_gp_orbit()
    # a computation internal to TATC (see get_observation_events)
    times, events = (
        yield from orbit._get_observation_events(  # pylint: disable=protected-access
            point, start, end, min_elevation_angle
        )
    )

    # build the observation periods
    window_start = pd.Timestamp(start.astimezone(tz=timezone.utc))
    window_end = pd.Timestamp(end.astimezone(tz=timezone.utc))
    obs_periods = []
    if len(events) == 0:
        # no rise, culminate, or set event was captured in [start, end]. This
        # means the elevation angle never crossed min_elevation_angle and had
        # no interior local maximum in this window -- which happens both
        # when the point is never visible, and when [start, end] falls
        # entirely within a longer visible pass (no rise/set inside the
        # window, and the window is too narrow, or off-center, to contain
        # the pass's culmination). Disambiguate by sampling the true
        # elevation angle at the window's midpoint.
        (elevation_angle,) = yield from _get_elevation_angles(
            point, satellite, [start + (end - start) / 2]
        )
        if elevation_angle > min_elevation_angle:
            # continuously visible for the entire window
            obs_periods += [pd.Interval(left=window_start, right=window_end)]
        return pd.Series(obs_periods, dtype="interval")
    # Skyfield's find_events can miss the rise or set of a grazing pass,
    # whose maximum elevation barely exceeds min_elevation_angle, so a set
    # without a rise (or a rise without a set) does not by itself mean that
    # the point is visible at the start (or end) of the window: check the
    # elevation angle there, and otherwise bound the period by the maximum
    # access time (which, a period being refined afterward, need not be
    # tight; an unbounded period could span several passes)
    visible_at_start, visible_at_end = (
        angle > min_elevation_angle
        for angle in (yield from _get_elevation_angles(point, satellite, [start, end]))
    )
    utc = pd.DatetimeIndex(np.atleast_1d(times.utc_datetime()))
    if np.all(events == 1):
        # only culminations: the window lies within one pass (visible
        # throughout) or the culminations are of grazing passes
        if visible_at_start and visible_at_end:
            obs_periods += [pd.Interval(left=window_start, right=window_end)]
        else:
            for culmination in utc:
                period = pd.Interval(
                    left=max(window_start, culmination - max_access_time),
                    right=min(window_end, culmination + max_access_time),
                )
                if obs_periods and period.left <= obs_periods[-1].right:
                    # merge overlapping periods
                    period = pd.Interval(
                        left=obs_periods.pop().left, right=period.right
                    )
                obs_periods.append(period)
        return pd.Series(obs_periods, dtype="interval")
    # otherwise, match rise/set events (converted to UTC timestamps once,
    # as converting Skyfield times is costly)
    rises = utc[events == 0]
    sets = utc[events == 2]
    if len(sets) > 0 and (len(rises) == 0 or sets[0] < rises[0]) and start < sets[0]:
        # if first event is a set, create a period from the start (or, if
        # its rise was missed, from the maximum access time before the set)
        obs_periods += [
            pd.Interval(
                left=(
                    window_start
                    if visible_at_start
                    else max(window_start, sets[0] - max_access_time)
                ),
                right=sets[0],
            )
        ]
    # create an observation period to match with each rise event if
    # there is a following set event within twice the maximum access time
    # (the events are in time order, so the first set after each rise)
    next_sets = sets.searchsorted(rises, side="right")
    obs_periods += [
        pd.Interval(left=rise, right=sets[i])
        for rise, i in zip(rises, next_sets)
        if i < len(sets) and sets[i] < rise + 2 * max_access_time
    ]
    if len(rises) > 0 and (len(sets) == 0 or rises[-1] > sets[-1]) and rises[-1] < end:
        # if last event is a rise, create a period to the end (or, if its
        # set was missed, to the maximum access time after the rise)
        obs_periods += [
            pd.Interval(
                left=rises[-1],
                right=(
                    window_end
                    if visible_at_end
                    else min(window_end, rises[-1] + max_access_time)
                ),
            )
        ]
    return pd.Series(obs_periods, dtype="interval")


def compute_access_periods(
    point: Point | geo.Point,
    satellite: Satellite,
    start: datetime,
    end: datetime,
    min_elevation_angle: float = 0,
) -> pd.Series:
    """
    Compute the periods when a satellite is in view of a point: when its
    elevation angle, seen from the point, is at least a minimum. The rise
    and set times are refined to a millisecond (see
    `GeneralPerturbationsOrbit.get_observation_events`), and periods in
    progress at the start or end of the analysis period are truncated to it.
    For the periods when an instrument observes the point, see
    `collect_observations`.

    Args:
        point (Point | shapely.geometry.Point): The point: a TAT-C point or a
                shapely point (longitude, latitude, and optional elevation in
                meters).
        satellite (Satellite): The satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        min_elevation_angle (float): The minimum elevation angle (degrees).

    Returns:
        pandas.Series: the access periods (`pandas.Interval` of UTC
            timestamps), in time order.
    """
    _warn_low_perigees([_check_satellite(satellite)], stacklevel=3)
    _check_time_window(start, end)
    return _run(
        _compute_access_periods(point, satellite, start, end, min_elevation_angle)
    )


def _compute_access_periods(
    point: Point | geo.Point,
    satellite: Satellite,
    start: datetime,
    end: datetime,
    min_elevation_angle: float = 0,
) -> TimeRequest:
    """
    Compute the periods when a satellite is in view of a point, as a
    computation (see `compute_access_periods` and
    `tatc.utils.computation.TimeRequest`).
    """
    _, _, elevation = _get_point_coordinates(point)
    # use the apogee altitude above the polar radius (and above the point,
    # if below the ellipsoid) as a conservative upper bound for pairing rise
    # and set events
    max_altitude = (
        compute_apoapsis_radius(
            satellite.orbit.get_semimajor_axis(), satellite.orbit.get_eccentricity()
        )
        - EARTH_POLAR_RADIUS
        - min(elevation, 0)
    )
    return (
        yield from _get_visible_interval_series(
            point, satellite, min_elevation_angle, max_altitude, start, end
        )
    )


_SCREEN_STEP = 60.0
"""
The maximum time step (seconds) of the coarse orbit track from which the
passes of many points are screened together (see `_collect_batch_observations`).
"""
_SCREEN_MARGIN = 1.0
"""
The angle (degrees) added to the radius about the sub-satellite point within
which points are screened, for the curvature of the sub-satellite path between
samples and the difference between the geodetic and geocentric verticals.
"""
_PART_DURATION = 600.0
"""
The maximum duration (seconds) of the parts into which a screened window is
divided, so that the residual has few extrema within each part.
"""
_PART_SAMPLES = 11
"""The number of samples of the residual in each part of a screened window."""
_SCREEN_BLOCK = 64
"""The number of time steps screened against all points at once."""
_BLOCK_POINT_DAYS = 10000
"""
The number of point-days (points times days of the analysis period) whose
observations are refined together, which bounds the memory of the vectorized
stages.
"""


def _get_screening_angle(
    max_radius: float, min_radius: float, half_angle: float
) -> float:
    """
    Get the Earth central angle (degrees) from the sub-satellite point within
    which a point can lie within a field of regard: on a sphere of the
    minimum radius (the smallest distance of a point from the Earth's
    center), seen from the maximum radius (of the satellite), at a nadir
    angle of at most the half angle, and above the horizon. The angle grows
    with the satellite's radius and decreases with the point's, so it bounds
    the angle of every point and satellite position.

    Args:
        max_radius (float): The maximum radius (meters) of the satellite.
        min_radius (float): The minimum radius (meters) of the points.
        half_angle (float): The half angle (degrees) of the field of regard.

    Returns:
        float: the Earth central angle (degrees).
    """
    horizon = np.arccos(min_radius / max_radius)
    if half_angle >= 90:
        return float(np.degrees(horizon))
    ratio = max_radius / min_radius * np.sin(np.radians(half_angle))
    if ratio >= 1:
        return float(np.degrees(horizon))
    return float(np.degrees(min(np.arcsin(ratio) - np.radians(half_angle), horizon)))


def _screen_windows(
    directions: np.ndarray, units: np.ndarray, radius: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Screen the windows in which points may be observed: the runs of
    consecutive time steps in which a point lies within an Earth central
    angle of the sub-satellite path. The path between two samples is bounded
    by the circle about their mean direction whose radius is half of the
    angle between them, so a point is screened in a step if its angle from
    the mean direction is at most the radius plus that half angle.

    Args:
        directions (numpy.ndarray): The unit directions of the satellite from
            the Earth's center at each sample, shape (3, M).
        units (numpy.ndarray): The unit directions of the points from the
            Earth's center, shape (3, N).
        radius (numpy.ndarray): The Earth central angle (radians) of each
            point, shape (N,).

    Returns:
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: the point, first
            step, and last step of each window (step `i` spans samples `i`
            and `i + 1`), sorted by point and first step.
    """
    middle = directions[:, :-1] + directions[:, 1:]
    middle = middle / np.linalg.norm(middle, axis=0)
    half_step = (
        np.arccos(
            np.clip(np.sum(directions[:, :-1] * directions[:, 1:], axis=0), -1, 1)
        )
        / 2
    )
    points, steps = [], []
    for block in range(0, middle.shape[1], _SCREEN_BLOCK):
        m = middle[:, block : block + _SCREEN_BLOCK]
        # element-wise, so that each point's result does not depend on the others
        cosine = (
            units[0][:, None] * m[0]
            + units[1][:, None] * m[1]
            + units[2][:, None] * m[2]
        )
        hits = cosine >= np.cos(
            np.minimum(
                np.pi, radius[:, None] + half_step[block : block + _SCREEN_BLOCK]
            )
        )
        point, step = np.nonzero(hits)
        points.append(point)
        steps.append(step + block)
    point, step = np.concatenate(points), np.concatenate(steps)
    order = np.lexsort((step, point))
    point, step = point[order], step[order]
    if len(point) == 0:
        return point, step, step
    first = np.ones(len(point), dtype=bool)
    first[1:] = (point[1:] != point[:-1]) | (step[1:] != step[:-1] + 1)
    last = np.ones(len(point), dtype=bool)
    last[:-1] = first[1:]
    return point[first], step[first], step[last]


def _find_hidden_minima(
    evaluate: Callable[[np.ndarray, np.ndarray], TimeRequest],
    samples: np.ndarray,
    values: np.ndarray,
    index: np.ndarray,
) -> TimeRequest:
    """
    Find the minima of a residual between its samples that may be negative
    although every sample near them is positive (as when a point grazes a
    narrow field of regard between samples), as a computation (see
    `tatc.utils.computation.TimeRequest`). Each local minimum of the samples
    of a part that is positive but no greater than the larger change of the
    residual to its neighboring samples (which bounds how far the residual
    can dip between them) is refined by golden-section search between its
    neighbors, stopping when a negative residual is found.

    Args:
        evaluate (Callable[[numpy.ndarray, numpy.ndarray], TimeRequest]): The
            computation of the residual at times (seconds) for the indices.
        samples (numpy.ndarray): The sample times (seconds) of each part,
            shape (P, K).
        values (numpy.ndarray): The residual at the samples, shape (P, K).
        index (numpy.ndarray): The index of the point of each part, shape (P,).

    Returns:
        TimeRequest: the computation of the part, time (seconds), and
            residual of each negative minimum found (tuple[numpy.ndarray,
            numpy.ndarray, numpy.ndarray]).
    """
    padded = np.pad(values, ((0, 0), (1, 1)), constant_values=np.inf)
    left, right = padded[:, :-2], padded[:, 2:]
    change = np.maximum(
        np.where(np.isfinite(left), left - values, 0),
        np.where(np.isfinite(right), right - values, 0),
    )
    part, k = np.nonzero(
        (values > 0) & (values <= left) & (values <= right) & (values <= change)
    )
    if len(part) == 0:
        empty = np.array([], dtype=float)
        return np.array([], dtype=int), empty, empty
    count = samples.shape[1]
    lower = samples[part, np.maximum(k - 1, 0)]
    upper = samples[part, np.minimum(k + 1, count - 1)]
    point = index[part]
    ratio = (np.sqrt(5) - 1) / 2
    c, d = upper - ratio * (upper - lower), lower + ratio * (upper - lower)
    f_c = yield from evaluate(c, point)
    f_d = yield from evaluate(d, point)
    found_time = np.where(f_c <= f_d, c, d)
    found_value = np.minimum(f_c, f_d)
    for _ in range(40):
        active = (found_value > 0) & (upper - lower > 1e-3)
        if not np.any(active):
            break
        keep_left = f_c <= f_d
        # shrink each bracket toward its smaller interior value
        upper = np.where(active & keep_left, d, upper)
        lower = np.where(active & ~keep_left, c, lower)
        x = np.where(
            keep_left, upper - ratio * (upper - lower), lower + ratio * (upper - lower)
        )
        f_x = np.full(len(x), np.inf)
        f_x[active] = yield from evaluate(x[active], point[active])
        c, f_c, d, f_d = (
            np.where(active & keep_left, x, np.where(active, d, c)),
            np.where(active & keep_left, f_x, np.where(active, f_d, f_c)),
            np.where(active & keep_left, c, np.where(active, x, d)),
            np.where(active & keep_left, f_c, np.where(active, f_x, f_d)),
        )
        better = f_x < found_value
        found_time = np.where(better, x, found_time)
        found_value = np.where(better, f_x, found_value)
    negative = found_value <= 0
    return part[negative], found_time[negative], found_value[negative]


def _merge_window_periods(
    parts: list[list[tuple[float, float]]], window: np.ndarray
) -> list[list[tuple[float, float]]]:
    """
    Merge the periods of the consecutive parts of each window that meet at a
    boundary between parts.

    Args:
        parts (list[list[tuple[float, float]]]): The periods (seconds) of each
            part, in order of window and time.
        window (numpy.ndarray): The window of each part.

    Returns:
        list[list[tuple[float, float]]]: the periods of each window.
    """
    windows: list[list[tuple[float, float]]] = [
        [] for _ in range(int(window.max()) + 1 if len(window) else 0)
    ]
    for periods, w in zip(parts, window):
        merged = windows[w]
        for left, right in periods:
            if merged and merged[-1][1] == left:
                merged[-1] = (merged[-1][0], right)
            else:
                merged.append((left, right))
    return windows


def _refine_window_periods(
    residual: Callable[[np.ndarray, np.ndarray], TimeRequest],
    window_point: np.ndarray,
    window_lower: np.ndarray,
    window_upper: np.ndarray,
) -> TimeRequest:
    """
    Refine the periods in which a residual (of a time and a point) is not
    positive within screened windows, as a computation (see
    `tatc.utils.computation.TimeRequest`). Each window is divided into parts
    of at most `_PART_DURATION`, each part is sampled at `_PART_SAMPLES`
    times, minima of the residual between samples are found (see
    `_find_hidden_minima`), each change of sign is refined to a millisecond
    (see `tatc.analysis.sampling._find_crossings`), and the periods of
    consecutive parts that meet are merged.

    Args:
        residual (Callable[[numpy.ndarray, numpy.ndarray], TimeRequest]): The
            computation of the residual at times (seconds) for the points.
        window_point (numpy.ndarray): The point of each window.
        window_lower (numpy.ndarray): The start of each window (seconds).
        window_upper (numpy.ndarray): The end of each window (seconds).

    Returns:
        TimeRequest: the computation of the point, start, and end (seconds)
            of each period (tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]).
    """
    counts = np.maximum(
        1, np.ceil((window_upper - window_lower) / _PART_DURATION).astype(int)
    )
    window = np.repeat(np.arange(len(counts)), counts)
    offset = np.arange(len(window)) - np.repeat(np.cumsum(counts) - counts, counts)
    width = (window_upper - window_lower)[window] / counts[window]
    # each part ends exactly where the next begins, so that their periods merge
    part_lower = window_lower[window] + offset * width
    part_upper = np.where(
        offset == counts[window] - 1,
        window_upper[window],
        window_lower[window] + (offset + 1) * width,
    )
    part_point = window_point[window]
    samples = part_lower[:, None] + (part_upper - part_lower)[:, None] * np.linspace(
        0, 1, _PART_SAMPLES
    )
    samples[:, -1] = part_upper
    values = np.reshape(
        (yield from residual(samples.ravel(), np.repeat(part_point, _PART_SAMPLES))),
        samples.shape,
    )
    part, dip_time, dip_value = yield from _find_hidden_minima(
        residual, samples, values, part_point
    )
    sample_list, value_list = list(samples), list(values)
    for p in np.unique(part):
        times = np.concatenate([sample_list[p], dip_time[part == p]])
        order = np.argsort(times, kind="stable")
        sample_list[p] = times[order]
        value_list[p] = np.concatenate([value_list[p], dip_value[part == p]])[order]
    brackets = [
        (i, j)
        for i, value in enumerate(value_list)
        for j in np.flatnonzero((value[:-1] <= 0) != (value[1:] <= 0))
    ]
    crossings = {}
    if len(brackets) > 0:
        bracket_point = part_point[[i for i, _ in brackets]]
        crossing, _ = yield from _find_crossings(
            lambda seconds, index: residual(seconds, bracket_point[index]),
            np.array([sample_list[i][j] for i, j in brackets]),
            np.array([sample_list[i][j + 1] for i, j in brackets]),
        )
        crossings = dict(zip(brackets, crossing))
    window_periods = _merge_window_periods(
        _assemble_periods(sample_list, value_list, crossings), window
    )
    period_point = np.array(
        [window_point[w] for w, periods in enumerate(window_periods) for _ in periods],
        dtype=int,
    )
    lower = np.array([left for periods in window_periods for left, _ in periods])
    upper = np.array([right for periods in window_periods for _, right in periods])
    return period_point, lower, upper


def _get_batch_epochs(
    orbit: GeneralPerturbationsOrbit,
    instrument: Instrument,
    xyz: np.ndarray,
    period_point: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    reference: pd.Timestamp,
) -> TimeRequest:
    """
    Get the observation epochs of periods, as a computation (see
    `tatc.utils.computation.TimeRequest`): the time a pointed instrument's
    view sweeps over the point (see `_get_batch_view_crossing_times`), the
    times the point crosses a conical instrument's cone (see
    `_get_batch_cone_crossing_times`), or otherwise the period's midpoint.

    Args:
        orbit (GeneralPerturbationsOrbit): The orbit.
        instrument (Instrument): The observing instrument.
        xyz (numpy.ndarray): The Earth-fixed coordinates (meters) of the
            points, shape (3, N).
        period_point (numpy.ndarray): The point of each period.
        lower (numpy.ndarray): The start of each period (seconds).
        upper (numpy.ndarray): The end of each period (seconds).
        reference (pandas.Timestamp): The reference time of the seconds.

    Returns:
        TimeRequest: the computation of the point, period, and time
            (seconds) of each epoch (tuple[numpy.ndarray, numpy.ndarray,
            numpy.ndarray]).
    """
    if isinstance(instrument, PointedInstrument):
        epochs = yield from _get_batch_view_crossing_times(
            orbit, instrument, xyz, period_point, lower, upper, reference
        )
        return period_point, np.arange(len(lower)), epochs
    if isinstance(instrument, ConicalInstrument):
        epoch_period, epochs = yield from _get_batch_cone_crossing_times(
            orbit, instrument, xyz, period_point, lower, upper, reference
        )
        return period_point[epoch_period], epoch_period, epochs
    return period_point, np.arange(len(lower)), (lower + upper) / 2


def _collect_batch_observations(
    points: list[Point | geo.Point],
    satellite: Satellite,
    start: datetime,
    end: datetime,
    instrument_index: int = 0,
    omit_solar: bool = True,
) -> TimeRequest:
    """
    Collect single satellite observations of many points together, as a
    computation (see `collect_observations` and
    `tatc.utils.computation.TimeRequest`). Each point's observations do not
    depend on the other points (to within a microsecond), but each stage is
    vectorized across the points:

    1. the satellite's orbit is propagated once, at steps of at most
       `_SCREEN_STEP`, and the windows in which each point may lie within the
       field of regard are screened from the angle between the point and
       the sub-satellite path (see `_screen_windows`);
    2. each window is divided into parts, each part is sampled, and the
       periods in which the point lies within the field of regard (and the
       satellite is above the minimum elevation angle) are refined from the
       changes of sign between samples and from minima between them (see
       `_find_hidden_minima`); and
    3. the epochs (view or cone crossings) and validity of the observations
       are computed with one target for each time.
    """
    instrument = satellite.instruments[instrument_index]
    orbit = satellite.orbit.to_gp_orbit()
    coordinates = np.array([_get_point_coordinates(point) for point in points])
    longitude, latitude, elevation = coordinates.T
    xyz = np.array(wgs84.latlon(latitude, longitude, elevation).itrs_xyz.m)
    up = np.array(
        [
            np.cos(np.radians(latitude)) * np.cos(np.radians(longitude)),
            np.cos(np.radians(latitude)) * np.sin(np.radians(longitude)),
            np.sin(np.radians(latitude)),
        ]
    )
    geometries = [geo.Point(x, y, z) for x, y, z in coordinates]
    # use the apogee radius and the polar radius (below a point below the
    # ellipsoid) as conservative bounds of the satellite's and points' radii
    max_radius = compute_apoapsis_radius(
        satellite.orbit.get_semimajor_axis(), satellite.orbit.get_eccentricity()
    )
    # each point's own bounds, so that its result does not depend on the others
    min_radius = EARTH_POLAR_RADIUS + np.minimum(elevation, 0)
    min_elevation_angle = np.array(
        [
            max(
                0.0,
                compute_min_elevation_angle(
                    max_radius - radius, instrument.field_of_regard
                )
                - 1.0,
            )
            for radius in min_radius
        ]
    )
    half_angle = instrument.field_of_regard / 2
    reference = pd.Timestamp(start.astimezone(tz=timezone.utc))
    duration = (end - start).total_seconds()
    _warn_if_propagating_with_drag(orbit, start, end)

    def residual(seconds: np.ndarray, index: np.ndarray) -> TimeRequest:
        # positive outside the field of regard or below the minimum elevation angle
        t = _to_time_from_offsets(reference, seconds)
        yield t
        p_m, nadir, _, _ = _compute_nadir_axes(
            orbit.get_orbit_track_at_time(t),
            nadir_reference=instrument.nadir_reference,
        )
        los = _compute_line_of_sight(xyz[:, index], p_m)
        distance = np.linalg.norm(los, axis=0)
        elevation_angle = np.degrees(
            np.arcsin(np.clip(-np.sum(los * up[:, index], axis=0) / distance, -1, 1))
        )
        excess = min_elevation_angle[index] - elevation_angle
        if half_angle >= 90:
            return excess
        cone = np.degrees(
            np.arccos(np.clip(np.sum(los * nadir, axis=0) / distance, -1, 1))
        )
        return np.maximum(cone - half_angle, excess)

    # 1. the coarse orbit track, from which the windows of each point are screened
    steps = max(1, int(np.ceil(duration / _SCREEN_STEP)))
    screen_times = np.linspace(0, duration, steps + 1)
    t = _to_time_from_offsets(reference, screen_times)
    yield t
    position = np.array(orbit.get_orbit_track_at_time(t).frame_xyz(itrs).m)
    directions = position / np.linalg.norm(position, axis=0)
    units = xyz / np.linalg.norm(xyz, axis=0)
    screening_radius = np.radians(
        np.array(
            [
                _get_screening_angle(max_radius, radius, half_angle)
                for radius in min_radius
            ]
        )
        + _SCREEN_MARGIN
    )

    def collect_block(block: np.ndarray) -> TimeRequest:
        """Collect the observations of a block of points (global indices)."""
        local_point, first_step, last_step = _screen_windows(
            directions, units[:, block], screening_radius[block]
        )
        if len(local_point) == 0:
            return _get_empty_coverage_frame(omit_solar)
        window_point = block[local_point]
        window_lower = screen_times[first_step]
        window_upper = screen_times[last_step + 1]

        # 2. refine the periods in the parts of each window
        period_point, lower, upper = yield from _refine_window_periods(
            residual, window_point, window_lower, window_upper
        )
        # 3. the epochs and validity of the observations
        epoch_point, epoch_period, epochs = yield from _get_batch_epochs(
            orbit, instrument, xyz, period_point, lower, upper, reference
        )
        keep = (
            upper[epoch_period] - lower[epoch_period]
            >= instrument.min_access_time.total_seconds()
        )
        epoch_point, epoch_period, epochs = (
            epoch_point[keep],
            epoch_period[keep],
            epochs[keep],
        )
        observations, orbit_track, valid_point = [], None, np.array([], dtype=int)
        if len(epochs) > 0:
            t = _to_time_from_offsets(reference, epochs)
            yield t
            orbit_track = orbit.get_orbit_track_at_time(t)
            target = wgs84.latlon(
                latitude[epoch_point], longitude[epoch_point], elevation[epoch_point]
            )
            valid = np.atleast_1d(instrument.is_valid_observation(orbit_track, target))
            if isinstance(instrument, (PointedInstrument, ConicalInstrument)):
                valid = valid & np.atleast_1d(
                    instrument.is_in_field_of_view(orbit_track, target)
                )
            valid_index = np.flatnonzero(valid)
            valid_point = epoch_point[valid_index]
            observations = [
                (
                    pd.Interval(
                        left=reference
                        + pd.Timedelta(seconds=float(lower[epoch_period[i]])),
                        right=reference
                        + pd.Timedelta(seconds=float(upper[epoch_period[i]])),
                    ),
                    reference + pd.Timedelta(seconds=float(epochs[i])),
                    tuple(coordinates[epoch_point[i]]),
                )
                for i in valid_index
            ]
            if len(valid_index) < len(epochs):
                orbit_track = _index_orbit_track(orbit_track, valid_index)
        observed = np.unique(valid_point)
        hashes = dict(
            zip(observed, _hash_geometries([geometries[p] for p in observed]))
        )
        return _build_observation_frame(
            observations,
            [hashes[p] for p in valid_point],
            [geometries[p] for p in valid_point],
            satellite,
            instrument,
            omit_solar,
            orbit_track,
        )

    # points in blocks of at most _BLOCK_POINT_DAYS point-days, which bound the
    # size of the vectorized stages
    size = max(1, int(_BLOCK_POINT_DAYS * 86400 / max(duration, 1.0)))
    frames = []
    for first in range(0, len(points), size):
        frame = yield from collect_block(
            np.arange(first, min(len(points), first + size))
        )
        if len(frame) > 0:
            frames.append(frame)
    if len(frames) == 0:
        return _get_empty_coverage_frame(omit_solar)
    return frames[0] if len(frames) == 1 else pd.concat(frames, ignore_index=True)


def _warn_if_propagating_with_drag(
    orbit: GeneralPerturbationsOrbit, start: datetime, end: datetime
) -> None:
    """
    Warn if an orbit is propagated with drag between its repeat tracks (as in
    `GeneralPerturbationsOrbit.get_observation_events`).

    Args:
        orbit (GeneralPerturbationsOrbit): The orbit.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
    """
    first, last = orbit._get_repeat_tracks()  # pylint: disable=protected-access
    direct_start = start if first is None else max(start, first.epoch)
    direct_end = end if last is None else min(end, last.epoch)
    if direct_start < direct_end:
        orbit._warn_if_propagating_with_drag(  # pylint: disable=protected-access
            _to_time([direct_start, direct_end]),
            np.array(orbit.get_closest_element_index([direct_start, direct_end])),
        )


def _get_batch_view_crossing_times(
    orbit: GeneralPerturbationsOrbit,
    instrument: PointedInstrument,
    xyz: np.ndarray,
    point: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    reference: pd.Timestamp,
) -> TimeRequest:
    """
    Get the time in each period when a pointed instrument's view sweeps over
    its point: when the point crosses the plane of the view's boresight and
    cross-track axis (where its along-track view angle is zero), or, for a
    `ViewGeometry.SCAN` view, the scan plane. The along-track component of
    the unit line of sight, which is defined even when the point is behind
    the view, decreases through zero as the satellite passes the point. If a
    period does not contain a crossing (for example, a period truncated by
    the analysis window), uses the period's end closest to the crossing.

    Args:
        orbit (GeneralPerturbationsOrbit): The orbit.
        instrument (PointedInstrument): The observing instrument.
        xyz (numpy.ndarray): The Earth-fixed coordinates (meters) of the
            points, shape (3, N).
        point (numpy.ndarray): The point of each period.
        lower (numpy.ndarray): The start of each period (seconds).
        upper (numpy.ndarray): The end of each period (seconds).
        reference (pandas.Timestamp): The reference time of the seconds.

    Returns:
        TimeRequest: the computation of the crossing time (seconds) in each
            period (numpy.ndarray).
    """
    if len(lower) == 0:
        return np.array([], dtype=float)

    def residual(seconds: np.ndarray, index: np.ndarray) -> TimeRequest:
        # along-track component of the unit line of sight in the view frame
        t = _to_time_from_offsets(reference, seconds)
        yield t
        orbit_track = orbit.get_orbit_track_at_time(t)
        if instrument.view_geometry == ViewGeometry.SCAN:
            # the target crosses the scan plane (tilted by the pitch angle)
            position, _, along, _ = _compute_view_frame(
                orbit_track,
                0,
                0,
                instrument.velocity_frame,
                instrument.nadir_reference,
                instrument.get_pitch_angle(orbit_track),
            )
        else:
            position, _, along, _ = _compute_view_frame(
                orbit_track,
                instrument.get_roll_angle(orbit_track),
                instrument.get_pitch_angle(orbit_track),
                instrument.velocity_frame,
                instrument.nadir_reference,
            )
        los = _compute_line_of_sight(xyz[:, point[index]], position)
        return np.sum(los * along, axis=0) / np.linalg.norm(los, axis=0)

    crossing, _ = yield from _find_crossings(residual, lower, upper)
    return crossing


def _get_batch_cone_crossing_times(
    orbit: GeneralPerturbationsOrbit,
    instrument: ConicalInstrument,
    xyz: np.ndarray,
    point: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    reference: pd.Timestamp,
) -> TimeRequest:
    """
    Get the times in each period when its point crosses a conical
    instrument's cone: when its angle from nadir equals the cone angle.
    During a pass, this angle decreases to a minimum near the closest
    approach and then increases, so there are up to two crossings, entering
    the cone (ahead of the satellite) and leaving it (behind).

    Args:
        orbit (GeneralPerturbationsOrbit): The orbit.
        instrument (ConicalInstrument): The observing instrument.
        xyz (numpy.ndarray): The Earth-fixed coordinates (meters) of the
            points, shape (3, N).
        point (numpy.ndarray): The point of each period.
        lower (numpy.ndarray): The start of each period (seconds).
        upper (numpy.ndarray): The end of each period (seconds).
        reference (pandas.Timestamp): The reference time of the seconds.

    Returns:
        TimeRequest: the computation of the period and time (seconds) of
            each crossing, in order of period (tuple[numpy.ndarray,
            numpy.ndarray]).
    """
    n = len(lower)
    if n == 0:
        return np.array([], dtype=int), np.array([], dtype=float)

    def cone_angle(seconds: np.ndarray, index: np.ndarray) -> TimeRequest:
        t = _to_time_from_offsets(reference, seconds)
        yield t
        p_m, nadir, _, _ = _compute_nadir_axes(
            orbit.get_orbit_track_at_time(t),
            instrument.velocity_frame,
            instrument.nadir_reference,
        )
        los = _compute_line_of_sight(xyz[:, point[index % n]], p_m)
        return np.degrees(
            np.arccos(
                np.clip(
                    np.sum(los * nadir, axis=0) / np.linalg.norm(los, axis=0), -1, 1
                )
            )
        )

    def residual(seconds: np.ndarray, index: np.ndarray) -> TimeRequest:
        return (yield from cone_angle(seconds, index)) - instrument.cone_angle

    # time of the minimum angle from nadir in each period, from sampled times
    samples = lower[:, None] + (upper - lower)[:, None] * np.linspace(0, 1, 21)
    closest = samples[
        np.arange(n),
        np.argmin(
            (
                yield from cone_angle(samples.ravel(), np.repeat(np.arange(n), 21))
            ).reshape(samples.shape),
            axis=1,
        ),
    ]
    crossing, bracketed = yield from _find_crossings(
        residual, np.concatenate([lower, closest]), np.concatenate([closest, upper])
    )
    # the crossings of each period, entering (k = 0) and leaving (k = 1) the cone
    period = np.concatenate([np.arange(n), np.arange(n)])
    order = np.lexsort((np.repeat([0, 1], n), period))
    keep = order[bracketed[order]]
    return period[keep], crossing[keep]


def collect_observations(
    points: Point | geo.Point | list[Point | geo.Point],
    satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
    instrument_index: int | None = 0,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect observations of one or more geodetic points of interest by one
    or more satellites. Each point is a TAT-C `Point` or a shapely `Point`
    (whose x, y, and optional z coordinates are its longitude, latitude, and
    elevation in meters); TAT-C points are expected to be replaced by
    shapely points in the future. Observations record the point (with its
    elevation as its z coordinate) as their geometry, identified by its
    hash, `target_hash` (see `tatc.utils.geometry.hash_geometry`); a TAT-C
    point's `id` is not recorded. For a region, see
    `tatc.analysis.region_sampling.collect_region_observations`.

    Each observation spans a period when the point lies within the
    instrument's field of regard (when its angle from nadir is at most half
    the field of regard). Its epoch is the period's midpoint or, for
    a `PointedInstrument`, the time when the instrument's view sweeps over
    the point (when the point's along-track view angle relative to the view
    center is zero), or, for a `ConicalInstrument`, a time when the point
    crosses the scanned cone (up to two per period, entering and leaving the
    cone); at that time, the point must lie within the field of view.

    If the orbit is propagated with a repeat cycle (see
    `GeneralPerturbationsOrbit.repeat_cycle`), it is modeled as maintained on
    its repeat ground track before its first and after its last element's
    epoch (see `GeneralPerturbationsOrbit.get_orbit_track_at_time`).

    The observations of all points by each satellite and instrument are
    computed together, with each stage vectorized across the points: the
    orbit is propagated once to screen the windows in which each point may
    be observed, and the periods, epochs, and validity of the observations
    are then refined to a millisecond in all windows at once. The
    observations of a point do not depend on the other points (to within a
    microsecond), so many points are much faster to collect together than
    separately.

    Args:
        points (Point | shapely.geometry.Point | list[Point | shapely.geometry.Point]):
                The ground point(s) of interest.
        satellites (Satellite | list[Satellite]): The observing satellite(s).
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        instrument_index (int | None): The index of the observing instrument
                in each satellite, or `None` for every instrument of each
                satellite.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations:
            for a single point, satellite, and instrument (an integer
            `instrument_index`), in time order; otherwise, those of each
            point, satellite, and instrument, concatenated and sorted by
            start time.
    """
    single = _is_single(points, satellites, instrument_index=instrument_index)
    satellites = _check_satellites(satellites)
    _check_time_window(start, end)
    points = list(points) if isinstance(points, (list, tuple)) else [points]
    for point in points:
        if not isinstance(point, (Point, geo.Point)):
            raise TypeError(
                "point must be a Point or shapely Point, not a "
                f"{type(point).__name__} (see collect_region_observations for a region)"
            )
    # each computation is vectorized across the points, so they are run in
    # turn (running them together would only enlarge each step)
    gdfs = [
        _run(
            _collect_batch_observations(
                points, satellite, start, end, index, omit_solar
            )
        )
        for satellite in satellites
        for index in _get_instrument_indices(satellite, instrument_index)
    ]
    return _combine_results(
        gdfs,
        single,
        "start",
        lambda: _get_empty_coverage_frame(omit_solar),
        stable=False,
    )


def collect_multi_observations(
    point: Point | geo.Point,
    satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect multiple satellite observations of a geodetic point of interest
    by every instrument on every satellite in `satellites`: equivalent to
    `collect_observations` with `instrument_index=None`, kept for backwards
    compatibility.

    Args:
        point (Point | shapely.geometry.Point): The ground point of interest
                (see `collect_observations`).
        satellites (Satellite | list[Satellite]): The observing satellite(s),
                each contributing an observation per instrument it carries.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with all recorded observations.
    """
    return collect_observations(point, satellites, start, end, None, omit_solar)
