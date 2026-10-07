"""
Methods to perform coverage analysis of points.

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
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ..constants import EARTH_POLAR_RADIUS, de421, timescale
from ..schemas import (
    ConicalInstrument,
    GeneralPerturbationsOrbit,
    Instrument,
    Point,
    PointedInstrument,
    Satellite,
)
from ..utils.geometry import _get_point_coordinates
from ..utils.observation import (
    compute_max_access_time,
    compute_min_elevation_angle,
)
from ..utils.orbital import compute_apoapsis_radius
from ..utils.projection import (
    ViewGeometry,
    _compute_view_frame,
    compute_cone_and_azimuth,
)
from .validation import _check_satellite, _check_satellites


def _get_visible_interval_series(
    point: Point | geo.Point,
    satellite: Satellite,
    min_elevation_angle: float,
    max_altitude: float,
    start: datetime,
    end: datetime,
) -> pd.Series:
    """
    Get the series of visible intervals based on altitude angle constraints.

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
        pandas.Series: Series of observation intervals.
    """
    # compute the maximum access time to filter bad data
    max_access_time = timedelta(
        seconds=compute_max_access_time(max_altitude, min_elevation_angle)
    )
    # find the set of observation events
    times, events = satellite.orbit.to_gp_orbit().get_observation_events(
        point, start, end, min_elevation_angle
    )

    # build the observation periods
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
        mid = start + (end - start) / 2
        longitude, latitude, elevation = _get_point_coordinates(point)
        topos = wgs84.latlon(latitude, longitude, elevation)
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track([mid])[0]
        elevation_angle = (
            (orbit_track - topos.at(timescale.from_datetime(mid))).altaz()[0].degrees
        )
        if elevation_angle > min_elevation_angle:  # type: ignore
            # continuously visible for the entire window
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(start.astimezone(tz=timezone.utc)),
                    right=pd.Timestamp(end.astimezone(tz=timezone.utc)),
                )
            ]
    elif np.all(events == 1):
        # if all events are type 1 (culminate), create a period from start to end
        obs_periods += [
            pd.Interval(
                left=pd.Timestamp(start.astimezone(tz=timezone.utc)),
                right=pd.Timestamp(end.astimezone(tz=timezone.utc)),
            )
        ]
    else:
        # otherwise, match rise/set events
        rises = times[events == 0]
        sets = times[events == 2]
        if (
            len(sets) > 0
            and (len(rises) == 0 or sets[0].utc_datetime() < rises[0].utc_datetime())
            and start < sets[0].utc_datetime()
        ):
            # if first event is a set, create a period from the start
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(start.astimezone(tz=timezone.utc)),
                    right=pd.Timestamp(sets[0].utc_datetime()),
                )
            ]
        # create an observation period to match with each rise event if
        # there is a following set event within twice the maximum access time
        obs_periods += [
            pd.Interval(
                left=pd.Timestamp(rise.utc_datetime()),
                right=pd.Timestamp(
                    sets[
                        np.logical_and(
                            rise.utc_datetime() < sets.utc_datetime(),
                            sets.utc_datetime()
                            < rise.utc_datetime() + 2 * max_access_time,
                        )
                    ][0].utc_datetime()
                ),
            )
            for rise in rises
            if np.any(
                np.logical_and(
                    rise.utc_datetime() < sets.utc_datetime(),
                    sets.utc_datetime() < rise.utc_datetime() + 2 * max_access_time,
                )
            )
        ]
        if (
            len(rises) > 0
            and (len(sets) == 0 or rises[-1].utc_datetime() > sets[-1].utc_datetime())
            and rises[-1].utc_datetime() < end
        ):
            # if last event is a rise, create a period to the end
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(rises[-1].utc_datetime()),
                    right=pd.Timestamp(end.astimezone(tz=timezone.utc)),
                )
            ]
    return pd.Series(obs_periods, dtype="interval")


def _get_empty_coverage_frame(omit_solar: bool) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for coverage analysis results.

    Args:
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "sat_alt": pd.Series(dtype="float"),
        "sat_az": pd.Series(dtype="float"),
    }
    if not omit_solar:
        columns = {
            **columns,
            "sat_sunlit": pd.Series(dtype="bool"),
            "solar_alt": pd.Series(dtype="float"),
            "solar_az": pd.Series(dtype="float"),
            "solar_time": pd.Series(dtype="float"),
        }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _build_observation_frame(
    observations: list[tuple[pd.Interval, pd.Timestamp, tuple[float, float, float]]],
    point_id: int,
    geometry: geo.base.BaseGeometry,
    satellite: Satellite,
    instrument: Instrument,
    omit_solar: bool,
) -> gpd.GeoDataFrame:
    """
    Builds the data frame of observations of a point or region, with the
    satellite (and solar) angles of each observed point at its epoch.

    Args:
        observations (list[tuple[pandas.Interval, pandas.Timestamp, tuple[float, float, float]]]):
                Each observation's period, epoch, and observed point's longitude
                (degrees), latitude (degrees), and elevation (meters).
        point_id (int): The identifier recorded with each observation.
        geometry (shapely.geometry.base.BaseGeometry): The geometry recorded with
                each observation.
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations.
    """
    if len(observations) == 0:
        return _get_empty_coverage_frame(omit_solar)
    gdf = gpd.GeoDataFrame(
        [
            {
                "point_id": point_id,
                "geometry": geometry,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "start": (
                    period.left
                    if not instrument.access_time_fixed
                    else epoch - instrument.min_access_time / 2
                ),
                "end": (
                    period.right
                    if not instrument.access_time_fixed
                    else epoch + instrument.min_access_time / 2
                ),
                "epoch": epoch,
            }
            for period, epoch, _ in observations
        ],
        crs="EPSG:4326",
    )
    # observed point of each observation
    longitude, latitude, elevation = np.array(
        [target for _, _, target in observations]
    ).T
    topos = wgs84.latlon(latitude, longitude, elevation)
    ts = timescale.from_datetimes(gdf.epoch)
    orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(gdf.epoch.tolist())
    # append satellite altitude/azimuth columns
    sat_altaz = (orbit_track - topos.at(ts)).altaz()
    gdf["sat_alt"] = sat_altaz[0].degrees  # type: ignore
    gdf["sat_az"] = sat_altaz[1].degrees  # type: ignore
    if not omit_solar:
        # append satellite sunlit column
        gdf["sat_sunlit"] = orbit_track.is_sunlit(de421)
        # append solar altitude/azimuth columns
        sun_altaz = (
            (de421["earth"] + topos).at(ts).observe(de421["sun"]).apparent().altaz()
        )
        gdf["solar_alt"] = sun_altaz[0].degrees
        gdf["solar_az"] = sun_altaz[1].degrees
        # append local solar time column
        gdf["solar_time"] = (de421["earth"] + topos).at(ts).observe(
            de421["sun"]
        ).apparent().hadec()[0].hours + 12
    return gdf


def _find_crossings(
    residual: Callable[[np.ndarray, np.ndarray], np.ndarray],
    lower: np.ndarray,
    upper: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Find a zero of a residual function within each of a set of intervals by
    the Illinois variant of the regula falsi method, vectorized across the
    intervals.

    Args:
        residual (Callable[[numpy.ndarray, numpy.ndarray], numpy.ndarray]):
                The residual function, evaluated at an array of times
                (seconds) within the intervals with the given indices.
        lower (numpy.ndarray): The lower ends of the intervals (seconds).
        upper (numpy.ndarray): The upper ends of the intervals (seconds).

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: The zero in each interval (or,
            if the residual does not change sign, the end with the smaller
            residual) and whether the interval brackets a zero.
    """
    lower, upper = np.array(lower, dtype=float), np.array(upper, dtype=float)
    index = np.arange(len(lower))
    f_lower, f_upper = residual(lower, index), residual(upper, index)
    # without a sign change, use the end with the smaller residual
    crossing = np.where(np.abs(f_lower) <= np.abs(f_upper), lower, upper)
    bracketed = np.sign(f_lower) * np.sign(f_upper) < 0
    # side of the bracket replaced in the previous iteration (-1: lower, 1: upper)
    side = np.zeros(len(lower))
    for _ in range(50):
        active = bracketed & (upper - lower > 1e-3)
        if not np.any(active):
            break
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(
                active, (lower * f_upper - upper * f_lower) / (f_upper - f_lower), lower
            )
        f_x = np.zeros(len(lower))
        f_x[active] = residual(x[active], index[active])
        replace_lower = active & (np.sign(f_x) == np.sign(f_lower))
        replace_upper = active & ~replace_lower
        # Illinois modification: halve the residual of an end retained twice
        f_upper = np.where(replace_lower & (side == -1), f_upper / 2, f_upper)
        f_lower = np.where(replace_upper & (side == 1), f_lower / 2, f_lower)
        lower = np.where(replace_lower, x, lower)
        f_lower = np.where(replace_lower, f_x, f_lower)
        upper = np.where(replace_upper, x, upper)
        f_upper = np.where(replace_upper, f_x, f_upper)
        side = np.where(replace_lower, -1, np.where(replace_upper, 1, side))
        crossing = np.where(active, x, crossing)
    return crossing, bracketed


def _refine_access_periods(
    residual: Callable[[Geocentric], np.ndarray],
    orbit: GeneralPerturbationsOrbit,
    periods: list[pd.Interval],
    max_step: timedelta | None = None,
) -> list[pd.Interval]:
    """
    Refine visible periods to the times when a residual function of the
    orbit track is not positive: for example, a target's angle from nadir
    less half an instrument's field of regard. The visible periods, from a
    conservative condition, bracket these times. Each period is sampled at
    21 times (or more, if needed to sample at least every `max_step`), and
    each change of sign of the residual between samples is refined, so that
    a period may be divided into several parts (for example, as a
    satellite passes over separate parts of a region). Periods in which the
    residual is positive at every sample are removed; period ends at which
    it is not positive (for example, at the ends of the analysis period) are
    kept.

    Args:
        residual (Callable[[skyfield.positionlib.Geocentric], numpy.ndarray]):
                The residual function of an orbit track (at one or more times).
        orbit (GeneralPerturbationsOrbit): The orbit.
        periods (list[pandas.Interval]): The visible periods.
        max_step (datetime.timedelta | None): The maximum time between samples.

    Returns:
        list[pandas.Interval]: The refined periods.
    """
    if len(periods) == 0:
        return periods
    reference = periods[0].left

    def evaluate(seconds: np.ndarray, _index: np.ndarray) -> np.ndarray:
        return np.reshape(
            residual(
                orbit.get_orbit_track(
                    [reference + pd.Timedelta(seconds=float(x)) for x in seconds]
                )
            ),
            -1,
        )

    lower = np.array([(period.left - reference).total_seconds() for period in periods])
    upper = np.array([(period.right - reference).total_seconds() for period in periods])
    counts = np.full(len(periods), 21)
    if max_step is not None:
        counts = np.maximum(
            counts, np.ceil((upper - lower) / max_step.total_seconds()).astype(int) + 1
        )
    samples = [
        np.linspace(lo, hi, count) for lo, hi, count in zip(lower, upper, counts)
    ]
    values = np.split(
        evaluate(np.concatenate(samples), np.array([])), np.cumsum(counts)[:-1]
    )
    # brackets of each change of sign between samples
    brackets = [
        (i, j)
        for i, value in enumerate(values)
        for j in np.flatnonzero((value[:-1] <= 0) != (value[1:] <= 0))
    ]
    crossings = {}
    if len(brackets) > 0:
        crossing, _ = _find_crossings(
            evaluate,
            np.array([samples[i][j] for i, j in brackets]),
            np.array([samples[i][j + 1] for i, j in brackets]),
        )
        crossings = dict(zip(brackets, crossing))
    refined = []
    for i, (sample, value) in enumerate(zip(samples, values)):
        inside = value <= 0
        left = sample[0] if inside[0] else None
        for j in range(len(sample) - 1):
            if inside[j] == inside[j + 1]:
                continue
            if inside[j + 1]:
                left = crossings[(i, j)]
            else:
                refined.append((left, crossings[(i, j)]))
                left = None
        if left is not None:
            refined.append((left, sample[-1]))
    return [
        pd.Interval(
            left=reference + pd.Timedelta(seconds=float(left)),
            right=reference + pd.Timedelta(seconds=float(right)),
        )
        for left, right in refined
    ]


def _get_view_crossing_times(
    target: GeographicPosition,
    satellite: Satellite,
    instrument: PointedInstrument,
    periods: list[pd.Interval],
) -> list[pd.Timestamp]:
    """
    Get the time in each visible period when a pointed instrument's view
    sweeps over a target: when the target crosses the plane of the view's
    boresight and cross-track axis (where its along-track view angle is
    zero). The along-track component of the unit line of sight to the
    target, which is defined even when the target is behind the view,
    decreases through zero as the satellite passes the target. If a period
    does not contain a crossing (for example, a period truncated by the
    analysis window), uses the period's end closest to the crossing.

    Args:
        target (skyfield.toposlib.GeographicPosition): The target position.
        satellite (Satellite): The observing satellite.
        instrument (PointedInstrument): The observing instrument.
        periods (list[pandas.Interval]): The visible periods.

    Returns:
        list[pandas.Timestamp]: The view crossing time in each period.
    """
    if len(periods) == 0:
        return []
    orbit = satellite.orbit.to_gp_orbit()
    reference = periods[0].left

    def residual(seconds: np.ndarray, _index: np.ndarray) -> np.ndarray:
        # along-track component of the unit line of sight in the view frame
        orbit_track = orbit.get_orbit_track(
            [reference + pd.Timedelta(seconds=float(x)) for x in seconds]
        )
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
            los = np.reshape(np.array(target.itrs_xyz.m), (3, 1)) - position
            return np.sum(los * along, axis=0) / np.linalg.norm(los, axis=0)
        position, _, along, _ = _compute_view_frame(
            orbit_track,
            instrument.get_roll_angle(orbit_track),
            instrument.get_pitch_angle(orbit_track),
            instrument.velocity_frame,
            instrument.nadir_reference,
        )
        los = np.reshape(np.array(target.itrs_xyz.m), (3, 1)) - position
        return np.sum(los * along, axis=0) / np.linalg.norm(los, axis=0)

    crossing, _ = _find_crossings(
        residual,
        np.array([(period.left - reference).total_seconds() for period in periods]),
        np.array([(period.right - reference).total_seconds() for period in periods]),
    )
    return [reference + pd.Timedelta(seconds=float(x)) for x in crossing]


def _get_cone_crossing_times(
    target: GeographicPosition,
    satellite: Satellite,
    instrument: ConicalInstrument,
    periods: list[pd.Interval],
) -> list[list[pd.Timestamp]]:
    """
    Get the times in each visible period when a target crosses a conical
    instrument's cone: when the target's angle from nadir equals the cone
    angle. During a pass, this angle decreases to a minimum near the closest
    approach and then increases, so there are up to two crossings: entering
    the cone (ahead of the satellite) and leaving it (behind).

    Args:
        target (skyfield.toposlib.GeographicPosition): The target position.
        satellite (Satellite): The observing satellite.
        instrument (ConicalInstrument): The observing instrument.
        periods (list[pandas.Interval]): The visible periods.

    Returns:
        list[list[pandas.Timestamp]]: The cone crossing times in each period.
    """
    if len(periods) == 0:
        return []
    orbit = satellite.orbit.to_gp_orbit()
    reference = periods[0].left
    n = len(periods)

    def cone_angle(seconds: np.ndarray, _index: np.ndarray) -> np.ndarray:
        cone, _ = compute_cone_and_azimuth(
            orbit.get_orbit_track(
                [reference + pd.Timedelta(seconds=float(x)) for x in seconds]
            ),
            target,
            instrument.velocity_frame,
            instrument.nadir_reference,
        )
        return np.reshape(cone, -1)

    lower = np.array([(period.left - reference).total_seconds() for period in periods])
    upper = np.array([(period.right - reference).total_seconds() for period in periods])
    # time of the minimum angle from nadir in each period, from sampled times
    samples = lower[:, None] + (upper - lower)[:, None] * np.linspace(0, 1, 21)
    closest = samples[
        np.arange(n),
        np.argmin(
            cone_angle(samples.ravel(), np.repeat(np.arange(n), 21)).reshape(
                samples.shape
            ),
            axis=1,
        ),
    ]
    crossing, bracketed = _find_crossings(
        lambda seconds, index: cone_angle(seconds, index) - instrument.cone_angle,
        np.concatenate([lower, closest]),
        np.concatenate([closest, upper]),
    )
    return [
        [
            reference + pd.Timedelta(seconds=float(crossing[i + k * n]))
            for k in range(2)
            if bracketed[i + k * n]
        ]
        for i in range(n)
    ]


def collect_observations(
    point: Point | geo.Point,
    satellite: Satellite,
    start: datetime,
    end: datetime,
    instrument_index: int = 0,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect single satellite observations of a geodetic point of interest:
    a TAT-C `Point` or a shapely `Point` (whose x, y, and optional z
    coordinates are its longitude, latitude, and elevation in meters);
    TAT-C points are expected to be replaced by shapely points in the
    future. Observations record a TAT-C point's `id` as their `point_id`, or
    0 for a shapely point. For a region, see
    `tatc.analysis.region_coverage.collect_region_observations`.

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

    Args:
        point (Point | shapely.geometry.Point): The ground point of interest.
        satellite (Satellite): The observing satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        instrument_index (int): The index of the observing instrument in satellite.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations.
    """
    _check_satellite(satellite)
    if not isinstance(point, (Point, geo.Point)):
        raise TypeError(
            "point must be a Point or shapely Point, not a "
            f"{type(point).__name__} (see collect_region_observations for a region)"
        )
    point_id = point.id if isinstance(point, Point) else 0
    longitude, latitude, elevation = _get_point_coordinates(point)
    instrument = satellite.instruments[instrument_index]
    orbit = satellite.orbit.to_gp_orbit()
    # use the apogee altitude above the polar radius (and above the point,
    # if below the ellipsoid) as a conservative upper bound for computing
    # access times, which are then refined to the field of regard
    max_altitude = (
        compute_apoapsis_radius(
            satellite.orbit.get_semimajor_axis(), satellite.orbit.get_eccentricity()
        )
        - EARTH_POLAR_RADIUS
        - min(elevation, 0)
    )
    # compute the minimum altitude angle required for observation, less a
    # margin for the spherical approximation (but not below the horizon)
    min_elevation_angle = max(
        0.0,
        compute_min_elevation_angle(max_altitude, instrument.field_of_regard) - 1.0,
    )
    target = wgs84.latlon(latitude, longitude, elevation)
    periods = list(
        _get_visible_interval_series(
            point, satellite, min_elevation_angle, max_altitude, start, end
        )
    )
    # refine the periods to the field of regard: when the point's angle
    # from nadir is at most half the field of regard
    half_angle = instrument.field_of_regard / 2
    if half_angle < 90:
        periods = _refine_access_periods(
            lambda orbit_track: compute_cone_and_azimuth(
                orbit_track, target, nadir_reference=instrument.nadir_reference
            )[0]
            - half_angle,
            orbit,
            periods,
        )
    # observation epochs: the time a pointed instrument's view sweeps over
    # the point, the times the point crosses a conical instrument's cone, or
    # otherwise the midpoint of each visible period
    if isinstance(instrument, PointedInstrument):
        epochs = [
            [epoch]
            for epoch in _get_view_crossing_times(
                target, satellite, instrument, periods
            )
        ]
    elif isinstance(instrument, ConicalInstrument):
        epochs = _get_cone_crossing_times(target, satellite, instrument, periods)
    else:
        epochs = [[period.mid] for period in periods]
    observations = []
    for period, period_epochs in zip(periods, epochs):
        for epoch in period_epochs:
            # instrument validity (illumination, field of view) is only
            # checked at each period's epoch, as an approximation of the
            # whole interval; a more general approach would refine the exact
            # observation period boundaries with Skyfield's find_discrete
            # using the instrument's own validity condition, but that is out
            # of scope for now
            orbit_track = orbit.get_orbit_track([epoch])
            if (
                instrument.min_access_time <= period.right - period.left
                and instrument.is_valid_observation(orbit_track, target).all()
                and (
                    not isinstance(instrument, (PointedInstrument, ConicalInstrument))
                    or instrument.is_in_field_of_view(orbit_track, target).all()
                )
            ):
                observations.append((period, epoch, (longitude, latitude, elevation)))
    return _build_observation_frame(
        observations,
        point_id,
        geo.Point(longitude, latitude, elevation),
        satellite,
        instrument,
        omit_solar,
    )


def collect_multi_observations(
    point: Point | geo.Point,
    satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect multiple satellite observations of a geodetic point of interest:
    calls `collect_observations` for every instrument on every satellite in
    `satellites`, and concatenates the results into one data frame.

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
    gdfs = [
        collect_observations(point, satellite, start, end, instrument_index, omit_solar)
        for satellite in _check_satellites(satellites)
        for instrument_index in range(len(satellite.instruments))
    ]
    if len(gdfs) == 0:
        # an empty `satellites` list leaves nothing to concatenate
        return _get_empty_coverage_frame(omit_solar)
    # concatenate into one data frame, sort by start time, and re-index
    return pd.concat(gdfs).sort_values("start").reset_index(drop=True)


def _get_target_keys(gdf: gpd.GeoDataFrame) -> list[pd.Series]:
    """
    Gets the keys that identify the target (point or region) of each
    observation: its `point_id` and its geometry (as well-known binary), so
    that targets with distinct geometries are kept apart even if they share
    a `point_id` (for example, shapely points with the default identifier
    of 0).

    Args:
        gdf (geopandas.GeoDataFrame): The observations.

    Returns:
        list[pandas.Series]: the keys
    """
    return [gdf["point_id"], gdf.geometry.to_wkb().rename("geometry_key")]


def _get_empty_aggregate_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for aggregated coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def aggregate_observations(observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Aggregate constellation observations. Interleaves observations by multiple
    satellites to compute aggregate performance metrics including access
    (observation duration) and revisit (duration between observations).
    Overlapping (including fully nested) observations of the same target
    (the same `point_id` and geometry, see `_get_target_keys`),
    possibly from different satellites/instruments, are merged into a single
    continuous coverage period; `satellite`/`instrument` record every
    contributor to that period, comma-separated. `epoch` is reassigned to
    the midpoint of the merged period's `start`/`end` (a representative
    instant), not the mean of the constituent observations' own epochs.
    Per-observation columns that lose their meaning once merged across
    satellites and over a potentially much longer period -- e.g. `sat_alt`,
    `sat_az`, `sat_sunlit`, `solar_alt`, `solar_az`, `solar_time` -- are
    intentionally dropped, even if present on `observations`.

    Args:
        observations (geopandas.GeoDataFrame): The collected observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with aggregated observations.
    """
    if observations.empty:
        return _get_empty_aggregate_frame()
    gdfs = []
    # split into constituent data frames for each target
    for _, gdf in observations.groupby(_get_target_keys(observations)):
        # sort the values by start datetime
        gdf = gdf.sort_values("start")
        # assign the observation group number based on overlapping start/end times
        gdf["obs"] = (gdf["start"] > gdf["end"].shift().cummax()).cumsum()
        # perform the aggregation to group overlapping observations
        gdf = gdf.dissolve(
            "obs",
            aggfunc={
                "point_id": "first",
                "satellite": ", ".join,  # type: ignore
                "instrument": ", ".join,  # type: ignore
                "start": "min",
                "end": "max",
            },
        )
        # reassign epoch to the midpoint of the merged period, as a single
        # representative instant, rather than the mean of the constituent
        # observations' own (pre-merge) epochs
        gdf["epoch"] = gdf["start"] + (gdf["end"] - gdf["start"]) / 2
        # compute access and revisit metrics
        gdf["access"] = gdf["end"] - gdf["start"]
        gdf["revisit"] = gdf["start"] - gdf["end"].shift()
        # append to the list of data frames
        gdfs.append(gdf)
    # return a concatenated data frame and re-index
    return pd.concat(gdfs).reset_index(drop=True)


def _get_empty_reduce_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for reduced coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "access": pd.Series([], dtype="timedelta64[ns]"),
        "revisit": pd.Series([], dtype="timedelta64[ns]"),
        "samples": pd.Series([], dtype="int"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def reduce_observations(aggregated_observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reduce constellation observations: for each unique target (`point_id`
    and geometry, see `_get_target_keys`) in `aggregated_observations`, computes the mean access period, the mean
    revisit period, and the total number of samples (aggregated periods)
    over the analysis period. The first sample's revisit is undefined (no
    prior observation to measure a gap from) and is excluded from the mean
    rather than counted as zero, which would otherwise bias the mean
    downward; a point with only one sample accordingly has an undefined
    (NaT) mean revisit.

    Args:
        aggregated_observations (geopandas.GeoDataFrame): The aggregated observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced observations.
    """
    if aggregated_observations.empty:
        return _get_empty_reduce_frame()
    # operate on a copy of the data frame
    gdf = aggregated_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"].dt.total_seconds()
    gdf["revisit"] = gdf["revisit"].dt.total_seconds()
    # assign each record to one observation
    gdf["samples"] = 1
    # perform the aggregation operation for each target
    gdf = (
        gdf.dissolve(
            _get_target_keys(gdf),
            aggfunc={
                "access": "mean",
                "revisit": "mean",
                "samples": "sum",
            },
        )
        .reset_index()
        .drop(columns="geometry_key")
    )
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = pd.to_timedelta(gdf["access"], unit="s")
    gdf["revisit"] = pd.to_timedelta(gdf["revisit"], unit="s")
    return gdf


def grid_observations(
    reduced_observations: gpd.GeoDataFrame, cells: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Grid reduced observations to cells: for every cell, sums the number of
    samples across every point it contains, and combines those points'
    access/revisit statistics into a single representative value per cell.
    Both access (a per-event duration) and revisit (a time-between-events
    duration, i.e. the reciprocal of a sampling rate) use a sample-weighted
    mean -- arithmetic for access, harmonic for revisit, since revisit
    needs to be averaged as a rate to stay a representative statistic.

    Args:
        reduced_observations (geopandas.GeoDataFrame): The reduced observations.
        cells (geopandas.GeoDataFrame): The cell specification.

    Returns:
        geopandas.GeoDataFrame: The data frame with gridded observations.
    """
    if reduced_observations.empty:
        gdf = cells.copy()
        gdf["samples"] = 0
        gdf["access"] = None
        gdf["revisit"] = None
        return gdf
    # operate on a copy of the data frame
    gdf = reduced_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"].dt.total_seconds()
    gdf["revisit"] = gdf["revisit"].dt.total_seconds()
    # pre-transform so the means below reduce to plain sums: groupby().agg()
    # with a dict of {column: function} only ever hands a custom callable
    # its own column's Series, never a sibling column like "samples" needed
    # to compute a weighted statistic within the callable. The weighted
    # harmonic mean of revisit is sum(samples) / sum(samples/revisit).
    gdf["access_x_samples"] = gdf["access"] * gdf["samples"]
    gdf["samples_over_revisit"] = gdf["samples"] / gdf["revisit"]
    gdf = (
        cells.sjoin(gdf, how="inner", predicate="contains")
        .dissolve(
            by="cell_id",
            aggfunc={
                "samples": "sum",
                "access_x_samples": "sum",
                "samples_over_revisit": "sum",
            },
        )
        .reset_index()
    )
    # finish the aggregation: sample-weighted arithmetic mean for access,
    # sample-weighted harmonic mean for revisit
    gdf["access"] = gdf["access_x_samples"] / gdf["samples"]
    gdf["revisit"] = gdf["samples"] / gdf["samples_over_revisit"]
    gdf = gdf.drop(columns=["access_x_samples", "samples_over_revisit"])
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = pd.to_timedelta(gdf["access"], unit="s")
    gdf["revisit"] = pd.to_timedelta(gdf["revisit"], unit="s")
    return gdf
