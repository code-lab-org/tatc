"""
Methods to perform coverage analysis.

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
from ..utils.observation import (
    compute_max_access_time,
    compute_min_elevation_angle,
)
from ..utils.orbital import compute_apoapsis_radius
from ..utils.projection import _compute_view_frame, compute_cone_and_azimuth


def _get_visible_interval_series(
    point: Point,
    satellite: Satellite,
    min_elevation_angle: float,
    max_altitude: float,
    start: datetime,
    end: datetime,
) -> pd.Series:
    """
    Get the series of visible intervals based on altitude angle constraints.

    Args:
        point (Point): Point to observe.
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
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(mid)
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


def _get_orbit_track(
    orbit: GeneralPerturbationsOrbit,
    times: list[datetime],
    shifts: list[timedelta] | None = None,
) -> Geocentric:
    """
    Get the orbit track at a list of times. With nonzero shifts, the
    satellite's Earth-fixed position and velocity at each time are those
    propagated to the time minus the shift, a whole number of repeat cycles,
    as for an orbit maintained on its repeat ground track; they are
    expressed in the inertial (GCRS) frame at the time itself, so that
    quantities that depend on the time (such as the Sun's position) are
    unaffected.

    Args:
        orbit (GeneralPerturbationsOrbit): The orbit.
        times (list[datetime.datetime]): The times.
        shifts (list[datetime.timedelta] | None): The shift at each time.

    Returns:
        skyfield.positionlib.Geocentric: The orbit track.
    """
    if shifts is None or not any(shifts):
        return orbit.get_orbit_track(times)
    shifted = orbit.get_orbit_track([t - shift for t, shift in zip(times, shifts)])
    position, velocity = shifted.frame_xyz_and_velocity(itrs)
    t = timescale.from_datetimes(times)
    # invert the Earth-fixed transformation (r_itrs = R r, v_itrs = R v + V r_itrs)
    rotation = itrs.rotation_at(t)
    rate = itrs._dRdt_times_RT_at(t)  # pylint: disable=protected-access
    v_rotating = velocity.au_per_d - np.einsum("ij...,j...->i...", rate, position.au)
    return Geocentric(
        np.einsum("ji...,j...->i...", rotation, position.au),
        np.einsum("ji...,j...->i...", rotation, v_rotating),
        t,
        center=399,
    )


def _get_repeat_shifts(
    orbit: GeneralPerturbationsOrbit,
    start: datetime,
    end: datetime,
    periods: list[pd.Interval],
) -> list[timedelta]:
    """
    Get the shift of each visible period, a whole number of repeat cycles,
    if `get_observation_events` repeats the events of the first cycle to
    cover the analysis period (otherwise zero).

    Args:
        orbit (GeneralPerturbationsOrbit): The orbit.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        periods (list[pandas.Interval]): The visible periods.

    Returns:
        list[datetime.timedelta]: The shift of each period.
    """
    repeat_cycle = orbit.get_observation_repeat_cycle(start, end)
    if repeat_cycle is None:
        return [timedelta(0) for _ in periods]
    reference = pd.Timestamp(start.astimezone(tz=timezone.utc))
    return [
        int(np.floor((period.mid - reference) / repeat_cycle)) * repeat_cycle
        for period in periods
    ]


def _refine_access_periods(
    target: GeographicPosition,
    orbit: GeneralPerturbationsOrbit,
    instrument: Instrument,
    periods: list[pd.Interval],
    shifts: list[timedelta],
) -> tuple[list[pd.Interval], list[timedelta]]:
    """
    Refine visible periods to the times when a target lies within an
    instrument's field of regard: when its angle from nadir (using the
    instrument's nadir reference) is at most half the field of regard. The
    visible periods, from a conservative minimum elevation angle, bracket
    these times. Periods during which the target does not enter the field
    of regard are removed; period ends at which it is already inside (for
    example, at the ends of the analysis period) are kept.

    Args:
        target (skyfield.toposlib.GeographicPosition): The target position.
        orbit (GeneralPerturbationsOrbit): The orbit.
        instrument (Instrument): The observing instrument.
        periods (list[pandas.Interval]): The visible periods.
        shifts (list[datetime.timedelta]): The repeat-cycle shift of each period.

    Returns:
        tuple[list[pandas.Interval], list[datetime.timedelta]]: The refined
            periods and their shifts.
    """
    half_angle = instrument.field_of_regard / 2
    if len(periods) == 0 or half_angle >= 90:
        return periods, shifts
    reference = periods[0].left
    n = len(periods)

    def angle_from_nadir(seconds: np.ndarray, index: np.ndarray) -> np.ndarray:
        angle, _ = compute_cone_and_azimuth(
            _get_orbit_track(
                orbit,
                [reference + pd.Timedelta(seconds=float(x)) for x in seconds],
                [shifts[i % n] for i in index],
            ),
            target,
            nadir_reference=instrument.nadir_reference,
        )
        return np.reshape(angle, -1)

    lower = np.array([(period.left - reference).total_seconds() for period in periods])
    upper = np.array([(period.right - reference).total_seconds() for period in periods])
    # time of the minimum angle from nadir in each period, from sampled times
    samples = lower[:, None] + (upper - lower)[:, None] * np.linspace(0, 1, 21)
    angles = angle_from_nadir(samples.ravel(), np.repeat(np.arange(n), 21)).reshape(
        samples.shape
    )
    closest = samples[np.arange(n), np.argmin(angles, axis=1)]
    crossing, bracketed = _find_crossings(
        lambda seconds, index: angle_from_nadir(seconds, index) - half_angle,
        np.concatenate([lower, closest]),
        np.concatenate([closest, upper]),
    )
    refined, refined_shifts = [], []
    for i in range(n):
        if np.min(angles[i]) > half_angle:
            continue
        left = crossing[i] if bracketed[i] else lower[i]
        right = crossing[i + n] if bracketed[i + n] else upper[i]
        refined.append(
            pd.Interval(
                left=reference + pd.Timedelta(seconds=float(left)),
                right=reference + pd.Timedelta(seconds=float(right)),
            )
        )
        refined_shifts.append(shifts[i])
    return refined, refined_shifts


def _get_view_crossing_times(
    target: GeographicPosition,
    satellite: Satellite,
    instrument: PointedInstrument,
    periods: list[pd.Interval],
    shifts: list[timedelta] | None = None,
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
        shifts (list[datetime.timedelta] | None): The repeat-cycle shift of
                each period (see `_get_orbit_track`).

    Returns:
        list[pandas.Timestamp]: The view crossing time in each period.
    """
    if len(periods) == 0:
        return []
    orbit = satellite.orbit.to_gp_orbit()
    reference = periods[0].left
    shifts = shifts or [timedelta(0) for _ in periods]

    def residual(seconds: np.ndarray, index: np.ndarray) -> np.ndarray:
        # along-track component of the unit line of sight in the view frame
        orbit_track = _get_orbit_track(
            orbit,
            [reference + pd.Timedelta(seconds=float(x)) for x in seconds],
            [shifts[i] for i in index],
        )
        position, _, along, _ = _compute_view_frame(
            orbit_track,
            instrument.get_roll_angle(orbit_track),
            instrument.pitch_angle,
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
    shifts: list[timedelta] | None = None,
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
        shifts (list[datetime.timedelta] | None): The repeat-cycle shift of
                each period (see `_get_orbit_track`).

    Returns:
        list[list[pandas.Timestamp]]: The cone crossing times in each period.
    """
    if len(periods) == 0:
        return []
    orbit = satellite.orbit.to_gp_orbit()
    reference = periods[0].left
    shifts = shifts or [timedelta(0) for _ in periods]
    n = len(periods)

    def cone_angle(seconds: np.ndarray, index: np.ndarray) -> np.ndarray:
        cone, _ = compute_cone_and_azimuth(
            _get_orbit_track(
                orbit,
                [reference + pd.Timedelta(seconds=float(x)) for x in seconds],
                [shifts[i % n] for i in index],
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
    point: Point,
    satellite: Satellite,
    start: datetime,
    end: datetime,
    instrument_index: int = 0,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect single satellite observations of a geodetic point of interest.
    Each observation spans a period when the point lies within the
    instrument's field of regard (when its angle from nadir is at most half
    the field of regard). Its epoch is the period's midpoint or, for
    a `PointedInstrument`, the time when the instrument's view sweeps over
    the point (when the point's along-track view angle relative to the view
    center is zero), or, for a `ConicalInstrument`, a time when the point
    crosses the scanned cone (up to two per period, entering and leaving the
    cone); at that time, the point must lie within the field of view.

    If the orbit's observation events are repeated with its repeat cycle
    (see `GeneralPerturbationsOrbit.get_observation_repeat_cycle`), the
    orbit is modeled as maintained on its repeat ground track: within each
    repeated cycle, epochs and field of view are evaluated with the
    satellite's Earth-fixed position and velocity propagated to the
    corresponding time in the first cycle.

    Args:
        point (Point): The ground point of interest.
        satellite (Satellite): The observing satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        instrument_index (int): The index of the observing instrument in satellite.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations.
    """
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
        - min(point.elevation, 0)
    )
    # compute the minimum altitude angle required for observation, less a
    # margin for the spherical approximation (but not below the horizon)
    min_elevation_angle = max(
        0.0,
        compute_min_elevation_angle(max_altitude, instrument.field_of_regard) - 1.0,
    )
    target = wgs84.latlon(point.latitude, point.longitude, point.elevation)
    periods = list(
        _get_visible_interval_series(
            point, satellite, min_elevation_angle, max_altitude, start, end
        )
    )
    # whole repeat cycles by which repeated periods are shifted
    shifts = _get_repeat_shifts(orbit, start, end, periods)
    # refine the periods to the field of regard
    periods, shifts = _refine_access_periods(target, orbit, instrument, periods, shifts)
    # observation epochs: the time a pointed instrument's view sweeps over
    # the point, the times the point crosses a conical instrument's cone, or
    # otherwise the midpoint of each visible period
    if isinstance(instrument, PointedInstrument):
        epochs = [
            [epoch]
            for epoch in _get_view_crossing_times(
                target, satellite, instrument, periods, shifts
            )
        ]
    elif isinstance(instrument, ConicalInstrument):
        epochs = _get_cone_crossing_times(
            target, satellite, instrument, periods, shifts
        )
    else:
        epochs = [[period.mid] for period in periods]
    records, record_shifts = [], []
    for period, period_epochs, shift in zip(periods, epochs, shifts):
        for epoch in period_epochs:
            # instrument validity (illumination, field of view) is only
            # checked at each period's epoch, as an approximation of the
            # whole interval; a more general approach would refine the exact
            # observation period boundaries with Skyfield's find_discrete
            # using the instrument's own validity condition, but that is out
            # of scope for now
            orbit_track = _get_orbit_track(orbit, [epoch], [shift])
            if not (
                instrument.min_access_time <= period.right - period.left
                and instrument.is_valid_observation(orbit_track, target).all()
                and (
                    not isinstance(instrument, (PointedInstrument, ConicalInstrument))
                    or instrument.is_in_field_of_view(orbit_track, target).all()
                )
            ):
                continue
            records.append(
                {
                    "point_id": point.id,
                    "geometry": geo.Point(
                        point.longitude, point.latitude, point.elevation
                    ),
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
            )
            record_shifts.append(shift)

    # build the dataframe
    if len(records) > 0:
        gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        ts = timescale.from_datetimes(gdf.epoch)
        orbit_track = _get_orbit_track(orbit, gdf.epoch.tolist(), record_shifts)
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
    else:
        gdf = _get_empty_coverage_frame(omit_solar)
    return gdf


def collect_multi_observations(
    point: Point,
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
        point (Point): The ground point of interest.
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
        for satellite in (satellites if isinstance(satellites, list) else [satellites])
        for instrument_index in range(len(satellite.instruments))
    ]
    if len(gdfs) == 0:
        # an empty `satellites` list leaves nothing to concatenate
        return _get_empty_coverage_frame(omit_solar)
    # concatenate into one data frame, sort by start time, and re-index
    return pd.concat(gdfs).sort_values("start").reset_index(drop=True)


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
    Overlapping (including fully nested) observations for the same point,
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
    # split into constituent data frames based on point_id
    for _, gdf in observations.groupby("point_id"):
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
    Reduce constellation observations: for each unique point_id in
    `aggregated_observations`, computes the mean access period, the mean
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
    # perform the aggregation operation
    gdf = gdf.dissolve(
        "point_id",
        aggfunc={
            "access": "mean",
            "revisit": "mean",
            "samples": "sum",
        },
    ).reset_index()
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
