"""
Methods to perform coverage analysis of points (see `observations` for the
aggregation, reduction, and gridding of observations).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import geometry as geo
from skyfield.api import wgs84
from skyfield.toposlib import GeographicPosition

from ..constants import EARTH_POLAR_RADIUS, timescale
from ..schemas import ConicalInstrument, Point, PointedInstrument, Satellite
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
from .observations import (
    _build_observation_frame,
    _find_crossings,
    _get_empty_coverage_frame,
    _refine_access_periods,
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
    _check_satellite(satellite)
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
    return _get_visible_interval_series(
        point, satellite, min_elevation_angle, max_altitude, start, end
    )


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
        compute_access_periods(point, satellite, start, end, min_elevation_angle)
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
