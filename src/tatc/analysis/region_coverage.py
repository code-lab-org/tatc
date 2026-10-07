"""
Methods to perform coverage analysis of regions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely import geometry as geo
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric

from ..constants import (
    EARTH_ECCENTRICITY,
    EARTH_EQUATORIAL_RADIUS,
    EARTH_FLATTENING,
    EARTH_MU,
    EARTH_POLAR_RADIUS,
    EARTH_ROTATION_RATE,
    timescale,
)
from ..schemas import Satellite
from ..utils.geometry import (
    _get_angular_distance_to_arcs,
    _get_boundary_arcs,
    _get_geodetic_coordinates,
    _get_nearest_arc_points,
    _get_surface_directions,
    _get_surface_positions,
    project_polygon_to_elevation,
    split_polygon,
)
from ..utils.projection import NadirReference, VelocityFrame, _compute_view_frame
from .point_coverage import (
    _build_observation_frame,
    _get_empty_coverage_frame,
    _refine_access_periods,
)
from .validation import _check_satellite, _check_satellites


def _get_visible_polygon_interval_series(
    geometry: geo.Polygon | geo.MultiPolygon,
    satellite: Satellite,
    field_of_regard: float,
    start: datetime,
    end: datetime,
    elevation: float = 0,
    margin: float = 0.1,
    coarse_step: float = 10,
) -> pd.Series:
    """
    Get the series of periods when an instrument's field of regard (a cone
    about nadir) may observe any part of a region: a conservative superset
    of the observation periods, for culling, from which no observation is
    missed however brief.

    The region is observable when the angular distance `d(t)` from the
    satellite's geocentric direction to the region (zero within it) is at
    most the Earth central angle `lambda` from nadir to the edge of the field
    of regard (or the horizon), so the periods are those when

        g(t) = d(t) - lambda - margin

    is not positive. For a single point, this is equivalent to a minimum
    elevation angle condition (see
    `tatc.analysis.point_coverage._get_visible_interval_series`). The
    central angle `lambda` is evaluated at the orbit's apoapsis (its largest
    over the orbit) on a sphere of the region's smallest geocentric radius
    (larger than on the ellipsoid), and widened by the largest separation of
    the geodetic and geocentric nadir directions as seen on the ground, so
    that it is conservative over the whole orbit.

    `d(t)` changes no faster than the angular rate `omega` of the satellite's
    direction in the Earth-fixed frame (bounded by the orbital angular rate
    at periapsis plus the Earth's rotation), so no period lies within an
    interval `[t_a, t_b]` in which `g(t_a) + g(t_b) > omega * (t_b - t_a)`, and
    the whole interval lies within a period if `g(t_a) + g(t_b) < -omega *
    (t_b - t_a)`. Starting from samples spaced by `coarse_step` degrees of the
    satellite's motion, intervals that are neither are halved until they
    are, or until they are shorter than a millisecond (kept as part of a
    period).

    The region's boundary is treated as great circle arcs between points at
    most a degree apart along its edges (see
    `tatc.utils.geometry._get_boundary_arcs`). The region is split along the
    anti-meridian and poles (see `tatc.utils.geometry.split_polygon`), so
    that its parts can be tested in longitude and latitude.

    Args:
        geometry (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region (with EPSG:4326 coordinates) to observe.
        satellite (Satellite): Satellite doing the observation.
        field_of_regard (float): The instrument's field of regard (degrees).
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        elevation (float): Elevation (meters) of the region above the WGS 84 ellipsoid.
        margin (float): Additional central angle (degrees) to widen the field of regard.
        coarse_step (float): Angle (degrees) of the satellite's motion between initial samples.

    Returns:
        pandas.Series: Series of observation intervals.
    """
    geometry = split_polygon(geometry)
    shapely.prepare(geometry)
    arcs = _get_boundary_arcs(geometry, elevation)
    orbit = satellite.orbit.to_gp_orbit()
    shapes = [
        (element.get_semimajor_axis(), element.eccentricity)
        for element in orbit.elements
    ]
    # bounds on the orbit radius, with a margin for the short-period
    # perturbations of the osculating radius about the mean elements, and on
    # the angular rate of the satellite's direction in the Earth-fixed frame
    # (radians/day), with a 2 percent margin, over all elements
    apoapsis = max(a * (1 + e) for a, e in shapes) + 25e3
    rate = 86400 * (
        1.02
        * max(
            np.sqrt(EARTH_MU * (1 + e) / (a * (1 - e))) / (a * (1 - e))
            for a, e in shapes
        )
        + EARTH_ROTATION_RATE
    )
    # earth central angle from nadir to the edge of the field of regard, or
    # to the horizon, at apoapsis, on a sphere of the smallest (geocentric)
    # radius of the region
    _, min_latitude, _, max_latitude = geometry.bounds
    latitude = np.radians(max(abs(min_latitude), abs(max_latitude)))
    radius = elevation + np.sqrt(
        (
            (EARTH_EQUATORIAL_RADIUS**2 * np.cos(latitude)) ** 2
            + (EARTH_POLAR_RADIUS**2 * np.sin(latitude)) ** 2
        )
        / (
            (EARTH_EQUATORIAL_RADIUS * np.cos(latitude)) ** 2
            + (EARTH_POLAR_RADIUS * np.sin(latitude)) ** 2
        )
    )
    half_angle = np.radians(min(field_of_regard / 2, 90))
    ratio = apoapsis / radius * np.sin(half_angle)
    central_angle = (
        np.arcsin(ratio) - half_angle if ratio < 1 else np.arccos(radius / apoapsis)
    )
    # separation on the ground of the geodetic and geocentric nadir
    # directions (which differ by up to about the flattening, in radians)
    nadir_offset = (apoapsis - radius) / radius * EARTH_FLATTENING
    threshold = central_angle + nadir_offset + np.radians(margin)
    e2 = EARTH_ECCENTRICITY**2

    def excess(x: np.ndarray) -> np.ndarray:
        position = np.array(
            orbit.get_orbit_track_at_time(timescale.tt_jd(x)).frame_xyz(itrs).m
        ).reshape(3, -1)
        u = position / np.linalg.norm(position, axis=0)
        distance = _get_angular_distance_to_arcs(u, arcs)
        # geodetic latitude of the surface point in the satellite's direction
        longitude = np.degrees(np.arctan2(u[1], u[0]))
        latitude = np.degrees(np.arctan2(u[2], (1 - e2) * np.hypot(u[0], u[1])))
        distance[shapely.contains_xy(geometry, longitude, latitude)] = 0
        return distance - threshold

    t_start = timescale.from_datetime(start).tt
    t_end = max(timescale.from_datetime(end).tt, t_start)
    count = int(np.ceil((t_end - t_start) * rate / np.radians(coarse_step))) + 1
    x = np.linspace(t_start, t_end, max(count, 2))
    f = excess(x)
    lower, upper, f_lower, f_upper = x[:-1], x[1:], f[:-1], f[1:]
    tolerance = 1e-3 / 86400
    windows = []
    while len(lower) > 0:
        span = rate * (upper - lower)
        inside = f_lower + f_upper < -span
        unresolved = ~inside & (f_lower + f_upper <= span)
        # keep intervals within a period, and those too short to resolve
        small = unresolved & (upper - lower <= tolerance)
        windows.extend(zip(lower[inside | small], upper[inside | small]))
        halve = unresolved & ~small
        lower, upper = lower[halve], upper[halve]
        f_lower, f_upper = f_lower[halve], f_upper[halve]
        if len(lower) > 0:
            middle = (lower + upper) / 2
            f_middle = excess(middle)
            lower, upper = np.concatenate((lower, middle)), np.concatenate(
                (middle, upper)
            )
            f_lower = np.concatenate((f_lower, f_middle))
            f_upper = np.concatenate((f_middle, f_upper))
    # merge adjoining intervals into periods
    periods = []
    for left, right in sorted(windows):
        if len(periods) > 0 and left <= periods[-1][1]:
            periods[-1][1] = max(periods[-1][1], right)
        else:
            periods.append([left, right])
    if len(periods) == 0:
        return pd.Series([], dtype="interval")
    bounds = timescale.tt_jd(np.array(periods).flatten()).utc_datetime()
    first = pd.Timestamp(start.astimezone(tz=timezone.utc))
    last = pd.Timestamp(end.astimezone(tz=timezone.utc))
    # widen the periods by a millisecond to cover the precision of the
    # (Julian date) times at which they were found
    epsilon = pd.Timedelta(milliseconds=1)
    return pd.Series(
        [
            pd.Interval(
                left=max(pd.Timestamp(left) - epsilon, first),
                right=min(pd.Timestamp(right) + epsilon, last),
                closed="both",
            )
            for left, right in zip(bounds[::2], bounds[1::2])
        ],
        dtype="interval",
    )


def _get_region_elevation(region: geo.Polygon | geo.MultiPolygon) -> float:
    """
    Gets the elevation of a region: the mean of its z coordinates (meters
    above the WGS 84 ellipsoid), or zero if it has none.

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The region.

    Returns:
        float: the elevation (meters)
    """
    if not region.has_z:
        return 0.0
    return float(np.mean(shapely.get_coordinates(region, include_z=True)[:, 2]))


def _get_region_view(
    region: geo.Polygon | geo.MultiPolygon,
    arcs: tuple[np.ndarray, np.ndarray],
    orbit_track: Geocentric,
    nadir_reference: NadirReference,
    elevation: float = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Gets the view from a satellite of a region's point nearest to nadir:
    the nadir point itself, if the region contains it, or otherwise the
    point of its boundary at the least angular distance from the nadir
    point (see `tatc.utils.geometry._get_nearest_arc_points`).

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region, split along the anti-meridian and poles (see
                `tatc.utils.geometry.split_polygon`).
        arcs (tuple[numpy.ndarray, numpy.ndarray]): The region's boundary
                arcs (see `tatc.utils.geometry._get_boundary_arcs`).
        orbit_track (skyfield.positionlib.Geocentric): The satellite orbit track.
        nadir_reference (NadirReference): The definition of the nadir direction.
        elevation (float): The region's elevation (meters) above the WGS 84 ellipsoid.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:
            the point's angle from nadir (degrees), the satellite's elevation
            angle seen from the point (degrees), and the point's longitude
            and latitude (degrees)
    """
    position, nadir, _, _ = _compute_view_frame(
        orbit_track, 0, 0, VelocityFrame.EARTH_FIXED, nadir_reference
    )
    position, nadir = np.reshape(position, (3, -1)), np.reshape(nadir, (3, -1))
    # direction of the nadir point
    if nadir_reference == NadirReference.GEODETIC:
        subpoint = wgs84.geographic_position_of(orbit_track)
        origin = _get_surface_directions(
            np.reshape(subpoint.longitude.degrees, -1),
            np.reshape(subpoint.latitude.degrees, -1),
            elevation,
        )
    else:
        origin = position / np.linalg.norm(position, axis=0)
    longitude, latitude = _get_geodetic_coordinates(origin)
    inside = shapely.contains_xy(region, longitude, latitude)
    _, nearest = _get_nearest_arc_points(origin, arcs)
    nearest[:, inside] = origin[:, inside]
    longitude, latitude = _get_geodetic_coordinates(nearest)
    # line of sight from the satellite to the point
    los = _get_surface_positions(nearest, elevation) - position
    distance = np.linalg.norm(los, axis=0)
    angle = np.degrees(
        np.arccos(np.clip(np.sum(los * nadir, axis=0) / distance, -1, 1))
    )
    angle[inside] = 0
    # satellite elevation angle above the point's (geodetic) horizon
    lon, lat = np.radians(longitude), np.radians(latitude)
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    sat_elevation = np.degrees(
        np.arcsin(np.clip(-np.sum(los * up, axis=0) / distance, -1, 1))
    )
    return angle, sat_elevation, longitude, latitude


def collect_region_observations(
    region: geo.Polygon | geo.MultiPolygon,
    satellite: Satellite,
    start: datetime,
    end: datetime,
    instrument_index: int = 0,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect single satellite observations of a region of interest: a
    shapely `Polygon` or `MultiPolygon` in longitude and latitude (degrees),
    at the elevation of the mean of its z coordinates (meters above the
    WGS 84 ellipsoid), if any, or otherwise zero. The region is split along
    the anti-meridian and poles (see `tatc.utils.geometry.split_polygon`).
    For a point, see `tatc.analysis.point_coverage.collect_observations`.

    Each observation spans a period when any part of the region lies within
    the instrument's field of regard and above the satellite's horizon: when
    the angle from nadir of the region's point nearest to nadir (the nadir
    point itself, if within the region) is at most half the field of regard.
    Only the field of regard is considered, for any instrument: the view
    geometry of a pointed or conical instrument is not. Each observation's
    epoch is the period's midpoint, and its validity (illumination) and
    satellite (and solar) angles refer to the region's point nearest to
    nadir at the epoch.

    If the orbit is propagated with a repeat cycle (see
    `GeneralPerturbationsOrbit.repeat_cycle`), it is modeled as maintained on
    its repeat ground track before its first and after its last element's
    epoch (see `GeneralPerturbationsOrbit.get_orbit_track_at_time`).

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region of interest.
        satellite (Satellite): The observing satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        instrument_index (int): The index of the observing instrument in satellite.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations
            (with a `point_id` of 0).
    """
    _check_satellite(satellite)
    if not isinstance(region, (geo.Polygon, geo.MultiPolygon)):
        raise TypeError(
            "region must be a Polygon or MultiPolygon, not a "
            f"{type(region).__name__} (see collect_observations for a point)"
        )
    elevation = _get_region_elevation(region)
    # split along the anti-meridian and poles (which discards z coordinates),
    # restoring the region's elevation, if any
    geometry = split_polygon(region)
    if region.has_z:
        geometry = project_polygon_to_elevation(geometry, elevation)
    instrument = satellite.instruments[instrument_index]
    orbit = satellite.orbit.to_gp_orbit()
    shapely.prepare(geometry)
    arcs = _get_boundary_arcs(geometry, elevation)
    periods = list(
        _get_visible_polygon_interval_series(
            geometry, satellite, instrument.field_of_regard, start, end, elevation
        )
    )
    # refine the periods to the field of regard: when the angle from nadir
    # of the region's point nearest to nadir is at most half the field of
    # regard, and the satellite is above that point's horizon
    half_angle = instrument.field_of_regard / 2

    def residual(orbit_track: Geocentric) -> np.ndarray:
        angle, sat_elevation, _, _ = _get_region_view(
            geometry, arcs, orbit_track, instrument.nadir_reference, elevation
        )
        return np.maximum(angle - half_angle, -sat_elevation)

    periods = _refine_access_periods(
        residual, orbit, periods, max_step=timedelta(seconds=10)
    )
    observations = []
    for period in periods:
        # instrument validity (illumination) is only checked at each
        # period's epoch (its midpoint), for the region's point nearest to
        # nadir, as an approximation of the whole interval
        epoch = period.mid
        orbit_track = orbit.get_orbit_track([epoch])
        _, _, longitude, latitude = _get_region_view(
            geometry, arcs, orbit_track, instrument.nadir_reference, elevation
        )
        target = wgs84.latlon(latitude[0], longitude[0], elevation)
        if (
            instrument.min_access_time <= period.right - period.left
            and instrument.is_valid_observation(orbit_track, target).all()
        ):
            observations.append((period, epoch, (longitude[0], latitude[0], elevation)))
    return _build_observation_frame(
        observations,
        0,
        geometry,
        satellite,
        instrument,
        omit_solar,
    )


def collect_multi_region_observations(
    region: geo.Polygon | geo.MultiPolygon,
    satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect multiple satellite observations of a region of interest: calls
    `collect_region_observations` for every instrument on every satellite in
    `satellites`, and concatenates the results into one data frame.

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region of interest (see `collect_region_observations`).
        satellites (Satellite | list[Satellite]): The observing satellite(s),
                each contributing an observation per instrument it carries.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with all recorded observations.
    """
    gdfs = [
        collect_region_observations(
            region, satellite, start, end, instrument_index, omit_solar
        )
        for satellite in _check_satellites(satellites)
        for instrument_index in range(len(satellite.instruments))
    ]
    if len(gdfs) == 0:
        # an empty `satellites` list leaves nothing to concatenate
        return _get_empty_coverage_frame(omit_solar)
    # concatenate into one data frame, sort by start time, and re-index
    return pd.concat(gdfs).sort_values("start").reset_index(drop=True)
