# -*- coding: utf-8 -*-
"""
Utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""
import re
from typing import List, Union

import numpy as np
from numba import njit
import geopandas as gpd
from pyproj import Transformer
from sgp4.conveniences import sat_epoch_datetime
from sgp4.api import Satrec
from shapely import Geometry, make_valid, simplify
from shapely.geometry import (
    Point,
    Polygon,
    MultiPolygon,
    GeometryCollection,
    LineString,
)
from shapely.ops import split, transform
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.toposlib import GeographicPosition
from skyfield.positionlib import Geocentric
from spiceypy.spiceypy import edlimb, inelpl, nvp2pl, recgeo, surfpt
from spiceypy.utils.exceptions import NotFoundError

from . import constants
from . import config


@njit
def mean_anomaly_to_true_anomaly(mean_anomaly: float, eccentricity: float = 0) -> float:
    """
    Converts mean anomaly to true anomaly.

    Args:
        mean_anomaly (float): The mean anomaly (degrees).
        true_anomaly (float): The orbit eccentricity.

    Returns:
        float: The true anomaly (degrees).
    """
    mean_anomaly_rad = np.radians(mean_anomaly)
    true_anomaly_rad = (
        mean_anomaly_rad
        + (2 * eccentricity - (1 / 4) * eccentricity**3) * np.sin(mean_anomaly_rad)
        + (5 / 4) * eccentricity**2 * np.sin(2 * mean_anomaly_rad)
        + (13 / 12) * eccentricity**3 * np.sin(3 * mean_anomaly_rad)
    )
    return np.degrees(true_anomaly_rad)


@njit
def true_anomaly_to_mean_anomaly(true_anomaly: float, eccentricity: float = 0) -> float:
    """
    Converts true anomaly to mean anomaly.

    Args:
        true_anomaly (float): The true anomaly (degrees).
        eccentricity (float): The orbit eccentricity.

    Returns:
        float: The mean anomaly (degrees).
    """
    true_anomaly_rad = np.radians(true_anomaly)
    mean_anomaly_rad = (
        true_anomaly_rad
        - 2 * eccentricity * np.sin(true_anomaly_rad)
        + ((3 / 4) * eccentricity**2 + (1 / 8) * eccentricity**4)
        * np.sin(2 * true_anomaly_rad)
        - (1 / 3) * eccentricity**3 * np.sin(3 * true_anomaly_rad)
        + (5 / 32) * eccentricity**4 * np.sin(4 * true_anomaly_rad)
    )
    return np.degrees(mean_anomaly_rad)


@njit
def compute_number_samples(distance: float) -> int:
    """
    Compute the number of global samples required to achieve a typical
    sample distance (meters) assuming equal spacing.

    Args:
        distance (float): The typical distance between samples (meters).

    Returns:
        int: The number of global samples.
    """
    # compute the angular distance of each sample (assuming mean sphere)
    theta = distance / constants.EARTH_MEAN_RADIUS
    # compute the distance from the center of earth to conic plane (assuming sphere)
    radius = constants.EARTH_MEAN_RADIUS * np.cos(theta / 2)
    # compute the distance from the conic plane to the surface (assuming sphere)
    height = constants.EARTH_MEAN_RADIUS - radius
    # compute the sperical cap area covered by the sample (assuming sphere)
    # https://en.wikipedia.org/wiki/Spherical_cap
    sample_area = 2 * np.pi * constants.EARTH_MEAN_RADIUS * height
    # return the fraction of earth-to-sample area
    return int(constants.EARTH_SURFACE_AREA / sample_area)


@njit
def swath_width_to_field_of_regard(
    altitude: float, swath_width: float, elevation: float = 0
) -> float:
    """
    Fast conversion from swath width to field of regard.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        swath_width (float): Observation diameter (meters) at specified elevation.
        elevation (float): Elevation (meters) above WGS 84 datum to observe.

    Returns:
        float: The field of regard (degrees).
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = (constants.EARTH_MEAN_RADIUS + elevation) / (
        constants.EARTH_MEAN_RADIUS + altitude
    )
    # lambda is the Earth central angle
    sin_lambda = np.sin((swath_width / 2) / (constants.EARTH_MEAN_RADIUS + elevation))
    # eta is the angular radius of the region viewable by the satellite
    tan_eta = sin_rho * sin_lambda / (1 - sin_rho * np.cos(np.arcsin(sin_lambda)))
    return np.degrees(2 * np.arctan(tan_eta))


@njit
def swath_width_to_field_of_view(
    altitude: float, swath_width: float, look_angle: float = 0, elevation: float = 0
) -> float:
    """
    Fast conversion from swath width to field of view considering off-nadir pointing.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        swath_width (float): Observation diameter (meters) at specified elevation.
        look_angle (float): Off-nadir look angle (degrees) to observation center.
        elevation (float): Elevation (meters) above WGS 84 datum to observe.

    Returns:
        float: The field of view (degrees).
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = (constants.EARTH_MEAN_RADIUS + elevation) / (
        constants.EARTH_MEAN_RADIUS + altitude
    )
    # eta is the angular radius from sub-satellite point to center of view
    sin_eta = min(sin_rho, np.sin(np.radians(look_angle) / 2))
    # epsilon is the satellite elevation from the center of view
    cos_epsilon = sin_eta / sin_rho
    # lambda is the Earth central angle to the center of view
    _lambda = np.pi / 2 - np.arcsin(sin_eta) - np.arccos(cos_epsilon)
    sin_lambda_1 = np.sin(
        _lambda - (swath_width / 2) / (constants.EARTH_MEAN_RADIUS + elevation)
    )
    sin_lambda_2 = np.sin(
        _lambda + (swath_width / 2) / (constants.EARTH_MEAN_RADIUS + elevation)
    )
    # eta is the angular radius of the region viewable by the satellite
    tan_eta_1 = sin_rho * sin_lambda_1 / (1 - sin_rho * np.cos(np.arcsin(sin_lambda_1)))
    tan_eta_2 = sin_rho * sin_lambda_2 / (1 - sin_rho * np.cos(np.arcsin(sin_lambda_2)))
    return np.degrees(np.arctan(tan_eta_2) - np.arctan(tan_eta_1))


@njit
def field_of_regard_to_swath_width(
    altitude: float, field_of_regard: float, elevation: float = 0
) -> float:
    """
    Fast conversion from field of regard to swath width.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        field_of_regard (float): Angular width (degrees) of observation.
        elevation (float): Elevation (meters) above WGS 84 datum to observe.

    Returns:
        float: The observation diameter (meters) at the specified elevation.
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = (constants.EARTH_MEAN_RADIUS + elevation) / (
        constants.EARTH_MEAN_RADIUS + altitude
    )
    # eta is the angular radius of the region viewable by the satellite
    sin_eta = min(sin_rho, np.sin(np.radians(field_of_regard) / 2))
    # epsilon is the min satellite elevation for obs (grazing angle)
    cos_epsilon = sin_eta / sin_rho
    # lambda is the Earth central angle
    _lambda = np.pi / 2 - np.arcsin(sin_eta) - np.arccos(cos_epsilon)
    return 2 * (constants.EARTH_MEAN_RADIUS + elevation) * _lambda


@njit
def compute_field_of_regard(
    altitude: float, min_elevation_angle: float, elevation: float = 0
) -> float:
    """
    Fast computation of field of regard for observation with a minimum altitude angle.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        min_elevation_angle (float): The minimum elevation angle (degrees) for observation.
        elevation (float): Elevation (meters) above WGS 84 datum to observe.

    Returns:
        float: Angular width (degrees) of observation.
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = (constants.EARTH_MEAN_RADIUS + elevation) / (
        constants.EARTH_MEAN_RADIUS + altitude
    )
    # epsilon is the min satellite elevation for obs (grazing angle)
    cos_epsilon = np.cos(np.radians(min_elevation_angle))
    # eta is the angular radius of the region viewable by the satellite
    sin_eta = sin_rho * cos_epsilon
    return np.degrees(np.arcsin(sin_eta) * 2)


@njit
def compute_min_elevation_angle(
    altitude: float, field_of_regard: float, elevation: float = 0
) -> float:
    """
    Fast computation of minimum elevation angle required to observe a point.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        field_of_regard (float): Angular width (degrees) of observation.
        elevation (float): Elevation (meters) above WGS 84 datum to observe.

    Returns:
        float: The minimum elevation angle (degrees) for observation.
    """
    # eta is the angular radius of the region viewable by the satellite
    sin_eta = np.sin(np.radians(field_of_regard) / 2)
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = (constants.EARTH_MEAN_RADIUS + elevation) / (
        constants.EARTH_MEAN_RADIUS + altitude
    )
    # epsilon is the min satellite elevation for obs (grazing angle)
    cos_epsilon = sin_eta / sin_rho
    if cos_epsilon > 1:
        return 0
    return np.degrees(np.arccos(cos_epsilon))


@njit
def compute_orbit_period(altitude: float) -> float:
    """
    Fast computation of approximate orbital period.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.

    Returns:
        float: The orbital period (seconds).
    """
    semimajor_axis = constants.EARTH_MEAN_RADIUS + altitude
    mean_motion_rad_s = np.sqrt(constants.EARTH_MU / semimajor_axis**3)
    return 2 * np.pi / mean_motion_rad_s


@njit
def compute_max_access_time(altitude: float, min_elevation_angle: float) -> float:
    """
    Fast computation of maximum access time to observe a point.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        min_elevation_angle (float): Minimum elevation angle (degrees) for observation.

    Returns:
        float: The maximum access time (seconds) for observation.
    """
    orbital_distance = (constants.EARTH_MEAN_RADIUS + altitude) * (
        np.pi - 2 * np.radians(min_elevation_angle)
    )
    orbital_velocity = np.sqrt(
        constants.EARTH_MU / (constants.EARTH_MEAN_RADIUS + altitude)
    )
    return orbital_distance / orbital_velocity


@njit
def compute_ground_velocity(altitude: float, inclination: float) -> float:
    """
    Fast computation of mean ground velocity for a nadir-pointing instrument.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        inclination (float): Inclination (degrees) of the observing instrument orbit.

    Returns:
        float: The access time (seconds) for observation.
    """
    semimajor_axis = constants.EARTH_MEAN_RADIUS + altitude
    mean_motion_rad_s = np.sqrt(constants.EARTH_MU / semimajor_axis**3)
    return constants.EARTH_MEAN_RADIUS * (
        mean_motion_rad_s
        - (2 * np.pi * np.cos(np.degrees(inclination)) / constants.EARTH_SIDEREAL_DAY_S)
    )


@njit
def along_track_distance_to_access_time(
    altitude: float, inclination: float, along_track: float
) -> float:
    """
    Fast computation of mean access time for a specified along track distance.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        inclination (float): Inclination (degrees) of the observing instrument orbit.
        along_track (float): Along track distance (meters) observed during access.

    Returns:
        float: The access time (seconds) for observation.
    """
    semimajor_axis = constants.EARTH_MEAN_RADIUS + altitude
    mean_motion_rad_s = np.sqrt(constants.EARTH_MU / semimajor_axis**3)
    ground_velocity = constants.EARTH_MEAN_RADIUS * (
        mean_motion_rad_s
        - (2 * np.pi * np.cos(np.degrees(inclination)) / constants.EARTH_SIDEREAL_DAY_S)
    )
    return along_track / ground_velocity


@njit
def access_time_to_along_track_distance(
    altitude: float, inclination: float, access_time: float
) -> float:
    """
    Fast computation of along track distance for a specified access time.

    Args:
        altitude (float): Altitude (meters) above WGS 84 datum for the observing instrument.
        inclination (float): Inclination (degrees) of the observing instrument orbit.
        access_time (float): Access time (seconds) during observation.

    Returns:
        float: The observation along track distance (meters).
    """
    semimajor_axis = constants.EARTH_MEAN_RADIUS + altitude
    mean_motion_rad_s = np.sqrt(constants.EARTH_MU / semimajor_axis**3)
    ground_velocity = constants.EARTH_MEAN_RADIUS * (
        mean_motion_rad_s
        - (2 * np.pi * np.cos(np.degrees(inclination)) / constants.EARTH_SIDEREAL_DAY_S)
    )
    return ground_velocity * access_time


def compute_projected_ray_position(
    orbit_track: Geocentric,
    cross_track_field_of_view: float,
    along_track_field_of_view: float,
    roll_angle: float = 0,
    pitch_angle: float = 0,
    is_rectangular: bool = False,
    angle: float = 0,
    elevation: float = 0,
) -> GeographicPosition:
    """
    Get the location of a projected ray from an instrument.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        cross_track_field_of_view (float): the instrument cross-track
            (orthogonal to velocity vector) field of view (degrees).
        along_track_field_of_view (float): the instrument along-track
            (parallel to velocity vector) field of view (degrees).
        roll_angle (float): the instrument roll (right-hand about
            velocity vector) angle (degrees).
        pitch_angle (float): the instrument pitch (right-hand about
            orbit normal vector) angle (degrees).
        is_rectangular (bool): `True` if the instrument view has a rectangular
            shape (otherwise elliptical).
        angle (float): ray angle (degrees) counterclockwise from right-hand cross-track
            direction about the instrument field of view.
        elevation (float): The elevation (meters) at which project the footprint.

    Returns:
        (skyfield.toposlib.GeographicPosition): the geographic position of the projected ray
    """
    # extract earth-fixed position and velocity
    position, velocity = orbit_track.frame_xyz_and_velocity(itrs)
    # velocity unit vector
    v = np.divide(velocity.m_per_s, np.linalg.norm(velocity.m_per_s, axis=0))
    # binormal unit vector
    b = np.divide(-position.m, np.linalg.norm(position.m, axis=0))
    # normal unit vector
    if len(np.shape(position.m)) > 1:
        n = np.cross(v, b, 0, 0, -1).T
    else:
        n = np.cross(v, b)
    # construct projected ray
    if is_rectangular:
        # find orientation of rectangle corner
        theta = np.arctan(along_track_field_of_view / cross_track_field_of_view)
        # along track half width
        tan_a_2 = np.tan(np.radians(along_track_field_of_view / 2))
        # cross track half width
        tan_c_2 = np.tan(np.radians(cross_track_field_of_view / 2))
        # compose the ray with different equations for each side
        ray = (
            b
            + v * np.tan(np.radians(pitch_angle))
            + n * np.tan(np.radians(roll_angle))
            + v
            * (
                tan_a_2
                if theta <= angle <= np.pi - theta
                else (
                    -tan_a_2
                    if np.pi + theta <= angle <= 2 * np.pi - theta
                    else (
                        tan_c_2 * np.tan(angle)
                        if angle < theta
                        else (
                            tan_c_2 * np.tan(np.pi - angle)
                            if angle < np.pi
                            else (
                                -tan_c_2 * np.tan(angle - np.pi)
                                if angle < np.pi + theta
                                else -tan_c_2 * np.tan(2 * np.pi - angle)
                            )
                        )
                    )
                )
            )
            + n
            * (
                tan_c_2
                if angle <= theta or angle >= 2 * np.pi - theta
                else (
                    -tan_c_2
                    if np.pi - theta <= angle <= np.pi + theta
                    else (
                        tan_a_2 * np.tan(np.pi / 2 - angle)
                        if angle < np.pi / 2
                        else (
                            -tan_a_2 * np.tan(angle - np.pi / 2)
                            if angle < np.pi - theta
                            else (
                                -tan_a_2 * np.tan(3 * np.pi / 2 - angle)
                                if angle < 3 * np.pi / 2
                                else tan_a_2 * np.tan(angle - 3 * np.pi / 2)
                            )
                        )
                    )
                )
            )
        )
    else:
        ray = (
            b
            + v * np.tan(np.radians(pitch_angle))
            + n * np.tan(np.radians(roll_angle))
            + v * np.sin(angle) * np.tan(np.radians(along_track_field_of_view / 2))
            + n * np.cos(angle) * np.tan(np.radians(cross_track_field_of_view / 2))
        )
    geos = np.zeros_like(position.m)
    for i in range(np.size(geos, axis=1)) if geos.ndim > 1 else [-1]:
        _position = position.m[:, i].copy() if i >= 0 else position.m
        _ray = ray[:, i].copy() if i >= 0 else ray
        # find the intersection of the ray and the WGS 84 geoid
        try:
            pt = surfpt(
                _position,
                _ray,
                constants.EARTH_EQUATORIAL_RADIUS + elevation,
                constants.EARTH_EQUATORIAL_RADIUS + elevation,
                constants.EARTH_POLAR_RADIUS + elevation,
            )
            geo = recgeo(
                pt,
                constants.EARTH_EQUATORIAL_RADIUS,
                constants.EARTH_FLATTENING,
            )
            if i >= 0:
                geos[:, i] = geo
            else:
                geos[:] = geo
        except NotFoundError:
            # projected point does not fall on the WGS 84 geoid surface
            # compute the observable limb ellipse for WGS 84 geoid
            limb = edlimb(
                constants.EARTH_EQUATORIAL_RADIUS + elevation,
                constants.EARTH_EQUATORIAL_RADIUS + elevation,
                constants.EARTH_POLAR_RADIUS + elevation,
                _position,
            )
            # find the two intersection points between orthogonal plane and limb ellipse
            _v = v[:, i].copy() if i >= 0 else v
            _n = n[:, i].copy() if i >= 0 else n
            _, pt_1, pt_2 = inelpl(
                limb,
                nvp2pl(
                    _v * np.sin(np.pi / 2 + angle) + _n * np.cos(np.pi / 2 + angle),
                    _position,
                ),
            )
            # compute the angles between the ray and limb intersection points
            angle_1 = np.arccos(
                np.dot(_ray, pt_1) / np.linalg.norm(_ray) / np.linalg.norm(pt_1)
            )
            angle_2 = np.arccos(
                np.dot(_ray, pt_2) / np.linalg.norm(_ray) / np.linalg.norm(pt_2)
            )
            # use the limb intersection point closer to the ray
            limb_pt = pt_1 if angle_1 <= angle_2 else pt_2
            limb_geo = recgeo(
                limb_pt,
                constants.EARTH_EQUATORIAL_RADIUS,
                constants.EARTH_FLATTENING,
            )
            if i >= 0:
                geos[:, i] = limb_geo
            else:
                geos[:] = limb_geo
    # return resulting geographic position
    if len(np.shape(geos)) > 1:
        return wgs84.latlon(np.degrees(geos[1, :]), np.degrees(geos[0, :]), geos[2, :])
    return wgs84.latlon(np.degrees(geos[1]), np.degrees(geos[0]), geos[2])


def compute_footprint(
    orbit_track: Geocentric,
    cross_track_field_of_view: float,
    along_track_field_of_view: float,
    roll_angle: float = 0,
    pitch_angle: float = 0,
    is_rectangular: bool = False,
    number_points: int = None,
    elevation: float = 0,
) -> Union[Geometry, List[Geometry]]:
    """
    Compute the instanteous instrument footprint.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
        cross_track_field_of_view (float): The angular (degrees) view orthogonal to velocity.
        along_track_field_of_view (float): The angular (degrees) view in direction of velocity.
        pitch_angle (float): The fore/aft look angle (degrees); right-hand
            rotation about orbit normal vector.
        roll_angle (float): The left/right look angle (degrees); right-hand
            rotation about orbit velocity vector.
        is_rectangular (float): True, if this is a rectangular sensor.
        number_points (int): The required number of polygon points to generate.
        elevation (float): The elevation (meters) at which project the footprint.

    Returns:
        Union[shapely.Geometry, List[shapely.Geometry]: The instrument footprint(s).
    """
    if number_points is None:
        # default number of points
        if is_rectangular:
            number_points = config.rc.footprint_points_rectangular_side
        else:
            number_points = config.rc.footprint_points_elliptical
    if is_rectangular:
        theta = np.arctan(along_track_field_of_view / cross_track_field_of_view)
        angles = np.concatenate(
            (
                np.linspace(-theta, theta, number_points, endpoint=False),
                np.linspace(theta, np.pi - theta, number_points, endpoint=False),
                np.linspace(
                    np.pi - theta, np.pi + theta, number_points, endpoint=False
                ),
                np.linspace(
                    np.pi + theta, 2 * np.pi - theta, number_points, endpoint=False
                ),
            )
        )
    else:
        angles = np.linspace(0, 2 * np.pi, number_points)
    points = [
        compute_projected_ray_position(
            orbit_track,
            cross_track_field_of_view,
            along_track_field_of_view,
            roll_angle,
            pitch_angle,
            is_rectangular,
            angle,
            elevation,
        )
        for angle in angles
    ]
    if np.size(orbit_track.t) > 1:
        return [
            project_polygon_to_elevation(
                split_polygon(
                    Polygon(
                        [
                            (point.longitude.degrees[i], point.latitude.degrees[i])
                            for point in points
                        ]
                    )
                ),
                elevation,
            )
            for i in range(np.size(orbit_track.t))
        ]
    return project_polygon_to_elevation(
        split_polygon(
            Polygon(
                [(point.longitude.degrees, point.latitude.degrees) for point in points]
            )
        ),
        elevation,
    )


def compute_limb(
    orbit_track: Geocentric,
    number_points: int = 16,
    elevation: float = 0,
) -> Union[Geometry, List[Geometry]]:
    """
    Compute the instanteous limb.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
        number_points (int): The required number of polygon points to generate.
        elevation (float): The elevation (meters) at which project the limb.

    Returns:
        Union[shapely.Geometry, List[shapely.Geometry]: The limb(s).
    """
    position, _ = orbit_track.frame_xyz_and_velocity(itrs)
    polygons = [None] * np.size(position.m, axis=1) if position.m.ndim > 1 else None
    for i in range(len(polygons)) if position.m.ndim > 1 else [-1]:
        _position = position.m[:, i].copy() if i >= 0 else position.m
        limb = edlimb(
            constants.EARTH_EQUATORIAL_RADIUS + elevation,
            constants.EARTH_EQUATORIAL_RADIUS + elevation,
            constants.EARTH_POLAR_RADIUS + elevation,
            _position,
        )
        polygon = project_polygon_to_elevation(
            split_polygon(
                Polygon(
                    [
                        Point(np.degrees(g[0]), np.degrees(g[1]))
                        for p in [
                            limb.center
                            + np.cos(i) * limb.semi_major
                            + np.sin(i) * limb.semi_minor
                            for i in np.linspace(0, np.pi * 2, number_points)
                        ]
                        for g in [
                            recgeo(
                                p,
                                constants.EARTH_EQUATORIAL_RADIUS,
                                constants.EARTH_FLATTENING,
                            )
                        ]
                    ]
                )
            ),
            elevation,
        )
        if i >= 0:
            polygons[i] = polygon
        else:
            polygons = polygon
    return polygons


def buffer_footprint(
    geometry: Geometry,
    to_crs: Transformer,
    from_crs: Transformer,
    swath_width: float,
    elevation: float,
) -> Polygon:
    """
    Buffers a ground track point to create a footprint.

    Args:
        geometry (shapely.Geometry): The geometry to buffer.
        origin_crs (str): The origin coordinate reference system (CRS).
        buffer_crs (str): The buffering coordinate reference system (CRS).
        swath_width (float): The swath width (meters) to buffer.
        elevation (float): The elevation (meters) at which project the buffered polygon.

    Returns:
        Union[shapely.geometry.Polygon, shapely.geometry.MultiPolygon]: The buffered footprint.
    """
    # do the swath projection in the specified coordinate reference system
    # split polygons to wrap over the anti-meridian and poles
    # reproject to specified elevation (lost during buffer)
    return project_polygon_to_elevation(
        split_polygon(
            transform(
                from_crs.transform,
                transform(to_crs.transform, geometry).buffer(swath_width / 2),
            )
        ),
        elevation,
    )


def buffer_target(
    geometry: Geometry,
    altitude: float,
    inclination: float,
    field_of_regard: float,
    time_step: float,
    distance_crs: str = "EPSG:4087",
    distance_scaling: float = 1.0,
):
    """
    Buffers a target geometry to support culling operations. Selects a buffer distance
    equal to the distance traveled in one time step plus half of the field of regard
    swath width. Simplifies geometries to distance tolerances within 5% of the buffer distance.

    Args:
        geometry (shapely.Geometry): The target geometry (with EPSG:4326 coordinates) to buffer.
        altitude (float): The spacecraft orbit altitude (meters).
        inclination (float): The spacecraft orbit inclination (degrees).
        field_of_regard (float): The spacecraft instrument field of regard (degrees).
        time_step (float): The simulation time step (seconds).
        distance_crs (str): The coordinate reference system in which to perform
            distance calculations (default: EPSG:4087).
        distance_scaling (float): A multiplicative scaling factor to adjust the buffer
            distance (default: 1.0).

    Returns:
        Union[shapely.geometry.Polygon, shapely.geometry.MultiPolygon]: The buffered geometry.
    """
    to_crs = Transformer.from_crs("EPSG:4326", distance_crs, always_xy=True)
    from_crs = Transformer.from_crs(distance_crs, "EPSG:4326", always_xy=True)
    swath_width = field_of_regard_to_swath_width(altitude, field_of_regard)
    ground_distance = compute_ground_velocity(altitude, inclination) * time_step
    distance = (ground_distance + swath_width / 2) * distance_scaling
    return split_polygon(
        transform(
            from_crs.transform, transform(to_crs.transform, geometry).buffer(distance)
        )
    )


def project_polygon_to_elevation(
    polygon: Union[Polygon, MultiPolygon], elevation: float
) -> Union[Polygon, MultiPolygon]:
    """
    Projects a polygon to a specified elevation (z-coordinate).

    Args:
        polygon (Polygon or MultiPolygon): The polygon to project.
        elevation (float): The elevation (meters) above the WGS 84 geoid.

    Returns:
        Polygon or MultiPolygon: The projected polygon.
    """
    if isinstance(polygon, Polygon):
        return Polygon(
            [(p[0], p[1], elevation) for p in polygon.exterior.coords],
            [[(p[0], p[1], elevation) for p in i.coords] for i in polygon.interiors],
        )
    return MultiPolygon(
        [project_polygon_to_elevation(g, elevation) for g in polygon.geoms]
    )


def _wrap_polygon_over_north_pole(
    polygon: Union[Polygon, MultiPolygon],
) -> Union[Polygon, MultiPolygon]:
    """
    Wraps polygon coordinates over the North pole. Due to buffering and projection,
    sometimes latitudes exceed 90 degrees. This method wraps them to the correct
    latitude between -90 and 90 degrees and adjusts the longitude by 180 degrees.
    This method requires a polygon above 90 degrees latitude to be only on one
    side of the prime meridian.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to wrap.

    Returns:
       Polygon, or MultiPolygon: The wrapped polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[1] <= 90 for c in polygon.exterior.coords):
            # no wrapping necessary
            return polygon
        # map latitudes from [90, 180) to [90, -90), adjusting longitude by 180 degrees
        lat_shift = 180 if all(c[0] <= 0 for c in polygon.exterior.coords) else -180
        pgon = Polygon(
            [
                [
                    c[0] + lat_shift if c[1] >= 90 else c[0],
                    180 - c[1] if c[1] >= 90 else c[1],
                ]
                for c in polygon.exterior.coords
            ],
            [
                [
                    [
                        c[0] + lat_shift if c[1] >= 90 else c[0],
                        180 - c[1] if c[1] >= 90 else c[1],
                    ]
                    for c in i.coords
                ]
                for i in polygon.interiors
            ],
        )
        # give up and return original polygon if invalid
        if not pgon.is_valid:
            return polygon
        return pgon
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        polygons = [_wrap_polygon_over_north_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in polygons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _split_polygon_north_pole(
    polygon: Union[Polygon, MultiPolygon],
) -> Union[Polygon, MultiPolygon]:
    """
    Splits a Polygon into a MultiPolygon if it crosses north pole.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to split.

    Returns:
       Polygon, or MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[1] <= 90 for c in polygon.exterior.coords):
            # no splitting necessary
            return polygon
        # split polygon along north pole
        parts = split(polygon, LineString([(-360, 90), (360, 90)]))
        # check and split part over prime meridian if necessary
        for part in parts.geoms:
            if part.crosses(LineString([(0, 90), (0, 180)])):
                parts = GeometryCollection(
                    [g for g in parts.geoms if g != part]
                    + [g for g in split(part, LineString([(0, 90), (0, 180)])).geoms]
                )
        # convert to a multi polygon
        if isinstance(parts, GeometryCollection):
            parts = _convert_collection_to_polygon(parts)
        # return polygon with components wrapped over north pole
        return _wrap_polygon_over_north_pole(parts)
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        pgons = [_split_polygon_north_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in pgons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _wrap_polygon_over_south_pole(
    polygon: Union[Polygon, MultiPolygon],
) -> Union[Polygon, MultiPolygon]:
    """
    Wraps polygon coordinates over the South pole. Due to buffering and projection,
    sometimes latitudes exceed -90 degrees. This method wraps them to the correct
    latitude between -90 and 90 degrees and adjusts the longitude by 180 degrees.
    This method requires a polygon above 90 degrees latitude to be only on one
    side of the prime meridian.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to wrap.

    Returns:
       Polygon, or MultiPolygon: The wrapped polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[1] >= -90 for c in polygon.exterior.coords):
            # no splitting necessary
            return polygon
        # map latitudes from [-90, -180) to [-90, 90), adjusting longitude by 180 degrees
        lat_shift = 180 if all(c[0] <= 0 for c in polygon.exterior.coords) else -180
        pgon = Polygon(
            [
                [
                    c[0] + lat_shift if c[1] <= -90 else c[0],
                    -180 - c[1] if c[1] <= -90 else c[1],
                ]
                for c in polygon.exterior.coords
            ],
            [
                [
                    [
                        (c[0] + lat_shift if c[1] <= -90 else c[0],),
                        -180 - c[1] if c[1] <= -90 else c[1],
                    ]
                    for c in i.coords
                ]
                for i in polygon.interiors
            ],
        )
        # give up and return original polygon if invalid
        if not pgon.is_valid:
            return polygon
        return pgon
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        polygons = [_wrap_polygon_over_south_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in polygons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _split_polygon_south_pole(
    polygon: Union[Polygon, MultiPolygon],
) -> Union[Polygon, MultiPolygon]:
    """
    Splits a Polygon into a MultiPolygon if it crosses south pole.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to split.

    Returns:
       Polygon, or MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        lat = np.array([c[1] for c in polygon.exterior.coords])
        if np.all(lat >= -90):
            return polygon
        # split polygon along south pole
        parts = split(polygon, LineString([(-360, -90), (360, -90)]))
        # check and split part over prime meridian if necessary
        for part in parts.geoms:
            if part.crosses(LineString([(0, -90), (0, -180)])):
                parts = GeometryCollection(
                    [g for g in parts.geoms if g != part]
                    + [g for g in split(part, LineString([(0, -90), (0, -180)])).geoms]
                )
        # convert to a multi polygon
        if isinstance(parts, GeometryCollection):
            parts = _convert_collection_to_polygon(parts)
        # return polygon with components wrapped over south pole
        return _wrap_polygon_over_south_pole(parts)
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        pgons = [_split_polygon_south_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in pgons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _wrap_polygon_over_antimeridian(
    polygon: Union[Polygon, MultiPolygon],
) -> Union[Polygon, MultiPolygon]:
    """
    Wraps polygon coordinates over the antimeridian. Due to buffering and projection,
    sometimes longitudes exceed 180 degrees. This method wraps them to
    the correct longitude between -180 and 180 degrees.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to wrap.

    Returns:
       Polygon, or MultiPolygon: The wrapped polygon.
    """
    if isinstance(polygon, Polygon):
        if all(c[0] >= -180 and c[0] <= 180 for c in polygon.exterior.coords):
            # no wrapping necessary
            return polygon
        if all(c[0] <= -180 for c in polygon.exterior.coords):
            # map longitudes from (-540, -180] to (-180, 180]
            pgon = Polygon(
                [[c[0] + 360, c[1]] for c in polygon.exterior.coords],
                [[[c[0] + 360, c[1]] for c in i.coords] for i in polygon.interiors],
            )
        if all(c[0] >= 180 for c in polygon.exterior.coords):
            # map longitudes from [180, 540) to [-180, 180)
            pgon = Polygon(
                [[c[0] - 360, c[1]] for c in polygon.exterior.coords],
                [[[c[0] - 360, c[1]] for c in i.coords] for i in polygon.interiors],
            )
        # give up and return original polygon if invalid
        if not pgon.is_valid:
            return polygon
        return pgon
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        pgons = [_wrap_polygon_over_antimeridian(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in pgons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def _convert_collection_to_polygon(
    collection: GeometryCollection,
) -> Union[Polygon, MultiPolygon]:
    """
    Converts a GeometryCollection to a Polygon or MultiPolygon. Quick clipping
    can create dirty results with points or lines on boundaries. This method
    drops and lines or points from a GeometryCollection to return only the
    Polygon or MultiPolygon geometry.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to convert.

    Returns:
       Polygon, or MultiPolygon: The converted polygon.
    """
    pgons = [p for p in collection.geoms if isinstance(p, Polygon)] + [
        p for g in collection.geoms if isinstance(g, MultiPolygon) for p in g.geoms
    ]
    if len(pgons) == 1:
        return pgons[0]
    return MultiPolygon(pgons)


def _split_polygon_antimeridian(
    polygon: Union[Polygon, MultiPolygon],
) -> Union[Polygon, MultiPolygon]:
    """
    Splits a Polygon into a MultiPolygon if it crosses the anti-meridian after
    wrapping its coordinates using `wrap_coordinates_antimeridian`. Note: this
    function only supports polygons that span LESS than 360 degrees longitude.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to split.

    Returns:
       Polygon, or MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        lon = np.array([c[0] for c in polygon.exterior.coords])
        # check if any longitudes cross the anti-meridian
        # (adjacent coordinate longitude differs by more than 180 degrees)
        if all(np.abs(np.diff(lon)) < 180):
            return polygon
        # check if this polygon contains a pole
        if Polygon(zip(np.cos(np.radians(lon)), np.sin(np.radians(lon)))).contains(
            Point(0, 0)
        ):
            # extract and sort coords by longitude
            coords = polygon.exterior.coords[0:-1]
            coords.sort(key=lambda r: r[0])
            # determine if contains north or south pole based on sign of mean latitude
            n_s = 1 if np.array(coords)[:, 1].mean() > 0 else -1
            # interpolate latitude at antimeridian
            lat = np.interp(
                180, [coords[-1][0], coords[0][0] + 180], [coords[-1][1], coords[0][1]]
            )
            # reconstruct polygon (ccw) with added coords on antimeridian
            # TODO potential problem if provided polygon has z-dimension
            pgon = Polygon(
                [(-180, 90 * n_s), (-180, lat)]
                + coords
                + [(180, lat), (180, 90 * n_s), (-180, 90 * n_s)],
                polygon.interiors,
            )
            # return polygon split down prime meridian to improve handling
            parts = split(pgon, LineString([(0, -180), (0, 180)]))
            # convert to multi polygon
            if isinstance(parts, GeometryCollection):
                parts = _convert_collection_to_polygon(parts)
            return parts
        # find anti-meridian crossings and calculate shift direction
        # coords from W -> E (shift < 0) will add 360 degrees to E component
        # coords from E -> W (shift > 0) will subtract 360 degrees from W component
        shift = np.insert(np.cumsum(np.around(np.diff(lon) / 360)), 0, 0)
        pgon = Polygon(
            [
                (c[0] - 360 * shift[i], c[1])
                for i, c in enumerate(polygon.exterior.coords)
            ],
            [
                [
                    (
                        ic[0]
                        - 360 * np.interp(ic[0], np.sort(lon), shift[np.argsort(lon)]),
                        ic[1],
                    )
                    for ic in i.coords
                ]
                for i in polygon.interiors
            ],
        )
        # split along the anti-meridian (-180 for shift > 0; 180 for shift < 0)
        shift_dir = -180 if shift.max() >= 1 else 180
        parts = split(pgon, LineString([(shift_dir, -180), (shift_dir, 180)]))
        # convert to multi polygon
        if isinstance(parts, GeometryCollection):
            parts = _convert_collection_to_polygon(parts)
        # return polygon with components wrapped over anti-meridian
        return _wrap_polygon_over_antimeridian(parts)
    if isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        pgons = [_split_polygon_antimeridian(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in pgons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    raise ValueError("Unknown geometry: " + str(type(polygon)))


def split_polygon(
    polygon: Union[Polygon, MultiPolygon],
) -> Union[Polygon, MultiPolygon]:
    """
    Splits a Polygon into a MultiPolygon if it crosses the anti-meridian
    (180 degrees longitude), exceeds the north pole (90 degrees latitude), or
    exceeds the south pole (-90 degrees latitude). Note: this function
    only supports polygons that span LESS than 360 degrees longitude.

    Args:
        polygon (Polygon or MultiPolygon): The polygon to split.

    Returns:
        Polygon, or MultiPolygon: The split polygon.
    """
    polygon = _split_polygon_north_pole(
        _split_polygon_south_pole(_split_polygon_antimeridian(polygon))
    )
    # invalid polygons can arise from narrow sensor geometries in polar regions
    if not polygon.is_valid:
        # try to fix geometry
        polygon = make_valid(polygon)
        if isinstance(polygon, GeometryCollection):
            polygon = _convert_collection_to_polygon(polygon)
    return polygon


def normalize_geometry(
    geometry: Union[Polygon, MultiPolygon, gpd.GeoDataFrame],
) -> gpd.GeoDataFrame:
    """
    Normalize geometry to a GeoDataFrame with antimeridian wrapping.

    Args:
        geometry (geopandas.GeoDataFrame, geopandas.GeoSeries, Polygon, or MultiPolygon): The geometry to normalize.

    Returns:
        geopandas.GeoDataFrame: The normalized geometry.
    """
    if isinstance(geometry, (Polygon, MultiPolygon)):
        if not geometry.is_valid:
            raise ValueError("Geometry is not a valid Polygon or MultiPolygon.")
        geometry = gpd.GeoDataFrame(geometry=gpd.GeoSeries([geometry]), crs="EPSG:4326")
    elif isinstance(geometry, gpd.GeoSeries):
        geometry = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")
    if isinstance(geometry, gpd.GeoDataFrame):
        geometry["geometry"] = geometry.apply(
            lambda r: split_polygon(r.geometry),
            axis=1,
        )
    return geometry


def zero_pad(object_name: str, max_number: int, current_number: int) -> str:
    """
    Uses length of max_number to zero pad allowing for alphanumeric sorting.

    Args:
        object_name (str): Object name, to be concatenated with zero padded number.
        max_number (int): Maximum number, utilized as reference for zero padding.
        current_number (int): Index number to be zero padded.

    Returns:
        str: The object name with zero padded number appended.
    """
    max_length = len(str(max_number))
    return object_name + " " + str(current_number).zfill(max_length)


def is_even_length_list(v: List) -> List:
    """
    Validator to check for even-length lists.
    """
    if len(v) % 2 == 1:
        raise ValueError(f"{v} does not have even length")
    return v


def is_valid_tle(v: List[str]) -> List[str]:
    """
    Validate the two line element set.
    """
    for i in range(0, len(v), 2):
        # based on orekit's TLE.isFormatOK function
        if len(v[i]) != 69:
            raise ValueError(f"Invalid tle: line {i+1} incorrect length.")
        if len(v[i + 1]) != 69:
            raise ValueError(f"Invalid tle: line {i+2} incorrect length.")

        line_1_pattern = (
            r"1 [ 0-9A-HJ-NP-Z][ 0-9]{4}[A-Z] [ 0-9]{5}[ A-Z]{3} "
            + r"[ 0-9]{5}[.][ 0-9]{8} (?:(?:[ 0+-][.][ 0-9]{8})|(?: "
            + r"[ +-][.][ 0-9]{7})) [ +-][ 0-9]{5}[+-][ 0-9] "
            + r"[ +-][ 0-9]{5}[+-][ 0-9] [ 0-9] [ 0-9]{4}[ 0-9]"
        )
        if re.match(line_1_pattern, v[i]) is None:
            raise ValueError(f"Invalid tle: line {i+1} does not match pattern.")
        line_2_pattern = (
            r"2 [ 0-9A-HJ-NP-Z][ 0-9]{4} [ 0-9]{3}[.][ 0-9]{4} "
            + r"[ 0-9]{3}[.][ 0-9]{4} [ 0-9]{7} [ 0-9]{3}[.][ 0-9]{4} "
            + r"[ 0-9]{3}[.][ 0-9]{4} [ 0-9]{2}[.][ 0-9]{13}[ 0-9]"
        )
        if re.match(line_2_pattern, v[i + 1]) is None:
            raise ValueError(f"Invalid tle: line {i+2} does not match pattern.")

        def checksum(line):
            the_sum = 0
            for j in range(68):
                if line[j].isdigit():
                    the_sum += int(line[j])
                elif line[j] == "-":
                    the_sum += 1
            return the_sum % 10

        if int(v[i][68]) != checksum(v[i]):
            raise ValueError(f"Invalid tle: line {i+1} checksum failed.")
        if int(v[i + 1][68]) != checksum(v[i + 1]):
            raise ValueError(f"Invalid tle: line {1+2} checksum failed.")
    return v


def is_chronological_tle(v: List[str]) -> List[str]:
    """
    Validator to check for chronological TLEs.
    """
    epochs = np.array(
        [
            sat_epoch_datetime(Satrec.twoline2rv(v[i], v[i + 1]))
            for i in range(0, len(v), 2)
        ]
    )
    if not all(epochs[:-1] <= epochs[1:]):
        raise ValueError(f"Invalid multi-tle: not in chronological order.")
    return v
