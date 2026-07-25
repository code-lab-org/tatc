"""
Projection utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""
from __future__ import annotations

import numpy as np
from pyproj import Transformer
from shapely import Geometry
from shapely.geometry import (
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.ops import transform
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition
from spiceypy.spiceypy import edlimb, inelpl, nvp2pl, recgeo, surfpt
from spiceypy.utils.exceptions import NotFoundError

from .. import config, constants
from .geometry import project_polygon_to_elevation, split_polygon
from .observation import field_of_regard_to_swath_width
from .orbital import compute_ground_surface_velocity


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
    number_points: int | None = None,
    elevation: float = 0,
) -> Geometry | list[Geometry]:
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
        shapely.Geometry | list[shapely.Geometry]: The instrument footprint(s).
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
) -> Geometry | list[Geometry]:
    """
    Compute the instanteous limb.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
        number_points (int): The required number of polygon points to generate.
        elevation (float): The elevation (meters) at which project the limb.

    Returns:
        shapely.Geometry | list[shapely.Geometry]: The limb(s).
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
        shapely.geometry.Polygon: The buffered footprint.
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
) -> Polygon | MultiPolygon:
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
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The buffered geometry.
    """
    to_crs = Transformer.from_crs("EPSG:4326", distance_crs, always_xy=True)
    from_crs = Transformer.from_crs(distance_crs, "EPSG:4326", always_xy=True)
    swath_width = field_of_regard_to_swath_width(altitude, field_of_regard)
    ground_distance = compute_ground_surface_velocity(altitude, 0, inclination) * time_step
    distance = (ground_distance + swath_width / 2) * distance_scaling
    return split_polygon(
        transform(
            from_crs.transform, transform(to_crs.transform, geometry).buffer(distance)
        )
    )
