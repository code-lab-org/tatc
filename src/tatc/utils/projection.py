"""
Projection utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Iterable

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


def compute_projected_ray_position(  # pylint: disable=too-many-branches,too-many-statements
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
    Get the location of a projected ray from an instrument. The ray is cast
    from the satellite position toward the WGS 84 geoid at the specified
    elevation; if it misses the geoid entirely (e.g. an off-nadir angle
    pointing past the horizon), the projected position instead falls back
    to the nearest point on the visible Earth limb. Zero roll, pitch, and
    field of view center on the geodetic nadir (the WGS 84 ellipsoid
    surface normal through the satellite), matching Skyfield's
    `wgs84.subpoint_of`/`wgs84.geographic_position_of`.

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
    # convert to radians for internal use
    angle = np.radians(angle)
    # extract earth-fixed position and velocity
    position, velocity = orbit_track.frame_xyz_and_velocity(itrs)
    v_m_per_s = np.array(velocity.m_per_s)
    p_m = np.array(position.m)
    # velocity unit vector
    v = np.divide(v_m_per_s, np.linalg.norm(v_m_per_s, axis=0))
    # nadir unit vector: the geodetic vertical, i.e. the WGS 84
    # ellipsoid surface normal at the sub-satellite point, pointed inward
    # (toward the Earth). This is NOT simply the geocentric direction
    # (-position, toward the Earth's center): the two coincide only at the
    # equator and poles, and otherwise differ by up to the WGS 84
    # geodetic/geocentric latitude discrepancy (~0.19 degrees), which can
    # shift a projected nadir point by kilometers at typical LEO altitudes.
    subpoint = wgs84.geographic_position_of(orbit_track)
    lat = np.array(subpoint.latitude.radians)
    lon = np.array(subpoint.longitude.radians)
    n = -np.array(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
    )
    # whether orbit_track represents a single time or a vector of times
    is_vectorized = len(np.shape(p_m)) > 1
    # cross-track unit vector
    if is_vectorized:
        c = np.cross(v, n, 0, 0, -1).T
    else:
        c = np.cross(v, n)
    # ray pointed at the field of view center (before adding the field of view extent)
    base_ray = n + v * np.tan(np.radians(pitch_angle)) + c * np.tan(np.radians(roll_angle))
    # construct projected ray
    if is_rectangular:
        # find orientation of rectangle corner
        theta = np.arctan(along_track_field_of_view / cross_track_field_of_view)
        # along track half width
        tan_a_2 = np.tan(np.radians(along_track_field_of_view / 2))
        # cross track half width
        tan_c_2 = np.tan(np.radians(cross_track_field_of_view / 2))
        # compose the ray by walking around the rectangle boundary: the (v, c)
        # coefficients are determined together, one segment at a time, rather
        # than by two independently re-derived branch chains. Corners are at
        # theta, pi - theta, pi + theta, and 2*pi - theta; the pi/2, pi, and
        # 3*pi/2 splits are just internal subdivisions of a single flat edge
        # (each formula is continuous across them) chosen to keep every
        # tan() argument close to zero.
        if angle <= theta:
            # right edge, upper half
            v_coef, c_coef = tan_c_2 * np.tan(angle), tan_c_2
        elif angle < np.pi / 2:
            # top edge, right half
            v_coef, c_coef = tan_a_2, tan_a_2 * np.tan(np.pi / 2 - angle)
        elif angle <= np.pi - theta:
            # top edge, left half
            v_coef, c_coef = tan_a_2, -tan_a_2 * np.tan(angle - np.pi / 2)
        elif angle < np.pi:
            # left edge, upper half
            v_coef, c_coef = tan_c_2 * np.tan(np.pi - angle), -tan_c_2
        elif angle < np.pi + theta:
            # left edge, lower half
            v_coef, c_coef = -tan_c_2 * np.tan(angle - np.pi), -tan_c_2
        elif angle < 3 * np.pi / 2:
            # bottom edge, right half
            v_coef, c_coef = -tan_a_2, -tan_a_2 * np.tan(3 * np.pi / 2 - angle)
        elif angle <= 2 * np.pi - theta:
            # bottom edge, left half
            v_coef, c_coef = -tan_a_2, tan_a_2 * np.tan(angle - 3 * np.pi / 2)
        else:
            # right edge, lower half
            v_coef, c_coef = -tan_c_2 * np.tan(2 * np.pi - angle), tan_c_2
        ray = base_ray + v * v_coef + c * c_coef
    else:
        ray = (
            base_ray
            + v * np.sin(angle) * np.tan(np.radians(along_track_field_of_view / 2))
            + c * np.cos(angle) * np.tan(np.radians(cross_track_field_of_view / 2))
        )
    geos = np.zeros_like(p_m)
    for i in range(np.size(geos, axis=1)) if is_vectorized else [-1]:
        _position = p_m[:, i].copy() if is_vectorized else p_m
        _ray = ray[:, i].copy() if is_vectorized else ray
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
            if is_vectorized:
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
            _v = v[:, i].copy() if is_vectorized else v
            _c = c[:, i].copy() if is_vectorized else c
            _, pt_1, pt_2 = inelpl(
                limb,
                nvp2pl(
                    _v * np.sin(np.pi / 2 + angle) + _c * np.cos(np.pi / 2 + angle),
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
            if is_vectorized:
                geos[:, i] = limb_geo
            else:
                geos[:] = limb_geo
    # return resulting geographic position
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
) -> list[Polygon | MultiPolygon]:
    """
    Compute the instantaneous instrument footprint. Supports both a scalar
    (single-time) and vectorized (multi-time) `orbit_track`; the result is
    always a list, with one footprint per time (length 1 for a scalar
    `orbit_track`).

    Args:
        orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
        cross_track_field_of_view (float): The angular (degrees) view orthogonal to velocity.
        along_track_field_of_view (float): The angular (degrees) view in direction of velocity.
        roll_angle (float): The left/right look angle (degrees); right-hand
            rotation about orbit velocity vector.
        pitch_angle (float): The fore/aft look angle (degrees); right-hand
            rotation about orbit normal vector.
        is_rectangular (bool): True, if this is a rectangular sensor.
        number_points (int | None): The required number of polygon points to
            generate: per side for a rectangular sensor, or total for an
            elliptical sensor. Defaults to the runtime configuration.
        elevation (float): The elevation (meters) at which project the footprint.

    Returns:
        list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]: The instrument footprint(s).
    """
    if number_points is None:
        # default number of points
        if is_rectangular:
            number_points = config.get_rc().footprint_points_rectangular_side
        else:
            number_points = config.get_rc().footprint_points_elliptical
    if is_rectangular:
        theta = np.degrees(
            np.arctan(along_track_field_of_view / cross_track_field_of_view)
        )
        angles = np.concatenate(
            (
                np.linspace(-theta, theta, number_points, endpoint=False),
                np.linspace(theta, 180 - theta, number_points, endpoint=False),
                np.linspace(
                    180 - theta, 180 + theta, number_points, endpoint=False
                ),
                np.linspace(
                    180 + theta, 360 - theta, number_points, endpoint=False
                ),
            )
        )
    else:
        angles = np.linspace(0, 360, number_points)
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
    is_vectorized = len(np.shape(orbit_track.t)) > 0  # type: ignore
    return [
        project_polygon_to_elevation(
            split_polygon(
                Polygon(
                    [
                        (
                            point.longitude.degrees[i]
                            if is_vectorized
                            else point.longitude.degrees,
                            point.latitude.degrees[i]
                            if is_vectorized
                            else point.latitude.degrees,
                        )
                        for point in points
                    ]
                )
            ),
            elevation,
        )
        for i in (range(np.size(orbit_track.t)) if is_vectorized else [None])  # type: ignore
    ]


def _compute_limb_for_position(
    position: Iterable[float],
    number_points: int = 16,
    elevation: float = 0,
) -> Polygon | MultiPolygon:
    """
    Compute the visible Earth limb (the outline of the Earth's disk, at the
    specified elevation, as seen from a single observer position) as a
    polygon sampled at evenly spaced angles around the limb ellipse
    returned by SPICE's `edlimb`.

    Args:
        position (Iterable[float]): The observer position (meters), in an
            Earth-fixed frame, from which the limb is visible.
        number_points (int): The required number of polygon points to generate.
        elevation (float): The elevation (meters) at which to project the limb.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The limb.
    """
    limb = edlimb(
        constants.EARTH_EQUATORIAL_RADIUS + elevation,
        constants.EARTH_EQUATORIAL_RADIUS + elevation,
        constants.EARTH_POLAR_RADIUS + elevation,
        position,
    )
    return project_polygon_to_elevation(
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


def compute_limb(
    orbit_track: Geocentric,
    number_points: int = 16,
    elevation: float = 0,
) -> list[Polygon | MultiPolygon]:
    """
    Compute the instantaneous limb (the outline of the visible Earth disk).
    Supports both a scalar (single-time) and vectorized (multi-time)
    `orbit_track`; the result is always a list, with one limb per time
    (length 1 for a scalar `orbit_track`).

    Args:
        orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
        number_points (int): The required number of polygon points to generate.
        elevation (float): The elevation (meters) at which project the limb.

    Returns:
        list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]: The limb(s).
    """
    position, _ = orbit_track.frame_xyz_and_velocity(itrs)
    p_m = np.array(position.m)
    if len(np.shape(p_m)) > 1:
        return [
            _compute_limb_for_position(p_m[:, i].copy(), number_points, elevation)
            for i in range(np.size(p_m, axis=1))
        ]
    return [_compute_limb_for_position(p_m, number_points, elevation)]


def buffer_footprint(
    geometry: Geometry,
    to_crs: Transformer,
    from_crs: Transformer,
    swath_width: float,
    elevation: float,
) -> Polygon | MultiPolygon:
    """
    Buffers a ground track point (in EPSG:4326 coordinates) to create a
    footprint, by reprojecting to a distance-preserving CRS, buffering by
    half the swath width, and reprojecting back.

    Args:
        geometry (shapely.Geometry): The geometry (with EPSG:4326 coordinates) to buffer.
        to_crs (pyproj.Transformer): Transformer from EPSG:4326 to the CRS
            in which to perform the buffer (must use meters).
        from_crs (pyproj.Transformer): Transformer back from the buffering
            CRS to EPSG:4326.
        swath_width (float): The swath width (meters) to buffer.
        elevation (float): The elevation (meters) at which project the buffered polygon.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The buffered footprint.
    """
    # do the swath projection in the specified coordinate reference system
    # split polygons to wrap over the anti-meridian and poles
    # reproject to specified elevation (lost during buffer)
    return project_polygon_to_elevation(
        split_polygon(
            transform(
                from_crs.transform,
                transform(to_crs.transform, geometry).buffer(swath_width / 2),  # type: ignore
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
    Buffers a target geometry to support culling operations. Selects a
    buffer distance equal to the distance traveled in one time step plus
    half of the field of regard swath width, so that no point still
    reachable within the time step is incorrectly excluded. The ground
    distance traveled uses the ground velocity at the equator, which is
    the fastest (most conservative, largest-distance) point for any given
    inclination.

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
    ground_distance = (
        compute_ground_surface_velocity(altitude, 0, inclination) * time_step
    )
    distance = (ground_distance + swath_width / 2) * distance_scaling
    return split_polygon(
        transform(
            from_crs.transform, transform(to_crs.transform, geometry).buffer(distance)  # type: ignore
        )
    )
