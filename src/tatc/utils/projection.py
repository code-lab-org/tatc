"""
Projection utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable
from enum import Enum

import numpy as np
import numpy.typing as npt
from pyproj import Transformer
import shapely
from shapely import Geometry
from shapely.geometry import MultiPolygon, Point, Polygon
from shapely.ops import transform
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition
from spiceypy.spiceypy import edlimb, recgeo

from .. import config, constants
from .geometry import project_polygon_to_elevation, split_polygon
from .observation import field_of_regard_to_swath_width
from .orbital import compute_ground_surface_velocity


class VelocityFrame(str, Enum):
    """
    Enumeration of reference frames for the velocity vector that defines an
    instrument's along-track direction (and hence its cross-track direction,
    orthogonal to velocity and nadir).
    """

    EARTH_FIXED = "earth_fixed"
    """
    Velocity relative to the rotating Earth: the field of view is aligned
    with the ground track, as for a yaw-steered spacecraft (for example,
    zero-Doppler steering for synthetic aperture radar).
    """
    INERTIAL = "inertial"
    """
    Inertial (orbital) velocity: the field of view is aligned with the orbit
    plane, as for a spacecraft without yaw steering. Near the equator, the
    field of view is rotated by up to about 4 degrees (in low Earth orbit)
    relative to the ground track.
    """


class ViewGeometry(str, Enum):
    """
    Enumeration of the geometries in which a pointed instrument's view (its
    fields of view and pixels) is defined.
    """

    FRAME = "frame"
    """
    In a plane perpendicular to the boresight (as for the focal plane of a
    framing camera or a pushbroom array): a rectangle with half widths
    `tan(field_of_view / 2)`, or an ellipse inscribed in it. The angular
    extent along track narrows away from the view center across track (in
    proportion to the cosine of the cross-track angle).
    """
    SCAN = "scan"
    """
    In angles (as for a cross-track scanner or a radar): the pitch angle tilts
    the scan plane (a rotation about the cross-track axis, as for a scanner
    that tilts fore or aft), and each direction is the tilted nadir rotated
    by a cross-track (scan) angle about the tilted along-track axis and then
    by an along-track angle about the rotated cross-track axis, with the
    fields of view as the ranges of these angles about the roll angle and
    the scan plane. The angular extent along track is constant across the
    scan, so the footprint widens along track away from nadir.
    """


class NadirReference(str, Enum):
    """
    Enumeration of definitions of the nadir direction, from which an
    instrument's view is rotated. The two coincide over the equator and the
    poles and otherwise differ by up to about 0.19 degrees (at middle
    latitudes), which shifts a projected view by up to about 2.7 km from
    800 km altitude.
    """

    GEODETIC = "geodetic"
    """
    The geodetic vertical: the WGS 84 ellipsoid surface normal through the
    satellite, toward the Earth (matching Skyfield's sub-satellite point).
    """
    GEOCENTRIC = "geocentric"
    """
    The geocentric direction: from the satellite toward the Earth's center.
    """


def _compute_instrument_frame(
    orbit_track: Geocentric,
    velocity_frame: VelocityFrame = VelocityFrame.EARTH_FIXED,
    nadir_reference: NadirReference = NadirReference.GEODETIC,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Compute the Earth-fixed position and the unit vectors that orient an
    instrument's view: nadir, velocity (along-track), and cross-track
    (velocity cross nadir, which points to the left of the direction of
    motion).

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        velocity_frame (VelocityFrame): The reference frame of the velocity
            vector that defines the along-track direction.
        nadir_reference (NadirReference): The definition of the nadir
            direction.

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
            the position (meters) and the nadir, velocity, and cross-track unit
            vectors, each with shape (3,) or (3, N) in Earth-fixed coordinates.
    """
    # extract earth-fixed position and velocity
    position, velocity = orbit_track.frame_xyz_and_velocity(itrs)
    v_m_per_s = np.array(velocity.m_per_s)
    p_m = np.array(position.m)
    if velocity_frame == VelocityFrame.INERTIAL:
        # add the velocity of the rotating frame (omega x r, with omega along
        # the Earth-fixed z-axis) to express the inertial velocity in
        # Earth-fixed coordinates
        v_m_per_s = v_m_per_s + constants.EARTH_ROTATION_RATE * np.array(
            [-p_m[1], p_m[0], np.zeros_like(p_m[2])]
        )
    # velocity unit vector
    v = np.divide(v_m_per_s, np.linalg.norm(v_m_per_s, axis=0))
    if nadir_reference == NadirReference.GEOCENTRIC:
        # nadir unit vector: toward the Earth's center
        n = -np.divide(p_m, np.linalg.norm(p_m, axis=0))
    else:
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
    # cross-track unit vector
    c = np.cross(v, n, axis=0)
    return p_m, v, n, c


def _compute_view_frame(
    orbit_track: Geocentric,
    roll_angle: float = 0,
    pitch_angle: float = 0,
    velocity_frame: VelocityFrame = VelocityFrame.EARTH_FIXED,
    nadir_reference: NadirReference = NadirReference.GEODETIC,
    tilt_angle: float = 0,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Compute the Earth-fixed position and the orthonormal unit vectors that
    orient a pointed instrument's view: the boresight (view center), the
    along-track axis, and the cross-track axis (to the left of the direction
    of motion). The view is rotated rigidly from nadir: first by the tilt
    angle about the cross-track axis (positive forward), then by the roll
    angle about the (tilted) along-track axis (positive to the left), then
    by the pitch angle about the rolled cross-track axis (positive forward).

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        roll_angle (float): the instrument roll angle (degrees).
        pitch_angle (float): the instrument pitch angle (degrees).
        velocity_frame (VelocityFrame): The reference frame of the velocity
            vector that defines the along-track direction.
        nadir_reference (NadirReference): The definition of the nadir
            direction from which the view is rotated.
        tilt_angle (float): the tilt angle (degrees) of the plane in which
            the roll angle is measured, as for the scan plane of a scanner
            that tilts fore or aft.

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
            the position (meters) and the boresight, along-track, and
            cross-track unit vectors, each with shape (3,) or (3, N) in
            Earth-fixed coordinates.
    """
    p_m, _, n, c = _compute_instrument_frame(
        orbit_track, velocity_frame, nadir_reference
    )
    # orthonormal cross-track (left) and along-track (forward) axes at nadir
    c = c / np.linalg.norm(c, axis=0)
    a = np.cross(n, c, axis=0)
    # tilt about the cross-track axis
    tilt = np.radians(tilt_angle)
    n, a = np.cos(tilt) * n + np.sin(tilt) * a, np.cos(tilt) * a - np.sin(tilt) * n
    # roll about the along-track axis
    roll, pitch = np.radians(roll_angle), np.radians(pitch_angle)
    boresight = np.cos(roll) * n + np.sin(roll) * c
    cross = np.cos(roll) * c - np.sin(roll) * n
    # pitch about the rolled cross-track axis
    boresight, along = (
        np.cos(pitch) * boresight + np.sin(pitch) * a,
        np.cos(pitch) * a - np.sin(pitch) * boresight,
    )
    return p_m, boresight, along, cross


def compute_view_tangents(
    orbit_track: Geocentric,
    target: GeographicPosition,
    velocity_frame: VelocityFrame = VelocityFrame.EARTH_FIXED,
    roll_angle: float = 0,
    pitch_angle: float = 0,
    nadir_reference: NadirReference = NadirReference.GEODETIC,
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Compute the tangents of the along-track and cross-track view angles of a
    target relative to the center of an instrument's view, using the same
    frame as `compute_projected_ray_position`: the line of sight to the
    target is parallel to `boresight + along * along_axis + cross *
    cross_axis`, where the boresight and axes are rotated from nadir by the
    roll and pitch angles. With zero roll and pitch, a target with tangents
    `(along, cross)` lies at the center of a view with roll angle
    `arctan(cross)` and pitch angle `arctan(along / sqrt(1 + cross**2))`.
    Does not check whether the target is above the satellite's horizon.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        target (skyfield.toposlib.GeographicPosition): the target position.
        velocity_frame (VelocityFrame): The reference frame of the velocity
            vector that defines the along-track direction.
        roll_angle (float): the instrument roll angle (degrees).
        pitch_angle (float): the instrument pitch angle (degrees).
        nadir_reference (NadirReference): The definition of the nadir
            direction from which the view is rotated.

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the along-track and
            cross-track (positive left) view angle tangents (`NaN` if the
            target is not in front of the view).
    """
    p_m, boresight, along, cross = _compute_view_frame(
        orbit_track, roll_angle, pitch_angle, velocity_frame, nadir_reference
    )
    # line of sight from the satellite to the target
    los = np.reshape(np.array(target.itrs_xyz.m), (3,) + (1,) * (p_m.ndim - 1)) - p_m
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.sum(los * boresight, axis=0)
        scale = np.where(scale > 0, scale, np.nan)
        return np.sum(los * along, axis=0) / scale, np.sum(los * cross, axis=0) / scale


def compute_view_angles(
    orbit_track: Geocentric,
    target: GeographicPosition,
    velocity_frame: VelocityFrame = VelocityFrame.EARTH_FIXED,
    nadir_reference: NadirReference = NadirReference.GEODETIC,
    tilt_angle: float = 0,
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Compute the roll and pitch angles (degrees) that point a view's center at
    a target, following the rigid rotation of `compute_projected_ray_position`
    (after the tilt angle about the cross-track axis, first by the roll angle
    about the along-track axis, then by the pitch angle about the rolled
    cross-track axis): the angles of the target in a `ViewGeometry.SCAN` view,
    whose scan plane is tilted by the tilt angle (the target's scan angle and
    its along-track angle from the scan plane). Does not check whether the
    target is above the satellite's horizon.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        target (skyfield.toposlib.GeographicPosition): the target position.
        velocity_frame (VelocityFrame): The reference frame of the velocity
            vector that defines the along-track direction.
        nadir_reference (NadirReference): The definition of the nadir
            direction from which the view is rotated.
        tilt_angle (float): the tilt angle (degrees) of the scan plane about
            the cross-track axis (positive forward).

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the roll (positive
            left) and pitch (positive forward) angles (degrees).
    """
    p_m, nadir, along, cross = _compute_view_frame(
        orbit_track, 0, 0, velocity_frame, nadir_reference, tilt_angle
    )
    los = np.reshape(np.array(target.itrs_xyz.m), (3,) + (1,) * (p_m.ndim - 1)) - p_m
    los = los / np.linalg.norm(los, axis=0)
    roll = np.degrees(
        np.arctan2(np.sum(los * cross, axis=0), np.sum(los * nadir, axis=0))
    )
    pitch = np.degrees(np.arcsin(np.clip(np.sum(los * along, axis=0), -1, 1)))
    return roll, pitch


def compute_cone_and_azimuth(
    orbit_track: Geocentric,
    target: GeographicPosition,
    velocity_frame: VelocityFrame = VelocityFrame.EARTH_FIXED,
    nadir_reference: NadirReference = NadirReference.GEODETIC,
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Compute the cone angle (from nadir) and the scan azimuth
    (about the nadir, from the along-track direction, positive to the left of
    the direction of motion) of the line of sight from a satellite to a
    target, as used to describe a conically scanning instrument. Does not
    check whether the target is above the satellite's horizon.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        target (skyfield.toposlib.GeographicPosition): the target position.
        velocity_frame (VelocityFrame): The reference frame of the velocity
            vector that defines the along-track direction.
        nadir_reference (NadirReference): The definition of the nadir
            direction.

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the cone angle
            (degrees, from 0 to 180) and scan azimuth (degrees, from -180 to
            180).
    """
    p_m, nadir, along, cross = _compute_view_frame(
        orbit_track, 0, 0, velocity_frame, nadir_reference
    )
    # line of sight from the satellite to the target
    los = np.reshape(np.array(target.itrs_xyz.m), (3,) + (1,) * (p_m.ndim - 1)) - p_m
    down = np.sum(los * nadir, axis=0)
    forward = np.sum(los * along, axis=0)
    left = np.sum(los * cross, axis=0)
    return (
        np.degrees(np.arctan2(np.hypot(forward, left), down)),
        np.degrees(np.arctan2(left, forward)),
    )


def _intersect_ellipsoid(
    position: npt.NDArray, direction: npt.NDArray, elevation: float = 0
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Computes the first intersections of rays with the ellipsoid of the
    WGS 84 semi-axes extended by an elevation, as SPICE's `surfpt`,
    vectorized across rays.

    Args:
        position (numpy.typing.NDArray): The rays' origins (meters, Earth-fixed, shape (3, N)).
        direction (numpy.typing.NDArray): The rays' directions (Earth-fixed, shape (3, N)).
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid.

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the intersection
            points (meters, shape (3, N); undefined for rays that miss) and
            whether each ray intersects the ellipsoid
    """
    axes = np.array(
        [
            constants.EARTH_EQUATORIAL_RADIUS + elevation,
            constants.EARTH_EQUATORIAL_RADIUS + elevation,
            constants.EARTH_POLAR_RADIUS + elevation,
        ]
    )[:, np.newaxis]
    # in coordinates scaled by the semi-axes, the ellipsoid is the unit sphere
    p, d = position / axes, direction / axes
    a = np.sum(d * d, axis=0)
    b = np.sum(p * d, axis=0)
    c = np.sum(p * p, axis=0) - 1
    discriminant = b**2 - a * c
    root = np.sqrt(np.maximum(discriminant, 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        # the nearer root, from outside (c > 0, toward the ellipsoid: b < 0),
        # in the form that avoids cancellation; or the exit, from inside
        t = np.where(c >= 0, c / (root - b), (root - b) / a)
    found = (discriminant >= 0) & ((c < 0) | (b < 0))
    return position + np.where(found, t, 0) * direction, found


def _rectangular_to_geodetic(
    position: npt.NDArray,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Converts Earth-fixed positions to WGS 84 geodetic coordinates, as SPICE's
    `recgeo`, vectorized across positions, by fixed-point iteration on the
    geodetic latitude (which converges to well below a micrometer within a
    few iterations for positions near the Earth's surface).

    Args:
        position (numpy.typing.NDArray): The positions (meters, shape (3, N)).

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
            the longitudes (radians), latitudes (radians), and altitudes (meters)
    """
    x, y, z = position
    p = np.hypot(x, y)
    e2 = constants.EARTH_ECCENTRICITY**2
    latitude = np.arctan2(z, p * (1 - e2))
    for _ in range(6):
        sin_latitude = np.sin(latitude)
        n = constants.EARTH_EQUATORIAL_RADIUS / np.sqrt(1 - e2 * sin_latitude**2)
        latitude = np.arctan2(z + e2 * n * sin_latitude, p)
    sin_latitude = np.sin(latitude)
    altitude = (
        p * np.cos(latitude)
        + z * sin_latitude
        - constants.EARTH_EQUATORIAL_RADIUS * np.sqrt(1 - e2 * sin_latitude**2)
    )
    return np.arctan2(y, x), latitude, altitude


def _intersect_limb(
    position: npt.NDArray,
    direction: npt.NDArray,
    normal: npt.NDArray,
    elevation: float = 0,
) -> npt.NDArray:
    """
    Computes points of the limb of the ellipsoid of the WGS 84 semi-axes
    extended by an elevation (its points whose lines of sight from a position
    are tangent to it), as for rays that miss it: the intersection of the limb
    with the plane through the position with a normal, choosing of its two
    points the one whose direction from the Earth's center is closer to the
    ray's. Equivalent to SPICE's `edlimb`, `nvp2pl`, and `inelpl`, vectorized
    across rays.

    In coordinates scaled by the semi-axes, the ellipsoid is the unit sphere
    and the limb seen from a scaled position `p` is the circle of points `x`
    on it with `x . p = 1`: centered at `p / |p|^2`, with radius
    `sqrt(1 - 1 / |p|^2)`, in the plane normal to `p`.

    Args:
        position (numpy.typing.NDArray): The rays' origins (meters, Earth-fixed, shape (3, N)).
        direction (numpy.typing.NDArray): The rays' directions (Earth-fixed, shape (3, N)).
        normal (numpy.typing.NDArray): The planes' normals (Earth-fixed, shape (3, N)).
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid.

    Returns:
        numpy.typing.NDArray: the limb points (meters, shape (3, N))
    """
    axes = np.array(
        [
            constants.EARTH_EQUATORIAL_RADIUS + elevation,
            constants.EARTH_EQUATORIAL_RADIUS + elevation,
            constants.EARTH_POLAR_RADIUS + elevation,
        ]
    )[:, np.newaxis]
    p = position / axes
    q = np.sum(p * p, axis=0)
    center = p / q
    radius = np.sqrt(1 - 1 / q)
    # orthonormal basis of the plane of the limb circle
    u = p / np.sqrt(q)
    helper = np.where(np.abs(u[2]) < 0.9, [[0.0], [0.0], [1.0]], [[1.0], [0.0], [0.0]])
    e_1 = np.cross(u, helper, axis=0)
    e_1 /= np.linalg.norm(e_1, axis=0)
    e_2 = np.cross(u, e_1, axis=0)
    # the plane through the position, in scaled coordinates: n . x = d
    n = normal * axes
    d = np.sum(normal * position, axis=0)
    # points of the circle on the plane, at angles phi +/- delta from e_1
    a_1, a_2 = np.sum(n * e_1, axis=0), np.sum(n * e_2, axis=0)
    phi = np.arctan2(a_2, a_1)
    delta = np.arccos(
        np.clip((d - np.sum(n * center, axis=0)) / (radius * np.hypot(a_1, a_2)), -1, 1)
    )
    points = [
        (center + radius * (np.cos(theta) * e_1 + np.sin(theta) * e_2)) * axes
        for theta in (phi + delta, phi - delta)
    ]
    # the point whose direction from the Earth's center is closer to the ray's
    closeness = [
        np.sum(direction * x, axis=0) / np.linalg.norm(x, axis=0) for x in points
    ]
    return np.where(closeness[0] >= closeness[1], points[0], points[1])


def compute_projected_ray_position(  # pylint: disable=too-many-branches,too-many-statements
    orbit_track: Geocentric,
    cross_track_field_of_view: float,
    along_track_field_of_view: float,
    roll_angle: float = 0,
    pitch_angle: float = 0,
    is_rectangular: bool = False,
    angle: float = 0,
    elevation: float = 0,
    velocity_frame: VelocityFrame = VelocityFrame.EARTH_FIXED,
    nadir_reference: NadirReference = NadirReference.GEODETIC,
    tilt_angle: float = 0,
) -> GeographicPosition:
    """
    Get the location of a projected ray from an instrument. The ray is cast
    from the satellite position toward the WGS 84 geoid at the specified
    elevation; if it misses the geoid entirely (e.g. an off-nadir angle
    pointing past the horizon), the projected position instead falls back
    to the nearest point on the visible Earth limb. Zero roll, pitch, and
    field of view center on nadir: by default, the geodetic nadir (the WGS 84
    ellipsoid surface normal through the satellite), matching Skyfield's
    `wgs84.subpoint_of`/`wgs84.geographic_position_of`.

    The view is rotated rigidly from nadir, first by the roll angle about
    the along-track axis, then by the pitch angle about the rolled
    cross-track axis. Its shape is defined in the plane perpendicular to the
    boresight at unit distance (as for the focal plane of a camera): a
    rectangle with half widths `tan(field_of_view / 2)`, or an ellipse
    inscribed in it.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        cross_track_field_of_view (float): the instrument cross-track
            (orthogonal to velocity vector) field of view (degrees).
        along_track_field_of_view (float): the instrument along-track
            (parallel to velocity vector) field of view (degrees).
        roll_angle (float): the instrument roll angle (degrees), a rotation
            about the along-track axis, positive to the left of the
            direction of motion (right-hand about the velocity vector).
        pitch_angle (float): the instrument pitch angle (degrees), a rotation
            about the rolled cross-track axis, positive forward (in the
            direction of motion).
        is_rectangular (bool): `True` if the instrument view has a rectangular
            shape (otherwise elliptical).
        angle (float): ray angle (degrees) about the view center, measured
            from the cross-track axis on the left of the direction of motion
            (0 degrees) toward the along-track axis forward (90 degrees).
        elevation (float): The elevation (meters) at which project the footprint.
        velocity_frame (VelocityFrame): The reference frame of the velocity
            vector that defines the along-track direction.
        nadir_reference (NadirReference): The definition of the nadir
            direction from which the view is rotated.
        tilt_angle (float): a rotation (degrees) about the cross-track axis,
            positive forward, applied before the roll angle (tilting the
            plane in which the roll angle is measured, as for the scan plane
            of a `ViewGeometry.SCAN` view).

    Returns:
        (skyfield.toposlib.GeographicPosition): the geographic position of the projected ray
    """
    # convert to radians for internal use
    angle = np.radians(angle)
    # ray pointed at the field of view center (before adding the field of
    # view extent) and the along-track (v) and cross-track (c) axes of the view
    p_m, base_ray, v, c = _compute_view_frame(
        orbit_track,
        roll_angle,
        pitch_angle,
        velocity_frame,
        nadir_reference,
        tilt_angle,
    )
    # whether orbit_track represents a single time or a vector of times
    is_vectorized = len(np.shape(p_m)) > 1
    # construct projected ray
    if is_rectangular:
        # along track half width
        tan_a_2 = np.tan(np.radians(along_track_field_of_view / 2))
        # cross track half width
        tan_c_2 = np.tan(np.radians(cross_track_field_of_view / 2))
        # find orientation of rectangle corner (the ratio of the half-width
        # tangents, not of the fields of view, which differ for wide views)
        theta = np.arctan2(tan_a_2, tan_c_2)
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
    # find the intersections of the rays and the WGS 84 geoid (at the
    # elevation), vectorized across times
    position = np.reshape(p_m, (3, -1))
    direction = np.reshape(ray, (3, -1))
    points, found = _intersect_ellipsoid(position, direction, elevation)
    geos = np.array(_rectangular_to_geodetic(points))
    miss = ~found
    if np.any(miss):
        # projected points do not fall on the WGS 84 geoid surface: use the
        # points of its limb on the plane of each ray, normal to the ray's
        # angle about the view center
        normal = np.reshape(v, (3, -1)) * np.sin(np.pi / 2 + angle) + np.reshape(
            c, (3, -1)
        ) * np.cos(np.pi / 2 + angle)
        geos[:, miss] = np.array(
            _rectangular_to_geodetic(
                _intersect_limb(
                    position[:, miss], direction[:, miss], normal[:, miss], elevation
                )
            )
        )
    if not is_vectorized:
        geos = geos[:, 0]
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
    velocity_frame: VelocityFrame = VelocityFrame.EARTH_FIXED,
    nadir_reference: NadirReference = NadirReference.GEODETIC,
    view_geometry: ViewGeometry = ViewGeometry.FRAME,
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
        roll_angle (float): The left/right look angle (degrees), a rotation
            about the along-track axis, positive to the left of the direction
            of motion.
        pitch_angle (float): The fore/aft look angle (degrees), a rotation
            about the rolled cross-track axis, positive forward (in
            `ViewGeometry.SCAN`, the tilt of the scan plane: a rotation about
            the cross-track axis applied before the roll angle).
        is_rectangular (bool): True, if this is a rectangular sensor.
        number_points (int | None): The required number of polygon points to
            generate: per side for a rectangular sensor, or total for an
            elliptical sensor. Defaults to the runtime configuration.
        elevation (float): The elevation (meters) at which project the footprint.
        velocity_frame (VelocityFrame): The reference frame of the velocity
            vector that defines the along-track direction.
        nadir_reference (NadirReference): The definition of the nadir
            direction from which the view is rotated.
        view_geometry (ViewGeometry): The geometry in which the fields of view
            are defined.

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
        # space points evenly along each side of the rectangle (in the plane
        # one unit from the instrument), starting from the lower right corner;
        # evenly spaced clock angles would instead cluster points near the
        # middle of the long sides of an elongated view
        tan_a_2 = np.tan(np.radians(along_track_field_of_view / 2))
        tan_c_2 = np.tan(np.radians(cross_track_field_of_view / 2))
        s = np.linspace(-1, 1, number_points, endpoint=False)
        theta = np.degrees(np.arctan2(tan_a_2, tan_c_2))
        angles = np.degrees(
            np.concatenate(
                (
                    np.arctan2(s * tan_a_2, tan_c_2),  # right side
                    np.arctan2(tan_a_2, -s * tan_c_2),  # top side
                    np.arctan2(-s * tan_a_2, -tan_c_2),  # left side
                    np.arctan2(-tan_a_2, s * tan_c_2),  # bottom side
                )
            )
        )
        # wrap to the range [-theta, 360 - theta) expected for clock angles
        angles = np.where(angles < -theta, angles + 360, angles)
    else:
        angles = np.linspace(0, 360, number_points)
    if view_geometry == ViewGeometry.SCAN:
        # perimeter in (cross-track, along-track) angles, counterclockwise from
        # the left side as for the clock angles above, with points evenly
        # spaced in angle along each side
        half_c, half_a = cross_track_field_of_view / 2, along_track_field_of_view / 2
        if is_rectangular:
            s = np.linspace(-1, 1, number_points, endpoint=False)
            offsets = np.concatenate(
                (
                    np.stack([np.full_like(s, half_c), s * half_a], axis=1),
                    np.stack([-s * half_c, np.full_like(s, half_a)], axis=1),
                    np.stack([np.full_like(s, -half_c), -s * half_a], axis=1),
                    np.stack([s * half_c, np.full_like(s, -half_a)], axis=1),
                )
            )
        else:
            theta = np.radians(np.linspace(0, 360, number_points))
            offsets = np.stack([half_c * np.cos(theta), half_a * np.sin(theta)], axis=1)
        # rays at the angular offsets from the roll angle and the scan plane,
        # which the pitch angle tilts about the cross-track axis
        points = [
            compute_projected_ray_position(
                orbit_track,
                0,
                0,
                roll_angle + cross_offset,
                along_offset,
                False,
                0,
                elevation,
                velocity_frame,
                nadir_reference,
                pitch_angle,
            )
            for cross_offset, along_offset in offsets
        ]
    else:
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
                velocity_frame,
                nadir_reference,
            )
            for angle in angles
        ]
    # footprint polygons (one per time), built at once from their points
    longitude = np.reshape(
        [point.longitude.degrees for point in points], (len(points), -1)
    ).T
    latitude = np.reshape(
        [point.latitude.degrees for point in points], (len(points), -1)
    ).T
    footprints = shapely.polygons(np.stack([longitude, latitude], axis=-1))
    # split footprints that cross the anti-meridian or exceed the poles, and
    # repair invalid ones (see split_polygon), checking all footprints at once
    ring = np.concatenate([longitude, longitude[:, :1]], axis=1)
    planar = (
        np.all(np.abs(longitude) <= 180, axis=1)
        & np.all(np.abs(latitude) <= 90, axis=1)
        & np.all(np.abs(np.diff(ring, axis=1)) <= 180, axis=1)
    )
    for i in np.flatnonzero(~(planar & shapely.is_valid(footprints))):
        footprints[i] = split_polygon(footprints[i])
    # project the footprints to the elevation
    return list(shapely.force_3d(footprints, elevation))


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
    distance_crs: str | None = None,
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

    Deprecated: buffering in a planar (equidistant cylindrical) projection
    does not conservatively bound observations near the poles or across the
    anti-meridian. The analysis functions instead cull to the periods when
    an instrument's field of regard may observe a region (see
    `tatc.analysis.region_coverage._get_visible_polygon_interval_series`). This
    function will be removed in a future release.

    Args:
        geometry (shapely.Geometry): The target geometry (with EPSG:4326 coordinates) to buffer.
        altitude (float): The spacecraft orbit altitude (meters).
        inclination (float): The spacecraft orbit inclination (degrees).
        field_of_regard (float): The spacecraft instrument field of regard (degrees).
        time_step (float): The simulation time step (seconds).
        distance_crs (str | None): The coordinate reference system in which to
            perform distance calculations. Defaults to `None`, which builds an
            equidistant cylindrical projection whose true-scale parallel
            (`lat_ts`) is set to `geometry`'s own most poleward latitude. This
            keeps the buffer conservative (true ground distance >= the
            requested distance) at every latitude within `geometry`.
        distance_scaling (float): A multiplicative scaling factor to adjust the buffer
            distance (default: 1.0).

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The buffered geometry.
    """
    warnings.warn(
        "buffer_target is deprecated and will be removed in a future release",
        DeprecationWarning,
        stacklevel=2,
    )
    if distance_crs is None:
        _, min_lat, _, max_lat = geometry.bounds
        lat_ts = min(max(abs(min_lat), abs(max_lat)), 89.9)
        distance_crs = f"+proj=eqc +lat_ts={lat_ts} +datum=WGS84 +units=m"
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
