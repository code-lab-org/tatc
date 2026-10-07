"""
Geometry of the WGS 84 reference ellipsoid: conversions between Earth-fixed
(rectangular) and geodetic coordinates, intersections of rays with the
ellipsoid and its limb, and tangent points (the points of minimum geodetic
altitude) of lines.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from skyfield.framelib import itrs
from skyfield.timelib import Time

from ..constants import EARTH_ECCENTRICITY, EARTH_EQUATORIAL_RADIUS, EARTH_POLAR_RADIUS


def _get_axes(elevation: float = 0) -> npt.NDArray:
    """
    Gets the semi-axes (meters, shape (3, 1)) of the WGS 84 ellipsoid
    extended by an elevation, which approximates the surface at that
    elevation above the ellipsoid.
    """
    return np.array(
        [
            EARTH_EQUATORIAL_RADIUS + elevation,
            EARTH_EQUATORIAL_RADIUS + elevation,
            EARTH_POLAR_RADIUS + elevation,
        ]
    )[:, np.newaxis]


def geodetic_to_rectangular(
    longitude: npt.ArrayLike, latitude: npt.ArrayLike, elevation: npt.ArrayLike = 0
) -> npt.NDArray:
    """
    Converts WGS 84 geodetic coordinates to Earth-fixed (ITRS) rectangular
    coordinates, vectorized across positions.

    Args:
        longitude (numpy.typing.ArrayLike): The geodetic longitudes (degrees).
        latitude (numpy.typing.ArrayLike): The geodetic latitudes (degrees).
        elevation (numpy.typing.ArrayLike): The elevations (meters) above the
            WGS 84 ellipsoid.

    Returns:
        numpy.typing.NDArray: the positions (meters, shape (3, ...))
    """
    lon = np.radians(np.asarray(longitude, dtype=float))
    lat = np.radians(np.asarray(latitude, dtype=float))
    elevation = np.asarray(elevation, dtype=float)
    e2 = EARTH_ECCENTRICITY**2
    n = EARTH_EQUATORIAL_RADIUS / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.array(
        [
            (n + elevation) * np.cos(lat) * np.cos(lon),
            (n + elevation) * np.cos(lat) * np.sin(lon),
            (n * (1 - e2) + elevation) * np.sin(lat),
        ]
    )


def rectangular_to_geodetic(
    position: npt.ArrayLike,
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Converts Earth-fixed (ITRS) rectangular coordinates to WGS 84 geodetic
    coordinates, vectorized across positions, by fixed-point iteration on
    the geodetic latitude (stable at the poles, and converged to well below
    a micrometer for positions near the Earth's surface).

    Args:
        position (numpy.typing.ArrayLike): The positions (meters, shape (3, ...)).

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
            the geodetic longitudes (degrees), latitudes (degrees), and
            altitudes (meters)
    """
    x, y, z = np.asarray(position, dtype=float)
    p = np.hypot(x, y)
    e2 = EARTH_ECCENTRICITY**2
    latitude = np.arctan2(z, p * (1 - e2))
    for _ in range(6):
        sin_latitude = np.sin(latitude)
        n = EARTH_EQUATORIAL_RADIUS / np.sqrt(1 - e2 * sin_latitude**2)
        latitude = np.arctan2(z + e2 * n * sin_latitude, p)
    sin_latitude = np.sin(latitude)
    altitude = (
        p * np.cos(latitude)
        + z * sin_latitude
        - EARTH_EQUATORIAL_RADIUS * np.sqrt(1 - e2 * sin_latitude**2)
    )
    return np.degrees(np.arctan2(y, x)), np.degrees(latitude), altitude


def compute_ellipsoid_intersection(
    position: npt.ArrayLike, direction: npt.ArrayLike, elevation: float = 0
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Computes the first intersections of rays with the WGS 84 ellipsoid, with
    its semi-axes extended by an elevation (which approximates the surface at
    that elevation above the ellipsoid), vectorized across rays.

    Args:
        position (numpy.typing.ArrayLike): The rays' origins (meters,
            Earth-fixed, shape (3, N)).
        direction (numpy.typing.ArrayLike): The rays' directions (Earth-fixed,
            shape (3, N)).
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid.

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the intersection
            points (meters, shape (3, N); the origins for rays that miss) and
            whether each ray intersects the ellipsoid
    """
    position = np.asarray(position, dtype=float)
    direction = np.asarray(direction, dtype=float)
    axes = _get_axes(elevation)
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


def compute_tangent_point(
    position: npt.ArrayLike, direction: npt.ArrayLike, iterations: int = 3
) -> npt.NDArray:
    """
    Computes the tangent point of each line `position + s * direction`: its
    point of minimum WGS 84 geodetic altitude, where the line is tangent to
    the surface of constant geodetic altitude (as for the line of sight of a
    radio occultation or a limb sounder). The tangent point differs from the
    line's closest approach to the Earth's center by up to about 20 km
    horizontally at middle latitudes, but by only meters in altitude.

    Starts from the closest approach to the Earth's center and refines the
    tangent altitude `iterations` times (see `_find_tangent_distance`); three
    refinements converge to well below a millimeter in altitude.

    Args:
        position (numpy.typing.ArrayLike): Points on the lines (meters,
            Earth-fixed, shape (3, N)).
        direction (numpy.typing.ArrayLike): The lines' directions (Earth-fixed,
            shape (3, N)).
        iterations (int): The number of refinements of the tangent altitude.

    Returns:
        numpy.typing.NDArray: the tangent points (meters, Earth-fixed, shape (3, N))
    """
    position = np.asarray(position, dtype=float)
    direction = np.asarray(direction, dtype=float)
    return (
        position + _find_tangent_distance(position, direction, iterations) * direction
    )


def _itrs_rotation(t: Time) -> np.ndarray:
    """
    Gets the GCRS -> ITRS rotation matrices at Skyfield time(s), with a
    trailing axis of length one for a scalar time, so that they broadcast
    against position arrays of shape (3, N) via
    `np.einsum("ij...,j...->i...", rotation, position)`.
    """
    rotation = itrs.rotation_at(t)
    if rotation.ndim == 2:
        rotation = rotation[..., np.newaxis]
    return rotation


def _ellipsoidal_tangent_distance(
    sat_p_itrs: np.ndarray, d_itrs: np.ndarray, target_elevations: np.ndarray
) -> np.ndarray:
    """
    Computes the distance (in units of `d_itrs`) along each ray
    (`sat_p_itrs + s * d_itrs`, Earth-fixed, shape (3, N)) to its
    ellipsoidal tangent point at a known tangent altitude: the point of
    minimum WGS 84 geodetic altitude, where the ray is tangent to the surface
    of constant geodetic altitude rather than to a sphere.

    Uses the closed-form closest approach to the origin after scaling the
    polar (z) axis by `(a + h) / (b + h)`, which maps the ellipsoid of
    semi-axes `(a + h, b + h)` -- closely approximating the surface of
    constant geodetic altitude `h` -- onto a sphere. For tangent altitudes
    between about -250 km and +100 km, this agrees with the exact
    minimum-altitude point to within a fraction of a millimeter in
    altitude.
    """
    k = (EARTH_EQUATORIAL_RADIUS + target_elevations) / (
        EARTH_POLAR_RADIUS + target_elevations
    )
    scale = np.stack([np.ones_like(k), np.ones_like(k), k])
    p_s = sat_p_itrs * scale
    d_s = d_itrs * scale
    return -np.einsum("ij,ij->j", p_s, d_s) / np.einsum("ij,ij->j", d_s, d_s)


def _find_tangent_distance(
    p_itrs: np.ndarray, d_itrs: np.ndarray, iterations: int = 3
) -> np.ndarray:
    """
    Computes the distance (in units of `d_itrs`) along each line
    (`p_itrs + s * d_itrs`, Earth-fixed, shape (3, N)) to its ellipsoidal
    tangent point when the tangent altitude is not known in advance: starts
    from the line's closest approach to the origin and refines the
    altitude-dependent scale factor of `_ellipsoidal_tangent_distance`
    `iterations` times.
    """
    s = -np.einsum("ij,ij->j", p_itrs, d_itrs) / np.einsum("ij,ij->j", d_itrs, d_itrs)
    for _ in range(iterations):
        _, _, altitude = rectangular_to_geodetic(p_itrs + s * d_itrs)
        s = _ellipsoidal_tangent_distance(p_itrs, d_itrs, altitude)
    return s


def _ellipsoidal_tangent_point(
    p: np.ndarray, d: np.ndarray, t: Time, iterations: int = 3
) -> np.ndarray:
    """
    Computes the tangent point (see `compute_tangent_point`) of each line
    `p + s * d` given in the inertial (GCRS) frame (shape (3, N)) at
    Skyfield time(s) `t` (e.g. for a radio occultation's
    receiver-transmitter line), in the same frame.

    Returns:
        numpy.ndarray: tangent point positions (m, shape (3, N), GCRS).
    """
    rotation = _itrs_rotation(t)
    p_itrs = np.einsum("ij...,j...->i...", rotation, p)
    d_itrs = np.einsum("ij...,j...->i...", rotation, d)
    return p + _find_tangent_distance(p_itrs, d_itrs, iterations) * d


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
    axes = _get_axes(elevation)
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


def _get_surface_directions(
    longitude: npt.ArrayLike, latitude: npt.ArrayLike, elevation: float = 0
) -> npt.NDArray[np.float64]:
    """
    Gets the geocentric unit vectors (Earth-fixed, shape (3, N)) toward
    geodetic positions at an elevation above the WGS 84 ellipsoid.

    Args:
        longitude (numpy.typing.ArrayLike): The geodetic longitudes (degrees).
        latitude (numpy.typing.ArrayLike): The geodetic latitudes (degrees).
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid.

    Returns:
        numpy.typing.NDArray[numpy.float64]: the unit vectors
    """
    position = geodetic_to_rectangular(longitude, latitude, elevation)
    return position / np.linalg.norm(position, axis=0)


def _get_surface_positions(
    directions: npt.NDArray[np.float64], elevation: float = 0
) -> npt.NDArray[np.float64]:
    """
    Gets the Earth-fixed positions (meters, shape (3, N)) in geocentric
    directions (unit vectors, shape (3, N)) on the surface of an ellipsoid
    of the WGS 84 semi-axes extended by an elevation (which approximates
    the surface at that elevation above the WGS 84 ellipsoid).

    Args:
        directions (numpy.typing.NDArray[numpy.float64]): The unit vectors (shape (3, N)).
        elevation (float): The elevation (meters) above the WGS 84 ellipsoid.

    Returns:
        numpy.typing.NDArray[numpy.float64]: the positions (meters)
    """
    return directions / np.linalg.norm(directions / _get_axes(elevation), axis=0)


def _get_geodetic_coordinates(
    directions: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """
    Gets the geodetic longitude and latitude (degrees) of the points on the
    WGS 84 ellipsoid in geocentric directions (unit vectors, shape (3, N)).

    Args:
        directions (numpy.typing.NDArray[numpy.float64]): The unit vectors (shape (3, N)).

    Returns:
        tuple[numpy.typing.NDArray[numpy.float64], numpy.typing.NDArray[numpy.float64]]:
            the longitudes and latitudes (degrees)
    """
    longitude = np.degrees(np.arctan2(directions[1], directions[0]))
    latitude = np.degrees(
        np.arctan2(
            directions[2],
            (1 - EARTH_ECCENTRICITY**2) * np.hypot(directions[0], directions[1]),
        )
    )
    return longitude, latitude
