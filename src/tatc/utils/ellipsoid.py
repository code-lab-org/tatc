"""
Tangent point geometry relative to the WGS 84 reference ellipsoid, shared by
the limb sounding and radio occultation (RO) coverage analyses.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
from skyfield.framelib import itrs
from skyfield.timelib import Time

from ..constants import EARTH_ECCENTRICITY, EARTH_EQUATORIAL_RADIUS, EARTH_POLAR_RADIUS


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


def _geodetic_altitude(position_itrs: np.ndarray) -> np.ndarray:
    """
    Computes the WGS 84 geodetic altitude (m) of Earth-fixed (ITRS)
    positions (m, shape (3, N)), by fixed-point iteration on geodetic
    latitude (stable at the poles, and converged to well below a
    millimeter within a few iterations for near-Earth positions).
    """
    x, y, z = position_itrs
    p = np.hypot(x, y)
    e2 = EARTH_ECCENTRICITY**2
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(5):
        sin_lat = np.sin(lat)
        n = EARTH_EQUATORIAL_RADIUS / np.sqrt(1 - e2 * sin_lat**2)
        lat = np.arctan2(z + e2 * n * sin_lat, p)
    sin_lat = np.sin(lat)
    return (
        p * np.cos(lat)
        + z * sin_lat
        - EARTH_EQUATORIAL_RADIUS * np.sqrt(1 - e2 * sin_lat**2)
    )


def _ellipsoidal_tangent_distance(
    sat_p_itrs: np.ndarray, d_itrs: np.ndarray, target_elevations: np.ndarray
) -> np.ndarray:
    """
    Computes the distance (in units of `d_itrs`) along each ray
    (`sat_p_itrs + s * d_itrs`, Earth-fixed, shape (3, N)) to its
    ellipsoidal tangent point: the point of minimum WGS 84 geodetic
    altitude, where the ray is tangent to the surface of constant geodetic
    altitude rather than to a sphere.

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


def _ellipsoidal_tangent_point(
    p: np.ndarray, d: np.ndarray, t: Time, iterations: int = 3
) -> np.ndarray:
    """
    Computes the ellipsoidal tangent point (the point of minimum WGS 84
    geodetic altitude) of each line `p + s * d` (inertial GCRS, shape
    (3, N)) at Skyfield time(s) `t`, when the tangent altitude is not known
    in advance (e.g. for a radio occultation's receiver-transmitter line).

    Starts from the line's geocentric closest approach to the origin and
    refines the altitude-dependent scale factor of
    `_ellipsoidal_tangent_distance` `iterations` times; three refinements
    converge to well below a millimeter in altitude.

    Returns:
        numpy.ndarray: tangent point positions (m, shape (3, N), GCRS).
    """
    rotation = _itrs_rotation(t)
    p_itrs = np.einsum("ij...,j...->i...", rotation, p)
    d_itrs = np.einsum("ij...,j...->i...", rotation, d)
    s = -np.einsum("ij,ij->j", p, d) / np.einsum("ij,ij->j", d, d)
    for _ in range(iterations):
        altitude = _geodetic_altitude(p_itrs + s * d_itrs)
        s = _ellipsoidal_tangent_distance(p_itrs, d_itrs, altitude)
    return p + s * d
