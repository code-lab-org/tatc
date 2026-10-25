"""
Methods to compute the geometry of sun glint: the specular reflection of the
Sun on the WGS 84 ellipsoid (as on a calm water surface) toward a satellite.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from skyfield.timelib import Time
from skyfield.toposlib import GeographicPosition

from ..constants import de421
from .ellipsoid import geodetic_to_rectangular, rectangular_to_geodetic


def _unit(vector: npt.NDArray) -> npt.NDArray:
    """
    Normalizes vectors (stacked on the first axis).
    """
    return vector / np.linalg.norm(vector, axis=0)


def _get_sun_position(t: Time) -> npt.NDArray:
    """
    Gets the Earth-fixed (ITRS) position (meters, shape (3, N)) of the Sun
    relative to the Earth's center at times.
    """
    sun = de421["earth"].at(t).observe(de421["sun"]).frame_xyz(itrs).m
    return np.reshape(np.array(sun), (3, -1))


def _get_surface_axes(
    longitude: npt.ArrayLike, latitude: npt.ArrayLike
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Gets the geodetic up (normal), east, and north unit vectors (Earth-fixed,
    shape (3, ...)) at geodetic coordinates (degrees).
    """
    lon, lat = np.radians(longitude), np.radians(latitude)
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), np.zeros(np.shape(lon))])
    north = np.array(
        [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)]
    )
    return up, east, north


def _get_specular_residual(
    longitude: npt.NDArray,
    latitude: npt.NDArray,
    satellite: npt.NDArray,
    sun: npt.NDArray,
) -> npt.NDArray:
    """
    Gets the components (east and north) of the sum of the unit vectors
    toward the satellite and the Sun in the surface's tangent plane, which
    vanish at the specular point (where the surface normal bisects them).
    """
    position = geodetic_to_rectangular(longitude, latitude)
    bisector = _unit(satellite - position) + _unit(sun - position)
    _, east, north = _get_surface_axes(longitude, latitude)
    return np.array([np.sum(bisector * east, axis=0), np.sum(bisector * north, axis=0)])


def compute_specular_point(
    orbit_track: Geocentric, max_iterations: int = 30
) -> GeographicPosition:
    """
    Compute the specular point (the center of the sun glint) seen from a
    satellite: the point on the WGS 84 ellipsoid where the surface normal
    bisects the directions toward the satellite and toward the Sun, so that
    sunlight reflected by a level surface there (for example, calm water)
    reaches the satellite. At the specular point, the Sun and the satellite
    are at the same elevation angle, in opposite azimuths. The point is found
    by Newton's method on the geodetic coordinates, from the point beneath
    the bisector of the directions toward the satellite and the Sun from the
    Earth's center. It is undefined (`NaN`) where the Sun is below the
    horizon at the point, as when the satellite is over the night side of
    the Earth.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        max_iterations (int): The maximum number of iterations.

    Returns:
        skyfield.toposlib.GeographicPosition: the specular point (one for each
            time of the orbit track).
    """
    satellite = np.reshape(np.array(orbit_track.frame_xyz(itrs).m), (3, -1))
    sun = _get_sun_position(orbit_track.t)
    longitude, latitude, _ = rectangular_to_geodetic(
        _unit(_unit(satellite) + _unit(sun)) * 6371e3
    )
    step_size = 1e-6  # degrees, for the finite-difference derivatives
    for _ in range(max_iterations):
        residual = _get_specular_residual(longitude, latitude, satellite, sun)
        d_lon = (
            _get_specular_residual(longitude + step_size, latitude, satellite, sun)
            - residual
        ) / step_size
        d_lat = (
            _get_specular_residual(longitude, latitude + step_size, satellite, sun)
            - residual
        ) / step_size
        # solve the 2 x 2 system [d_lon d_lat] step = -residual at each time
        determinant = d_lon[0] * d_lat[1] - d_lat[0] * d_lon[1]
        with np.errstate(divide="ignore", invalid="ignore"):
            step_lon = (-residual[0] * d_lat[1] + residual[1] * d_lat[0]) / determinant
            step_lat = (-residual[1] * d_lon[0] + residual[0] * d_lon[1]) / determinant
        # limit each step (degrees) for stability far from the solution
        scale = np.minimum(1, 5 / np.maximum(np.hypot(step_lon, step_lat), 1e-12))
        step_lon = np.nan_to_num(step_lon * scale)
        step_lat = np.nan_to_num(step_lat * scale)
        longitude, latitude = longitude + step_lon, np.clip(latitude + step_lat, -90, 90)
        if np.all(np.hypot(step_lon, step_lat) < 1e-9):
            break
    position = geodetic_to_rectangular(longitude, latitude)
    up, _, _ = _get_surface_axes(longitude, latitude)
    # the Sun and the satellite above the horizon (the residual also vanishes
    # where both are below it)
    sunlit = np.sum(up * _unit(sun - position), axis=0) > 0
    visible = np.sum(up * _unit(satellite - position), axis=0) > 0
    converged = np.max(
        np.abs(_get_specular_residual(longitude, latitude, satellite, sun)), axis=0
    ) < 1e-6
    valid = sunlit & visible & converged
    longitude = np.where(valid, (longitude + 180) % 360 - 180, np.nan)
    latitude = np.where(valid, latitude, np.nan)
    if np.ndim(orbit_track.t.tt) == 0:
        return wgs84.latlon(latitude[0], longitude[0])
    return wgs84.latlon(latitude, longitude)


def compute_glint_angle(
    orbit_track: Geocentric, target: GeographicPosition
) -> npt.NDArray:
    """
    Compute the glint angle (degrees) at a target: the angle between the
    direction from the target toward the satellite and the direction of the
    Sun's specular reflection at the target (the direction toward the Sun
    mirrored about the WGS 84 surface normal). It is zero at the specular
    point (see `compute_specular_point`), where a level surface reflects
    sunlight toward the satellite, and grows away from it. Sun glint
    contaminates observations of water (for example, of ocean color) at
    small glint angles, typically below about 20 to 30 degrees, whereas
    glint observations (for example, of greenhouse gases over the ocean)
    seek them.

    Args:
        orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
        target (skyfield.toposlib.GeographicPosition): the target position
            (a single target, or one target for each time of the orbit track).

    Returns:
        numpy.typing.NDArray: the glint angle (degrees).
    """
    satellite = np.reshape(np.array(orbit_track.frame_xyz(itrs).m), (3, -1))
    sun = _get_sun_position(orbit_track.t)
    position = np.reshape(np.array(target.itrs_xyz.m), (3, -1))
    up, _, _ = _get_surface_axes(
        np.reshape(target.longitude.degrees, -1), np.reshape(target.latitude.degrees, -1)
    )
    to_sun = _unit(sun - position)
    reflected = 2 * np.sum(up * to_sun, axis=0) * up - to_sun
    angle = np.degrees(
        np.arccos(
            np.clip(np.sum(_unit(satellite - position) * reflected, axis=0), -1, 1)
        )
    )
    if np.ndim(orbit_track.t.tt) == 0 and np.ndim(target.latitude.degrees) == 0:
        return angle[0]
    return angle
