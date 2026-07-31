"""
Orbital utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
from numba import njit

from .. import constants


@njit
def compute_orbit_inertial_velocity(mean_altitude: float) -> float:
    """
    Fast computation of orbit inertial velocity (the orbital speed relative
    to a non-rotating, Earth-centered frame) assuming a circular orbit
    around a spherical Earth (using the mean Earth radius).

    Args:
        mean_altitude (float): Orbit mean altitude (meters) above the mean Earth radius.

    Returns:
        float: Inertial orbit velocity (meters/second).
    """
    return np.sqrt(constants.EARTH_MU / (constants.EARTH_MEAN_RADIUS + mean_altitude))


@njit
def compute_ground_inertial_velocity(
    mean_altitude: float, elevation: float = 0
) -> float:
    """
    Fast computation of the inertial velocity (relative to a non-rotating,
    Earth-centered frame) of the sub-satellite point projected onto a
    sphere at the specified elevation, assuming a circular orbit around a
    spherical Earth (using the mean Earth radius). Since the projected
    point shares the satellite's angular velocity, its linear velocity
    scales with its (smaller, or larger if elevation exceeds mean_altitude)
    radius relative to the orbit's radius.

    Args:
        mean_altitude (float): Orbit mean altitude (meters) above the mean Earth radius.
        elevation (float): Surface elevation (meters) above the mean Earth
            radius at which to project the ground point.

    Returns:
        float: Ground point inertial velocity (meters/second).
    """
    v_orbital = compute_orbit_inertial_velocity(mean_altitude)
    return (
        v_orbital
        * (constants.EARTH_MEAN_RADIUS + elevation)
        / (constants.EARTH_MEAN_RADIUS + mean_altitude)
    )


@njit
def compute_ground_surface_velocity(
    mean_altitude: float,
    elevation: float = 0,
    inclination: float = 0,
    latitude: float = 0,
) -> float:
    """
    Fast computation of ground surface velocity (relative to the rotating
    Earth's surface) assuming a circular orbit around a spherical Earth
    (using the mean Earth radius).

    Args:
        mean_altitude (float): Orbit mean altitude (meters) above the mean Earth radius.
        elevation (float): Surface elevation (meters) above the mean Earth
            radius at which to project the ground point.
        inclination (float): Orbit inclination (degrees).
        latitude (float): Surface latitude (degrees) at which to evaluate the ground velocity.

    Returns:
        float: Ground surface velocity (meters/second), relative to the rotating Earth.
    """
    v_inertial = compute_ground_inertial_velocity(mean_altitude, elevation)
    # compute flight path angle beta relative to due North
    sin_beta = np.cos(np.deg2rad(inclination)) / np.cos(np.deg2rad(latitude))
    sin_beta = max(-1.0, min(1.0, sin_beta))
    beta = np.arcsin(sin_beta)
    v_surface_north = v_inertial * np.cos(beta)
    v_surface_east = v_inertial * np.sin(beta)
    v_earth_west = (
        2
        * np.pi
        / constants.EARTH_SIDEREAL_DAY_S
        * constants.EARTH_MEAN_RADIUS
        * np.cos(np.deg2rad(latitude))
    )
    return np.sqrt((v_surface_east - v_earth_west) ** 2 + v_surface_north**2)


@njit
def semimajor_axis_to_mean_motion(semimajor_axis: float) -> float:
    """
    Fast computation of mean motion (average angular rate) from Kepler's
    third law, assuming a circular orbit around a spherical Earth.

    Args:
        semimajor_axis (float): Orbit semimajor axis (meters).

    Returns:
        float: Orbit mean motion (degrees/second).
    """
    return np.degrees(np.sqrt(constants.EARTH_MU / semimajor_axis**3))


@njit
def mean_motion_to_orbit_period(mean_motion: float) -> float:
    """
    Fast computation of orbital period: the time (360 degrees of mean
    motion) to complete one revolution. This function is its own inverse:
    calling it again on a period (seconds) recovers the mean motion
    (degrees/second), since both are 360 divided by the other.

    Args:
        mean_motion (float): Orbit mean motion (degrees/second).

    Returns:
        float: Orbital period (seconds).
    """
    return 360 / mean_motion


@njit
def mean_motion_to_semimajor_axis(mean_motion: float) -> float:
    """
    Fast computation of semimajor axis from Kepler's third law, assuming a
    circular orbit around a spherical Earth. This is the inverse of
    `semimajor_axis_to_mean_motion`.

    Args:
        mean_motion (float): Orbit mean motion (degrees/second).

    Returns:
        float: The semimajor axis (meters).
    """
    return np.cbrt(constants.EARTH_MU / (np.radians(mean_motion) ** 2))


@njit
def semimajor_axis_to_orbit_period(semimajor_axis: float) -> float:
    """
    Fast computation of orbital period from Kepler's third law, assuming a
    circular orbit around a spherical Earth.

    Args:
        semimajor_axis (float): Orbit semimajor axis (meters).

    Returns:
        float: Orbital period (seconds).
    """
    return mean_motion_to_orbit_period(semimajor_axis_to_mean_motion(semimajor_axis))


@njit
def mean_anomaly_to_true_anomaly(mean_anomaly: float, eccentricity: float = 0) -> float:
    """
    Approximates orbit true anomaly using a third-order series expansion.

    Args:
        mean_anomaly (float): Orbit mean anomaly (degrees).
        eccentricity (float): Orbit eccentricity.

    Returns:
        float: Orbit true anomaly (degrees).
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
    Approximates orbit mean anomaly using a third-order series expansion.

    Args:
        true_anomaly (float): Orbit true anomaly (degrees).
        eccentricity (float): Orbit eccentricity.

    Returns:
        float: Orbit mean anomaly (degrees).
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
def compute_j2_raan_rate(
    semimajor_axis: float, inclination: float, eccentricity: float
) -> float:
    """
    Fast computation of the secular precession rate of the right ascension
    of the ascending node (nodal regression) due to Earth's J2 oblateness
    perturbation. The rate is negative (westward regression) for prograde
    orbits (inclination < 90 degrees), zero for polar orbits (cos(90) = 0),
    and positive for retrograde orbits. Sun-synchronous orbits are defined
    by choosing an inclination that makes this rate equal to the Earth's
    mean motion around the Sun (about 360/365.2422 degrees/day).

    Args:
        semimajor_axis (float): Orbit semimajor axis (meters).
        inclination (float): Orbit inclination (degrees).
        eccentricity (float): Orbit eccentricity.

    Returns:
        float: The right ascension of ascending node precession rate (degrees/second).
    """
    return (
        -3
        / 2
        * constants.EARTH_J2
        * semimajor_axis_to_mean_motion(semimajor_axis)
        * (constants.EARTH_MEAN_RADIUS / (semimajor_axis * (1 - eccentricity**2))) ** 2
        * np.cos(np.radians(inclination))
    )


@njit
def compute_j2_aop_rate(
    semimajor_axis: float, inclination: float, eccentricity: float
) -> float:
    """
    Fast computation of the secular precession rate of the argument of
    periapsis due to Earth's J2 oblateness perturbation. The rate is
    positive below the critical inclination (~63.43 degrees, where
    5*cos^2(inclination) - 1 = 0), negative above it, and exactly zero at
    the critical inclination itself. Molniya-type orbits use this
    critical inclination specifically so their periapsis (and apoapsis)
    location does not drift over time.

    Args:
        semimajor_axis (float): Orbit semimajor axis (meters).
        inclination (float): Orbit inclination (degrees).
        eccentricity (float): Orbit eccentricity.

    Returns:
        float: The argument of periapsis rate (degrees/second).
    """
    return (
        3
        / 4
        * constants.EARTH_J2
        * semimajor_axis_to_mean_motion(semimajor_axis)
        * (constants.EARTH_MEAN_RADIUS / (semimajor_axis * (1 - eccentricity**2))) ** 2
        * (5 * np.cos(np.radians(inclination)) ** 2 - 1)
    )
