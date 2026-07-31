"""
Observation utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
from numba import njit

from .. import constants
from .orbital import compute_ground_surface_velocity, semimajor_axis_to_mean_motion


@njit
def swath_width_to_field_of_regard(
    altitude: float, swath_width: float, elevation: float = 0
) -> float:
    """
    Fast computation of the field of regard (the total angular width an
    instrument must be able to point or scan across) required to observe a
    specified ground swath width, assuming a spherical Earth (using the
    mean Earth radius) and a circular orbit at constant altitude.

    Args:
        altitude (float): Altitude (meters) above the mean Earth radius for the
            observing instrument.
        swath_width (float): Ground swath width (meters): the cross-track distance,
            measured along the Earth's surface, at the specified elevation.
        elevation (float): Elevation (meters) above the mean Earth radius of the observed swath.

    Returns:
        float: The field of regard (degrees): the full angular width, centered on
        nadir, that the instrument must be able to point across to observe the swath.
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
    Fast computation of the field of view (the angular extent, as seen
    from the satellite, spanning the near to far edge of a swath) required
    to observe a specified ground swath width centered on a given off-nadir
    look angle, assuming a spherical Earth (using the mean Earth radius)
    and a circular orbit at constant altitude. Unlike
    `swath_width_to_field_of_regard`, the swath is not centered on nadir,
    so the field of view spans asymmetrically between the near and far
    edges of the swath.

    Args:
        altitude (float): Altitude (meters) above the mean Earth radius for the
            observing instrument.
        swath_width (float): Ground swath width (meters): the cross-track distance,
            measured along the Earth's surface, centered on the look angle direction.
        look_angle (float): Off-nadir look angle (degrees), measured at the satellite,
            to the center of the swath. Saturates at the horizon-limited maximum.
        elevation (float): Elevation (meters) above the mean Earth radius of the observed swath.

    Returns:
        float: The field of view (degrees): the angular extent, as seen from the
        satellite, spanning the near to far edge of the swath.
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = (constants.EARTH_MEAN_RADIUS + elevation) / (
        constants.EARTH_MEAN_RADIUS + altitude
    )
    # eta is the angular radius from sub-satellite point to center of view
    sin_eta = min(sin_rho, np.sin(np.radians(look_angle)))
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
    Fast computation of the ground swath width observable for a specified
    field of regard, assuming a spherical Earth (using the mean Earth
    radius) and a circular orbit at constant altitude. This is the inverse
    of `swath_width_to_field_of_regard`. A `field_of_regard` at or beyond
    the horizon-limited maximum (i.e. pointing to the horizon) saturates
    to the maximum observable swath width rather than growing without bound.

    Args:
        altitude (float): Altitude (meters) above the mean Earth radius for the
            observing instrument.
        field_of_regard (float): The full angular width (degrees), centered on
            nadir, that the instrument points or scans across.
        elevation (float): Elevation (meters) above the mean Earth radius of the observed swath.

    Returns:
        float: The ground swath width (meters): the cross-track distance,
        measured along the Earth's surface, at the specified elevation.
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
def compute_max_access_time(mean_altitude: float, min_elevation_angle: float) -> float:
    """
    Fast computation of maximum access time to observe a point.

    Args:
        mean_altitude (float): Orbit mean altitude (meters).
        min_elevation_angle (float): Minimum elevation angle (degrees) for observation.

    Returns:
        float: The maximum access time (seconds) for observation.
    """
    # angular distance from sub-satellite point to edge of viewable region
    earth_angle = np.degrees(
        np.arccos(
            constants.EARTH_MEAN_RADIUS
            / (constants.EARTH_MEAN_RADIUS + mean_altitude)
            * np.cos(np.radians(min_elevation_angle))
        )
        - np.radians(min_elevation_angle)
    )
    # max access time is twice the earth central angle divided by the mean motion of the orbit
    return (
        2
        * earth_angle
        / semimajor_axis_to_mean_motion(constants.EARTH_MEAN_RADIUS + mean_altitude)
    )


@njit
def compute_max_transit_time(
    mean_altitude: float, inclination: float, along_track: float
) -> float:
    """
    Fast computation of the maximum (conservative, worst-case) transit time
    to cover a specified along-track distance, using the slowest ground
    track velocity attained over the orbit, which occurs at the orbit's
    extreme latitude (min(inclination, 180 - inclination)), not the equator.

    Args:
        mean_altitude (float): The mean orbit altitude (meters) above WGS 84 datum.
        inclination (float): The orbit inclination (degrees).
        along_track (float): The along track distance (meters) observed during access.

    Returns:
        float: The maximum access time (seconds) to traverse the along track distance.
    """
    # slowest velocity occurs at the orbit's extreme latitude
    extreme_latitude = min(inclination, 180 - inclination)
    v_slowest = compute_ground_surface_velocity(
        mean_altitude, 0, inclination, extreme_latitude
    )
    return along_track / v_slowest


@njit
def compute_min_along_track_distance(
    mean_altitude: float, inclination: float, access_time: float
) -> float:
    """
    Fast computation of the minimum (conservative, worst-case) along-track
    distance observed in a specified access time, using the slowest ground
    track velocity attained over the orbit, which occurs at the orbit's
    extreme latitude (min(inclination, 180 - inclination)), not the
    equator. This is the inverse of `compute_max_transit_time`.

    Args:
        mean_altitude (float): The mean orbit altitude (meters) above WGS 84 datum.
        inclination (float): The orbit inclination (degrees).
        access_time (float): The access time (seconds) during observation.

    Returns:
        float: The minimum along track distance (meters) observed during the access time.
    """
    # slowest velocity occurs at the orbit's extreme latitude
    extreme_latitude = min(inclination, 180 - inclination)
    v_slowest = compute_ground_surface_velocity(
        mean_altitude, 0, inclination, extreme_latitude
    )
    return access_time * v_slowest
