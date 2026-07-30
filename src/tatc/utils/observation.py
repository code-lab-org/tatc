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
    Fast computation of maximum transit time to cover a specified along track distance.

    Args:
        mean_altitude (float): The mean orbit altitude (meters) above WGS 84 datum.
        inclination (float): The orbit inclination (degrees).
        along_track (float): The along track distance (meters) observed during access.

    Returns:
        float: The access time (seconds) for observation.
    """
    # slowest velocity occurs at the equator due to Earth rotation
    v_slowest = compute_ground_surface_velocity(mean_altitude, 0, inclination, 0)
    return along_track / v_slowest


@njit
def compute_min_along_track_distance(
    mean_altitude: float, inclination: float, access_time: float
) -> float:
    """
    Fast computation of minimum along track distance for a specified access time.

    Args:
        mean_altitude (float): The mean orbit altitude (meters) above WGS 84 datum.
        inclination (float): The orbit inclination (degrees).
        access_time (float): The access time (seconds) during observation.

    Returns:
        float: The observation along track distance (meters).
    """
    # slowest velocity occurs at the equator due to Earth rotation
    v_slowest = compute_ground_surface_velocity(mean_altitude, 0, inclination, 0)
    return access_time * v_slowest
