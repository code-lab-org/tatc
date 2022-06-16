# -*- coding: utf-8 -*-
"""
Utility functions.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

import numpy as np
from numba import njit
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon

from . import constants


@njit
def compute_number_samples(distance):
    """
    Compute the number of global samples required to achieve a typical
    sample distance (meters) assuming equal spacing.

    Args:
        distance (float): The typical distance between samples (meters).

    Returns:
        int: The number of global samples.
    """
    # compute the angular disance of each sample (assuming mean sphere)
    theta = distance / constants.earth_mean_radius
    # compute the distance from the center of earth to conic plane (assuming sphere)
    r = constants.earth_mean_radius * np.cos(theta / 2)
    # compute the distance from the conic plane to the surface (assuming sphere)
    h = constants.earth_mean_radius - r
    # compute the conic section area covered by the sample (assuming sphere)
    sample_area = 2 * np.pi * constants.earth_mean_radius * h
    # return the fraction of earth-to-sample area
    return int(constants.earth_surface_area / sample_area)


@njit
def swath_width_to_field_of_regard(height, swath_width):
    """
    Fast conversion from swath width to field of regard.

    Args:
        height (float): Height (meters) above surface of the observation.
        swath_width (float): Observation diameter (meters).

    Returns:
        float: The field of regard (degrees).
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = constants.earth_mean_radius / (constants.earth_mean_radius + height)
    # lambda is the Earth central angle
    sin_lambda = np.sin((swath_width / 2) / constants.earth_mean_radius)
    # eta is the angular radius of the region viewable by the satellite
    tan_eta = sin_rho * sin_lambda / (1 - sin_rho * np.cos(np.arcsin(sin_lambda)))
    return np.degrees(2 * np.arctan(tan_eta))


@njit
def field_of_regard_to_swath_width(height, field_of_regard):
    """
    Fast conversion from field of regard to swath width.

    Args:
        height (float): Height (meters) above surface of the observation.
        field_of_regard (float): Angular width (degrees) of observation.

    Returns:
        float: The observation diameter (meters).
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = constants.earth_mean_radius / (constants.earth_mean_radius + height)
    # eta is the angular radius of the region viewable by the satellite
    sin_eta = min(sin_rho, np.sin(np.radians(field_of_regard) / 2))
    # epsilon is the min satellite elevation for obs (grazing angle)
    cos_epsilon = sin_eta / sin_rho
    # lambda is the Earth central angle
    _lambda = np.pi / 2 - np.arcsin(sin_eta) - np.arccos(cos_epsilon)
    return 2 * constants.earth_mean_radius * _lambda


@njit
def compute_field_of_regard(height, min_altitude):
    """
    Fast computation of field of regard for observation with a minimum altitude angle.

    Args:
        height (float): Height (meters) above surface of the observation.
        min_altitude (float): The minimum altitude angle (degrees) for observation.

    Returns:
        float: Angular width (degrees) of observation.
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = constants.earth_mean_radius / (constants.earth_mean_radius + height)
    # epsilon is the min satellite elevation for obs (grazing angle)
    cos_epsilon = np.cos(np.radians(min_altitude))
    # eta is the angular radius of the region viewable by the satellite
    sin_eta = sin_rho * cos_epsilon
    return np.degrees(np.arcsin(sin_eta) * 2)


@njit
def compute_min_altitude(height, field_of_regard):
    """
    Fast computation of minimum altitude angle required to observe a point.

    Args:
        height (float): Height (meters) above surface of the observation.
        field_of_regard (float): Angular width (degrees) of observation.

    Returns:
        float: The minimum altitude angle (degrees) for observation.
    """
    # eta is the angular radius of the region viewable by the satellite
    sin_eta = np.sin(np.radians(field_of_regard) / 2)
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = constants.earth_mean_radius / (constants.earth_mean_radius + height)
    # epsilon is the min satellite elevation for obs (grazing angle)
    cos_epsilon = sin_eta / sin_rho
    if cos_epsilon > 1:
        return 0
    return np.degrees(np.arccos(cos_epsilon))


@njit
def compute_orbit_period(height):
    """
    Fast computation of approximate orbital period.

    Args:
        height (float): Height (meters) above surface of the observation.

    Returns:
        float: The orbital period (seconds).
    """
    return (
        2
        * np.pi
        * np.sqrt(
            np.power(constants.earth_mean_radius + height, 3) / constants.earth_mu
        )
    )


@njit
def compute_max_access_time(height, min_altitude):
    """
    Fast computation of maximum access time to observe a point.

    Args:
        height (float): Height (meters) above surface of the observation.
        min_altitude (float): Minimum altitude angle (degrees) for observation.

    Returns:
        float: The maximum access time (seconds) for observation.
    """
    orbital_distance = height * (np.pi - 2 * np.radians(min_altitude))
    orbital_velocity = np.sqrt(
        constants.earth_mu / (constants.earth_mean_radius + height)
    )
    return orbital_distance / orbital_velocity


def wrap_coordinates_antimeridian(coords):
    """
    Wraps negative longitudes around the antimeridian to positive values.

    Args:
        coords (list(tuple(float))): the raw coordinates (longitude, latitude).

    Returns:
        list(tuple(float)): The updated coordinates (longitude, latitude).
    """
    lon = np.array([c[0] for c in coords])
    lat = np.array([c[1] for c in coords])
    # wrap negative coords crossing antimeridian to positive values
    if np.any(lon < 0) and np.any(lon > 0) and 180 < np.ptp(lon) < 360:
        lon[lon < 0] += 360
    return list(zip(lon, lat))


def normalize_geometry(geometry):
    """
    Normalize geometry to a GeoDataFrame with antimeridian wrapping.

    Args:
        geometry (:obj:`GeoDataFrame`, :obj:`GeoSeries`, :obj:`Polygon`,
              or :obj:`MultiPolygon`, optional): The geometry to normalize.

    Returns:
        :obj:`GeoDataFrame`: The normalized geometry.
    """
    if isinstance(geometry, Polygon) or isinstance(geometry, MultiPolygon):
        geometry = gpd.GeoDataFrame(geometry=gpd.GeoSeries([geometry]), crs="EPSG:4326")
    elif isinstance(geometry, gpd.GeoSeries):
        geometry = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")
    if isinstance(geometry, gpd.GeoDataFrame):
        geometry["geometry"] = geometry.apply(
            lambda r: Polygon(wrap_coordinates_antimeridian(r.geometry.exterior.coords))
            if isinstance(r.geometry, Polygon)
            else r.geometry,
            axis=1,
        )
    return geometry
