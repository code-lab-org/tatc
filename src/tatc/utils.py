# -*- coding: utf-8 -*-
"""
Utility functions.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

import numpy as np
from numba import njit
import geopandas as gpd
from typing import Union
from shapely.geometry import Polygon, MultiPolygon, GeometryCollection
from shapely.ops import clip_by_rect

from . import constants


@staticmethod
@njit
def mean_anomaly_to_true_anomaly(mean_anomaly, eccentricity=0):
    """
    Converts mean anomaly to true anomaly.

    Args:
        mean_anomaly (float): The mean anomaly (degrees).
        true_anomaly (float): The orbit eccentricity.

    Returns:
        float: The true anomaly (degrees).
    """
    M = np.radians(mean_anomaly)
    e = eccentricity
    nu = (
        M
        + (2 * e - (1 / 4) * e**3) * np.sin(M)
        + (5 / 4) * e**2 * np.sin(2 * M)
        + (13 / 12) * e**3 * np.sin(3 * M)
    )
    return np.degrees(nu)


@staticmethod
@njit
def true_anomaly_to_mean_anomaly(true_anomaly, eccentricity=0):
    """
    Converts true anomaly to mean anomaly.

    Args:
        true_anomaly (float): The true anomaly (degrees).
        eccentricity (float): The orbit eccentricity.

    Returns:
        float: The mean anomaly (degrees).
    """
    e = eccentricity
    nu = np.radians(true_anomaly)
    M = (
        nu
        - 2 * e * np.sin(nu)
        + ((3 / 4) * e**2 + (1 / 8) * e**4) * np.sin(2 * nu)
        - (1 / 3) * e**3 * np.sin(3 * nu)
        + (5 / 32) * e**4 * np.sin(4 * nu)
    )
    return np.degrees(M)


@staticmethod
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
    # compute the angular distance of each sample (assuming mean sphere)
    theta = distance / constants.earth_mean_radius
    # compute the distance from the center of earth to conic plane (assuming sphere)
    r = constants.earth_mean_radius * np.cos(theta / 2)
    # compute the distance from the conic plane to the surface (assuming sphere)
    h = constants.earth_mean_radius - r
    # compute the sperical cap area covered by the sample (assuming sphere)
    # https://en.wikipedia.org/wiki/Spherical_cap
    sample_area = 2 * np.pi * constants.earth_mean_radius * h
    # return the fraction of earth-to-sample area
    return int(constants.earth_surface_area / sample_area)


@staticmethod
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


@staticmethod
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


@staticmethod
@njit
def compute_field_of_regard(height, min_elevation_angle):
    """
    Fast computation of field of regard for observation with a minimum altitude angle.

    Args:
        height (float): Height (meters) above surface of the observation.
        min_elevation_angle (float): The minimum elevation angle (degrees) for observation.

    Returns:
        float: Angular width (degrees) of observation.
    """
    # rho is the angular radius of the earth viewed by the satellite
    sin_rho = constants.earth_mean_radius / (constants.earth_mean_radius + height)
    # epsilon is the min satellite elevation for obs (grazing angle)
    cos_epsilon = np.cos(np.radians(min_elevation_angle))
    # eta is the angular radius of the region viewable by the satellite
    sin_eta = sin_rho * cos_epsilon
    return np.degrees(np.arcsin(sin_eta) * 2)


@staticmethod
@njit
def compute_min_elevation_angle(height, field_of_regard):
    """
    Fast computation of minimum elevation angle required to observe a point.

    Args:
        height (float): Height (meters) above surface of the observation.
        field_of_regard (float): Angular width (degrees) of observation.

    Returns:
        float: The minimum elevation angle (degrees) for observation.
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


@staticmethod
@njit
def compute_orbit_period(height):
    """
    Fast computation of approximate orbital period.

    Args:
        height (float): Height (meters) above surface of the observation.

    Returns:
        float: The orbital period (seconds).
    """
    semimajor_axis = constants.earth_mean_radius + height
    mean_motion_rad_s = np.sqrt(constants.earth_mu / semimajor_axis**3)
    return 2 * np.pi / mean_motion_rad_s


@staticmethod
@njit
def compute_max_access_time(height, min_elevation_angle):
    """
    Fast computation of maximum access time to observe a point.

    Args:
        height (float): Height (meters) above surface of the observation.
        min_elevation_angle (float): Minimum elevation angle (degrees) for observation.

    Returns:
        float: The maximum access time (seconds) for observation.
    """
    orbital_distance = height * (np.pi - 2 * np.radians(min_elevation_angle))
    orbital_velocity = np.sqrt(
        constants.earth_mu / (constants.earth_mean_radius + height)
    )
    return orbital_distance / orbital_velocity


@staticmethod
def _wrap_polygon_over_north_pole(polygon: Polygon) -> Polygon:
    """
    Wraps polygon coordinates over the North pole. Due to buffering and projection,
    sometimes latitudes exceed 90 degrees. This method wraps them to the correct
    latitude between -90 and 90 degrees and adjusts the longitude by 180 degrees.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (Polygon): The polygon to wrap.

    Returns:
       Polygon: The wrapped polygon.
    """
    # map latitudes from [90, 180) to [90, -90), adjusting longitude by 180 degrees
    polygon = Polygon(
        [
            [
                (c[0] + 180 if c[0] <= 0 else c[0] - 180) if c[1] >= 90 else c[0],
                180 - c[1] if c[1] >= 90 else c[1],
            ]
            for c in polygon.exterior.coords
        ],
        [
            [
                [
                    (c[0] + 180 if c[0] <= 0 else c[0] - 180) if c[1] >= 90 else c[0],
                    180 - c[1] if c[1] >= 90 else c[1],
                ]
                for c in i.coords
            ]
            for i in polygon.interiors
        ],
    )
    # attempt to fix invalid polygon using a zero buffer
    if not polygon.is_valid:
        polygon = polygon.buffer(0)
    return polygon


@staticmethod
def _split_polygon_north_pole(
    polygon: Union[Polygon, MultiPolygon]
) -> Union[Polygon, MultiPolygon]:
    """
    Splits a Polygon into a MultiPolygon if it crosses north pole.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to split.

    Returns:
       Polygon, or MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        lat = np.array([c[1] for c in polygon.exterior.coords])
        if np.any(lat > 90):
            # find the "bottom" portion of the polygon
            bottom = clip_by_rect(polygon, -180, -180, 180, 90)
            # check for dirty clip and fix if necessary
            if isinstance(bottom, GeometryCollection):
                bottom = _convert_collection_to_polygon(top)
            # find the "top" portion of the polygon
            top = clip_by_rect(polygon, -180, 90, 180, 180)
            # check for dirty clip and fix if necessary
            if isinstance(top, GeometryCollection):
                top = _convert_collection_to_polygon(top)
            # wrap polygon over the north pole
            if isinstance(top, Polygon):
                top = _wrap_polygon_over_north_pole(top)
            elif isinstance(top, MultiPolygon):
                top = MultiPolygon(
                    [_wrap_polygon_over_north_pole(p) for p in top.geoms]
                )
            else:
                raise ValueError("Geometry not implemented: " + str(type(top)))
            # return the composite multi polygon
            return MultiPolygon(
                [bottom, top]
                if (isinstance(bottom, Polygon) and isinstance(top, Polygon))
                else [bottom] + list(top.geoms)
                if isinstance(bottom, Polygon)
                else list(bottom.geoms) + [top]
                if isinstance(top, Polygon)
                else list(bottom.geoms) + list(top.geoms)
            )
        else:
            return polygon
    elif isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        polygons = [_split_polygon_north_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in polygons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    else:
        raise ValueError("Unknown geometry: " + str(type(polygon)))


@staticmethod
def _wrap_polygon_over_south_pole(polygon: Polygon) -> Polygon:
    """
    Wraps polygon coordinates over the South pole. Due to buffering and projection,
    sometimes latitudes exceed -90 degrees. This method wraps them to the correct
    latitude between -90 and 90 degrees and adjusts the longitude by 180 degrees.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (Polygon): The polygon to wrap.

    Returns:
       Polygon: The wrapped polygon.
    """
    # map latitudes from [-90, -180) to [-90, 90), adjusting longitude by 180 degrees
    polygon = Polygon(
        [
            [
                (c[0] + 180 if c[0] <= 0 else c[0] - 180) if c[1] <= -90 else c[0],
                -180 - c[1] if c[1] <= -90 else c[1],
            ]
            for c in polygon.exterior.coords
        ],
        [
            [
                [
                    (c[0] + 180 if c[0] <= 0 else c[0] - 180) if c[1] <= -90 else c[0],
                    -180 - c[1] if c[1] <= -90 else c[1],
                ]
                for c in i.coords
            ]
            for i in polygon.interiors
        ],
    )
    # attempt to fix invalid polygon using a zero buffer
    if not polygon.is_valid:
        polygon = polygon.buffer(0)
    return polygon


@staticmethod
def _split_polygon_south_pole(
    polygon: Union[Polygon, MultiPolygon]
) -> Union[Polygon, MultiPolygon]:
    """
    Splits a Polygon into a MultiPolygon if it crosses south pole.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to split.

    Returns:
       Polygon, or MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        lat = np.array([c[1] for c in polygon.exterior.coords])
        if np.any(lat < -90):
            # find the "top" portion of the polygon
            top = clip_by_rect(polygon, -180, -90, 180, 180)
            # check for dirty clip and fix if necessary
            if isinstance(top, GeometryCollection):
                top = _convert_collection_to_polygon(top)
            # find the "bottom" portion of the polygon
            bottom = clip_by_rect(polygon, -180, -180, 180, -90)
            # check for dirty clip and fix if necessary
            if isinstance(bottom, GeometryCollection):
                bottom = _convert_collection_to_polygon(bottom)
            # wrap polygon over the south pole
            if isinstance(bottom, Polygon):
                bottom = _wrap_polygon_over_south_pole(bottom)
            elif isinstance(bottom, MultiPolygon):
                bottom = MultiPolygon(
                    [_wrap_polygon_over_south_pole(p) for p in bottom.geoms]
                )
            else:
                raise ValueError("Geometry not implemented: " + str(type(bottom)))
            # return the composite multi polygon
            return MultiPolygon(
                [top, bottom]
                if (isinstance(top, Polygon) and isinstance(bottom, Polygon))
                else [top] + list(bottom.geoms)
                if isinstance(top, Polygon)
                else list(top.geoms) + [bottom]
                if isinstance(bottom, Polygon)
                else list(top.geoms) + list(bottom.geoms)
            )
        else:
            return polygon
    elif isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        polygons = [_split_polygon_south_pole(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in polygons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    else:
        raise ValueError("Unknown geometry: " + str(type(polygon)))


@staticmethod
def _wrap_polygon_over_antimeridian(polygon: Polygon) -> Polygon:
    """
    Wraps polygon coordinates over the antimeridian. Due to buffering and projection,
    sometimes longitudes exceed 180 degrees. This method wraps them to
    the correct longitude between -180 and 180 degrees.

    Note: this method only changes coordinates: it does not create a MultiPolygon.

    Args:
       polygon (Polygon): The polygon to wrap.

    Returns:
       Polygon: The wrapped polygon.
    """
    # map longitudes from [180, 360) to [-180, 0)
    polygon = Polygon(
        [[c[0] - 360 if c[0] >= 180 else c[0], c[1]] for c in polygon.exterior.coords],
        [
            [[c[0] - 360 if c[0] >= 180 else c[0], c[1]] for c in i.coords]
            for i in polygon.interiors
        ],
    )
    # attempt to fix invalid polygon using a zero buffer
    if not polygon.is_valid:
        polygon = polygon.buffer(0)
    return polygon


@staticmethod
def _convert_collection_to_polygon(
    collection: GeometryCollection,
) -> Union[Polygon, MultiPolygon]:
    """
    Converts a GeometryCollection to a Polygon or MultiPolygon. Quick clipping
    can create dirty results with points or lines on boundaries. This method
    drops and lines or points from a GeometryCollection to return only the
    Polygon or MultiPolygon geometry.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to convert.

    Returns:
       Polygon, or MultiPolygon: The converted polygon.
    """
    pgons = [p for p in collection.geoms if isinstance(p, Polygon)] + [
        p.geoms
        for g in collection.geoms
        for p in g.geoms
        if isinstance(p, MultiPolygon)
    ]
    if len(pgons) == 1:
        return pgons[0]
    else:
        return MultiPolygon(pgons)


@staticmethod
def _split_polygon_antimeridian(
    polygon: Union[Polygon, MultiPolygon]
) -> Union[Polygon, MultiPolygon]:
    """
    Splits a Polygon into a MultiPolygon if it crosses the anti-meridian after
    wrapping its coordinates using `wrap_coordinates_antimeridian`.

    Args:
       polygon (Polygon or MultiPolygon): The polygon to split.

    Returns:
       Polygon, or MultiPolygon: The split polygon.
    """
    if isinstance(polygon, Polygon):
        lon = np.array([c[0] for c in polygon.exterior.coords])
        # check if any longitudes cross the anti-meridian
        if np.any(np.abs(np.diff(lon)) > 100):
            # map longitudes from [-180, 45) to [225, 405)
            polygon = Polygon(
                [
                    (c[0] + 360 if c[0] < 45 else c[0], c[1])
                    for c in polygon.exterior.coords
                ],
                [
                    [(c[0] + 360 if c[0] < 45 else c[0], c[1]) for c in i.coords]
                    for i in polygon.interiors
                ],
            )
            # attempt to fix invalid polygon using a zero buffer
            if not polygon.is_valid:
                polygon = polygon.buffer(0)
            # find the "left" portion of the polygon
            left = clip_by_rect(polygon, -180, -180, 180, 180)
            # check for dirty clip and fix if necessary
            if isinstance(left, GeometryCollection):
                left = _convert_collection_to_polygon(left)
            # find the "right" portion of the polygon
            right = clip_by_rect(polygon, 180, -180, 405, 180)
            # check for dirty clip and fix if necessary
            if isinstance(right, GeometryCollection):
                right = _convert_collection_to_polygon(right)
            # wrap polygon over the anti-meridian
            if isinstance(right, Polygon):
                right = _wrap_polygon_over_antimeridian(right)
            elif isinstance(right, MultiPolygon):
                right = MultiPolygon(
                    [_wrap_polygon_over_antimeridian(p) for p in right.geoms]
                )
            else:
                raise ValueError("Geometry not implemented: " + str(type(right)))
            # return the composite multi polygon
            return MultiPolygon(
                [left, right]
                if (isinstance(left, Polygon) and isinstance(right, Polygon))
                else [left] + list(right.geoms)
                if isinstance(left, Polygon)
                else list(left.geoms) + [right]
                if isinstance(right, Polygon)
                else list(left.geoms) + list(right.geoms)
            )
        else:
            return polygon
    elif isinstance(polygon, MultiPolygon):
        # recursive call for each polygon
        polygons = [_split_polygon_antimeridian(p) for p in polygon.geoms]
        return MultiPolygon(
            [
                g
                for p in polygons
                for g in (p.geoms if isinstance(p, MultiPolygon) else [p])
            ]
        )
    else:
        raise ValueError("Unknown geometry: " + str(type(polygon)))


@staticmethod
def split_polygon(
    polygon: Union[Polygon, MultiPolygon]
) -> Union[Polygon, MultiPolygon]:
    """
    Splits a Polygon into a MultiPolygon if it crosses the anti-meridian
    (180 degrees longitude), exceeds the north pole (90 degrees latitude), or
    exceeds the south pole (-90 degrees latitude).

    Args:
        polygon (Polygon or MultiPolygon): The polygon to split.

    Returns:
        Polygon, or MultiPolygon: The split polygon.
    """
    polygon = _split_polygon_north_pole(
        _split_polygon_south_pole(_split_polygon_antimeridian(polygon))
    )
    # attempt to fix invalid polygon using a zero buffer
    if not polygon.is_valid:
        polygon = polygon.buffer(0)
    return polygon


@staticmethod
def normalize_geometry(
    geometry: Union[Polygon, MultiPolygon, gpd.GeoDataFrame]
) -> gpd.GeoDataFrame:
    """
    Normalize geometry to a GeoDataFrame with antimeridian wrapping.

    Args:
        geometry (geopandas.GeoDataFrame, geopandas.GeoSeries, Polygon, or MultiPolygon): The geometry to normalize.

    Returns:
        geopandas.GeoDataFrame: The normalized geometry.
    """
    if isinstance(geometry, Polygon) or isinstance(geometry, MultiPolygon):
        if not geometry.is_valid:
            raise ValueError("Geometry is not a valid Polygon or MultiPolygon.")
        geometry = gpd.GeoDataFrame(geometry=gpd.GeoSeries([geometry]), crs="EPSG:4326")
    elif isinstance(geometry, gpd.GeoSeries):
        geometry = gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326")
    if isinstance(geometry, gpd.GeoDataFrame):
        geometry["geometry"] = geometry.apply(
            lambda r: split_polygon(r.geometry),
            axis=1,
        )
    return geometry
