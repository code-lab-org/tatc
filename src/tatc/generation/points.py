"""
Methods to generate geospatial points to sample data.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon

from ..constants import EARTH_MEAN_RADIUS
from ..utils.geometry import _hash_geometries
from ..utils.surface import compute_number_samples
from ._grid import clip_to_mask, generate_indices_uniform_spacing


def generate_points_fibonacci_lattice(
    distance: float,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic points following a Fibonacci lattice.

    See: Gonzalez (2010). "Measurement of areas on a sphere using Fibonacci
    and latitude-longitude lattices", Mathematical Geosciences 42(49).
    doi: 10.1007/s11004-009-9257-x

    Note: this implementation differs slightly from Gonzalez. Gonzalez
    requires an odd number of points. This implementation allows any number
    of points. In agreement with Gonzalez, no points are placed at poles.

    Args:
        distance (float): The typical surface distance (meters) between points.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain points
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.

    Returns:
        geopandas.GeoDataFrame: the data frame of generated points, each
            identified by the hash of its geometry (`point_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """

    # determine the number of global samples to achieve average sample distance
    samples = compute_number_samples(distance)
    index = np.arange(samples)
    # compute latitudes, starting from the southern hemisphere and placing
    # neither first nor last points at poles
    latitudes = np.degrees(np.arcsin(2 * (index + 1) / (samples + 2) - 1))
    # compute longitudes on the interval [-180, 180]
    phi = (1 + np.sqrt(5)) / 2  # golden ratio
    longitudes = np.mod(360 * index / phi, 360)
    longitudes = np.where(longitudes > 180, longitudes - 360, longitudes)
    if mask is not None:
        if isinstance(mask, (Polygon, MultiPolygon)):
            if not mask.is_valid:
                raise ValueError("Mask is not a valid Polygon or MultiPolygon.")
            total_bounds = mask.bounds
        else:
            total_bounds = [-180, -90, 180, 90]
        # use the total_bounds to filter relevant points
        min_longitude = total_bounds[0]
        min_latitude = total_bounds[1]
        max_longitude = 180 if total_bounds[2] == -180 else total_bounds[2]
        max_latitude = total_bounds[3]
        # shift longitudes to the interval [0, 360] for masks beyond 180
        if max_longitude > 180:
            longitudes = np.where(longitudes < 0, longitudes + 360, longitudes)
        in_bounds = (
            (min_latitude <= latitudes)
            & (latitudes <= max_latitude)
            & (min_longitude <= longitudes)
            & (longitudes <= max_longitude)
        )
        longitudes = longitudes[in_bounds]
        latitudes = latitudes[in_bounds]
    # create a geodataframe in the WGS84 coordinate reference system (EPSG:4326)
    gdf = gpd.GeoDataFrame(
        geometry=shapely.points(
            longitudes, latitudes, np.full_like(longitudes, elevation)
        ),
        crs="EPSG:4326",
    )
    # clip the geodataframe to the supplied mask, if required
    if mask is not None:
        gdf = clip_to_mask(gdf, mask)
    # identify each point by the hash of its geometry
    gdf.insert(0, "point_id", _hash_geometries(gdf.geometry))
    # return the final geodataframe
    return gdf


def generate_points_uniform_spacing(
    distance: float,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic points at the centroid of regular equally spaced grid
    cells.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    Args:
        distance (float):  The typical surface distance (meters) between points.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain points
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.

    Returns:
        geopandas.GeoDataFrame: the data frame of generated points, each
            identified by the hash of its geometry (`point_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """
    # compute the angular disance of each sample (assuming sphere)
    theta_longitude = np.degrees(distance / EARTH_MEAN_RADIUS)
    theta_latitude = np.degrees(distance / EARTH_MEAN_RADIUS)
    return generate_points_uniform_angular_distance(
        theta_longitude, theta_latitude, elevation, mask
    )


def generate_points_uniform_angular_distance(
    theta_longitude: float,
    theta_latitude: float,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic cells following regular equally spaced grid.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    Args:
        theta_longitude (float): The angular difference in longitude (degrees)
            between points.
        theta_latitude (float): The angular difference in latitude (degrees)
            between points.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain points
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.

    Returns:
        geopandas.GeoDataFrame: the data frame of generated points, each
            identified by the hash of its geometry (`point_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """

    # generate grid cells over the filtered region
    indices = generate_indices_uniform_spacing(
        theta_longitude,
        theta_latitude,
        mask,
    )
    longitudes = -180 + (indices[:, 0] + 0.5) * theta_longitude
    latitudes = -90 + (indices[:, 1] + 0.5) * theta_latitude
    # create a geodataframe in the WGS84 reference frame
    gdf = gpd.GeoDataFrame(
        geometry=shapely.points(
            longitudes, latitudes, np.full_like(longitudes, elevation)
        ),
        crs="EPSG:4326",
    )
    # clip the geodataframe to the supplied mask, if required
    if mask is not None:
        gdf = clip_to_mask(gdf, mask)
    # identify each point by the hash of its geometry
    gdf.insert(0, "point_id", _hash_geometries(gdf.geometry))
    # return the final geodataframe
    return gdf
