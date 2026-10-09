"""
Methods to generate geospatial cells to aggregate data.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon

from ..constants import EARTH_MEAN_RADIUS
from ..utils.geometry import _hash_geometries, get_planar_bounds
from ._grid import clip_to_mask, generate_indices_uniform_spacing


def generate_cells_uniform_spacing(
    distance: float,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
    strips: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic polygons over a regular equally spaced grid.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    Args:
        distance (float):  The typical surface distance (meters) between points.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain cells
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        strips (str | None): Option to generate strip-cells along latitude (`"lat"`),
            longitude (`"lon"`), or none (`None`).

    Returns:
        geopandas.GeoDataFrame: the data frame of generated cells, each
            identified by the hash of its geometry (`cell_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """
    # compute the angular disance of each sample (assuming sphere)
    theta_longitude = np.degrees(distance / EARTH_MEAN_RADIUS)
    theta_latitude = np.degrees(distance / EARTH_MEAN_RADIUS)
    return generate_cells_uniform_angular_spacing(
        theta_longitude, theta_latitude, elevation, mask, strips
    )


def generate_cells_uniform_angular_spacing(
    theta_longitude: float,
    theta_latitude: float,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
    strips: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic polygons over a regular equally spaced grid.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    Args:
        theta_longitude (float): The angular difference in longitude (degrees)
            between cell centroids.
        theta_latitude (float): The angular difference in latitude (degrees)
            between cell centroids.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain cells
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        strips (str | None): Option to generate strip-cells along latitude (`"lat"`),
            longitude (`"lon"`), or none (`None`).

    Returns:
        geopandas.GeoDataFrame: the data frame of generated cells, each
            identified by the hash of its geometry (`cell_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """

    # generate indices of grid cells over the filtered region
    indices = generate_indices_uniform_spacing(
        theta_longitude,
        theta_latitude,
        mask,
        strips,
    )
    i = indices[:, 0]
    j = indices[:, 1]
    # get the bounds of the mask
    min_longitude, min_latitude, max_longitude, max_latitude = get_planar_bounds(mask)
    # compute the bounds of each cell, spanning the mask for strips
    west = (
        np.full(len(i), min_longitude)
        if strips == "lat"
        else -180 + i * theta_longitude
    )
    east = (
        np.full(len(i), max_longitude)
        if strips == "lat"
        else -180 + (i + 1) * theta_longitude
    )
    south = (
        np.full(len(j), min_latitude) if strips == "lon" else -90 + j * theta_latitude
    )
    north = (
        np.full(len(j), max_latitude)
        if strips == "lon"
        else -90 + (j + 1) * theta_latitude
    )
    # trace each cell clockwise from its south-east corner
    longitudes = np.column_stack((east, east, west, west, east))
    latitudes = np.column_stack((south, north, north, south, south))
    elevations = np.full_like(longitudes, elevation)
    # create a geodataframe in the WGS84 reference frame
    gdf = gpd.GeoDataFrame(
        geometry=shapely.polygons(
            np.stack((longitudes, latitudes, elevations), axis=-1)
        ),
        crs="EPSG:4326",
    )
    # clip the geodataframe to the supplied mask, if required
    if mask is not None:
        gdf = clip_to_mask(gdf, mask)
        # convert each cell to a convex hull to simplify presentation
        gdf.geometry = gdf.geometry.convex_hull
    # identify each cell by the hash of its geometry
    gdf.insert(0, "cell_id", _hash_geometries(gdf.geometry))
    # return the final geodataframe
    return gdf
