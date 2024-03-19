# -*- coding: utf-8 -*-
"""
Methods to generate geospatial cells to aggregate data.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import Optional, Union

import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon

from .points import (
    _compute_cubed_sphere_point_id,
    _generate_cubed_sphere_indices,
    _get_bounds,
)

from ..constants import EARTH_MEAN_RADIUS


def generate_cubed_sphere_cells(
    distance: float,
    elevation: float = 0,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    strips: str = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic polygons over a regular cubed-sphere grid.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    Args:
        distance (float):  The typical surface distance (meters) between points.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (Polygon or MultiPolygon):  An optional mask to constrain cells
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        strips (str): Option to generate strip-cells along latitude (`"lat"`),
            longitude (`"lon"`), or none (`None`).

    Returns:
        geopandas.GeoDataFrame: the data frame of generated cells
    """
    # compute the angular disance of each sample (assuming sphere)
    theta_longitude = np.degrees(distance / EARTH_MEAN_RADIUS)
    theta_latitude = np.degrees(distance / EARTH_MEAN_RADIUS)
    return _generate_cubed_sphere_cells(
        theta_longitude, theta_latitude, elevation, mask, strips
    )


def _generate_cubed_sphere_cells(
    theta_longitude: float,
    theta_latitude: float,
    elevation: float = 0,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    strips: str = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic polygons over a regular cubed-sphere grid.

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
        mask (Polygon or MultiPolygon):  An optional mask to constrain cells
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        strips (str): Option to generate strip-cells along latitude (`"lat"`),
            longitude (`"lon"`), or none (`None`).

    Returns:
        geopandas.GeoDataFrame: the data frame of generated cells
    """

    # generate indices of grid cells over the filtered region
    indices = _generate_cubed_sphere_indices(
        theta_longitude,
        theta_latitude,
        mask,
        strips,
    )
    # get the bounds of the mask
    min_longitude, min_latitude, max_longitude, max_latitude = _get_bounds(mask)
    # create a geodataframe in the WGS84 reference frame
    gdf = gpd.GeoDataFrame(
        {
            "cell_id": [
                _compute_cubed_sphere_point_id(i, j, theta_longitude, theta_latitude)
                for (i, j) in indices
            ],
            "geometry": [
                Polygon(
                    [
                        (
                            max_longitude
                            if strips == "lat"
                            else -180 + (i + 1) * theta_longitude,
                            min_latitude
                            if strips == "lon"
                            else -90 + j * theta_latitude,
                            elevation,
                        ),
                        (
                            max_longitude
                            if strips == "lat"
                            else -180 + (i + 1) * theta_longitude,
                            max_latitude
                            if strips == "lon"
                            else -90 + (j + 1) * theta_latitude,
                            elevation,
                        ),
                        (
                            min_longitude
                            if strips == "lat"
                            else -180 + i * theta_longitude,
                            max_latitude
                            if strips == "lon"
                            else -90 + (j + 1) * theta_latitude,
                            elevation,
                        ),
                        (
                            min_longitude
                            if strips == "lat"
                            else -180 + i * theta_longitude,
                            min_latitude
                            if strips == "lon"
                            else -90 + j * theta_latitude,
                            elevation,
                        ),
                        (
                            max_longitude
                            if strips == "lat"
                            else -180 + (i + 1) * theta_longitude,
                            min_latitude
                            if strips == "lon"
                            else -90 + j * theta_latitude,
                            elevation,
                        ),
                    ]
                )
                for (i, j) in indices
            ],
        },
        crs="EPSG:4326",
    )
    # clip the geodataframe to the supplied mask, if required
    if mask is not None:
        gdf = gpd.clip(gdf, mask).reset_index(drop=True)
        # convert each cell to a convex hull to simplify presentation
        gdf.geometry = gdf.geometry.convex_hull
    # return the final geodataframe
    return gdf
