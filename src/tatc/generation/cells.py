# -*- coding: utf-8 -*-
"""
Methods to generate geospatial cells to aggregate data.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

import numpy as np
import geopandas as gpd
from numba import njit
from shapely.geometry import box, Polygon, MultiPolygon
from shapely.errors import TopologicalError
from typing import Optional, Union

from ..constants import earth_mean_radius


@staticmethod
def generate_cubed_sphere_cells(
    distance: float,
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
        mask (Polygon or MultiPolygon):  An optional mask to constrain cells
                using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        strips (str): Option to generate strip-cells along latitude (`"lat"`), longitude (`"lon"`), or none (`None`).

    Returns:
        geopandas.GeoDataFrame: the data frame of generated cells
    """
    # compute the angular disance of each sample (assuming sphere)
    theta_longitude = np.degrees(distance / earth_mean_radius)
    theta_latitude = np.degrees(distance / earth_mean_radius)
    return _generate_cubed_sphere_cells(theta_longitude, theta_latitude, mask, strips)


@staticmethod
def _generate_cubed_sphere_cells(
    theta_longitude: float,
    theta_latitude: float,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    strips: str = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic polygons over a regular cubed-sphere grid.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    Args:
        theta_longitude (float): The angular difference in longitude (degrees) between cell centroids.
        theta_latitude (float): The angular difference in latitude (degrees) between cell centroids.
        mask (Polygon or MultiPolygon):  An optional mask to constrain cells
                using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        strips (str): Option to generate strip-cells along latitude (`"lat"`), longitude (`"lon"`), or none (`None`).

    Returns:
        geopandas.GeoDataFrame: the data frame of generated cells
    """

    @njit
    def _compute_id(i, j, theta_i, theta_j):
        """
        Fast method to compute the flattened id for a cubed sphere grid cell.
        Indices increment west-to-east followed by south-to-north with a first
        point at -180 degrees latitude and close to -90 degrees latitude.

        Args:
            i (int): The zero-based longitude index.
            j (int): The zero-based latitude index.
            theta_i (float): The angular step in longitude (degrees).
            theta_j (float): The angular step in latitude (degrees).

        Returns:
            int: The id of this cell.
        """
        return int(j * int(360 / theta_j) + np.mod(i, int(360 / theta_i)))

    if isinstance(mask, Polygon) or isinstance(mask, MultiPolygon):
        if not mask.is_valid:
            raise ValueError("Mask is not a valid Polygon or MultiPolygon.")
        total_bounds = mask.bounds
    else:
        total_bounds = [-180, -90, 180, 90]
    min_longitude = total_bounds[0]
    min_latitude = total_bounds[1]
    max_longitude = 180 if total_bounds[2] == -180 else total_bounds[2]
    max_latitude = total_bounds[3]

    if strips == "lat":
        # if latitude strips, only generate grid cells for variable latitude
        indices = [
            (0, j)
            for j in range(
                int(np.round((min_latitude + 90) / theta_latitude)),
                int(np.round((max_latitude + 90) / theta_latitude)),
            )
        ]
    elif strips == "lon":
        # if longitude strips, only generate grid cells for variable longitude
        indices = [
            (i, 0)
            for i in range(
                int(np.round((min_longitude + 180) / theta_longitude)),
                int(np.round((max_longitude + 180) / theta_longitude)),
            )
        ]
    else:
        # generate grid cells over the filtered latitude/longitude range
        indices = [
            (i, j)
            for j in range(
                int(np.round((min_latitude + 90) / theta_latitude)),
                int(np.round((max_latitude + 90) / theta_latitude)),
            )
            for i in range(
                int(np.round((min_longitude + 180) / theta_longitude)),
                int(np.round((max_longitude + 180) / theta_longitude)),
            )
        ]
    # create a geodataframe in the WGS84 reference frame
    gdf = gpd.GeoDataFrame(
        {
            "cell_id": [
                _compute_id(i, j, theta_longitude, theta_latitude) for (i, j) in indices
            ],
            "geometry": [
                box(
                    min_longitude if strips == "lat" else -180 + i * theta_longitude,
                    min_latitude if strips == "lon" else -90 + j * theta_latitude,
                    max_longitude
                    if strips == "lat"
                    else -180 + (i + 1) * theta_longitude,
                    max_latitude if strips == "lon" else -90 + (j + 1) * theta_latitude,
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
