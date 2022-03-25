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


def generate_cubed_sphere_cells(
    distance: float,
    mask: Optional[Union[Polygon, MultiPolygon]] = None,
    strips: str = None,
):
    """
    Generates geodetic polygons over a regular cubed-sphere grid.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    :param distance: The typical surface distance (meters) between points.
    :type distance: float
    :param strips: An optional mask to constrain generated cells, defaults to None
    :type mask: :class:`shapely.geometry.Polygon` or :class:`shapely.geometry.MultiPolygon`, optional
    :param strips: An optional argument to generate strips of
        constant latitude (`lat`) or strips of constant longitude (`lon`), defaults to None
    :type strips: str, optional
    :return: An instance of :class:`geopandas.GeoDataFrame` specifying the cells.
    :rtype: :class:`geopandas.GeoDataFrame`
    """

    @njit
    def _compute_id(i, j, theta):
        """
        Fast method to compute the flattened id for a cubed sphere grid point.
        Indices increment west-to-east followed by south-to-north with a first
        point at -180 degrees latitude and close to -90 degrees latitude.

        :param i: The zero-based longitude index
        :type i: int
        :param j: The zero-based latitude index
        :type j: int
        :return: The latitude (degrees) of this point.
        :rtype: float
        """
        return int(j * int(360 / theta) + np.mod(i, int(360 / theta)))

    # compute the angular disance of each sample (assuming sphere)
    theta = np.degrees(distance / earth_mean_radius)
    if isinstance(mask, Polygon) or isinstance(mask, MultiPolygon):
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
                int(np.round((min_latitude + 90) / theta)),
                int(np.round((max_latitude + 90) / theta)),
            )
        ]
    elif strips == "lon":
        # if longitude strips, only generate grid cells for variable longitude
        indices = [
            (i, 0)
            for i in range(
                int(np.round((min_longitude + 180) / theta)),
                int(np.round((max_longitude + 180) / theta)),
            )
        ]
    else:
        # generate grid cells over the filtered latitude/longitude range
        indices = [
            (i, j)
            for j in range(
                int(np.round((min_latitude + 90) / theta)),
                int(np.round((max_latitude + 90) / theta)),
            )
            for i in range(
                int(np.round((min_longitude + 180) / theta)),
                int(np.round((max_longitude + 180) / theta)),
            )
        ]
    # create a geodataframe in the WGS84 reference frame
    gdf = gpd.GeoDataFrame(
        {
            "cell_id": [_compute_id(i, j, theta) for (i, j) in indices],
            "geometry": [
                box(
                    min_longitude if strips == "lat" else -180 + i * theta,
                    min_latitude if strips == "lon" else -90 + j * theta,
                    max_longitude if strips == "lat" else -180 + (i + 1) * theta,
                    max_latitude if strips == "lon" else -90 + (j + 1) * theta,
                )
                for (i, j) in indices
            ],
        },
        crs="EPSG:4326",
    )
    # clip the geodataframe to the supplied mask, if required
    if mask is not None:
        try:
            gdf = gpd.clip(gdf, mask).reset_index(drop=True)
            # convert each cell to a convex hull
            gdf.geometry = gdf.geometry.convex_hull
        except TopologicalError:
            pass
    # return the final geodataframe
    return gdf
