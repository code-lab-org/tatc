"""
Internal helpers for regular equally-spaced latitude/longitude grids, shared
by the points and cells generation modules.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import pandas as pd
import shapely
from shapely.geometry import MultiPolygon, Polygon

from ..utils.geometry import get_planar_bounds


def generate_indices_uniform_spacing(
    theta_longitude: float,
    theta_latitude: float,
    mask: Polygon | MultiPolygon | None = None,
    strips: str | None = None,
) -> npt.NDArray[np.int64]:
    """
    Generates the indices for an equally spaced grid.

    Args:
        theta_longitude (float): The angular difference in longitude (degrees)
            between points.
        theta_latitude (float): The angular difference in latitude (degrees)
            between points.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain points
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        strips (str | None): Option to generate one-dimensional strips along latitude
            (`"lat"`), longitude (`"lon"`), or none (`None`).

    Returns:
        numpy.ndarray: array (n, 2) of (longitude, latitude) indices, ordered
            by latitude then longitude index
    """
    # get the bounds of the mask
    min_longitude, min_latitude, max_longitude, max_latitude = get_planar_bounds(mask)

    # if latitude strips, only generate indices for variable latitude
    i_range = (
        np.arange(1)
        if strips == "lat"
        else np.arange(
            int(np.round((min_longitude + 180) / theta_longitude)),
            int(np.round((max_longitude + 180) / theta_longitude)),
        )
    )
    # if longitude strips, only generate indices for variable longitude
    j_range = (
        np.arange(1)
        if strips == "lon"
        else np.arange(
            int(np.round((min_latitude + 90) / theta_latitude)),
            int(np.round((max_latitude + 90) / theta_latitude)),
        )
    )
    # generate indices over the two-dimensional latitude/longitude range
    i, j = np.meshgrid(i_range, j_range)
    return np.column_stack((i.ravel(), j.ravel())).astype(np.int64)


def clip_to_mask(
    gdf: gpd.GeoDataFrame, mask: Polygon | MultiPolygon
) -> gpd.GeoDataFrame:
    """
    Clips a data frame of geometries to a mask, retaining row order. Only
    geometries not properly contained by the mask are clipped (see
    `geopandas.clip`), as the others are unchanged by clipping.

    Args:
        gdf (geopandas.GeoDataFrame): The data frame of geometries.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The mask
            using WGS84 (EPSG:4326) geodetic coordinates.

    Returns:
        geopandas.GeoDataFrame: the clipped data frame, with a new index
    """
    shapely.prepare(mask)
    inside = shapely.contains_properly(mask, gdf.geometry.values)
    if inside.all():
        return gdf.reset_index(drop=True)
    return (
        pd.concat([gdf[inside], gpd.clip(gdf[~inside], mask)])
        .sort_index()
        .reset_index(drop=True)
    )
