"""
Helpers shared by the exporters: the places and times of observations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from ..analysis.check import _get_observation_times


def get_places(
    observations: gpd.GeoDataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Gets the longitude and latitude (degrees) of each observation (of its
    point, or the representative point of its region), and its height
    (meters, the z-coordinate of a point, or `NaN`).
    """
    geometry = observations.geometry.to_numpy()
    is_point = shapely.get_type_id(geometry) == shapely.GeometryType.POINT
    points = np.where(is_point, geometry, shapely.point_on_surface(geometry))
    height = np.where(
        is_point, shapely.get_z(np.where(is_point, geometry, None)), np.nan
    ).astype(float)
    return (
        shapely.get_x(points).astype(float),
        shapely.get_y(points).astype(float),
        height,
    )


def get_times(observations: pd.DataFrame, time_column: str | None) -> pd.DatetimeIndex:
    """Gets the time (UTC) of each observation (see `_get_observation_times`)."""
    return pd.DatetimeIndex(_get_observation_times(observations, time_column))
