"""
Internal helpers for regular equally-spaced latitude/longitude grids, shared
by the points and cells generation modules.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from ..utils.geometry import get_planar_bounds


def generate_indices_uniform_spacing(
    theta_longitude: float,
    theta_latitude: float,
    mask: Polygon | MultiPolygon | None = None,
    strips: str | None = None,
) -> list:
    """
    Generates a list of indices for an equally spaced grid.

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
        list: list of indices
    """
    # get the bounds of the mask
    min_longitude, min_latitude, max_longitude, max_latitude = get_planar_bounds(mask)

    if strips == "lat":
        # if latitude strips, only generate indices for variable latitude
        return [
            (0, j)
            for j in range(
                int(np.round((min_latitude + 90) / theta_latitude)),
                int(np.round((max_latitude + 90) / theta_latitude)),
            )
        ]
    if strips == "lon":
        # if longitude strips, only generate indices for variable longitude
        return [
            (i, 0)
            for i in range(
                int(np.round((min_longitude + 180) / theta_longitude)),
                int(np.round((max_longitude + 180) / theta_longitude)),
            )
        ]
    # generate indices over the two-dimensional latitude/longitude range
    return [
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
