"""
Methods to convert heights above mean sea level (the geoid) to heights above
the WGS 84 ellipsoid, the vertical reference of TAT-C.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import os
from functools import cache

import numpy as np
import numpy.typing as npt

from ..constants import resources_dir

_GEOID_FILE = os.path.join(resources_dir, "egm2008_15.npz")
"""
EGM2008 geoid heights at the nodes of a 15-minute (0.25-degree) grid,
subsampled from the 2.5-minute grid of the National Geospatial-Intelligence
Agency (public domain), as distributed by PROJ (`us_nga_egm08_25.tif`); see
`validation/geoid/build_geoid_grid.py` of the repository.
"""


@cache
def _load_geoid() -> tuple[npt.NDArray, float]:
    """
    Loads the grid of geoid heights.

    Returns:
        tuple[numpy.typing.NDArray, float]: the geoid heights (meters, rows
            from 90 to -90 degrees latitude and columns from -180 to 180
            degrees longitude, inclusive), and the grid spacing (degrees)
    """
    with np.load(_GEOID_FILE) as data:
        heights = data["height_cm"].astype(float) / 100
        spacing = float(data["spacing_deg"])
    # repeat the column at -180 degrees longitude at 180 degrees
    return np.concatenate([heights, heights[:, :1]], axis=1), spacing


def compute_geoid_height(
    longitude: float | npt.ArrayLike, latitude: float | npt.ArrayLike
) -> float | npt.NDArray:
    """
    Computes the height (meters) of the EGM2008 geoid (mean sea level) above
    the WGS 84 ellipsoid, interpolated bilinearly from a 0.25-degree grid
    (with errors of about 0.08 m root mean square, 0.33 m at the 99th
    percentile, and at most 2.7 m, in steep terrain, compared with the
    2.5-minute grid of EGM2008). A height above mean sea level plus the
    geoid height is a height above the ellipsoid, the vertical reference of
    TAT-C's altitudes and elevations: for example, the heights of the
    Copernicus DEM (see `tatc.preprocess.sample_dem_elevation`), flight
    plans, and flight recordings are above mean sea level. The geoid height
    ranges from about -107 to 86 meters.

    Args:
        longitude (float | numpy.typing.ArrayLike): The longitude(s)
            (degrees).
        latitude (float | numpy.typing.ArrayLike): The geodetic latitude(s)
            (degrees, from -90 to 90).

    Returns:
        float | numpy.typing.NDArray: The geoid height(s) (meters), a float
            for scalar coordinates.
    """
    heights, spacing = _load_geoid()
    lon, lat = np.broadcast_arrays(
        np.asarray(longitude, dtype=float), np.asarray(latitude, dtype=float)
    )
    if np.any(np.abs(lat) > 90):
        raise ValueError("latitudes must be from -90 to 90 degrees")
    row = (90 - lat) / spacing
    column = np.mod(lon + 180, 360) / spacing
    i = np.clip(np.floor(row).astype(int), 0, heights.shape[0] - 2)
    j = np.clip(np.floor(column).astype(int), 0, heights.shape[1] - 2)
    u, v = row - i, column - j
    result = (
        (1 - u) * (1 - v) * heights[i, j]
        + (1 - u) * v * heights[i, j + 1]
        + u * (1 - v) * heights[i + 1, j]
        + u * v * heights[i + 1, j + 1]
    )
    return float(result) if result.ndim == 0 else result
