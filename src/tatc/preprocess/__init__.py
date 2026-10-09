"""
Preprocessing utilities that derive TAT-C model inputs from external
geospatial datasets (e.g. a digital elevation model). These utilities
require additional dependencies not installed by default; install them via
`pip install tatc[preprocess]`.
"""

from .raster import read_raster_weights
from .terrain import (
    compute_terrain_mask,
    compute_terrain_mask_for_station,
    get_copernicus_dem_tile_urls,
    sample_dem_elevation,
)

__all__ = [
    "compute_terrain_mask",
    "compute_terrain_mask_for_station",
    "get_copernicus_dem_tile_urls",
    "read_raster_weights",
    "sample_dem_elevation",
]
