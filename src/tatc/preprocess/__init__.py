"""
Preprocessing utilities that derive TAT-C model inputs from external
geospatial datasets (e.g. a digital elevation model). These utilities
require additional dependencies not installed by default; install them via
`pip install tatc[preprocess]`.
"""

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
    "sample_dem_elevation",
]
