"""
Raster preprocessing utilities: read gridded weights (for example,
population counts or densities) from a GeoTIFF for weighted random point
generation (`tatc.generation.generate_points_random`).

Requires the optional `rasterio` dependency: `pip install tatc[preprocess]`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

try:
    import rasterio
    from rasterio.windows import Window
except ImportError as error:
    raise ImportError(
        "tatc.preprocess.raster requires the optional 'rasterio' "
        "dependency; install it via `pip install tatc[preprocess]`."
    ) from error

from ..utils.geometry import get_planar_bounds
from .terrain import _normalize_dem_path

# maximum number of raster cells read at once
_MAX_STRIP_CELLS = 2**24


def _get_aligned_range(start: int, stop: int, size: int, aggregate: int) -> range:
    """
    Expands a range of raster indices to a whole number of aggregate blocks,
    shifted to stay within the raster where possible.

    Args:
        start (int): The first index.
        stop (int): The index after the last.
        size (int): The number of indices in the raster.
        aggregate (int): The number of indices in each block.

    Returns:
        range: The expanded range, which may extend beyond the raster only
            if it is smaller than a whole number of blocks.
    """
    stop = start + math.ceil((stop - start) / aggregate) * aggregate
    shift = max(min(stop - size, start), 0)
    return range(start - shift, stop - shift)


def read_raster_weights(
    path: str,
    mask: Polygon | MultiPolygon | None = None,
    aggregate: int = 1,
    density: bool = False,
    band: int = 1,
) -> tuple[np.ndarray, tuple[float, float, float, float]]:
    """
    Reads a raster of weights (for example, population counts or densities)
    from a GeoTIFF in WGS84 (EPSG:4326) geodetic coordinates, for
    `tatc.generation.generate_points_random`. Only the part of the raster
    overlapping the bounds of the mask is read. Missing (nodata or
    not-a-number) values have zero weight.

    Large rasters (for example, a global population count at 30 arc-seconds)
    can be aggregated to coarser cells of `aggregate` by `aggregate` raster
    cells: summed for totals per cell (counts) or averaged for values per
    unit area (`density`).

    Args:
        path (str): A local file path or HTTP(S) URL to a raster.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        aggregate (int): The number of raster cells (in each direction) to
            combine into each weight.
        density (bool): True, if values are per unit area (averaged when
            aggregated) rather than totals per cell (summed when aggregated).
        band (int): The raster band to read.

    Returns:
        tuple[numpy.ndarray, tuple[float, float, float, float]]: The weights,
            with rows from north to south and columns from west to east, and
            their bounds (west, south, east, north) in degrees.
    """
    if aggregate < 1:
        raise ValueError(f"Aggregate must be at least 1, got {aggregate}.")
    min_longitude, min_latitude, max_longitude, max_latitude = get_planar_bounds(mask)
    with rasterio.open(_normalize_dem_path(path)) as dataset:
        if dataset.crs is None or dataset.crs.to_epsg() != 4326:
            raise ValueError(
                f"Raster must use WGS84 (EPSG:4326) coordinates, not {dataset.crs}; "
                "reproject it first (for example, `gdalwarp -t_srs EPSG:4326`)."
            )
        transform = dataset.transform
        if transform.b != 0 or transform.d != 0 or transform.e >= 0:
            raise ValueError("Raster must be north-up without rotation.")
        left, bottom, right, top = dataset.bounds
        if (
            min(left, min_longitude) < -180 or max(right, max_longitude) > 180
        ) and not left <= min_longitude <= max_longitude <= right:
            # read all longitudes for masks or rasters beyond 180 degrees,
            # which may overlap at longitudes a whole turn apart
            min_longitude, max_longitude = left, right
        if (
            max_longitude <= left
            or min_longitude >= right
            or max_latitude <= bottom
            or min_latitude >= top
        ):
            raise ValueError("Raster does not overlap the mask.")
        width = transform.a
        height = -transform.e
        # find the raster indices overlapping the bounds
        columns = _get_aligned_range(
            max(math.floor((min_longitude - left) / width), 0),
            min(math.ceil((max_longitude - left) / width), dataset.width),
            dataset.width,
            aggregate,
        )
        rows = _get_aligned_range(
            max(math.floor((top - max_latitude) / height), 0),
            min(math.ceil((top - min_latitude) / height), dataset.height),
            dataset.height,
            aggregate,
        )
        west = left + columns.start * width
        east = left + columns.stop * width
        north = top - rows.start * height
        south = top - rows.stop * height
        if south < -90 or north > 90:
            raise ValueError(
                f"Aggregate {aggregate} extends the raster beyond the poles; "
                "choose an aggregate that divides its number of rows."
            )
        # read strips of whole blocks to limit memory
        strip = aggregate * max(1, _MAX_STRIP_CELLS // (len(columns) * aggregate))
        weights = []
        for row in range(rows.start, rows.stop, strip):
            # read the strip, padding any part beyond the raster
            read_rows = range(max(row, 0), min(row + strip, rows.stop, dataset.height))
            read_columns = range(
                max(columns.start, 0), min(columns.stop, dataset.width)
            )
            values = np.zeros((min(strip, rows.stop - row), len(columns)))
            valid = np.zeros(values.shape)
            data = dataset.read(
                band,
                window=Window(
                    read_columns.start,
                    read_rows.start,
                    len(read_columns),
                    len(read_rows),
                ),
                masked=True,
            )
            target = (
                slice(read_rows.start - row, read_rows.stop - row),
                slice(
                    read_columns.start - columns.start,
                    read_columns.stop - columns.start,
                ),
            )
            values[target] = np.nan_to_num(data.astype(np.float64).filled(0.0))
            valid[target] = 1
            # combine blocks of cells
            shape = (
                values.shape[0] // aggregate,
                aggregate,
                values.shape[1] // aggregate,
                aggregate,
            )
            total = values.reshape(shape).sum(axis=(1, 3))
            if density:
                # average over the cells within the raster
                count = valid.reshape(shape).sum(axis=(1, 3))
                total = np.divide(
                    total, count, out=np.zeros_like(total), where=count > 0
                )
            weights.append(total)
    return np.concatenate(weights), (west, south, east, north)
