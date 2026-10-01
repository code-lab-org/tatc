"""
Terrain-based preprocessing utilities: derive a `TerrainMask` for a
ground-based radar station by sampling a digital elevation model (DEM) for
terrain obstructions.

A convenient, public, global, 30-meter-resolution DEM requiring no
authentication is the Copernicus DEM GLO-30 dataset, hosted on AWS Open
Data (https://registry.opendata.aws/copernicus-dem/) as Cloud-Optimized
GeoTIFFs (COGs) tiled on a 1x1-degree grid; `get_copernicus_dem_tile_urls`
builds the tile URL(s) covering a region of interest. Any other DEM
readable by `rasterio` (a local GeoTIFF, another public COG dataset, etc.)
also works.

Requires the optional `rasterio` dependency: `pip install tatc[preprocess]`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math

import numpy as np

try:
    import rasterio
    import rasterio.windows
except ImportError as error:
    raise ImportError(
        "tatc.preprocess.terrain requires the optional 'rasterio' "
        "dependency; install it via `pip install tatc[preprocess]`."
    ) from error

from .. import constants
from ..schemas.surface.radar import RadarStation, TerrainMask
from ..utils.geometry import geodesic_destination
from ..utils.radar import compute_terrain_elevation_angle

_COPERNICUS_DEM_URL_TEMPLATE = (
    "https://copernicus-dem-30m.s3.amazonaws.com/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/"
    "Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif"
)


def get_copernicus_dem_tile_urls(
    longitude: float, latitude: float, radius: float = 0
) -> list[str]:
    """
    Builds the Copernicus DEM GLO-30 tile URL(s) covering a circular region
    of a specified radius around a longitude/latitude, using the dataset's
    1x1-degree tile grid. The dataset is public and requires no
    authentication (see
    https://registry.opendata.aws/copernicus-dem/); each URL is a direct
    HTTPS link to a Cloud-Optimized GeoTIFF (COG), readable by `rasterio`
    (including via `/vsicurl/` for partial, range-request reads without
    downloading the entire file -- `compute_terrain_mask` and
    `sample_dem_elevation` do this automatically).

    Args:
        longitude (float): Longitude (decimal degrees) of the center point.
        latitude (float): Latitude (decimal degrees) of the center point.
        radius (float): Radius (meters) around the center point to cover;
            `0` (the default) returns only the single tile containing the
            center point.

    Returns:
        list[str]: HTTPS URLs to the covering Copernicus DEM GLO-30 tiles.
    """
    lat_buffer = math.degrees(radius / constants.EARTH_MEAN_RADIUS)
    lon_buffer = lat_buffer / max(math.cos(math.radians(latitude)), 0.01)
    min_lat, max_lat = latitude - lat_buffer, latitude + lat_buffer
    min_lon, max_lon = longitude - lon_buffer, longitude + lon_buffer
    urls = []
    for tile_lat in range(math.floor(min_lat), math.floor(max_lat) + 1):
        for tile_lon in range(math.floor(min_lon), math.floor(max_lon) + 1):
            urls.append(
                _COPERNICUS_DEM_URL_TEMPLATE.format(
                    ns="N" if tile_lat >= 0 else "S",
                    lat=abs(tile_lat),
                    ew="E" if tile_lon >= 0 else "W",
                    lon=abs(tile_lon),
                )
            )
    return urls


def _normalize_dem_path(path: str) -> str:
    """
    Normalizes a DEM source path for `rasterio`/GDAL: a plain HTTP(S) URL
    is prefixed with GDAL's `/vsicurl/` virtual file system handler,
    enabling partial, range-request reads of a remote Cloud-Optimized
    GeoTIFF (COG) without downloading the entire file. A local file path
    or a path already using a GDAL virtual file system prefix is returned
    unchanged.

    Args:
        path (str): A local file path, HTTP(S) URL, or GDAL virtual file
            system path.

    Returns:
        str: The normalized path.
    """
    if path.startswith(("http://", "https://")) and not path.startswith("/vsicurl/"):
        return f"/vsicurl/{path}"
    return path


def sample_dem_elevation(
    dem_paths: str | list[str], longitude: float, latitude: float
) -> float:
    """
    Samples terrain elevation at a specified longitude/latitude, from the
    first of one or more DEM raster sources whose extent contains the
    point and reports valid (non-nodata) data. Elevation is reported in
    whatever vertical datum the DEM uses (for most public global DEMs,
    including Copernicus DEM GLO-30, this is close to the EGM2008 geoid,
    a reasonable proxy for the WGS 84 datum used elsewhere in TAT-C).

    Args:
        dem_paths (str | list[str]): One or more DEM raster sources, each
            either a local file path or an HTTP(S) URL to a Cloud-
            Optimized GeoTIFF (COG), tried in order.
        longitude (float): Longitude (decimal degrees) of the sample point.
        latitude (float): Latitude (decimal degrees) of the sample point.

    Returns:
        float: The sampled terrain elevation (meters).

    Raises:
        ValueError: If none of the provided DEM sources cover the
            specified point with valid data.
    """
    paths = [dem_paths] if isinstance(dem_paths, str) else list(dem_paths)
    for path in paths:
        with rasterio.open(_normalize_dem_path(path)) as dataset:
            left, bottom, right, top = dataset.bounds
            if not (left <= longitude <= right and bottom <= latitude <= top):
                continue
            value = next(dataset.sample([(longitude, latitude)]))[0]
            if dataset.nodata is None or value != dataset.nodata:
                return float(value)
    raise ValueError(
        f"No DEM data available at ({longitude}, {latitude}) from the provided source(s)."
    )


def _read_windowed_tile(
    path: str, min_lon: float, min_lat: float, max_lon: float, max_lat: float
) -> tuple[np.ndarray, rasterio.Affine, float | None] | None:
    """
    Opens a DEM source and reads the windowed region overlapping a
    specified longitude/latitude bounding box, or `None` if the source
    does not overlap that box at all. Reading only the needed window (a
    single request for a remote Cloud-Optimized GeoTIFF) avoids
    downloading the entire source, and avoids a separate network request
    per sample point.

    Args:
        path (str): A local file path or HTTP(S) URL to a DEM raster.
        min_lon (float): Minimum longitude (decimal degrees) of the box.
        min_lat (float): Minimum latitude (decimal degrees) of the box.
        max_lon (float): Maximum longitude (decimal degrees) of the box.
        max_lat (float): Maximum latitude (decimal degrees) of the box.

    Returns:
        tuple[numpy.ndarray, rasterio.Affine, float | None] | None: The
        windowed elevation array, its affine transform, and the source's
        nodata value; or `None` if there is no overlap.
    """
    with rasterio.open(_normalize_dem_path(path)) as dataset:
        left, bottom, right, top = dataset.bounds
        if max_lon < left or min_lon > right or max_lat < bottom or min_lat > top:
            return None
        window = (
            rasterio.windows.from_bounds(
                max(min_lon, left),
                max(min_lat, bottom),
                min(max_lon, right),
                min(max_lat, top),
                transform=dataset.transform,
            )
            .round_lengths()
            .round_offsets()
        )
        array = dataset.read(1, window=window)
        transform = rasterio.windows.transform(window, dataset.transform)
        nodata = dataset.nodata
    return array, transform, nodata


def _lookup_elevation(
    tiles: list[tuple[np.ndarray, rasterio.Affine, float | None]],
    longitude: float,
    latitude: float,
) -> float | None:
    """
    Looks up terrain elevation at a specified longitude/latitude from the
    first windowed tile array (as returned by `_read_windowed_tile`) whose
    extent contains it with valid (non-nodata) data.

    Args:
        tiles (list[tuple[numpy.ndarray, rasterio.Affine, float | None]]):
            Windowed elevation arrays, as returned by `_read_windowed_tile`.
        longitude (float): Longitude (decimal degrees) of the sample point.
        latitude (float): Latitude (decimal degrees) of the sample point.

    Returns:
        float | None: The sampled elevation (meters), or `None` if no tile
        covers the point with valid data.
    """
    for array, transform, nodata in tiles:
        col, row = ~transform * (longitude, latitude)
        row_index, col_index = int(row), int(col)
        if 0 <= row_index < array.shape[0] and 0 <= col_index < array.shape[1]:
            value = array[row_index, col_index]
            if nodata is None or value != nodata:
                return float(value)
    return None


def compute_terrain_mask(
    dem_paths: str | list[str],
    longitude: float,
    latitude: float,
    station_elevation: float,
    search_radius: float = 100000,
    min_range: float = 1000,
    number_azimuths: int = 360,
    number_range_samples: int = 100,
) -> TerrainMask:
    """
    Computes a `TerrainMask` for a ground-based radar station by ray-
    tracing a digital elevation model (DEM): at each sampled azimuth,
    terrain elevation is sampled at a series of ranges out to
    `search_radius`, converted to a curvature-corrected elevation angle as
    seen from the station (`tatc.utils.compute_terrain_elevation_angle`),
    and the maximum (highest obstruction) angle along that ray is taken as
    the azimuth's blocking angle. This is the standard approach used by
    radio/radar line-of-sight terrain-masking tools, under the simplifying
    assumption that a single nearest/highest ridge along each azimuth
    dominates (see `tatc.schemas.surface.radar.TerrainMask`).

    `search_radius` defaults to a much shorter distance than a typical
    radar's full hardware range: because the curvature-drop correction
    grows with the square of distance, only unrealistically tall terrain
    far away can ever produce the same blocking angle as modest terrain
    nearby, so nearby terrain dominates in practice. This also keeps the
    DEM read volume (and, for a remote source, the download size) modest.
    Override it if distant, very tall terrain is a specific concern.

    Args:
        dem_paths (str | list[str]): One or more DEM raster sources, each
            either a local file path or an HTTP(S) URL to a Cloud-
            Optimized GeoTIFF (COG); see `get_copernicus_dem_tile_urls`
            for a public, global, no-authentication option. Only sources
            that overlap the search region are read; provide enough
            sources to fully cover it, or some azimuths/ranges may be
            skipped for lack of data (see return value notes).
        longitude (float): Longitude (decimal degrees) of the radar station.
        latitude (float): Latitude (decimal degrees) of the radar station.
        station_elevation (float): Elevation (meters) of the radar
            station's antenna, above the same vertical datum as the DEM
            (typically close to the WGS 84 datum for most public global
            DEMs). This should include antenna tower height, not just
            bare-ground elevation; see `sample_dem_elevation` to determine
            the ground elevation at the station's location.
        search_radius (float): The maximum ground distance (meters) to
            examine for obstructions.
        min_range (float): The minimum ground distance (meters) to
            examine, excluding the immediate vicinity of the antenna itself.
        number_azimuths (int): The number of evenly-spaced azimuth samples
            (one terrain profile per azimuth) spanning a full revolution.
        number_range_samples (int): The number of range samples per
            azimuthal profile, evenly spaced between `min_range` and
            `search_radius`.

    Returns:
        TerrainMask: The resulting azimuthal terrain blockage mask. An
        azimuth with no valid DEM data at any sampled range (e.g. the
        provided `dem_paths` do not cover that direction) reports a
        minimum elevation angle of `-90` degrees (no blockage assumed).

    Raises:
        ValueError: If none of the provided DEM sources overlap the
            search region at all.
    """
    paths = [dem_paths] if isinstance(dem_paths, str) else list(dem_paths)
    lat_buffer = math.degrees(search_radius / constants.EARTH_MEAN_RADIUS)
    lon_buffer = lat_buffer / max(math.cos(math.radians(latitude)), 0.01)
    tiles = [
        tile
        for path in paths
        if (
            tile := _read_windowed_tile(
                path,
                longitude - lon_buffer,
                latitude - lat_buffer,
                longitude + lon_buffer,
                latitude + lat_buffer,
            )
        )
        is not None
    ]
    if not tiles:
        raise ValueError(
            "None of the provided DEM sources overlap the search region "
            f"around ({longitude}, {latitude})."
        )
    azimuths = np.linspace(0, 360, number_azimuths, endpoint=False)
    ranges = np.linspace(min_range, search_radius, number_range_samples)
    min_elevation_angles = []
    for azimuth in azimuths:
        max_angle = -90.0
        for distance in ranges:
            sample_longitude, sample_latitude = geodesic_destination(
                longitude, latitude, float(azimuth), float(distance)
            )
            terrain_elevation = _lookup_elevation(
                tiles, sample_longitude, sample_latitude
            )
            if terrain_elevation is None:
                continue
            angle = compute_terrain_elevation_angle(
                distance, terrain_elevation, station_elevation
            )
            max_angle = max(max_angle, angle)
        min_elevation_angles.append(max_angle)
    return TerrainMask(
        azimuth=[float(a) for a in azimuths], min_elevation_angle=min_elevation_angles
    )


def compute_terrain_mask_for_station(
    dem_paths: str | list[str],
    station: RadarStation,
    search_radius: float = 100000,
    min_range: float = 1000,
    number_azimuths: int = 360,
    number_range_samples: int = 100,
) -> TerrainMask:
    """
    Convenience wrapper around `compute_terrain_mask` that reads station
    location and elevation directly from a `RadarStation`.

    Args:
        dem_paths (str | list[str]): See `compute_terrain_mask`.
        station (RadarStation): The radar station (its `longitude`,
            `latitude`, and `elevation` are used).
        search_radius (float): See `compute_terrain_mask`.
        min_range (float): See `compute_terrain_mask`.
        number_azimuths (int): See `compute_terrain_mask`.
        number_range_samples (int): See `compute_terrain_mask`.

    Returns:
        TerrainMask: The resulting azimuthal terrain blockage mask.
    """
    return compute_terrain_mask(
        dem_paths,
        station.longitude,
        station.latitude,
        station.elevation,
        search_radius=search_radius,
        min_range=min_range,
        number_azimuths=number_azimuths,
        number_range_samples=number_range_samples,
    )
