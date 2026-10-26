"""
Terrain-based preprocessing utilities: derive a `TerrainMask` for a
ground-based radar station or a fixed surface platform by sampling a
digital elevation model (DEM) for terrain obstructions, and the terrain
slope, aspect, and local incidence angle at observed targets.

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

import geopandas as gpd
import numpy as np
import pandas as pd

try:
    import rasterio
    import rasterio.windows
except ImportError as error:
    raise ImportError(
        "tatc.preprocess.terrain requires the optional 'rasterio' "
        "dependency; install it via `pip install tatc[preprocess]`."
    ) from error

from .. import constants
from ..schemas.surface import RadarStation, SurfacePlatform, TerrainMask
from ..utils.geometry import geodesic_destination
from ..utils.observation import compute_local_incidence_angle
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
    whatever vertical datum the DEM uses: for most public global DEMs,
    including Copernicus DEM GLO-30, this is mean sea level (the EGM2008
    geoid), which differs from the WGS 84 ellipsoid by the local geoid
    undulation (up to about 100 m in magnitude).

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
    path: str,
    min_lon: float,
    min_lat: float,
    max_lon: float,
    max_lat: float,
    padding: int = 0,
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
        padding (int): The number of cells by which to extend the window
            on each side (within the source), as for interpolation between
            the cells around points at the edge of the box.

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
        if padding:
            window = rasterio.windows.Window(
                window.col_off - padding,
                window.row_off - padding,
                window.width + 2 * padding,
                window.height + 2 * padding,
            ).intersection(rasterio.windows.Window(0, 0, dataset.width, dataset.height))
        array = dataset.read(1, window=window)
        transform = rasterio.windows.transform(window, dataset.transform)
        nodata = dataset.nodata
    return array, transform, nodata


def _to_pixel_coordinates(
    transform: rasterio.Affine, longitudes: np.ndarray, latitudes: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Converts longitudes/latitudes to fractional (column, row) pixel
    coordinates of a raster with an affine transform, applying the inverse
    transform's coefficients directly (rather than its operators, which
    differ between versions of `affine`).

    Args:
        transform (rasterio.Affine): The raster's affine transform.
        longitudes (numpy.ndarray): Longitudes (decimal degrees).
        latitudes (numpy.ndarray): Latitudes (decimal degrees).

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: The column and row coordinates.
    """
    inverse = ~transform
    longitudes, latitudes = np.asarray(longitudes), np.asarray(latitudes)
    return (
        inverse.a * longitudes + inverse.b * latitudes + inverse.c,
        inverse.d * longitudes + inverse.e * latitudes + inverse.f,
    )


def _lookup_elevations(
    tiles: list[tuple[np.ndarray, rasterio.Affine, float | None]],
    longitudes: np.ndarray,
    latitudes: np.ndarray,
) -> np.ndarray:
    """
    Looks up terrain elevations at specified longitudes/latitudes, each
    from the first windowed tile array (as returned by
    `_read_windowed_tile`) whose extent contains it with valid
    (non-nodata) data.

    Args:
        tiles (list[tuple[numpy.ndarray, rasterio.Affine, float | None]]):
            Windowed elevation arrays, as returned by `_read_windowed_tile`.
        longitudes (numpy.ndarray): Longitudes (decimal degrees) of the
            sample points.
        latitudes (numpy.ndarray): Latitudes (decimal degrees) of the
            sample points, with the same shape as `longitudes`.

    Returns:
        numpy.ndarray: The sampled elevations (meters), with the same
        shape as `longitudes`; `numpy.nan` where no tile covers a point
        with valid data.
    """
    elevations = np.full(np.shape(longitudes), np.nan)
    for array, transform, nodata in tiles:
        cols, rows = _to_pixel_coordinates(transform, longitudes, latitudes)
        rows, cols = np.floor(rows).astype(int), np.floor(cols).astype(int)
        candidates = (
            np.isnan(elevations)
            & (rows >= 0)
            & (rows < array.shape[0])
            & (cols >= 0)
            & (cols < array.shape[1])
        )
        values = array[rows[candidates], cols[candidates]].astype(float)
        if nodata is not None:
            values[values == nodata] = np.nan
        elevations[candidates] = values
    return elevations


def compute_terrain_mask(
    dem_paths: str | list[str],
    longitude: float,
    latitude: float,
    station_elevation: float,
    search_radius: float = 100000,
    min_range: float = 1000,
    number_azimuths: int = 360,
    number_range_samples: int = 400,
    range_spacing: str = "geometric",
    effective_earth_radius_factor: float = constants.EFFECTIVE_EARTH_RADIUS_FACTOR,
) -> TerrainMask:
    """
    Computes a `TerrainMask` for a ground-based radar station (or other
    fixed antenna or instrument) by ray-tracing a digital elevation model
    (DEM): at each sampled azimuth,
    terrain elevation is sampled at a series of ranges out to
    `search_radius`, converted to a curvature-corrected elevation angle as
    seen from the station (`tatc.utils.compute_terrain_elevation_angle`),
    and the maximum (highest obstruction) angle along that ray is taken as
    the azimuth's blocking angle. This is the standard approach used by
    radio/radar line-of-sight terrain-masking tools, under the simplifying
    assumption that a single nearest/highest ridge along each azimuth
    dominates (see `tatc.schemas.surface.radar.TerrainMask`).

    `search_radius` defaults to a shorter distance than a typical radar's
    full hardware range, which keeps the DEM read volume (and, for a
    remote source, the download size) modest. The curvature-drop
    correction grows with the square of distance (about 0.3 degrees of
    angle at 75 km), so distant terrain must be much taller than nearby
    terrain to block the same angle; in mountainous regions, however,
    ranges 50 to 100 km away commonly set the blocking angle (validation
    against NEXRAD radars in Arizona found differences of up to about
    1.4 degrees between 50 and 100 km search radii). Choose a radius that
    reaches the significant relief around the station, and a range
    sampling (`number_range_samples`) fine enough not to skip narrow
    ridge crests close to the station, where small sampling errors
    produce the largest angle errors.

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
            (mean sea level for Copernicus DEM GLO-30 and most public
            global DEMs, not the WGS 84 ellipsoid; NOAA's NEXRAD station
            tables report antenna heights above mean sea level). This
            should include antenna tower height, not just
            bare-ground elevation; see `sample_dem_elevation` to determine
            the ground elevation at the station's location.
        search_radius (float): The maximum ground distance (meters) to
            examine for obstructions.
        min_range (float): The minimum ground distance (meters) to
            examine, excluding the immediate vicinity of the antenna itself.
        number_azimuths (int): The number of evenly-spaced azimuth samples
            (one terrain profile per azimuth) spanning a full revolution.
        number_range_samples (int): The number of range samples per
            azimuthal profile, between `min_range` and `search_radius`.
        range_spacing (str): The spacing of range samples: `"geometric"`
            (the default) spaces samples in proportion to their range,
            densely near the station and sparsely far away, so a missed
            ridge crest causes a similar angle error at every range;
            `"linear"` spaces samples evenly.
        effective_earth_radius_factor (float): The ratio of the effective
            Earth radius to the mean Earth radius, accounting for refraction
            (see `tatc.utils.compute_terrain_elevation_angle`): 4/3 for radio
            waves in a standard atmosphere (the default), or 1 for straight
            lines of sight (as for the footprints of a `SurfacePlatform`).

    Returns:
        TerrainMask: The resulting azimuthal terrain blockage mask. An
        azimuth with no valid DEM data at any sampled range (e.g. the
        provided `dem_paths` do not cover that direction) reports a
        minimum elevation angle of `-90` degrees (no blockage assumed).

    Raises:
        ValueError: If `range_spacing` is not `"geometric"` or `"linear"`,
            if `min_range` is not positive for geometric spacing, or if
            none of the provided DEM sources overlap the search region.
    """
    if range_spacing == "geometric":
        if min_range <= 0:
            raise ValueError("min_range must be positive for geometric range spacing.")
        ranges = np.geomspace(min_range, search_radius, number_range_samples)
    elif range_spacing == "linear":
        ranges = np.linspace(min_range, search_radius, number_range_samples)
    else:
        raise ValueError(
            f"range_spacing must be 'geometric' or 'linear', not '{range_spacing}'."
        )
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
    azimuth_grid, range_grid = np.meshgrid(azimuths, ranges, indexing="ij")
    sample_longitudes, sample_latitudes = geodesic_destination(
        np.full(azimuth_grid.shape, longitude),
        np.full(azimuth_grid.shape, latitude),
        azimuth_grid,
        range_grid,
    )
    terrain_elevations = _lookup_elevations(tiles, sample_longitudes, sample_latitudes)
    angles = compute_terrain_elevation_angle(
        range_grid, terrain_elevations, station_elevation, effective_earth_radius_factor
    )
    # the highest obstruction along each azimuth (-90 if no valid DEM data)
    min_elevation_angles = np.max(np.where(np.isnan(angles), -90.0, angles), axis=1)
    return TerrainMask(
        azimuth=[float(a) for a in azimuths],
        min_elevation_angle=[float(a) for a in min_elevation_angles],
    )


def compute_terrain_mask_for_station(
    dem_paths: str | list[str],
    station: RadarStation,
    search_radius: float = 100000,
    min_range: float = 1000,
    number_azimuths: int = 360,
    number_range_samples: int = 400,
    range_spacing: str = "geometric",
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
        range_spacing (str): See `compute_terrain_mask`.

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
        range_spacing=range_spacing,
    )


def compute_terrain_mask_for_platform(
    dem_paths: str | list[str],
    platform: SurfacePlatform,
    search_radius: float = 100000,
    min_range: float = 1000,
    number_azimuths: int = 360,
    number_range_samples: int = 400,
    range_spacing: str = "geometric",
    effective_earth_radius_factor: float = 1,
) -> TerrainMask:
    """
    Convenience wrapper around `compute_terrain_mask` for the instruments of
    a `SurfacePlatform`, a `height` above the ground at its location. The
    ground elevation is sampled from the DEM (see `sample_dem_elevation`),
    in the DEM's own vertical datum, rather than taken from the platform's
    `elevation` (above the WGS 84 ellipsoid), so the instruments' height
    above the surrounding terrain does not depend on the geoid undulation.

    Args:
        dem_paths (str | list[str]): See `compute_terrain_mask`.
        platform (SurfacePlatform): The platform (its `longitude`,
            `latitude`, and `height` are used).
        search_radius (float): See `compute_terrain_mask`.
        min_range (float): See `compute_terrain_mask`.
        number_azimuths (int): See `compute_terrain_mask`.
        number_range_samples (int): See `compute_terrain_mask`.
        range_spacing (str): See `compute_terrain_mask`.
        effective_earth_radius_factor (float): See `compute_terrain_mask`.
            Defaults to 1 (straight lines of sight, without refraction), as
            for the platform's footprints.

    Returns:
        TerrainMask: The resulting azimuthal terrain blockage mask.
    """
    ground = sample_dem_elevation(dem_paths, platform.longitude, platform.latitude)
    return compute_terrain_mask(
        dem_paths,
        platform.longitude,
        platform.latitude,
        ground + platform.height,
        search_radius=search_radius,
        min_range=min_range,
        number_azimuths=number_azimuths,
        number_range_samples=number_range_samples,
        range_spacing=range_spacing,
        effective_earth_radius_factor=effective_earth_radius_factor,
    )


def _interpolate_elevations(
    tiles: list[tuple[np.ndarray, rasterio.Affine, float | None]],
    longitudes: np.ndarray,
    latitudes: np.ndarray,
) -> np.ndarray:
    """
    Interpolates terrain elevations bilinearly between the centers of DEM
    cells at specified longitudes/latitudes, each on the grid of the first
    windowed tile array (as returned by `_read_windowed_tile`) within a cell
    of it whose surrounding cell centers have valid (non-nodata) data, in
    any tile (so that points between the cells of adjacent, aligned tiles,
    as of the Copernicus DEM, are interpolated across the seam). Unlike the
    nearest-cell lookup of `_lookup_elevations`, elevations vary
    continuously, so finite differences over distances near the cell size
    give smooth slopes.

    Args:
        tiles (list[tuple[numpy.ndarray, rasterio.Affine, float | None]]):
            Windowed elevation arrays, as returned by `_read_windowed_tile`.
        longitudes (numpy.ndarray): Longitudes (decimal degrees) of the
            sample points.
        latitudes (numpy.ndarray): Latitudes (decimal degrees) of the
            sample points, with the same shape as `longitudes`.

    Returns:
        numpy.ndarray: The interpolated elevations (meters), with the same
        shape as `longitudes`; `numpy.nan` where no tile surrounds a point
        with valid data.
    """
    longitudes, latitudes = np.asarray(longitudes), np.asarray(latitudes)
    elevations = np.full(longitudes.shape, np.nan)
    for array, transform, _ in tiles:
        cols, rows = _to_pixel_coordinates(transform, longitudes, latitudes)
        # fractional indices relative to the cell centers
        cols, rows = cols - 0.5, rows - 0.5
        col0, row0 = np.floor(cols), np.floor(rows)
        candidates = (
            np.isnan(elevations)
            & (row0 >= -1)
            & (row0 < array.shape[0])
            & (col0 >= -1)
            & (col0 < array.shape[1])
        )
        c, r = col0[candidates], row0[candidates]
        # the elevations at the surrounding cell centers, from any tile
        corners = [
            _lookup_elevations(
                tiles,
                transform.a * (c + dc + 0.5)
                + transform.b * (r + dr + 0.5)
                + transform.c,
                transform.d * (c + dc + 0.5)
                + transform.e * (r + dr + 0.5)
                + transform.f,
            )
            for dr, dc in [(0, 0), (0, 1), (1, 0), (1, 1)]
        ]
        fc, fr = cols[candidates] - c, rows[candidates] - r
        elevations[candidates] = (1 - fr) * (
            (1 - fc) * corners[0] + fc * corners[1]
        ) + fr * ((1 - fc) * corners[2] + fc * corners[3])
    return elevations


def compute_terrain_slope_aspect(
    dem_paths: str | list[str],
    longitudes: float | np.ndarray,
    latitudes: float | np.ndarray,
    baseline: float = 30,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes the terrain slope and aspect at specified longitudes/latitudes
    from a digital elevation model (DEM), by central differences of the
    elevations (interpolated bilinearly) a `baseline` distance east, west,
    north, and south of each point. The baseline sets the scale of the
    terrain the slope represents: the DEM's cell size (30 meters for
    Copernicus DEM GLO-30, the default) for the steepest local slopes, or a
    larger distance, such as an instrument's ground sample size, for the
    slope averaged over its samples.

    Args:
        dem_paths (str | list[str]): One or more DEM raster sources (see
            `compute_terrain_mask`), covering the points and their baselines.
        longitudes (float | numpy.ndarray): Longitudes (decimal degrees) of
            the points.
        latitudes (float | numpy.ndarray): Latitudes (decimal degrees) of
            the points, with the same shape as `longitudes`.
        baseline (float): The distance (meters) from each point at which
            elevations are differenced.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: The slope (degrees from
        horizontal) and aspect (the azimuth of the downslope direction, in
        degrees clockwise from north, or 0 on flat terrain) at each point,
        with the shape of `longitudes`; `numpy.nan` where the DEM sources do
        not cover a point with valid data.

    Raises:
        ValueError: If `baseline` is not positive, or if none of the
            provided DEM sources overlap the points.
    """
    if baseline <= 0:
        raise ValueError("baseline must be positive.")
    longitudes, latitudes = np.broadcast_arrays(
        np.asarray(longitudes, dtype=float), np.asarray(latitudes, dtype=float)
    )
    # sample points east, west, north, and south of each point
    shape = (4,) + longitudes.shape
    azimuths = np.array([90.0, 270.0, 0.0, 180.0]).reshape(
        (4,) + (1,) * longitudes.ndim
    )
    sample_longitudes, sample_latitudes = (
        np.asarray(values)
        for values in geodesic_destination(
            np.broadcast_to(longitudes, shape),
            np.broadcast_to(latitudes, shape),
            np.broadcast_to(azimuths, shape),
            np.full(shape, float(baseline)),
        )
    )
    paths = [dem_paths] if isinstance(dem_paths, str) else list(dem_paths)
    tiles = [
        tile
        for path in paths
        if (
            tile := _read_windowed_tile(
                path,
                float(np.min(sample_longitudes)),
                float(np.min(sample_latitudes)),
                float(np.max(sample_longitudes)),
                float(np.max(sample_latitudes)),
                padding=2,
            )
        )
        is not None
    ]
    if not tiles:
        raise ValueError("None of the provided DEM sources overlap the points.")
    east, west, north, south = _interpolate_elevations(
        tiles, sample_longitudes, sample_latitudes
    )
    # the elevation gradient (rise per meter) toward the east and north
    gradient_east = (east - west) / (2 * baseline)
    gradient_north = (north - south) / (2 * baseline)
    slope = np.degrees(np.arctan(np.hypot(gradient_east, gradient_north)))
    # the downslope direction is opposite the gradient
    aspect = np.where(
        (gradient_east == 0) & (gradient_north == 0),
        0.0,
        np.degrees(np.arctan2(-gradient_east, -gradient_north)) % 360,
    )
    return slope, aspect


def add_local_incidence_angles(
    observations: gpd.GeoDataFrame,
    dem_paths: str | list[str],
    baseline: float = 30,
) -> gpd.GeoDataFrame:
    """
    Adds the terrain slope and aspect (see `compute_terrain_slope_aspect`)
    and the local incidence angle (see
    `tatc.utils.compute_local_incidence_angle`) at each target of a data
    frame of observations with point geometries and satellite angles, such
    as the ground pixels of `tatc.analysis.collect_ground_pixels` (with
    `sat_altaz=True`) or the observations of
    `tatc.analysis.collect_observations`. Observations may then be
    constrained by their local incidence angle, for example, for
    side-looking instruments over mountainous terrain (an angle above 90
    degrees means the target is hidden behind its own slope).

    Args:
        observations (geopandas.GeoDataFrame): The observations, with point
            geometries (longitude and latitude, in degrees) and the
            satellite's elevation angle and azimuth (`sat_alt` and
            `sat_az` columns, in degrees) at each.
        dem_paths (str | list[str]): One or more DEM raster sources (see
            `compute_terrain_mask`), covering the points.
        baseline (float): See `compute_terrain_slope_aspect`.

    Returns:
        geopandas.GeoDataFrame: A copy of the observations, with added
        columns `terrain_slope`, `terrain_aspect`, and
        `local_incidence_angle` (degrees; `numpy.nan` where the DEM sources
        do not cover a target with valid data).

    Raises:
        ValueError: If the observations lack the `sat_alt` or `sat_az`
            columns, or have geometries other than points.
    """
    missing = [c for c in ["sat_alt", "sat_az"] if c not in observations.columns]
    if missing:
        raise ValueError(
            f"observations lack the column(s) {missing}; collect ground pixels "
            "with sat_altaz=True."
        )
    result = observations.copy()
    columns = ["terrain_slope", "terrain_aspect", "local_incidence_angle"]
    if result.empty:
        for column in columns:
            result[column] = pd.Series(dtype=float)
        return result
    if not (result.geometry.geom_type == "Point").all():
        raise ValueError("observations must have point geometries.")
    slope, aspect = compute_terrain_slope_aspect(
        dem_paths, result.geometry.x.values, result.geometry.y.values, baseline
    )
    result["terrain_slope"] = slope
    result["terrain_aspect"] = aspect
    result["local_incidence_angle"] = compute_local_incidence_angle(
        result.sat_alt.values, result.sat_az.values, slope, aspect
    )
    return result
