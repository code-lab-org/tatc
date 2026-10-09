"""
Methods to generate geospatial points to sample data.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import numpy.typing as npt
import shapely
from shapely.geometry import MultiPolygon, Polygon

from ..constants import EARTH_MEAN_RADIUS
from ..utils.geometry import _hash_geometries, get_planar_bounds
from ..utils.surface import compute_number_samples
from ._grid import clip_to_mask, generate_indices_uniform_spacing

# limits on the rounds and size of random draws to sample points in a mask
_MAX_SAMPLE_ROUNDS = 100
_MAX_SAMPLE_SIZE = 10_000_000


def generate_points_fibonacci_lattice(
    distance: float,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic points following a Fibonacci lattice.

    See: Gonzalez (2010). "Measurement of areas on a sphere using Fibonacci
    and latitude-longitude lattices", Mathematical Geosciences 42(49).
    doi: 10.1007/s11004-009-9257-x

    Note: this implementation differs slightly from Gonzalez. Gonzalez
    requires an odd number of points. This implementation allows any number
    of points. In agreement with Gonzalez, no points are placed at poles.

    Args:
        distance (float): The typical surface distance (meters) between points.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain points
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.

    Returns:
        geopandas.GeoDataFrame: the data frame of generated points, each
            identified by the hash of its geometry (`point_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """

    # determine the number of global samples to achieve average sample distance
    samples = compute_number_samples(distance)
    index = np.arange(samples)
    # compute latitudes, starting from the southern hemisphere and placing
    # neither first nor last points at poles
    latitudes = np.degrees(np.arcsin(2 * (index + 1) / (samples + 2) - 1))
    # compute longitudes on the interval [-180, 180]
    phi = (1 + np.sqrt(5)) / 2  # golden ratio
    longitudes = np.mod(360 * index / phi, 360)
    longitudes = np.where(longitudes > 180, longitudes - 360, longitudes)
    if mask is not None:
        if isinstance(mask, (Polygon, MultiPolygon)):
            if not mask.is_valid:
                raise ValueError("Mask is not a valid Polygon or MultiPolygon.")
            total_bounds = mask.bounds
        else:
            total_bounds = [-180, -90, 180, 90]
        # use the total_bounds to filter relevant points
        min_longitude = total_bounds[0]
        min_latitude = total_bounds[1]
        max_longitude = 180 if total_bounds[2] == -180 else total_bounds[2]
        max_latitude = total_bounds[3]
        # shift longitudes to the interval [0, 360] for masks beyond 180
        if max_longitude > 180:
            longitudes = np.where(longitudes < 0, longitudes + 360, longitudes)
        in_bounds = (
            (min_latitude <= latitudes)
            & (latitudes <= max_latitude)
            & (min_longitude <= longitudes)
            & (longitudes <= max_longitude)
        )
        longitudes = longitudes[in_bounds]
        latitudes = latitudes[in_bounds]
    # create a geodataframe in the WGS84 coordinate reference system (EPSG:4326)
    gdf = gpd.GeoDataFrame(
        geometry=shapely.points(
            longitudes, latitudes, np.full_like(longitudes, elevation)
        ),
        crs="EPSG:4326",
    )
    # clip the geodataframe to the supplied mask, if required
    if mask is not None:
        gdf = clip_to_mask(gdf, mask)
    # identify each point by the hash of its geometry
    gdf.insert(0, "point_id", _hash_geometries(gdf.geometry))
    # return the final geodataframe
    return gdf


def generate_points_uniform_spacing(
    distance: float,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic points at the centroid of regular equally spaced grid
    cells.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    Args:
        distance (float):  The typical surface distance (meters) between points.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain points
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.

    Returns:
        geopandas.GeoDataFrame: the data frame of generated points, each
            identified by the hash of its geometry (`point_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """
    # compute the angular disance of each sample (assuming sphere)
    theta_longitude = np.degrees(distance / EARTH_MEAN_RADIUS)
    theta_latitude = np.degrees(distance / EARTH_MEAN_RADIUS)
    return generate_points_uniform_angular_distance(
        theta_longitude, theta_latitude, elevation, mask
    )


def generate_points_uniform_angular_distance(
    theta_longitude: float,
    theta_latitude: float,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic cells following regular equally spaced grid.

    See: Putman and Lin (2007). "Finite-volume transport on various
    cubed-sphere grids", Journal of Computational Physics, 227(1).
    doi: 10.1016/j.jcp.2007.07.022

    Args:
        theta_longitude (float): The angular difference in longitude (degrees)
            between points.
        theta_latitude (float): The angular difference in latitude (degrees)
            between points.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain points
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.

    Returns:
        geopandas.GeoDataFrame: the data frame of generated points, each
            identified by the hash of its geometry (`point_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """

    # generate grid cells over the filtered region
    indices = generate_indices_uniform_spacing(
        theta_longitude,
        theta_latitude,
        mask,
    )
    longitudes = -180 + (indices[:, 0] + 0.5) * theta_longitude
    latitudes = -90 + (indices[:, 1] + 0.5) * theta_latitude
    # create a geodataframe in the WGS84 reference frame
    gdf = gpd.GeoDataFrame(
        geometry=shapely.points(
            longitudes, latitudes, np.full_like(longitudes, elevation)
        ),
        crs="EPSG:4326",
    )
    # clip the geodataframe to the supplied mask, if required
    if mask is not None:
        gdf = clip_to_mask(gdf, mask)
    # identify each point by the hash of its geometry
    gdf.insert(0, "point_id", _hash_geometries(gdf.geometry))
    # return the final geodataframe
    return gdf


def _get_weighted_cells(
    weights: npt.ArrayLike,
    weights_bounds: tuple[float, float, float, float],
    density: bool,
    bounds: tuple[float, float, float, float],
) -> tuple[npt.NDArray[np.float64], ...]:
    """
    Gets the cells of a weights raster with positive weight that overlap a
    bounding box, and the (unnormalized) probability of sampling each.

    Args:
        weights (numpy.typing.ArrayLike): The two-dimensional weights raster,
            with rows from north to south and columns from west to east.
        weights_bounds (tuple[float, float, float, float]): The bounds (west,
            south, east, north) of the raster in degrees.
        density (bool): True, if weights are per unit area.
        bounds (tuple[float, float, float, float]): The bounding box (min
            longitude, min latitude, max longitude, max latitude) in degrees.

    Returns:
        tuple[numpy.ndarray, ...]: The west, south, east, and north bounds
            (degrees) and the probability of each cell.
    """
    values = np.asarray(weights, dtype=np.float64)
    if values.ndim != 2:
        raise ValueError("Weights must be a two-dimensional array.")
    if np.isinf(values).any():
        raise ValueError("Weights must be finite.")
    values = np.nan_to_num(values, nan=0.0)
    if (values < 0).any():
        raise ValueError("Weights must not be negative.")
    west, south, east, north = weights_bounds
    if not (west < east and -90 <= south < north <= 90):
        raise ValueError(f"Invalid weights bounds: {weights_bounds}.")
    min_longitude, min_latitude, max_longitude, max_latitude = bounds
    longitudes = np.linspace(west, east, values.shape[1] + 1)
    latitudes = np.linspace(north, south, values.shape[0] + 1)
    # shift each column by a whole turn, if required, to overlap the bounds
    shift = np.full(values.shape[1], np.nan)
    for turn in (0, 360, -360):
        overlap = (longitudes[:-1] + turn < max_longitude) & (
            longitudes[1:] + turn > min_longitude
        )
        shift[np.isnan(shift) & overlap] = turn
    columns = np.flatnonzero(~np.isnan(shift))
    rows = np.flatnonzero(
        (latitudes[1:] < max_latitude) & (latitudes[:-1] > min_latitude)
    )
    values = values[np.ix_(rows, columns)]
    # identify cells with positive weight
    row, column = np.nonzero(values)
    cell_west = longitudes[columns[column]] + shift[columns[column]]
    cell_east = longitudes[columns[column] + 1] + shift[columns[column]]
    cell_south = latitudes[rows[row] + 1]
    cell_north = latitudes[rows[row]]
    probability = values[row, column]
    if density:
        # scale by the cell area (on a unit sphere)
        probability = probability * (
            np.radians(cell_east - cell_west)
            * (np.sin(np.radians(cell_north)) - np.sin(np.radians(cell_south)))
        )
    return cell_west, cell_south, cell_east, cell_north, probability


def generate_points_random(
    count: int,
    elevation: float = 0,
    mask: Polygon | MultiPolygon | None = None,
    weights: npt.ArrayLike | None = None,
    weights_bounds: tuple[float, float, float, float] = (-180, -90, 180, 90),
    density: bool = False,
    seed: int | np.random.Generator | None = None,
) -> gpd.GeoDataFrame:
    """
    Generates geodetic points at random, optionally weighted by a raster
    (for example, of population). Without weights, points are distributed
    uniformly by surface area (assuming a sphere). With weights, a raster
    cell is chosen with probability proportional to its weight (or, for
    `density`, its weight times its area) and a point is placed uniformly by
    area within it. Points outside the mask are discarded and redrawn, so
    the distribution within the mask is unchanged.

    Weights can be read from a GeoTIFF with
    `tatc.preprocess.read_raster_weights`.

    Args:
        count (int): The number of points to generate.
        elevation (float): The elevation (meters) above the datum in the WGS 84
            coordinate system.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | None):  An optional mask to constrain points
            using WGS84 (EPSG:4326) geodetic coordinates in a Polygon or MultiPolygon.
        weights (numpy.typing.ArrayLike | None): An optional two-dimensional
            raster of non-negative weights in WGS84 (EPSG:4326) geodetic
            coordinates, with rows from north to south and columns from west
            to east (as in a GeoTIFF). Not-a-number values have zero weight.
        weights_bounds (tuple[float, float, float, float]): The bounds (west,
            south, east, north) of the weights raster in degrees.
        density (bool): True, if weights are per unit area (for example,
            people per square kilometer) rather than totals per cell (for
            example, people).
        seed (int | numpy.random.Generator | None): An optional seed (or
            generator) for reproducible random numbers.

    Returns:
        geopandas.GeoDataFrame: the data frame of generated points, each
            identified by the hash of its geometry (`point_id`, see
            `tatc.utils.geometry.hash_geometry`)
    """
    if count < 0:
        raise ValueError(f"Count must not be negative, got {count}.")
    rng = np.random.default_rng(seed)
    bounds = get_planar_bounds(mask)
    if weights is None:
        # sample the bounding box of the mask as a single cell
        cell_west, cell_south, cell_east, cell_north = (np.array([b]) for b in bounds)
        probability = np.ones(1)
    else:
        cell_west, cell_south, cell_east, cell_north, probability = _get_weighted_cells(
            weights, weights_bounds, density, bounds
        )
    if mask is not None:
        shapely.prepare(mask)
        if weights is not None:
            # discard cells outside the mask
            inside = shapely.intersects(
                mask, shapely.box(cell_west, cell_south, cell_east, cell_north)
            )
            cell_west, cell_south, cell_east, cell_north, probability = (
                cell_west[inside],
                cell_south[inside],
                cell_east[inside],
                cell_north[inside],
                probability[inside],
            )
    cumulative = np.cumsum(probability)
    if len(cumulative) == 0 or cumulative[-1] <= 0:
        raise ValueError("Weights must be positive somewhere within the mask.")
    sin_south = np.sin(np.radians(cell_south))
    sin_north = np.sin(np.radians(cell_north))
    longitudes = []
    latitudes = []
    remaining = count
    drawn = accepted = 0
    for _ in range(_MAX_SAMPLE_ROUNDS):
        if remaining == 0:
            break
        # draw enough points to replace those expected outside the mask
        size = min(
            int(np.ceil(1.1 * remaining * (drawn + 1) / (accepted + 1))) + 16,
            _MAX_SAMPLE_SIZE,
        )
        cell = np.minimum(
            np.searchsorted(cumulative, rng.random(size) * cumulative[-1], "right"),
            len(cumulative) - 1,
        )
        longitude = cell_west[cell] + rng.random(size) * (
            cell_east[cell] - cell_west[cell]
        )
        # distribute latitudes uniformly by area
        latitude = np.degrees(
            np.arcsin(
                sin_south[cell] + rng.random(size) * (sin_north[cell] - sin_south[cell])
            )
        )
        if mask is not None:
            keep = shapely.intersects_xy(mask, longitude, latitude)
            longitude = longitude[keep]
            latitude = latitude[keep]
        drawn += size
        accepted += len(longitude)
        longitudes.append(longitude[:remaining])
        latitudes.append(latitude[:remaining])
        remaining -= len(longitudes[-1])
    if remaining > 0:
        raise ValueError("Could not sample points within the mask.")
    longitudes = np.concatenate(longitudes) if longitudes else np.empty(0)
    latitudes = np.concatenate(latitudes) if latitudes else np.empty(0)
    if mask is None:
        # wrap longitudes of cells beyond 180 degrees to the interval [-180, 180]
        longitudes = np.where(longitudes > 180, longitudes - 360, longitudes)
        longitudes = np.where(longitudes < -180, longitudes + 360, longitudes)
    # create a geodataframe in the WGS84 reference frame
    gdf = gpd.GeoDataFrame(
        geometry=shapely.points(
            longitudes, latitudes, np.full_like(longitudes, elevation)
        ),
        crs="EPSG:4326",
    )
    # identify each point by the hash of its geometry
    gdf.insert(0, "point_id", _hash_geometries(gdf.geometry))
    # return the final geodataframe
    return gdf
