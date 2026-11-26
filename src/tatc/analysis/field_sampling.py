"""
Methods to sample fields (such as the variables of a nature run in an
observing system simulation experiment, OSSE) at the times and places of
observations: at their points, over the footprints of their samples, inside
their regions, at levels or by weights of levels, and as profiles at heights,
on regular, curvilinear, and unstructured grids.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from .. import constants
from ..utils._grids import FieldGrid
from .check import _get_observation_times


def _get_points(geometries: gpd.GeoSeries) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the longitudes and latitudes (degrees) of geometries: points, or the
    representative points of other geometries.
    """
    points = geometries.representative_point()
    return points.x.to_numpy(), points.y.to_numpy()


def _get_dimension(field: Any, names: Sequence[str], given: str | None) -> str:
    """Gets the name of a dimension of a field."""
    if given is not None:
        return given
    for name in names:
        if name in field.dims:
            return name
    raise ValueError(f"the field has no dimension named any of {list(names)}")


def _field_times(field: Any, dimension: str, times: pd.Series) -> np.ndarray:
    """
    Converts times to those of the time coordinate of a field: naive (UTC) if
    the field's are.
    """
    times = pd.to_datetime(times, utc=True)
    if getattr(field[dimension].dtype, "tz", None) is None:
        times = times.dt.tz_convert("UTC").dt.tz_localize(None)
    return times.to_numpy()


def _field_longitudes(field: Any, dimension: str, longitude: np.ndarray) -> np.ndarray:
    """
    Converts longitudes to the convention of a field: [0, 360) if its
    longitudes exceed 180 degrees, else [-180, 180).
    """
    if float(field[dimension].max()) > 180:
        return np.mod(longitude, 360)
    return (longitude + 180) % 360 - 180


def _coordinate_extent(values: np.ndarray) -> tuple[float, float]:
    """
    Gets the extent of a coordinate: half a spacing beyond its outer values
    (the cells of which they are the centers), or unbounded for one value.
    """
    values = np.unique(np.asarray(values, dtype=float))
    if len(values) < 2:
        return -np.inf, np.inf
    return (
        values[0] - (values[1] - values[0]) / 2,
        values[-1] + (values[-1] - values[-2]) / 2,
    )


def _outside_field(
    field: Any,
    lat_dim: str,
    lon_dim: str,
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> np.ndarray:
    """
    Finds the places outside the extent of a field (see `_coordinate_extent`),
    as of a regional nature run, in latitude and, unless the field's
    longitudes go around the globe, in longitude (in the field's convention,
    see `_field_longitudes`).
    """
    low, high = _coordinate_extent(field[lat_dim].to_numpy())
    outside = (latitude < low) | (latitude > high)
    low, high = _coordinate_extent(field[lon_dim].to_numpy())
    if high - low < 360 - 1e-9:
        outside |= (longitude < low) | (longitude > high)
    return outside


def _get_grid(field: Any, dimensions: Mapping[str, str], time_dim: str) -> FieldGrid:
    """
    Gets the horizontal grid of a field (see `tatc.utils._grids.FieldGrid`),
    by the names of its latitudes and longitudes (`lat` and `lon` of the
    dimensions, if given).
    """
    return FieldGrid.from_field(
        field, dimensions.get("lat"), dimensions.get("lon"), time_dim
    )


def _check_method(grid: FieldGrid, method: str) -> None:
    """
    Checks a method of sampling a field on a grid: linear interpolation
    requires a regular grid.
    """
    if method == "linear" and grid.kind != "regular":
        raise ValueError(
            "linear interpolation requires a field on a regular grid (with "
            + f"latitude and longitude dimensions), not a {grid.kind} grid of "
            + f"dimensions {grid.dims}: use method 'nearest'"
        )


def _to_nanoseconds(times: Any) -> np.ndarray:
    """Converts times (naive, UTC) to nanoseconds since 1970."""
    return np.asarray(times).astype("datetime64[ns]").astype(np.int64)


def _outside_period(field: Any, time_dim: str, times: np.ndarray) -> np.ndarray:
    """
    Finds the times outside the period of a field: by more than half a time
    step before its first time or after its last (see `_coordinate_extent`,
    as for the extent of its grid), or none for a field of one time.
    """
    low, high = _coordinate_extent(_to_nanoseconds(field[time_dim].to_numpy()))
    times = _to_nanoseconds(times)
    return (times < low) | (times > high)


def _nearest_times(field: Any, time_dim: str, times: np.ndarray) -> np.ndarray:
    """
    Gets the index of the nearest time of a field to each time (the later of
    two as near, as `xarray.DataArray.sel` does).
    """
    coordinate = _to_nanoseconds(field[time_dim].to_numpy())
    order = np.argsort(coordinate, kind="stable")
    ordered = coordinate[order]
    if len(ordered) == 1:
        return np.zeros(len(times), dtype=np.intp)
    times = _to_nanoseconds(times)
    after = np.clip(np.searchsorted(ordered, times), 1, len(ordered) - 1)
    nearest = np.where(
        times - ordered[after - 1] < ordered[after] - times, after - 1, after
    )
    return order[nearest]


def _select_places(
    field: Any,
    longitude: np.ndarray,
    latitude: np.ndarray,
    method: str,
    dimensions: Mapping[str, str],
    time_dim: str,
    dim: str,
) -> tuple[Any, np.ndarray]:
    """
    Selects a field at places, along a new dimension, keeping its other
    dimensions (such as time): at the nearest cell of its grid, or (of a
    regular grid) by linear interpolation; and whether each place is outside
    the extent of the field (for the nearest cells, see `_outside_field` and
    `tatc.utils._grids.FieldGrid.nearest`).

    Args:
        field (xarray.DataArray): The field.
        longitude (numpy.ndarray): The longitudes (degrees) of the places.
        latitude (numpy.ndarray): The latitudes (degrees) of the places.
        method (str): `"nearest"` or `"linear"`.
        dimensions (Mapping[str, str]): The names of the field's dimensions.
        time_dim (str): The name of the field's time dimension.
        dim (str): The name of the dimension of the places.

    Returns:
        tuple[xarray.DataArray, numpy.ndarray]: the field at the places, and
            whether each is outside its extent
    """
    # pylint: disable=import-outside-toplevel
    import xarray as xr

    grid = _get_grid(field, dimensions, time_dim)
    _check_method(grid, method)
    if grid.kind == "regular":
        lat_dim, lon_dim = grid.dims
        where = {
            lat_dim: xr.DataArray(latitude, dims=dim),
            lon_dim: xr.DataArray(
                _field_longitudes(field, lon_dim, longitude), dims=dim
            ),
        }
        if method == "linear":
            return field.interp(where), np.zeros(len(latitude), dtype=bool)
        outside = _outside_field(
            field, lat_dim, lon_dim, latitude, where[lon_dim].to_numpy()
        )
        return field.sel(where, method="nearest"), outside
    index, outside = grid.nearest(longitude, latitude)
    return (
        field.isel(
            {
                name: xr.DataArray(position, dims=dim)
                for name, position in zip(grid.dims, grid.unravel(index))
            }
        ),
        outside,
    )


def _sample_values(
    field: Any,
    times: pd.Series,
    longitude: np.ndarray,
    latitude: np.ndarray,
    method: str,
    dimensions: Mapping[str, str] | None,
) -> np.ndarray:
    """
    Samples a field at times and places (see `sample_field`): its nearest
    value (`NaN` outside the period or extent of the field) or (on a regular
    grid) its linear interpolation in time and space.

    Args:
        field (xarray.DataArray): The field.
        times (pandas.Series): The times of the samples.
        longitude (numpy.ndarray): The longitudes (degrees) of the samples.
        latitude (numpy.ndarray): The latitudes (degrees) of the samples.
        method (str): `"nearest"` or `"linear"`.
        dimensions (Mapping[str, str] | None): The names of the field's
            dimensions of `time`, `lat`, and `lon`, if not the defaults.

    Returns:
        numpy.ndarray: the values
    """
    # pylint: disable=import-outside-toplevel
    import xarray as xr

    if method not in ("nearest", "linear"):
        raise ValueError(f"method must be 'nearest' or 'linear', not {method!r}")
    dimensions = dict(dimensions or {})
    time_dim = _get_dimension(field, ["time"], dimensions.get("time"))
    grid = _get_grid(field, dimensions, time_dim)
    _check_method(grid, method)
    sample_times = _field_times(field, time_dim, pd.Series(times))
    if grid.kind == "regular":
        lat_dim, lon_dim = grid.dims
        where = {
            time_dim: xr.DataArray(sample_times, dims="sample"),
            lat_dim: xr.DataArray(latitude, dims="sample"),
            lon_dim: xr.DataArray(
                _field_longitudes(field, lon_dim, longitude), dims="sample"
            ),
        }
        if method == "linear":
            return field.interp(where).transpose("sample", ...).to_numpy().astype(float)
        values = (
            field.sel(where, method="nearest")
            .transpose("sample", ...)
            .to_numpy()
            .astype(float)
        )
        outside = _outside_field(
            field, lat_dim, lon_dim, latitude, where[lon_dim].to_numpy()
        )
    else:
        index, outside = grid.nearest(longitude, latitude)
        indexers = {
            time_dim: xr.DataArray(
                _nearest_times(field, time_dim, sample_times), dims="sample"
            )
        }
        for name, position in zip(grid.dims, grid.unravel(index)):
            indexers[name] = xr.DataArray(position, dims="sample")
        values = field.isel(indexers).transpose("sample", ...).to_numpy().astype(float)
    values[outside | _outside_period(field, time_dim, sample_times)] = np.nan
    return values


_METHODS = ("nearest", "linear", "footprint", "area")
"""The methods of sampling a field."""

_FOOTPRINTS = ("gaussian", "uniform")
"""The shapes of footprints over which to weight a field."""

_GAUSSIAN_REACH = 1.5
"""The distance (in the larger of a Gaussian footprint's full widths at half
maximum) within which cells are weighted: beyond it, weights are below about
0.2% of the center's."""

_SIZE_COLUMNS = ("sample_along", "sample_cross", "sample_azimuth")
"""The columns of the sizes and orientation of samples (see
`compute_sample_sizes`) by which footprints are weighted."""

_CHUNK = 10_000_000
"""The most cells to weight at once (over the samples of a chunk)."""


def _get_slice(field: Any, grid: FieldGrid, time_dim: str, index: int) -> np.ndarray:
    """
    Gets the values of a field at a time (by index), flattened in the order of
    the horizontal dimensions of its grid.
    """
    values = field.isel({time_dim: int(index)}).transpose(*grid.dims)
    return np.ravel(values.to_numpy()).astype(float)


def _check_dimensions(field: Any, grid: FieldGrid, time_dim: str) -> None:
    """
    Checks that a field has no dimensions besides time and those of its
    grid, as needed to weight its cells.
    """
    other = [dim for dim in field.dims if dim != time_dim and dim not in grid.dims]
    if other:
        raise ValueError(
            "the field must have no dimensions besides time and those of its "
            + f"grid {grid.dims} to weight its cells, not {other}: select them first"
        )


def _get_time_steps(
    field: Any, time_dim: str, times: np.ndarray, time_method: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Gets the time steps of a field at which to sample each time: the nearest
    step (`time_method="nearest"`), or the steps before and after it and the
    fraction of the interval between them (`"linear"`), and whether each time
    is outside the period of the field (by more than half a time step for
    the nearest step; at all for linear interpolation).

    Returns:
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: the first step
            (by index), the fraction toward the next step (zero for the
            nearest step), and whether each time is outside
    """
    if time_method == "nearest":
        return (
            _nearest_times(field, time_dim, times),
            np.zeros(len(times)),
            _outside_period(field, time_dim, times),
        )
    coordinate = _to_nanoseconds(field[time_dim].to_numpy())
    order = np.argsort(coordinate, kind="stable")
    ordered = coordinate[order]
    times = _to_nanoseconds(times)
    outside = (times < ordered[0]) | (times > ordered[-1])
    if len(ordered) == 1:
        return np.zeros(len(times), dtype=np.intp), np.zeros(len(times)), outside
    before = np.clip(
        np.searchsorted(ordered, times, side="right") - 1, 0, len(ordered) - 2
    )
    fraction = (times - ordered[before]) / (ordered[before + 1] - ordered[before])
    fraction = np.clip(fraction, 0, 1)
    # the next step, if any weight falls on it (as is not the case at the last time)
    return order[before], fraction, outside


def _get_footprints(
    observations: gpd.GeoDataFrame, footprint: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Gets the footprint of each observation's sample (see
    `compute_sample_sizes`): its along-track and cross-track sizes (meters,
    full widths at half maximum of a Gaussian footprint, or the sides of a
    uniform one), the azimuth of its along-track axis, the distance (meters)
    within which cells are weighted, and its area (square meters, the
    integral of its weights).
    """
    missing = [column for column in _SIZE_COLUMNS if column not in observations]
    if missing:
        raise ValueError(
            f"observations need the columns {missing} to weight a field over "
            + "their footprints: see tatc.analysis.compute_sample_sizes"
        )
    along = observations["sample_along"].to_numpy(dtype=float)
    cross = observations["sample_cross"].to_numpy(dtype=float)
    azimuth = observations["sample_azimuth"].to_numpy(dtype=float)
    if footprint == "gaussian":
        reach = _GAUSSIAN_REACH * np.fmax(along, cross)
        area = np.pi * along * cross / (4 * np.log(2))
    else:
        reach = np.hypot(along, cross) / 2
        area = along * cross
    return along, cross, azimuth, reach, area


def _local_offsets(
    longitude: np.ndarray,
    latitude: np.ndarray,
    cell_longitude: np.ndarray,
    cell_latitude: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the offsets (meters) east and north of cells from places, in the
    plane tangent to a sphere at each place (from Earth-centered unit
    vectors, so that they are continuous across the antimeridian).
    """
    lon, lat = np.radians(longitude)[:, None], np.radians(latitude)[:, None]
    cell_lon, cell_lat = np.radians(cell_longitude), np.radians(cell_latitude)
    x = np.cos(cell_lat) * np.cos(cell_lon)
    y = np.cos(cell_lat) * np.sin(cell_lon)
    z = np.sin(cell_lat)
    east = -np.sin(lon) * x + np.cos(lon) * y
    north = -np.sin(lat) * np.cos(lon) * x - np.sin(lat) * np.sin(lon) * y
    north = north + np.cos(lat) * z
    facing = np.cos(lat) * np.cos(lon) * x + np.cos(lat) * np.sin(lon) * y
    facing = facing + np.sin(lat) * z
    radius = constants.EARTH_MEAN_RADIUS
    # cells on the far side of the sphere are not near
    east = np.where(facing > 0, east * radius, np.nan)
    north = np.where(facing > 0, north * radius, np.nan)
    return east, north


class _Weights:
    """
    The weights of the cells of a grid near samples: over the footprint of
    each sample (`method="footprint"`), or inside the region of each
    observation (`method="area"`), with the area of each sample (the total of
    its weights, if the field had values everywhere), from which to find the
    share of a sample's weight on cells with values.
    """

    def __init__(
        self,
        observations: gpd.GeoDataFrame,
        grid: FieldGrid,
        areas: np.ndarray,
        method: str,
        footprint: str,
    ):
        self.grid = grid
        self.areas = areas
        self.method = method
        self.footprint = footprint
        centroid = shapely.centroid(observations.geometry.to_numpy())
        self.longitude = shapely.get_x(centroid)
        self.latitude = shapely.get_y(centroid)
        if method == "footprint":
            (
                self.along,
                self.cross,
                self.azimuth,
                self.reach,
                self.area,
            ) = _get_footprints(observations, footprint)
        else:
            self.geometry = observations.geometry.to_numpy()
            self.reach, self.area = _get_regions(self.geometry, centroid)
        self.nearest, self.outside = grid.nearest(self.longitude, self.latitude)

    def is_weighted(self) -> np.ndarray:
        """
        Gets whether each sample is weighted over cells: whether its
        footprint (or region) is known, and larger than its nearest cell
        (else it takes the value of its nearest cell).
        """
        with np.errstate(invalid="ignore"):
            return (
                np.isfinite(self.reach)
                & (self.reach > 0)
                & np.isfinite(self.area)
                & (self.area > self.areas[self.nearest])
            )

    def get(self, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Gets the cells (flat indices, `-1` for none) near samples (by
        position) and their weights (with their areas), each with shape
        (samples, cells).
        """
        index = self.grid.candidates(
            self.longitude[rows], self.latitude[rows], self.reach[rows]
        )
        cell_longitude, cell_latitude = self.grid.positions(index)
        if self.method == "footprint":
            east, north = _local_offsets(
                self.longitude[rows], self.latitude[rows], cell_longitude, cell_latitude
            )
            azimuth = np.radians(self.azimuth[rows])[:, None]
            along = east * np.sin(azimuth) + north * np.cos(azimuth)
            cross = east * np.cos(azimuth) - north * np.sin(azimuth)
            size_along = self.along[rows][:, None]
            size_cross = self.cross[rows][:, None]
            if self.footprint == "gaussian":
                weight = np.exp(
                    -4
                    * np.log(2)
                    * ((along / size_along) ** 2 + (cross / size_cross) ** 2)
                )
                weight = np.where(
                    np.hypot(along, cross) <= self.reach[rows][:, None], weight, 0
                )
            else:
                weight = (
                    (np.abs(along) <= size_along / 2)
                    & (np.abs(cross) <= size_cross / 2)
                ).astype(float)
        else:
            weight = shapely.contains_xy(
                self.geometry[rows][:, None],
                (cell_longitude + 180) % 360 - 180,
                cell_latitude,
            ).astype(float)
        weight = np.where(index >= 0, np.nan_to_num(weight), 0)
        return index, weight * self.areas[np.maximum(index, 0)]


def _get_regions(
    geometry: np.ndarray, centroid: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the distance (meters) from the centroid of each region within which
    to find its cells (to the farthest corner of its bounds), and its
    geodesic area (square meters, zero for points and lines).
    """
    # pylint: disable=import-outside-toplevel
    from pyproj import Geod

    geod = Geod(ellps="WGS84")
    bounds = shapely.bounds(geometry)
    longitude, latitude = shapely.get_x(centroid), shapely.get_y(centroid)
    reach = np.zeros(len(geometry))
    for corner in ((0, 1), (0, 3), (2, 1), (2, 3)):
        _, _, distance = geod.inv(
            longitude, latitude, bounds[:, corner[0]], bounds[:, corner[1]]
        )
        reach = np.fmax(reach, np.asarray(distance))
    area = np.array([abs(geod.geometry_area_perimeter(shape)[0]) for shape in geometry])
    return reach, area


def _sample_weighted(
    observations: gpd.GeoDataFrame,
    field: Any,
    times: pd.Series,
    method: str,
    footprint: str,
    time_method: str,
    min_valid_fraction: float,
    dimensions: Mapping[str, str],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Samples a field over the footprint or region of each observation (see
    `sample_field`): the mean of its cells, weighted by the footprint (or
    inside the region) and by their areas, at the nearest time step or
    interpolated linearly between two, and the share of the footprint's
    weight on cells with values (`NaN` below `min_valid_fraction`).

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: the values, and the share of
            each sample's weight on cells with values
    """
    time_dim = _get_dimension(field, ["time"], dimensions.get("time"))
    grid = _get_grid(field, dimensions, time_dim)
    _check_dimensions(field, grid, time_dim)
    weights = _Weights(observations, grid, grid.cell_areas(field), method, footprint)
    weighted = weights.is_weighted()
    step, fraction, outside = _get_time_steps(
        field, time_dim, _field_times(field, time_dim, times), time_method
    )
    shares = (1 - fraction, fraction)
    count = len(observations)
    total = np.zeros(count)
    weight_sum = np.zeros(count)
    nearest_value = np.full(count, np.nan)
    for first in np.unique(step[~outside]):
        group = np.flatnonzero((step == first) & ~outside)
        slices = [_get_slice(field, grid, time_dim, first)]
        if np.any(fraction[group] > 0):
            slices.append(_get_slice(field, grid, time_dim, first + 1))
        nearest = np.zeros(len(group))
        for values, share in zip(slices, shares):
            nearest += np.where(
                share[group] > 0, share[group] * values[weights.nearest[group]], 0
            )
        nearest_value[group] = nearest
        rows = group[weighted[group]]
        rows = rows[np.argsort(weights.reach[rows], kind="stable")]
        for chunk, index, weight in _chunked(weights, rows):
            combined = np.zeros(index.shape)
            known = np.ones(index.shape, dtype=bool)
            for values, share in zip(slices, shares):
                part = share[chunk][:, None]
                cell = np.where(index >= 0, values[np.maximum(index, 0)], np.nan)
                known &= np.isfinite(cell) | (part == 0)
                combined += np.where(part > 0, part * np.nan_to_num(cell), 0)
            valid = known & (weight > 0)
            total[chunk] = np.sum(np.where(valid, weight * combined, 0), axis=1)
            weight_sum[chunk] = np.sum(np.where(valid, weight, 0), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        value = np.where(weighted, total / weight_sum, nearest_value)
        share = np.where(
            weighted,
            np.minimum(weight_sum / weights.area, 1),
            np.isfinite(nearest_value).astype(float),
        )
    # samples at their nearest cells outside the grid have no value
    beyond = ~weighted & weights.outside
    value[beyond] = np.nan
    share[beyond] = 0.0
    value[outside] = np.nan
    share[outside] = 0.0
    if method == "footprint":
        # footprints without sizes: of an unknown share
        share[~np.isfinite(weights.reach) & ~outside & ~beyond] = np.nan
    with np.errstate(invalid="ignore"):
        value[(share < min_valid_fraction) | ~np.isfinite(value)] = np.nan
    return value, share


def _chunked(weights: _Weights, rows: np.ndarray):
    """
    Gets the cells near samples (see `_Weights.get`) in chunks of samples,
    each of at most about `_CHUNK` cells (halving a chunk of more), as the
    samples, their cells, and the cells' weights.
    """
    size = 1024
    start = 0
    while start < len(rows):
        chunk = rows[start : start + size]
        index, weight = weights.get(chunk)
        if index.size > _CHUNK and len(chunk) > 1:
            size = max(1, len(chunk) // 2)
            continue
        yield chunk, index, weight
        start += len(chunk)


def _get_vertical(
    field: Any, grid: FieldGrid, time_dim: str, dimensions: Mapping[str, str]
) -> str | None:
    """
    Gets the vertical dimension of a field: its dimension besides time and
    those of its grid (named `level` of the dimensions, if given), or `None`
    if it has none.
    """
    other = [dim for dim in field.dims if dim != time_dim and dim not in grid.dims]
    given = dimensions.get("level")
    if given is not None:
        if given not in other:
            raise ValueError(f"the field has no vertical dimension named {given!r}")
        other.remove(given)
    if other and (given is not None or len(other) > 1):
        raise ValueError(
            "the field must have no dimensions besides time, those of its grid "
            + f"{grid.dims}, and a vertical dimension, not {other}: select them first"
        )
    return given or (str(other[0]) if other else None)


def _select_level(field: Any, vertical: str, level: float) -> Any:
    """
    Selects a level of a field: the values at a value of its vertical
    coordinate, interpolated linearly between the levels around it (`NaN`
    beyond its levels).
    """
    coordinate = field[vertical].to_numpy()
    if np.any(coordinate == level):
        return field.sel({vertical: level})
    return field.interp({vertical: level})


def _get_weights(weights: Any, field: Any, vertical: str, count: int) -> np.ndarray:
    """
    Gets the weights of a field's levels: along its vertical dimension (by
    its coordinate, if the weights have one), for every observation or along
    a dimension `observation` (in the order of the observations), with shape
    (levels,) or (observations, levels).
    """
    other = set(weights.dims) - {vertical, "observation"}
    if vertical not in weights.dims or other:
        raise ValueError(
            f"weights must be over the field's vertical dimension {vertical!r} "
            + f"(and optionally 'observation'), not {tuple(weights.dims)}"
        )
    if vertical in weights.coords:
        weights = weights.sel({vertical: field[vertical].to_numpy()})
    elif weights.sizes[vertical] != field.sizes[vertical]:
        raise ValueError(
            f"weights must have the {field.sizes[vertical]} levels of the field"
        )
    if "observation" in weights.dims:
        if weights.sizes["observation"] != count:
            raise ValueError(
                f"weights must have the {count} observations along 'observation'"
            )
        return weights.transpose("observation", vertical).to_numpy().astype(float)
    return weights.to_numpy().astype(float)


def _weigh_field(field: Any, vertical: str, weights: np.ndarray) -> Any:
    """
    Sums the levels of a field weighted by their weights (see
    `_weigh_levels`), over the levels with values, as a field without its
    vertical dimension.
    """
    # pylint: disable=import-outside-toplevel
    import xarray as xr

    weight = xr.DataArray(weights, dims=vertical)
    total = field.fillna(0).dot(weight, dim=vertical)
    return total.where(field.notnull().any(vertical)).rename(field.name)


def _weigh_levels(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """
    Sums the values of levels (with shape (observations, levels)) weighted by
    their weights, over the levels with values (as of pressure levels below
    the ground), or `NaN` if none has a value.
    """
    weights = np.broadcast_to(weights, values.shape)
    finite = np.isfinite(values)
    total = np.sum(np.where(finite, weights * values, 0), axis=1)
    return np.where(np.any(finite, axis=1), total, np.nan)


def sample_field(
    observations: gpd.GeoDataFrame,
    field: Any,
    name: str | None = None,
    method: str = "nearest",
    time_column: str | None = None,
    dimensions: Mapping[str, str] | None = None,
    footprint: str = "gaussian",
    time_method: str = "nearest",
    min_valid_fraction: float = 0.5,
    level: float | None = None,
    weights: Any = None,
) -> gpd.GeoDataFrame:
    """
    Samples a field (such as a variable of a nature run) at the time and place
    of each observation: at its point (or the representative point of its
    region), by the nearest value or by linear interpolation in time and
    space; or over the footprint of its sample or inside its region, by the
    mean of the field's cells weighted by the footprint or the region and by
    their areas. Observations outside the period of the field (by more than
    half a time step before its first time or after its last), or outside its
    extent, have no value (`NaN`). Requires xarray.

    The field may be on a regular grid, with latitude and longitude
    dimensions; on a curvilinear grid, with two-dimensional latitude and
    longitude coordinates (as of a regional model on a projected grid, such
    as WRF's `XLAT` and `XLONG`, taken at the first time if they vary in
    time); or on an unstructured grid, with latitude and longitude
    coordinates of a dimension of cells (as of an MPAS mesh's `latCell` and
    `lonCell`, in degrees or, if their units say so, radians). A regular
    grid's extent reaches half a spacing beyond its outer latitudes and,
    unless they go around the globe, longitudes. Grids that are not regular
    are sampled at their nearest cells (found on the sphere, so they may
    cross the antimeridian), and their extent reaches about half the diagonal
    of a cell beyond their outer cells; they require scipy, and cannot be
    interpolated linearly.

    With `method="footprint"`, each observation's sample (see
    `compute_sample_sizes`, whose `sample_along`, `sample_cross`, and
    `sample_azimuth` it needs) is centered at its centroid, its along-track
    size lies in the direction of its azimuth on the ground, and its
    cross-track size perpendicular to it. Its footprint is `"gaussian"`
    (default), with the sizes as the full widths at half maximum, as of the
    response of a beam or pixel, weighted within 1.5 times the larger; or
    `"uniform"`, a rectangle of the sizes. With `method="area"`, the cells
    whose centers lie inside each observation's region (such as the viewed
    part of a feature, see `collect_feature_observations`) are weighted
    equally. A footprint or region smaller than the cell nearest its center
    (or a region without cells inside), and a footprint without sizes (as of
    an instrument without samples), take the value of that cell.

    The share of a sample's weight on cells with values (the weights of
    those cells, with their areas, relative to the area of the footprint or
    region, at most 1) is given as `{name}_valid_fraction`: it is less than 1
    where the footprint extends beyond the extent of the field (as of a
    regional nature run) or over cells without values (such as those of a
    field defined only over the ocean). Samples whose share is less than
    `min_valid_fraction` have no value; those without sizes have an unknown
    share (`NaN`). Footprints are found in the plane tangent to the Earth at
    their centers, which suits footprints much smaller than the Earth.

    A field with a vertical dimension (such as the levels of a nature run,
    besides time and those of its grid) is sampled at a `level` (a value of
    its vertical coordinate, such as a pressure, interpolated linearly between
    levels), or as the sum of its levels weighted by `weights` (an
    `xarray.DataArray` over its vertical dimension, as of a column amount's
    layer thicknesses or an instrument's weighting function or averaging
    kernel, which the user gives; levels without values, as below the
    ground, are left out of the sum). Weights over the vertical dimension
    alone are applied to the field before sampling it (by any method);
    weights that vary by observation (along a dimension `observation`, in
    the order of the observations, as of averaging kernels that depend on the
    scene) are applied to the nearest or interpolated values of each
    observation's levels. See `sample_profiles` for the values of each level.

    The time of a footprint or region is that of the nearest time step of the
    field (`time_method="nearest"`, as suits a field of means over intervals
    centered on its times), or its values are interpolated linearly between
    the steps before and after it (`"linear"`, with no value outside the
    field's period).

    Args:
        observations (geopandas.GeoDataFrame): The observations, with the time
            (`time_column`) of each sample (and, for footprints, the sizes and
            azimuth of each sample, see `compute_sample_sizes`).
        field (xarray.DataArray): The field, with a dimension of time and the
            latitudes and longitudes of its grid (by default named `time`,
            and any of `lat`, `latitude`, `XLAT`, `latCell`, or `nav_lat`
            and of `lon`, `longitude`, `XLONG`, `lonCell`, or `nav_lon`),
            with longitudes in [-180, 180) or [0, 360). The areas of its
            cells may be given by a coordinate named `cell_area`, `area`, or
            `areaCell` over the dimensions of its grid.
        name (str | None): The column of the sampled values, or `None` for the
            field's name (or `"value"`).
        method (str): `"nearest"`, `"linear"` (of a regular grid),
            `"footprint"`, or `"area"`.
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end`.
        dimensions (Mapping[str, str] | None): The names of the field's
            `time` dimension and of its latitudes (`lat`) and longitudes
            (`lon`): the dimensions of a regular grid, or the coordinates of
            another grid, if not the defaults; and of its vertical dimension
            (`level`), if it has others.
        footprint (str): The shape of footprints (`method="footprint"`):
            `"gaussian"` or `"uniform"`.
        time_method (str): The times of footprints and regions: `"nearest"`
            or `"linear"`.
        min_valid_fraction (float): The least share of a footprint's or
            region's weight on cells with values for it to have a value.
        level (float | None): The value of the field's vertical coordinate
            at which to sample it, if it has a vertical dimension.
        weights (xarray.DataArray | None): The weights of the field's levels
            (over its vertical dimension, and optionally `observation`), by
            which to sum them, if it has a vertical dimension.

    Returns:
        geopandas.GeoDataFrame: the observations, with the sampled values
            (and, of footprints and regions, the share of their weight on
            cells with values)
    """
    if method not in _METHODS:
        raise ValueError(f"method must be one of {list(_METHODS)}, not {method!r}")
    if footprint not in _FOOTPRINTS:
        raise ValueError(
            f"footprint must be one of {list(_FOOTPRINTS)}, not {footprint!r}"
        )
    if time_method not in ("nearest", "linear"):
        raise ValueError(
            f"time_method must be 'nearest' or 'linear', not {time_method!r}"
        )
    name = name or field.name or "value"
    result = observations.copy()
    weighted = method in ("footprint", "area")
    if len(observations) == 0:
        result[name] = pd.Series(dtype=float)
        if weighted:
            result[f"{name}_valid_fraction"] = pd.Series(dtype=float)
        return result
    dimensions = dict(dimensions or {})
    time_dim = _get_dimension(field, ["time"], dimensions.get("time"))
    vertical = _get_vertical(
        field, _get_grid(field, dimensions, time_dim), time_dim, dimensions
    )
    if vertical is None and (level is not None or weights is not None):
        raise ValueError("a level or weights need a field with a vertical dimension")
    if vertical is not None and level is not None:
        if weights is not None:
            raise ValueError("give a level or weights, not both")
        field, vertical = _select_level(field, vertical, level), None
    if vertical is not None and weights is None:
        raise ValueError(
            f"the field has a vertical dimension {vertical!r}: give a level or "
            + "weights to sample it (or see sample_profiles)"
        )
    level_weights = (
        None
        if vertical is None
        else _get_weights(weights, field, vertical, len(observations))
    )
    if level_weights is not None and level_weights.ndim == 1:
        # weights of the levels alone: the weighted column, sampled by any method
        field = _weigh_field(field, vertical, level_weights)
        level_weights, vertical = None, None
    if level_weights is not None and weighted:
        raise ValueError(
            "weights that vary by observation need method 'nearest' or 'linear'"
        )
    times = _get_observation_times(observations, time_column)
    if weighted:
        result[name], result[f"{name}_valid_fraction"] = _sample_weighted(
            observations,
            field,
            times,
            method,
            footprint,
            time_method,
            min_valid_fraction,
            dimensions,
        )
        return result
    longitude, latitude = _get_points(observations.geometry)
    values = _sample_values(field, times, longitude, latitude, method, dimensions)
    if level_weights is not None:
        values = _weigh_levels(values, level_weights)
    result[name] = values
    return result


def _interpolate_heights(
    heights: np.ndarray, values: np.ndarray, target: np.ndarray
) -> np.ndarray:
    """
    Interpolates values linearly in height at target heights, for each row of
    heights and values (with shape (rows, levels), in any order of height,
    as of levels from the top down), with `NaN` beyond each row's heights.
    """
    heights = np.where(np.isfinite(values), heights, np.nan)
    order = np.argsort(np.where(np.isfinite(heights), heights, np.inf), axis=1)
    heights = np.take_along_axis(heights, order, axis=1)
    values = np.take_along_axis(values, order, axis=1)
    count = np.sum(np.isfinite(heights), axis=1)
    target = target[:, None]
    above = np.sum(heights < target, axis=1)
    # a target at the lowest height is between it and the next
    above = np.where((above == 0) & (heights[:, 0] == target[:, 0]), 1, above)
    valid = (above > 0) & (above < np.maximum(count, 1)) | (
        (count == 1) & (heights[:, 0] == target[:, 0])
    )
    upper = np.clip(above, 1, heights.shape[1] - 1) if heights.shape[1] > 1 else 0
    lower = np.maximum(upper - 1, 0) if heights.shape[1] > 1 else 0
    rows = np.arange(len(heights))
    h_low, h_high = heights[rows, lower], heights[rows, upper]
    v_low, v_high = values[rows, lower], values[rows, upper]
    with np.errstate(invalid="ignore", divide="ignore"):
        share = np.where(h_high > h_low, (target[:, 0] - h_low) / (h_high - h_low), 0)
    result = v_low + share * (v_high - v_low)
    return np.where(valid, result, np.nan)


def sample_profiles(
    observations: gpd.GeoDataFrame,
    field: Any,
    name: str | None = None,
    heights: Any = None,
    method: str = "nearest",
    time_column: str | None = None,
    dimensions: Mapping[str, str] | None = None,
) -> gpd.GeoDataFrame:
    """
    Samples the levels of a field with a vertical dimension (such as a
    nature run on pressure or model levels) at the time and place of each
    observation (see `sample_field`, by the nearest value or by linear
    interpolation in time and space), in long format: a row for each level
    of each observation, with its `observation` (the index of the
    observation, whose other columns it repeats), and its `level` (the value
    of the field's vertical coordinate) or `height`, and the sampled value.

    The points of an observation at heights (the vertices of a point or
    multipoint with z-coordinates, as the tangent points of a limb scan, see
    `collect_limb_observations`, or the positions of an aircraft) are each
    sampled at their height (meters above the WGS 84 ellipsoid), interpolated
    linearly in height between the levels around it (`NaN` beyond them),
    each a row with its `height` (and no `level`). The heights of the levels
    are those of `heights` (such as the geopotential height of each level,
    sampled as the field, which may vary in time and space), or else the
    field's vertical coordinate (in meters). Nature runs usually give heights
    above the geoid (mean sea level) rather than the ellipsoid: they differ
    by up to about 100 meters (see `tatc.utils.compute_geoid_height`).

    Other observations (points without heights, and regions, at their
    representative points) are sampled at every level, each a row with its
    `level` (and its `height`, from `heights`, if given), as the profile of a
    sounding below a satellite. The geometry of each row is the point at
    which it is sampled.

    Args:
        observations (geopandas.GeoDataFrame): The observations, with the time
            (`time_column`) of each.
        field (xarray.DataArray): The field, with a dimension of time, the
            latitudes and longitudes of its grid (see `sample_field`), and a
            vertical dimension.
        name (str | None): The column of the sampled values, or `None` for the
            field's name (or `"value"`).
        heights (xarray.DataArray | None): The heights (meters) of the field's
            levels, over its vertical dimension (and optionally its time and
            grid), or `None` for the field's vertical coordinate.
        method (str): `"nearest"` or `"linear"` (of a regular grid), in time
            and space.
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end`.
        dimensions (Mapping[str, str] | None): The names of the field's
            dimensions (see `sample_field`), including its vertical dimension
            (`level`), if not the defaults.

    Returns:
        geopandas.GeoDataFrame: a row for each level (or height) of each
            observation
    """
    if method not in ("nearest", "linear"):
        raise ValueError(f"method must be 'nearest' or 'linear', not {method!r}")
    name = name or field.name or "value"
    dimensions = dict(dimensions or {})
    time_dim = _get_dimension(field, ["time"], dimensions.get("time"))
    vertical = _get_vertical(
        field, _get_grid(field, dimensions, time_dim), time_dim, dimensions
    )
    if vertical is None:
        raise ValueError("profiles need a field with a vertical dimension")
    levels = field[vertical].to_numpy()
    times = _get_observation_times(observations, time_column).reset_index(drop=True)
    geometry = observations.geometry.to_numpy()
    # the vertices of points and multipoints, at heights or not
    pointed = np.isin(
        shapely.get_type_id(geometry),
        [shapely.GeometryType.POINT, shapely.GeometryType.MULTIPOINT],
    )
    coordinates, owner = shapely.get_coordinates(
        np.where(pointed, geometry, None), include_z=True, return_index=True
    )
    at_height = np.isfinite(coordinates[:, 2]) if len(coordinates) else np.zeros(0)
    vertex_owner = owner[at_height.astype(bool)]
    vertices = coordinates[at_height.astype(bool)]
    # the other observations: at their representative points
    elevated = np.zeros(len(observations), dtype=bool)
    elevated[vertex_owner] = True
    flat = np.flatnonzero(~elevated)
    flat_longitude, flat_latitude = _get_points(observations.geometry.iloc[flat])
    frames = []
    columns = observations.drop(columns=observations.geometry.name)
    if len(vertices):
        values = _sample_values(
            field,
            times.iloc[vertex_owner],
            vertices[:, 0],
            vertices[:, 1],
            method,
            dimensions,
        )
        level_heights = _sample_heights(
            heights,
            field,
            vertical,
            times.iloc[vertex_owner],
            vertices,
            method,
            dimensions,
        )
        frame = columns.iloc[vertex_owner].reset_index(names="observation")
        frame["level"] = np.nan
        frame["height"] = vertices[:, 2]
        frame[name] = _interpolate_heights(level_heights, values, vertices[:, 2])
        frame["geometry"] = shapely.points(vertices)
        frames.append(frame)
    if len(flat):
        values = _sample_values(
            field,
            times.iloc[flat],
            flat_longitude,
            flat_latitude,
            method,
            dimensions,
        )
        frame = columns.iloc[np.repeat(flat, len(levels))].reset_index(
            names="observation"
        )
        frame["level"] = np.tile(levels, len(flat))
        if heights is None:
            frame["height"] = np.nan
        else:
            positions = np.column_stack(
                [flat_longitude, flat_latitude, np.zeros(len(flat))]
            )
            frame["height"] = _sample_heights(
                heights,
                field,
                vertical,
                times.iloc[flat],
                positions,
                method,
                dimensions,
            ).ravel()
        frame[name] = np.reshape(values, -1)
        frame["geometry"] = np.repeat(
            shapely.points(flat_longitude, flat_latitude), len(levels)
        )
        frames.append(frame)
    if not frames:
        frame = columns.iloc[:0].reset_index(names="observation")
        for column in ("level", "height", name):
            frame[column] = pd.Series(dtype=float)
        frame["geometry"] = pd.Series(dtype=object)
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True)
    order = np.argsort(
        pd.Index(observations.index).get_indexer(result["observation"]), kind="stable"
    )
    return gpd.GeoDataFrame(
        result.iloc[order].reset_index(drop=True),
        geometry="geometry",
        crs=observations.crs,
    )


def _sample_heights(
    heights: Any,
    field: Any,
    vertical: str,
    times: pd.Series,
    positions: np.ndarray,
    method: str,
    dimensions: Mapping[str, str],
) -> np.ndarray:
    """
    Gets the heights of a field's levels at times and positions (see
    `sample_profiles`): those of `heights`, over the vertical dimension alone
    or sampled as the field, or else the field's vertical coordinate.
    """
    count = len(positions)
    if heights is None:
        heights = field[vertical]
    if set(heights.dims) == {vertical}:
        values = heights.sel({vertical: field[vertical].to_numpy()}).to_numpy()
        return np.broadcast_to(values.astype(float), (count, len(values))).copy()
    sampled = _sample_values(
        heights, times, positions[:, 0], positions[:, 1], method, dimensions
    )
    return np.reshape(sampled, (count, -1))
