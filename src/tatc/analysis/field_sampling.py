"""
Methods to sample fields (such as the variables of a nature run in an
observing system simulation experiment, OSSE) at the times and places of
observations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd

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
            return field.interp(where).to_numpy().astype(float)
        values = field.sel(where, method="nearest").to_numpy().astype(float)
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
        values = field.isel(indexers).to_numpy().astype(float)
    values[outside | _outside_period(field, time_dim, sample_times)] = np.nan
    return values


def sample_field(
    observations: gpd.GeoDataFrame,
    field: Any,
    name: str | None = None,
    method: str = "nearest",
    time_column: str | None = None,
    dimensions: Mapping[str, str] | None = None,
) -> gpd.GeoDataFrame:
    """
    Samples a field (such as a variable of a nature run) at the time and place
    of each observation: at its point (or the representative point of its
    region), by the nearest value or by linear interpolation in time and
    space. Observations outside the period of the field (by more than half a
    time step before its first time or after its last), or outside its
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

    Args:
        observations (geopandas.GeoDataFrame): The observations, with the time
            (`time_column`) of each sample.
        field (xarray.DataArray): The field, with a dimension of time and the
            latitudes and longitudes of its grid (by default named `time`,
            and any of `lat`, `latitude`, `XLAT`, `latCell`, or `nav_lat`
            and of `lon`, `longitude`, `XLONG`, `lonCell`, or `nav_lon`),
            with longitudes in [-180, 180) or [0, 360).
        name (str | None): The column of the sampled values, or `None` for the
            field's name (or `"value"`).
        method (str): `"nearest"` or `"linear"`.
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end`.
        dimensions (Mapping[str, str] | None): The names of the field's
            `time` dimension and of its latitudes (`lat`) and longitudes
            (`lon`): the dimensions of a regular grid, or the coordinates of
            another grid, if not the defaults.

    Returns:
        geopandas.GeoDataFrame: the observations, with the sampled values
    """
    if method not in ("nearest", "linear"):
        raise ValueError(f"method must be 'nearest' or 'linear', not {method!r}")
    name = name or field.name or "value"
    result = observations.copy()
    if len(observations) == 0:
        result[name] = pd.Series(dtype=float)
        return result
    longitude, latitude = _get_points(observations.geometry)
    result[name] = _sample_values(
        field,
        _get_observation_times(observations, time_column),
        longitude,
        latitude,
        method,
        dimensions,
    )
    return result
