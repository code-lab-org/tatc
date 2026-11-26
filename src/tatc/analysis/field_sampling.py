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
    value (`NaN` outside the period or extent of the field) or its linear
    interpolation in time and space.

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
    lat_dim = _get_dimension(field, ["lat", "latitude"], dimensions.get("lat"))
    lon_dim = _get_dimension(field, ["lon", "longitude"], dimensions.get("lon"))
    where = {
        time_dim: xr.DataArray(
            _field_times(field, time_dim, pd.Series(times)), dims="sample"
        ),
        lat_dim: xr.DataArray(latitude, dims="sample"),
        lon_dim: xr.DataArray(
            _field_longitudes(field, lon_dim, longitude), dims="sample"
        ),
    }
    if method == "linear":
        return field.interp(where).to_numpy().astype(float)
    values = field.sel(where, method="nearest").to_numpy().astype(float)
    sample_times = where[time_dim].to_numpy()
    coordinate = field[time_dim].to_numpy()
    values[(sample_times < coordinate.min()) | (sample_times > coordinate.max())] = (
        np.nan
    )
    values[
        _outside_field(field, lat_dim, lon_dim, latitude, where[lon_dim].to_numpy())
    ] = np.nan
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
    space. Observations outside the period of the field, or outside its
    extent (as of a regional nature run: by more than half a grid spacing
    beyond its outer latitudes or, unless they go around the globe,
    longitudes), have no value (`NaN`). Requires xarray.

    Args:
        observations (geopandas.GeoDataFrame): The observations, with the time
            (`time_column`) of each sample.
        field (xarray.DataArray): The field, with dimensions of time,
            latitude, and longitude (by default named `time`, `lat` or
            `latitude`, and `lon` or `longitude`), in degrees, with longitudes
            in [-180, 180) or [0, 360).
        name (str | None): The column of the sampled values, or `None` for the
            field's name (or `"value"`).
        method (str): `"nearest"` or `"linear"`.
        time_column (str | None): The column of the times of the samples,
            or `None` (default) for the first of `epoch`, `time`, or
            `closest_time`, or else the midpoint of `start` and `end`.
        dimensions (Mapping[str, str] | None): The names of the field's
            dimensions of `time`, `lat`, and `lon`, if not the defaults.

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
