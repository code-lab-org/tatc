"""
Exports observations as CF-netCDF point files (discrete sampling geometries
of `featureType` point). Requires xarray.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from os import PathLike
from typing import Any

import geopandas as gpd
import numpy as np
import pandas as pd

from .. import __version__
from ._common import get_places, get_times


def to_cf_netcdf(
    observations: gpd.GeoDataFrame,
    path: str | PathLike | None = None,
    variables: Sequence[str] | None = None,
    attributes: Mapping[str, Mapping[str, Any]] | None = None,
    global_attributes: Mapping[str, Any] | None = None,
    time_column: str | None = None,
) -> Any:
    """
    Exports observations (such as those sampled from a nature run, see
    `tatc.analysis.sample_field`) as a CF-netCDF point file: a discrete
    sampling geometry of `featureType` point, with a dimension `obs` of the
    observations, their `time` (UTC), `lat` and `lon` (degrees, of each point,
    or a point on the surface of each region), `height` (meters above the
    WGS 84 ellipsoid, if any observation is a point with a z-coordinate), and
    the values of their columns. Missing values are `NaN`. Requires xarray
    (and a netCDF engine, such as netCDF4, to write a file).

    Args:
        observations (geopandas.GeoDataFrame): The observations.
        path (str | os.PathLike | None): The file to write, or `None` to
            only return the dataset.
        variables (Sequence[str] | None): The columns to export, or `None`
            for every numeric column.
        attributes (Mapping[str, Mapping[str, Any]] | None): The attributes
            (such as `units`, `long_name`, or `standard_name`) of each
            exported column, by name.
        global_attributes (Mapping[str, Any] | None): The global attributes
            of the file (besides `Conventions`, `featureType`, and `source`).
        time_column (str | None): The column of the times of the
            observations, or `None` (default) for the first of `epoch`,
            `time`, or `closest_time`, or else the midpoint of `start` and
            `end`.

    Returns:
        xarray.Dataset: the dataset (written to `path`, if given)
    """
    # pylint: disable=import-outside-toplevel
    import xarray as xr

    if variables is None:
        variables = [
            column
            for column in observations.columns
            if column != observations.geometry.name
            and pd.api.types.is_numeric_dtype(observations[column])
            and not pd.api.types.is_bool_dtype(observations[column])
        ]
    missing = [column for column in variables if column not in observations.columns]
    if missing:
        raise ValueError(f"the observations have no columns {missing}")
    reserved = [column for column in variables if column in ("time", "lat", "lon")]
    if reserved:
        raise ValueError(f"columns {reserved} are reserved for the coordinates")
    longitude, latitude, height = get_places(observations)
    times = get_times(observations, time_column).tz_convert("UTC").tz_localize(None)
    coordinates: dict[str, Any] = {
        "time": (
            "obs",
            times.to_numpy(),
            {"standard_name": "time", "axis": "T"},
        ),
        "lat": (
            "obs",
            latitude,
            {"standard_name": "latitude", "units": "degrees_north", "axis": "Y"},
        ),
        "lon": (
            "obs",
            longitude,
            {"standard_name": "longitude", "units": "degrees_east", "axis": "X"},
        ),
    }
    if np.any(np.isfinite(height)):
        coordinates["height"] = (
            "obs",
            height,
            {
                "standard_name": "height_above_reference_ellipsoid",
                "units": "m",
                "positive": "up",
                "axis": "Z",
            },
        )
    attributes = attributes or {}
    data = {
        column: (
            "obs",
            np.asarray(observations[column].to_numpy(), dtype=float),
            dict(attributes.get(column, {})),
        )
        for column in variables
    }
    dataset = xr.Dataset(
        data,
        coords=coordinates,
        attrs={
            "Conventions": "CF-1.8",
            "featureType": "point",
            "source": f"TAT-C {__version__}",
            **dict(global_attributes or {}),
        },
    )
    if path is not None:
        dataset.to_netcdf(
            path,
            encoding={
                "time": {
                    "units": "seconds since 1970-01-01T00:00:00Z",
                    "dtype": "float64",
                }
            },
        )
    return dataset
