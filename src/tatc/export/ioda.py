"""
Exports observations as the IODA files of the Joint Effort for Data
assimilation Integration (JEDI), in the layout written by JCSDA's
ioda-converters (`pyiodaconv/ioda_conv_engines.py`): a dimension `Location`
of the observations, and groups `MetaData` (places, times, and geometry),
`ObsValue` (values), `ObsError` (standard deviations of their errors), and
`PreQC` (quality flags). Requires netCDF4.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping
from os import PathLike
from typing import Any

import geopandas as gpd
import numpy as np

from .. import __version__
from ._common import get_places, get_times

FLOAT_FILL = np.float32(9.969209968386869e36)
"""The fill value of float32 variables (`get_default_fill_val` of
ioda-converters' `ioda_conv_engines.py`)."""

INT32_FILL = np.int32(-2147483647)
"""The fill value of int32 variables (as of ioda-converters)."""

INT64_FILL = np.int64(-9223372036854775806)
"""The fill value of int64 variables (as of ioda-converters)."""

DATETIME_UNITS = "seconds since 1970-01-01T00:00:00Z"
"""The units of `MetaData/dateTime` (int64 seconds, as of ioda-converters)."""

_GEOMETRY = {
    "sensorZenithAngle": ("platform_alt", lambda angle: 90 - angle),
    "sensorAzimuthAngle": ("platform_az", lambda angle: angle),
    "solarZenithAngle": ("solar_alt", lambda angle: 90 - angle),
    "solarAzimuthAngle": ("solar_az", lambda angle: angle),
}
"""The `MetaData` of the viewing geometry: each from a column of
observations (see `tatc.analysis.collect_observations`), in degrees."""


def _float(values: Any) -> np.ndarray:
    """Gets float32 values, with the fill value for missing values."""
    values = np.asarray(values, dtype=float)
    return np.where(np.isfinite(values), values, FLOAT_FILL).astype(np.float32)


def _error_column(observations: gpd.GeoDataFrame, column: str) -> str | None:
    """Gets the column of the standard deviations of the errors of a column of
    values (see `tatc.analysis.simulate_observation_errors`), if any."""
    base = column[: -len("_observed")] if column.endswith("_observed") else column
    for candidate in (f"{column}_error_std", f"{base}_error_std"):
        if candidate in observations.columns:
            return candidate
    return None


def to_ioda(
    observations: gpd.GeoDataFrame,
    path: str | PathLike,
    variables: Mapping[str, str],
    errors: Mapping[str, str] | None = None,
    qc: Mapping[str, str] | None = None,
    units: Mapping[str, str] | None = None,
    metadata: Mapping[str, str] | None = None,
    global_attributes: Mapping[str, Any] | None = None,
    time_column: str | None = None,
) -> None:
    """
    Exports observations (such as those sampled from a nature run with
    simulated errors, see `tatc.analysis.simulate_observation_errors`) as an
    IODA file of JEDI, as the observations of an observing system
    simulation experiment (OSSE) to assimilate: a dimension `Location` of the
    observations, and groups:

    - `MetaData`: `latitude` and `longitude` (degrees, of each point, or a
      point on the surface of each region), `dateTime` (int64 seconds since
      1970-01-01T00:00:00Z), `height` (meters above the WGS 84 ellipsoid, if
      any observation is a point with a z-coordinate), `sensorZenithAngle`,
      `sensorAzimuthAngle`, `solarZenithAngle`, and `solarAzimuthAngle`
      (degrees, from the columns `platform_alt`, `platform_az`, `solar_alt`,
      and `solar_az`, if given), `sequenceNumber` (the `observation` of
      profiles, see `tatc.analysis.sample_profiles`), and other numeric
      columns (`metadata`);
    - `ObsValue`: the values of each variable (float32), named as JEDI's
      UFO names them (such as `airTemperature`), which the user gives;
    - `ObsError`: the standard deviation of each variable's errors (by
      default, its `{column}_error_std` column, or that of the column without
      its `_observed` suffix, if any);
    - `PreQC`: quality flags of each variable (int32), from columns, if given.

    Missing values are the fill values of ioda-converters (`FLOAT_FILL`,
    `INT32_FILL`, and `INT64_FILL`). The variable `Location` holds the index of
    each observation; ioda-converters leave it to pyioda, whose values were not
    confirmed. Requires netCDF4.

    Args:
        observations (geopandas.GeoDataFrame): The observations.
        path (str | os.PathLike): The file to write.
        variables (Mapping[str, str]): The IODA name of the variable of each
            column of values (such as `{"temperature_observed":
            "airTemperature"}`).
        errors (Mapping[str, str] | None): The column of the standard
            deviations of the errors of each column of values, or `None` to
            find them (see above).
        qc (Mapping[str, str] | None): The column of the quality flags of
            each column of values, if any.
        units (Mapping[str, str] | None): The units of each column of
            values, if any.
        metadata (Mapping[str, str] | None): The IODA name of other numeric
            columns to write to `MetaData`.
        global_attributes (Mapping[str, Any] | None): The global attributes
            of the file (besides `converter`).
        time_column (str | None): The column of the times of the
            observations, or `None` (default) for the first of `epoch`,
            `time`, or `closest_time`, or else the midpoint of `start` and
            `end`.
    """
    # pylint: disable=import-outside-toplevel,too-many-locals
    import netCDF4

    errors = dict(errors or {})
    qc = dict(qc or {})
    units = dict(units or {})
    metadata = dict(metadata or {})
    needed = (
        list(variables) + list(errors.values()) + list(qc.values()) + list(metadata)
    )
    missing = [column for column in needed if column not in observations.columns]
    if missing:
        raise ValueError(f"the observations have no columns {missing}")
    longitude, latitude, height = get_places(observations)
    times = get_times(observations, time_column)
    seconds = np.round(
        (times.tz_convert("UTC").tz_localize(None) - np.datetime64("1970-01-01"))
        / np.timedelta64(1, "s")
    ).astype(np.int64)
    count = len(observations)
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("Location", count)
        location = dataset.createVariable("Location", "i4", ("Location",))
        location[:] = np.arange(count, dtype=np.int32)
        dataset.setncatts(
            {"converter": f"TAT-C {__version__}", **dict(global_attributes or {})}
        )
        groups = {
            name: dataset.createGroup(name)
            for name in ("MetaData", "ObsValue", "ObsError", "PreQC")
        }

        def write(group: str, name: str, values: np.ndarray, **attributes: Any):
            dtype = values.dtype
            fill = {
                np.dtype(np.float32): FLOAT_FILL,
                np.dtype(np.int32): INT32_FILL,
                np.dtype(np.int64): INT64_FILL,
            }[dtype]
            variable = groups[group].createVariable(
                name, dtype, ("Location",), fill_value=fill
            )
            variable[:] = values
            if attributes:
                variable.setncatts(attributes)

        write("MetaData", "latitude", _float(latitude), units="degrees_north")
        write("MetaData", "longitude", _float(longitude), units="degrees_east")
        write("MetaData", "dateTime", seconds, units=DATETIME_UNITS)
        if np.any(np.isfinite(height)):
            write("MetaData", "height", _float(height), units="m")
        for name, (column, convert) in _GEOMETRY.items():
            if column in observations.columns:
                values = convert(observations[column].to_numpy(dtype=float))
                write("MetaData", name, _float(values), units="degree")
        if "observation" in observations.columns:
            sequence = np.unique(observations["observation"], return_inverse=True)[1]
            write("MetaData", "sequenceNumber", sequence.astype(np.int32), units="1")
        for column, name in metadata.items():
            write("MetaData", name, _float(observations[column]))
        for column, name in variables.items():
            attributes = {"units": units[column]} if column in units else {}
            write("ObsValue", name, _float(observations[column]), **attributes)
            error = errors.get(column, _error_column(observations, column))
            if error is not None:
                write("ObsError", name, _float(observations[error]), **attributes)
            if column in qc:
                flags = observations[qc[column]].to_numpy(dtype=float)
                write(
                    "PreQC",
                    name,
                    np.where(np.isfinite(flags), flags, INT32_FILL).astype(np.int32),
                )
