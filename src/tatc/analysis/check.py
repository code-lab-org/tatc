"""
Input validation and result handling shared by the analysis functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Iterable
from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd

from ..schemas import Satellite, SurfacePlatform
from .. import constants
from ..schemas.space.base_constellation import BaseConstellation

_MIN_PERIGEE_ALTITUDE = 100e3
"""
Perigee altitude (meters) below which a satellite's orbit is assumed to be
specified in the wrong units (see `_warn_low_perigees`).
"""


def _check_satellite(satellite: object, name: str = "satellite") -> Satellite:
    """
    Checks that `satellite` is a `Satellite`, raising a `TypeError` otherwise.

    A constellation carries its lead member's name, orbit, and instruments,
    so without this check it would be silently analyzed as that single
    satellite; the error instead directs to its `generate_members` method.

    Args:
        satellite (object): the value to check.
        name (str): the argument name to report in the error message.

    Returns:
        Satellite: the checked satellite.
    """
    if isinstance(satellite, Satellite):
        return satellite
    if isinstance(satellite, BaseConstellation):
        raise TypeError(
            f"{name} must be a Satellite, not a {type(satellite).__name__}; "
            "use its generate_members() method to pass its member satellites"
        )
    raise TypeError(f"{name} must be a Satellite, not a {type(satellite).__name__}")


def _check_platforms(value: object, kind: type, name: str) -> list:
    """
    Checks that `value` is an instance of `kind` (for example, an
    `Aircraft`) or a list of them, raising a `TypeError` otherwise.

    Args:
        value (object): The value to check.
        kind (type): The type of each platform.
        name (str): The argument name to report in the error message.

    Returns:
        list: the checked platforms, as a list.
    """
    checked = list(value) if isinstance(value, (list, tuple)) else [value]
    for i, platform in enumerate(checked):
        if not isinstance(platform, kind):
            label = f"{name}[{i}]" if isinstance(value, (list, tuple)) else name
            article = "an" if kind.__name__[0] in "AEIOU" else "a"
            raise TypeError(
                f"{label} must be {article} {kind.__name__}, "
                f"not a {type(platform).__name__}"
            )
    return checked


def _check_satellites(
    satellites: object, name: str = "satellites", allow_single: bool = True
) -> list[Satellite]:
    """
    Checks that `satellites` is a list of `Satellite`s (or, if
    `allow_single`, a single `Satellite`), raising a `TypeError` otherwise
    (see `_check_satellite`).

    Args:
        satellites (object): the value to check.
        name (str): the argument name to report in the error message.
        allow_single (bool): `True`, to also accept a single `Satellite`.

    Returns:
        list[Satellite]: the checked satellites, as a list.
    """
    if isinstance(satellites, list):
        checked = [
            _check_satellite(satellite, f"{name}[{i}]")
            for i, satellite in enumerate(satellites)
        ]
        _warn_low_perigees(checked)
        return checked
    if allow_single:
        checked = [_check_satellite(satellites, name)]
        _warn_low_perigees(checked)
        return checked
    if isinstance(satellites, BaseConstellation):
        raise TypeError(
            f"{name} must be a list of Satellites, not a "
            f"{type(satellites).__name__}; use its generate_members() method"
        )
    raise TypeError(
        f"{name} must be a list of Satellites, not a {type(satellites).__name__}"
    )


def _warn_low_perigees(satellites: list[Satellite], stacklevel: int = 4) -> None:
    """
    Warns (once, for all of them) if any satellites have a perigee altitude
    below 100 km, as from an orbit altitude specified in kilometers rather
    than meters, which would otherwise silently give empty or invalid
    results. Warned here, once per analysis, rather than as each orbit is
    created, since orbits are also created internally (as for constellation
    members or conversions to general perturbations orbits).

    Args:
        satellites (list[Satellite]): the satellites to check.
        stacklevel (int): the stack level of the warning, so that it refers
            to the analysis function's caller.
    """
    low = {}
    for satellite in satellites:
        orbit = satellite.orbit
        perigee = (
            orbit.get_semimajor_axis() * (1 - orbit.get_eccentricity())
            - constants.EARTH_MEAN_RADIUS
        )
        if perigee < _MIN_PERIGEE_ALTITUDE:
            low[satellite.name] = min(perigee, low.get(satellite.name, perigee))
    if low:
        names = ", ".join(repr(name) for name in list(low)[:3])
        if len(low) > 3:
            names += f", and {len(low) - 3} more"
        warnings.warn(
            f"satellites ({names}) have a perigee altitude below "
            f"{_MIN_PERIGEE_ALTITUDE / 1e3:.0f} km (as low as "
            f"{min(low.values()) / 1e3:.1f} km): check that orbit altitudes and "
            "semimajor axes are in meters",
            stacklevel=stacklevel,
        )


def _check_aware(name: str, value: datetime) -> None:
    """
    Checks that a datetime is timezone-aware, raising a `ValueError` otherwise
    (a naive datetime would otherwise fail deep within Skyfield).

    Args:
        name (str): the name of the datetime, for the error message.
        value (datetime): the datetime.
    """
    if isinstance(value, datetime) and value.utcoffset() is None:
        raise ValueError(
            f"{name} must be a timezone-aware datetime, "
            "e.g. datetime(2025, 1, 1, tzinfo=timezone.utc)"
        )


def _check_times(times: datetime | Iterable[datetime] | None) -> None:
    """
    Checks that the times of an analysis are timezone-aware datetimes,
    raising a `ValueError` naming the first that is not.

    Args:
        times (datetime | Iterable[datetime] | None): the time(s), or `None`.
    """
    if times is None:
        return
    for value in [times] if isinstance(times, datetime) else times:
        if isinstance(value, datetime):
            _check_aware(f"time {value.isoformat()}", value)


def _check_time_window(start: datetime, end: datetime) -> None:
    """
    Checks that an analysis time window has timezone-aware `start` and `end`
    datetimes, with `end` no earlier than `start`, raising a `ValueError`
    otherwise. A naive datetime would otherwise fail deep within Skyfield,
    and a reversed window would silently give an empty result.

    Args:
        start (datetime): the start of the time window.
        end (datetime): the end of the time window.
    """
    _check_aware("start", start)
    _check_aware("end", end)
    if end < start:
        raise ValueError(
            f"end ({end.isoformat()}) is before start ({start.isoformat()})"
        )


def _check_instrument_index(
    satellite: Satellite | SurfacePlatform,
    instrument_index: int,
    kind: str = "satellite",
) -> int:
    """
    Checks that `instrument_index` indexes one of a satellite's (or other
    platform's) instruments, raising an `IndexError` that names it otherwise.

    Args:
        satellite (Satellite | SurfacePlatform): the satellite or platform.
        instrument_index (int): the instrument index.
        kind (str): the kind of platform, as named in the error.

    Returns:
        int: the checked instrument index.
    """
    count = len(satellite.instruments)
    if not -count <= instrument_index < count:
        raise IndexError(
            f"instrument_index {instrument_index} is out of range for {kind} "
            f"{satellite.name!r}, which has {count} "
            f"instrument{'' if count == 1 else 's'}"
        )
    return instrument_index


def _is_single(*values: object, instrument_index: int | None = 0) -> bool:
    """
    Checks whether an analysis of targets and satellites is of a single one
    of each (none of `values` is a list or tuple) and a single instrument
    (an integer `instrument_index`), so that its result is not combined
    (see `_combine_results`).

    Args:
        values (object): The targets, satellites, or other inputs, each a
            single value or a list (or tuple) of them.
        instrument_index (int | None): The instrument index (None for every
            instrument).

    Returns:
        bool: True, if every input is single
    """
    return instrument_index is not None and not any(
        isinstance(value, (list, tuple)) for value in values
    )


def _get_instrument_indices(
    satellite: Satellite | SurfacePlatform,
    instrument_index: int | None,
    kind: str = "satellite",
) -> range | list[int]:
    """
    Gets the indices of a satellite's (or other platform's) observing
    instruments: the given index (see `_check_instrument_index`), or every
    instrument's if None.

    Args:
        satellite (Satellite | SurfacePlatform): The satellite or platform.
        instrument_index (int | None): The instrument index (None for every
            instrument).
        kind (str): the kind of platform, as named in an error.

    Returns:
        range | list[int]: the instrument indices
    """
    if instrument_index is None:
        return range(len(satellite.instruments))
    return [_check_instrument_index(satellite, instrument_index, kind)]


def _combine_results(
    results: list[gpd.GeoDataFrame],
    single: bool,
    sort_by: str,
    empty: Callable[[], gpd.GeoDataFrame],
    stable: bool = True,
) -> gpd.GeoDataFrame:
    """
    Combines the results of an analysis of one or more targets, satellites,
    and instruments: the result of a single one (see `_is_single`), or
    otherwise all of them, concatenated, sorted, and re-indexed (or an
    empty result if there are none).

    Args:
        results (list[geopandas.GeoDataFrame]): The results.
        single (bool): True, if the analysis is of a single one of each.
        sort_by (str): The column by which to sort combined results.
        empty (Callable[[], geopandas.GeoDataFrame]): Gets an empty result.
        stable (bool): True, to keep the order of results that sort equally.

    Returns:
        geopandas.GeoDataFrame: the combined results
    """
    if single:
        return results[0]
    if len(results) == 0:
        return empty()
    return (
        pd.concat(results)
        .sort_values(sort_by, kind="stable" if stable else "quicksort")
        .reset_index(drop=True)
    )


_PLATFORM_COLUMNS = ("satellite", "aircraft", "surface_platform", "radar_station")
"""
Columns that name the observing platform of a result, one for each kind of
platform: satellites, aircraft, fixed surface platforms, and radar stations.
"""


def _get_platform_columns(frame: pd.DataFrame) -> list[str]:
    """
    Gets the columns of a result that name its observing platforms (see
    `_PLATFORM_COLUMNS`), in their standard order.

    Args:
        frame (pandas.DataFrame): The result.

    Returns:
        list[str]: the platform columns
    """
    return [column for column in _PLATFORM_COLUMNS if column in frame.columns]


def _get_platforms(frame: pd.DataFrame) -> np.ndarray:
    """
    Gets the observing platform of each row of a result, named by its
    `satellite`, `aircraft`, `surface_platform`, or `radar_station` column,
    as its kind and name (so that a satellite and an aircraft of the same
    name are different platforms).

    Args:
        frame (pandas.DataFrame): The result.

    Returns:
        numpy.ndarray: the platform (kind and name) of each row
    """
    platforms = pd.Series(None, index=frame.index, dtype=object)
    for column in _get_platform_columns(frame):
        named = platforms.isna() & frame[column].notna()
        platforms[named] = column + ": " + frame.loc[named, column].astype(str)
    if platforms.isna().any():
        raise ValueError(
            "each row must name its platform in a "
            + ", ".join(_PLATFORM_COLUMNS[:-1])
            + f", or {_PLATFORM_COLUMNS[-1]} column"
        )
    return platforms.to_numpy()


def _join_names(names: pd.Series) -> str | None:
    """
    Joins the names of the contributors to a combined result (such as the
    platforms of a coincident observation), comma-separated, omitting missing
    names (as of the platforms of another kind).

    Args:
        names (pandas.Series): The names.

    Returns:
        str | None: the joined names, or `None` if all are missing
    """
    present = names.dropna()
    return ", ".join(present) if len(present) > 0 else None


def _aggregate_names(
    grouped: pd.core.groupby.DataFrameGroupBy, column: str
) -> pd.Series:
    """
    Joins the names in a column of each group of a result (see
    `_join_names`), quickly if no name is missing.

    Args:
        grouped (pandas.core.groupby.DataFrameGroupBy): The grouped result.
        column (str): The column of names.

    Returns:
        pandas.Series: the joined names of each group
    """
    names = grouped[column]
    if names.obj.notna().all():
        return names.agg(", ".join)
    return names.agg(_join_names)


def _dissolve_footprints(
    track: gpd.GeoDataFrame, by: str | list[str] | None = None
) -> gpd.GeoDataFrame:
    """
    Dissolves the footprints of valid observations (`valid_obs`, if given)
    of a track into one geometry, or one for each group of `by` (see
    `tatc.analysis.footprints.dissolve_footprints`).

    Args:
        track (geopandas.GeoDataFrame): The track.
        by (str | list[str] | None): The column(s) by which to group the
            footprints, or `None` for one geometry.

    Returns:
        geopandas.GeoDataFrame: the dissolved footprints
    """
    columns = [] if by is None else [by] if isinstance(by, str) else list(by)
    crs = track.crs if track.crs is not None else "EPSG:4326"
    if "valid_obs" in track.columns:
        track = track[track["valid_obs"].eq(True).to_numpy()]
    if track.empty:
        return gpd.GeoDataFrame(
            {
                column: pd.Series(
                    [], dtype=track[column].dtype if column in track else "object"
                )
                for column in columns
            },
            geometry=gpd.GeoSeries([], crs=crs),
            crs=crs,
        )
    geometry = track.geometry.name
    if not columns:
        return track[[geometry]].dissolve().reset_index(drop=True)
    dissolved = track[[*columns, geometry]].dissolve(by=columns, dropna=False)
    return dissolved.reset_index()[[*columns, geometry]]


def _get_empty_observation_frame(
    platform_columns: list[str] | None = None,
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame of observations of targets (such as aggregated
    or coincident observations), with the columns that name their platforms.

    Args:
        platform_columns (list[str] | None): The columns that name the
            platforms (see `_PLATFORM_COLUMNS`), by default `satellite`.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "target_hash": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        **{
            column: pd.Series([], dtype="str")
            for column in platform_columns or ["satellite"]
        },
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")
