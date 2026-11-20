"""
Input validation and result handling shared by the analysis functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd

from .. import constants
from .._warnings import warn
from ..schemas import Aircraft, RadarStation, Satellite, SurfacePlatform
from ..schemas.space.base_constellation import BaseConstellation

_MIN_PERIGEE_ALTITUDE = 100e3
"""
Perigee altitude (meters) below which a satellite's orbit is assumed to be
specified in the wrong units (see `_warn_low_perigees`).
"""


@dataclass(frozen=True)
class _PlatformKind:
    """
    A kind of observing platform, and how analyses name it.

    Attributes:
        type (type): The schema of platforms of this kind.
        column (str): The column that names a platform in results.
        label (str): The kind, as named in errors.
        argument (str): The argument that takes platforms of this kind in
            the analyses of a single kind (such as `collect_ground_track`).
        fixed (bool): `True` if platforms of this kind do not move in the
            Earth-fixed frame, so that their footprints do not change with
            time (and may be collected without times).
        instrumented (bool): `True` if platforms of this kind carry
            instruments (`instruments`), `False` if a platform is itself the
            instrument (as a radar station is).
    """

    type: type
    column: str
    label: str
    argument: str
    fixed: bool = False
    instrumented: bool = True


SATELLITE = _PlatformKind(Satellite, "satellite", "satellite", "satellites")
AIRCRAFT = _PlatformKind(Aircraft, "aircraft", "aircraft", "aircraft")
SURFACE_PLATFORM = _PlatformKind(
    SurfacePlatform, "surface_platform", "surface platform", "platforms", fixed=True
)
RADAR_STATION = _PlatformKind(
    RadarStation,
    "radar_station",
    "radar station",
    "stations",
    fixed=True,
    instrumented=False,
)

_PLATFORM_KINDS = (SATELLITE, AIRCRAFT, SURFACE_PLATFORM, RADAR_STATION)
"""
The kinds of observing platforms, in the order their results are combined
and their columns are listed.
"""

_PLATFORM_COLUMNS = tuple(kind.column for kind in _PLATFORM_KINDS)
"""
Columns that name the observing platform of a result, one for each kind of
platform: satellites, aircraft, fixed surface platforms, and radar stations.
"""


def _get_platform_kind(platform: object) -> _PlatformKind | None:
    """
    Gets the kind of an observing platform (see `_PLATFORM_KINDS`).

    Args:
        platform (object): The platform.

    Returns:
        _PlatformKind | None: the kind, or `None` if it is not a platform
    """
    return next(
        (kind for kind in _PLATFORM_KINDS if isinstance(platform, kind.type)), None
    )


def _get_constellation_hint(value: object) -> str:
    """
    Gets a hint for an error about a value that is not a platform: to pass
    the members of a constellation, which carries its lead member's name,
    orbit, and instruments (so that, without the error, it would be silently
    analyzed as that single satellite).

    Args:
        value (object): The value.

    Returns:
        str: the hint (beginning with a semicolon), or an empty string
    """
    if isinstance(value, BaseConstellation):
        return "; use its generate_members() method to pass its member satellites"
    return ""


def _group_platforms(platforms: object) -> dict[_PlatformKind, list]:
    """
    Groups observing platforms of any kinds by kind (see `_PLATFORM_KINDS`),
    raising a `TypeError` for any other value.

    Args:
        platforms (object): A platform or a list of them.

    Returns:
        dict[_PlatformKind, list]: the platforms of each kind, in their order
    """
    is_list = isinstance(platforms, (list, tuple))
    groups: dict[_PlatformKind, list] = {kind: [] for kind in _PLATFORM_KINDS}
    for i, platform in enumerate(platforms if is_list else [platforms]):  # type: ignore
        kind = _get_platform_kind(platform)
        if kind is None:
            label = f"platforms[{i}]" if is_list else "platforms"
            names = [kind.type.__name__ for kind in _PLATFORM_KINDS]
            raise TypeError(
                f"{label} must be a {', '.join(names[:-1])}, or {names[-1]}, "
                f"not a {type(platform).__name__}{_get_constellation_hint(platform)}"
            )
        groups[kind].append(platform)
    for kind, members in groups.items():
        _check_distinct_names(members, kind, "platforms")
    return groups


def _check_platform(platform: object, kind: _PlatformKind, name: str) -> object:
    """
    Checks that `platform` is a platform of a kind (for example, an
    `Aircraft`), raising a `TypeError` otherwise (directing a constellation
    to its `generate_members` method, see `_get_constellation_hint`).

    Args:
        platform (object): The value to check.
        kind (_PlatformKind): The kind of platform.
        name (str): The argument name to report in the error message.

    Returns:
        object: the checked platform.
    """
    if isinstance(platform, kind.type):
        return platform
    article = "an" if kind.type.__name__[0] in "AEIOU" else "a"
    raise TypeError(
        f"{name} must be {article} {kind.type.__name__}, "
        f"not a {type(platform).__name__}{_get_constellation_hint(platform)}"
    )


def _check_satellite(satellite: object, name: str = "satellite") -> Satellite:
    """
    Checks that `satellite` is a `Satellite`, raising a `TypeError` otherwise
    (see `_check_platform`).

    Args:
        satellite (object): the value to check.
        name (str): the argument name to report in the error message.

    Returns:
        Satellite: the checked satellite.
    """
    return _check_platform(satellite, SATELLITE, name)  # type: ignore


def _check_distinct_names(platforms: list, kind: _PlatformKind, name: str) -> None:
    """
    Checks that platforms of a kind have distinct names, raising a
    `ValueError` otherwise: results name each platform in the column of its
    kind (see `_get_platforms`), so platforms of the same name would be
    merged (as in coincidences) or lost.

    Args:
        platforms (list): The platforms of the kind.
        kind (_PlatformKind): The kind of platform.
        name (str): The argument name to report in the error message.
    """
    seen = set()
    for platform in platforms:
        if platform.name in seen:
            raise ValueError(
                f"{name} must have distinct names: more than one {kind.label} "
                f"is named {platform.name!r}"
            )
        seen.add(platform.name)


def _check_platforms(
    value: object,
    kind: _PlatformKind,
    name: str | None = None,
    allow_single: bool = True,
) -> list:
    """
    Checks that `value` is a list (or tuple) of platforms of a kind (for
    example, `Aircraft`) or, if `allow_single`, a single one, raising a
    `TypeError` otherwise (see `_check_platform`), and warns of satellites
    with low perigees (see `_warn_low_perigees`).

    Args:
        value (object): The value to check.
        kind (_PlatformKind): The kind of platform.
        name (str | None): The argument name to report in the error
            message, by default that of the kind (see `_PlatformKind`).
        allow_single (bool): `True`, to also accept a single platform.

    Returns:
        list: the checked platforms, as a list.
    """
    name = kind.argument if name is None else name
    if isinstance(value, (list, tuple)):
        checked = [
            _check_platform(platform, kind, f"{name}[{i}]")
            for i, platform in enumerate(value)
        ]
    elif allow_single:
        checked = [_check_platform(value, kind, name)]
    else:
        raise TypeError(
            f"{name} must be a list of {kind.type.__name__}s, "
            f"not a {type(value).__name__}{_get_constellation_hint(value)}"
        )
    if kind is SATELLITE:
        _warn_low_perigees(checked)
    return checked


def _check_satellites(
    satellites: object,
    name: str = "satellites",
    allow_single: bool = True,
) -> list[Satellite]:
    """
    Checks that `satellites` is a list (or tuple) of `Satellite`s (or, if
    `allow_single`, a single `Satellite`), raising a `TypeError` otherwise,
    and warns of low perigees (see `_check_platforms`).

    Args:
        satellites (object): the value to check.
        name (str): the argument name to report in the error message.
        allow_single (bool): `True`, to also accept a single `Satellite`.

    Returns:
        list[Satellite]: the checked satellites, as a list.
    """
    return _check_platforms(satellites, SATELLITE, name, allow_single)


def _warn_low_perigees(satellites: list[Satellite]) -> None:
    """
    Warns (once, for all of them) if any satellites have a perigee altitude
    below 100 km, as from an orbit altitude specified in kilometers rather
    than meters, which would otherwise silently give empty or invalid
    results. Warned here, once per analysis, rather than as each orbit is
    created, since orbits are also created internally (as for constellation
    members or conversions to general perturbations orbits).

    Args:
        satellites (list[Satellite]): the satellites to check. The warning
            refers to the code that called the analysis (see
            `tatc._warnings.warn`).
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
        warn(
            f"satellites ({names}) have a perigee altitude below "
            f"{_MIN_PERIGEE_ALTITUDE / 1e3:.0f} km (as low as "
            f"{min(low.values()) / 1e3:.1f} km): check that orbit altitudes and "
            "semimajor axes are in meters"
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
    platform: Satellite | Aircraft | SurfacePlatform, instrument_index: int
) -> int:
    """
    Checks that `instrument_index` indexes one of a platform's instruments,
    raising an `IndexError` that names the platform otherwise.

    Args:
        platform (Satellite | Aircraft | SurfacePlatform): the platform.
        instrument_index (int): the instrument index.

    Returns:
        int: the checked instrument index.
    """
    count = len(platform.instruments)
    if not -count <= instrument_index < count:
        kind = _get_platform_kind(platform)
        raise IndexError(
            f"instrument_index {instrument_index} is out of range for "
            f"{kind.label if kind is not None else 'platform'} "
            f"{platform.name!r}, which has {count} "
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
    platform: Satellite | Aircraft | SurfacePlatform, instrument_index: int | None
) -> range | list[int]:
    """
    Gets the indices of a platform's observing instruments: the given index
    (see `_check_instrument_index`), or every instrument's if None.

    Args:
        platform (Satellite | Aircraft | SurfacePlatform): The platform.
        instrument_index (int | None): The instrument index (None for every
            instrument).

    Returns:
        range | list[int]: the instrument indices
    """
    if instrument_index is None:
        return range(len(platform.instruments))
    return [_check_instrument_index(platform, instrument_index)]


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


def _get_platform_kinds(frame: pd.DataFrame) -> np.ndarray:
    """
    Gets the kind of the observing platform of each row of a result: the
    first of its `satellite`, `aircraft`, `surface_platform`, or
    `radar_station` columns (see `_PLATFORM_COLUMNS`) that names it.

    Args:
        frame (pandas.DataFrame): The result.

    Returns:
        numpy.ndarray: the platform column of each row
    """
    kinds = pd.Series(None, index=frame.index, dtype=object)
    for column in _get_platform_columns(frame):
        kinds[kinds.isna() & frame[column].notna()] = column
    if kinds.isna().any():
        raise ValueError(
            "each row must name its platform in a "
            + ", ".join(_PLATFORM_COLUMNS[:-1])
            + f", or {_PLATFORM_COLUMNS[-1]} column"
        )
    return kinds.to_numpy()


def _get_platforms(frame: pd.DataFrame) -> np.ndarray:
    """
    Gets the observing platform of each row of a result, named by its
    `satellite`, `aircraft`, `surface_platform`, or `radar_station` column
    (see `_get_platform_kinds`), as its kind and name (so that a satellite
    and an aircraft of the same name are different platforms).

    Args:
        frame (pandas.DataFrame): The result.

    Returns:
        numpy.ndarray: the platform (kind and name) of each row
    """
    kinds = _get_platform_kinds(frame)
    names = np.empty(len(frame), dtype=object)
    for column in pd.unique(kinds):
        rows = kinds == column
        names[rows] = _get_platform_key(
            column, frame.loc[rows, column].astype(str)
        ).to_numpy()
    return names


def _get_platform_key(column: str, name: str | pd.Series) -> str | pd.Series:
    """
    Gets the key that identifies a platform (or each of a series of
    platforms) of a kind across kinds: its kind's column and its name (see
    `_get_platforms`).

    Args:
        column (str): The column of the platform's kind (see
            `_PLATFORM_COLUMNS`).
        name (str | pandas.Series): The platform's name, or a series of names.

    Returns:
        str | pandas.Series: the key(s)
    """
    return column + ": " + name


def _get_instruments(frame: pd.DataFrame) -> np.ndarray:
    """
    Gets the instrument of each row of a result: its `instrument` column, or
    `None` for every row if it has none (as for the footprints of radar
    stations alone).

    Args:
        frame (pandas.DataFrame): The result.

    Returns:
        numpy.ndarray: the instrument of each row
    """
    if "instrument" in frame.columns:
        return frame["instrument"].to_numpy()
    return np.full(len(frame), None, dtype=object)


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
