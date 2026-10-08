"""
Input validation and result handling shared by the analysis functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from datetime import datetime

import geopandas as gpd
import pandas as pd

from ..schemas import Satellite
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
    for name, value in (("start", start), ("end", end)):
        if isinstance(value, datetime) and value.utcoffset() is None:
            raise ValueError(
                f"{name} must be a timezone-aware datetime, "
                "e.g. datetime(2025, 1, 1, tzinfo=timezone.utc)"
            )
    if end < start:
        raise ValueError(
            f"end ({end.isoformat()}) is before start ({start.isoformat()})"
        )


def _check_instrument_index(satellite: Satellite, instrument_index: int) -> int:
    """
    Checks that `instrument_index` indexes one of a satellite's instruments,
    raising an `IndexError` that names the satellite otherwise.

    Args:
        satellite (Satellite): the satellite.
        instrument_index (int): the instrument index.

    Returns:
        int: the checked instrument index.
    """
    count = len(satellite.instruments)
    if not -count <= instrument_index < count:
        raise IndexError(
            f"instrument_index {instrument_index} is out of range for satellite "
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
    satellite: Satellite, instrument_index: int | None
) -> range | list[int]:
    """
    Gets the indices of a satellite's observing instruments: the given
    index (see `_check_instrument_index`), or every instrument's if None.

    Args:
        satellite (Satellite): The satellite.
        instrument_index (int | None): The instrument index (None for every
            instrument).

    Returns:
        range | list[int]: the instrument indices
    """
    if instrument_index is None:
        return range(len(satellite.instruments))
    return [_check_instrument_index(satellite, instrument_index)]


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
