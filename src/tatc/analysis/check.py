"""
Input validation and result handling shared by the analysis functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable

import geopandas as gpd
import pandas as pd

from ..schemas import Satellite
from ..schemas.space.base_constellation import BaseConstellation


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
        return [
            _check_satellite(satellite, f"{name}[{i}]")
            for i, satellite in enumerate(satellites)
        ]
    if allow_single:
        return [_check_satellite(satellites, name)]
    if isinstance(satellites, BaseConstellation):
        raise TypeError(
            f"{name} must be a list of Satellites, not a "
            f"{type(satellites).__name__}; use its generate_members() method"
        )
    raise TypeError(
        f"{name} must be a list of Satellites, not a {type(satellites).__name__}"
    )


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
    index, or every instrument's if None.

    Args:
        satellite (Satellite): The satellite.
        instrument_index (int | None): The instrument index (None for every
            instrument).

    Returns:
        range | list[int]: the instrument indices
    """
    if instrument_index is None:
        return range(len(satellite.instruments))
    return [instrument_index]


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
