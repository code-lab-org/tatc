"""
Input validation shared by the analysis functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

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
