"""
Events of the deployments of constellations: the production and launches of
satellites (see `Deployment` and `MissionHistory`), and values given for
each type of satellite.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

PerType = float | Mapping[Any, float]
"""A value for every type of satellite, or for each (by its name)."""


def _per_type(value: Any, kind: Any, name: str) -> Any:
    """Gets the value of a type of satellite."""
    if isinstance(value, Mapping):
        if kind not in value:
            raise ValueError(f"no {name} of satellites of kind {kind!r}")
        return value[kind]
    return value


@dataclass(frozen=True)
class Launch:
    """
    A launch of satellites of a plane (see `Deployment` and
    `MissionHistory`).
    """

    time: float
    """float: The time (years from the start of development)."""
    plane: int
    """int: The plane (its index, in the order of the layouts' planes)."""
    kind: Any
    """Any: The type of its satellites."""
    first: int
    """int: The position in its plane of its first satellite (-1 for the
    launches of restorations and spares)."""
    satellites: int
    """int: The number of satellites (including spares)."""
    generation: int = 0
    """int: The generation (0 for the first, others planned replacements)."""
    purpose: str = "deployment"
    """str: `"deployment"`, `"replacement"` (planned, at the end of a
    design life), `"restoration"` (of failed satellites), or `"spare"` (a
    spare in orbit replenished)."""
    failed: bool = False
    """bool: Whether the launch failed, losing its satellites."""


@dataclass(frozen=True)
class Production:
    """
    The production of satellites of a type (see `Deployment` and
    `MissionHistory`).
    """

    kind: Any
    """Any: The type of the satellites."""
    begin: float
    """float: The start of their production (years)."""
    end: float
    """float: The end of their production (years)."""
    satellites: int
    """int: The number of satellites."""
    purpose: str = "deployment"
    """str: `"deployment"`, `"replacement"`, `"restoration"`, or `"spare"`."""
