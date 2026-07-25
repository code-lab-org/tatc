"""
Object schemas for orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .circular import CircularOrbit
from .gp import GeneralPerturbationsOrbit
from .keplerian import KeplerianOrbit
from .molniya import MolniyaOrbit
from .sun_synchronous import SunSynchronousOrbit
from .tundra import TundraOrbit

AllOrbits = (
    CircularOrbit
    | GeneralPerturbationsOrbit
    | KeplerianOrbit
    | MolniyaOrbit
    | SunSynchronousOrbit
    | TundraOrbit
)

__all__ = [
    "AllOrbits",
    "CircularOrbit",
    "GeneralPerturbationsOrbit",
    "KeplerianOrbit",
    "MolniyaOrbit",
    "SunSynchronousOrbit",
    "TundraOrbit",
]
