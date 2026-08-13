"""
Object schemas for orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .circular import CircularOrbit
from .geosynchronous import GeosynchronousOrbit
from .gp import GeneralPerturbationsOrbit
from .keplerian import KeplerianOrbit
from .molniya import MolniyaOrbit
from .sun_synchronous import SunSynchronousOrbit
from .tundra import TundraOrbit

AllOrbits = (
    CircularOrbit
    | GeneralPerturbationsOrbit
    | GeosynchronousOrbit
    | KeplerianOrbit
    | MolniyaOrbit
    | SunSynchronousOrbit
    | TundraOrbit
)

__all__ = [
    "AllOrbits",
    "CircularOrbit",
    "GeneralPerturbationsOrbit",
    "GeosynchronousOrbit",
    "KeplerianOrbit",
    "MolniyaOrbit",
    "SunSynchronousOrbit",
    "TundraOrbit",
]
