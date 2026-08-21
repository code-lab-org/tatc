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
from .two_line_elements import TwoLineElements

AllOrbits = (
    CircularOrbit
    | GeneralPerturbationsOrbit
    | GeosynchronousOrbit
    | KeplerianOrbit
    | MolniyaOrbit
    | SunSynchronousOrbit
    | TundraOrbit
    | TwoLineElements
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
    "TwoLineElements",
]
