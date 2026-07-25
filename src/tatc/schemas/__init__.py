"""
Defines object schemas.
"""

from .architecture import Architecture
from .instrument import AllInstruments, Instrument, PointedInstrument
from .orbit import (
    AllOrbits,
    CircularOrbit,
    GeneralPerturbationsOrbit,
    KeplerianOrbit,
    MolniyaOrbit,
    SunSynchronousOrbit,
    TundraOrbit,
)
from .space import (
    AllSpaceObjects,
    MOGConstellation,
    Satellite,
    SOCConstellation,
    TrainConstellation,
    WalkerConfiguration,
    WalkerConstellation,
)
from .surface import AllSurfaceObjects, GroundStation, Point

__all__ = [
    "AllInstruments",
    "AllOrbits",
    "AllSpaceObjects",
    "AllSurfaceObjects",
    "Architecture",
    "CircularOrbit",
    "GeneralPerturbationsOrbit",
    "GroundStation",
    "Instrument",
    "KeplerianOrbit",
    "MOGConstellation",
    "MolniyaOrbit",
    "Point",
    "PointedInstrument",
    "SOCConstellation",
    "Satellite",
    "SunSynchronousOrbit",
    "TrainConstellation",
    "TundraOrbit",
    "WalkerConfiguration",
    "WalkerConstellation",
]