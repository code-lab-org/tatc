"""
Defines object schemas.
"""

from .architecture import Architecture
from .instrument import AllInstruments, ConicalInstrument, Instrument, PointedInstrument
from .orbit import (
    AllOrbits,
    CircularOrbit,
    GeneralPerturbationsOrbit,
    GeosynchronousOrbit,
    KeplerianOrbit,
    MolniyaOrbit,
    SunSynchronousOrbit,
    TundraOrbit,
    TwoLineElements,
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
from .surface import (
    AllSurfaceObjects,
    GroundStation,
    Point,
    RadarBand,
    RadarStation,
    TerrainMask,
)

__all__ = [
    "AllInstruments",
    "AllOrbits",
    "AllSpaceObjects",
    "AllSurfaceObjects",
    "Architecture",
    "CircularOrbit",
    "ConicalInstrument",
    "GeneralPerturbationsOrbit",
    "GeosynchronousOrbit",
    "GroundStation",
    "Instrument",
    "KeplerianOrbit",
    "MOGConstellation",
    "MolniyaOrbit",
    "Point",
    "PointedInstrument",
    "RadarBand",
    "RadarStation",
    "SOCConstellation",
    "Satellite",
    "SunSynchronousOrbit",
    "TerrainMask",
    "TrainConstellation",
    "TundraOrbit",
    "TwoLineElements",
    "WalkerConfiguration",
    "WalkerConstellation",
]
