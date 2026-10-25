"""
Defines object schemas.
"""

from .architecture import Architecture
from .attitude import Attitude, AttitudeSchedule, AttitudeSegment
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
    SurfacePlatform,
    TerrainMask,
)

__all__ = [
    "AllInstruments",
    "AllOrbits",
    "AllSpaceObjects",
    "AllSurfaceObjects",
    "Architecture",
    "Attitude",
    "AttitudeSchedule",
    "AttitudeSegment",
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
    "SurfacePlatform",
    "TerrainMask",
    "TrainConstellation",
    "TundraOrbit",
    "TwoLineElements",
    "WalkerConfiguration",
    "WalkerConstellation",
]
