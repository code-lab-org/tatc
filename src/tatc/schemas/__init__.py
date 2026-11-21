"""
Defines object schemas.
"""

from .air import Aircraft, Trajectory
from .architecture import Architecture
from .attitude import (
    Attitude,
    AttitudeSchedule,
    AttitudeSegment,
    GlintAttitude,
    TargetAttitude,
)
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

AllPlatforms = Satellite | Aircraft | SurfacePlatform | RadarStation
"""
Observing platforms: satellites, aircraft, fixed surface platforms, and
radar stations (see `tatc.analysis.collect_footprints`).
"""

__all__ = [
    "Aircraft",
    "AllInstruments",
    "AllOrbits",
    "AllPlatforms",
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
    "TargetAttitude",
    "GlintAttitude",
    "TerrainMask",
    "Trajectory",
    "TrainConstellation",
    "TundraOrbit",
    "TwoLineElements",
    "WalkerConfiguration",
    "WalkerConstellation",
]
