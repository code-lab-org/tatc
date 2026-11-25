"""
Defines object schemas.
"""

from .air import Aircraft, Trajectory
from .architecture import ActiveWindow, Architecture, Relay
from .attitude import (
    Attitude,
    AttitudeSchedule,
    AttitudeSegment,
    GlintAttitude,
    TargetAttitude,
)
from .communication import (
    AllReceivers,
    AllTransmitters,
    Band,
    DataFlow,
    DataGeneration,
    LinkMode,
    Receiver,
    RelayCommunications,
    ShannonModem,
    SpacecraftCommunications,
    StationCommunications,
    Transmitter,
)
from .instrument import (
    AllInstruments,
    ConicalInstrument,
    InSituInstrument,
    Instrument,
    PointedInstrument,
)
from .orbit import (
    AllOrbits,
    CircularOrbit,
    GeneralPerturbationsOrbit,
    GeosynchronousOrbit,
    KeplerianOrbit,
    MolniyaOrbit,
    RepeatGroundTrackOrbit,
    SunSynchronousOrbit,
    TundraOrbit,
    TwoLineElements,
)
from .resources import InstrumentPerformance, InstrumentResources, SpacecraftResources
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
    "AllReceivers",
    "AllTransmitters",
    "Band",
    "DataFlow",
    "DataGeneration",
    "LinkMode",
    "Receiver",
    "RelayCommunications",
    "ShannonModem",
    "SpacecraftCommunications",
    "StationCommunications",
    "Transmitter",
    "ActiveWindow",
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
    "InSituInstrument",
    "Instrument",
    "InstrumentPerformance",
    "InstrumentResources",
    "KeplerianOrbit",
    "MOGConstellation",
    "MolniyaOrbit",
    "RepeatGroundTrackOrbit",
    "Point",
    "PointedInstrument",
    "RadarBand",
    "RadarStation",
    "Relay",
    "SOCConstellation",
    "Satellite",
    "SpacecraftResources",
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
