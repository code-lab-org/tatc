"""
Defines object schemas.
"""

from .architecture import Architecture
from .instrument import Instrument, PointedInstrument
from .orbit import (
    TwoLineElements,
    CircularOrbit,
    SunSynchronousOrbit,
    KeplerianOrbit,
    MolniyaOrbit,
    TundraOrbit,
)
from .point import Point, GroundStation
from .satellite import (
    Satellite,
    TrainConstellation,
    WalkerConfiguration,
    WalkerConstellation,
    MOGConstellation,
    SOCConstellation,
)
