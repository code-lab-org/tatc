"""
Defines object schemas.
"""

from .architecture import Architecture
from .instrument import Instrument
from .orbit import TwoLineElements, CircularOrbit, SunSynchronousOrbit, KeplerianOrbit
from .point import Point, GroundStation
from .satellite import (
    Satellite,
    TrainConstellation,
    WalkerConfiguration,
    WalkerConstellation,
)
