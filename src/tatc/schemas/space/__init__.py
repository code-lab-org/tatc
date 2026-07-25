"""
Object schemas for space objects.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .mog import MOGConstellation
from .satellite import Satellite
from .soc import SOCConstellation
from .train import TrainConstellation
from .walker import WalkerConfiguration, WalkerConstellation

AllSpaceObjects = (
    MOGConstellation
    | Satellite
    | SOCConstellation
    | TrainConstellation
    | WalkerConstellation
)

__all__ = [
    "AllSpaceObjects",
    "MOGConstellation",
    "SOCConstellation",
    "Satellite",
    "TrainConstellation",
    "WalkerConfiguration",
    "WalkerConstellation",
]
