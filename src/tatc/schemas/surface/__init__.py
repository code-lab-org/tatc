"""
Object schemas for surface objects.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .point import Point
from .radar import RadarStation, TerrainMask
from .station import GroundStation

AllSurfaceObjects = Point | GroundStation | RadarStation

__all__ = [
    "AllSurfaceObjects",
    "GroundStation",
    "Point",
    "RadarStation",
    "TerrainMask",
]
