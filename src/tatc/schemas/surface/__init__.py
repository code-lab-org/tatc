"""
Object schemas for surface objects.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .point import Point
from .radar import RadarBand, RadarStation, TerrainMask
from .station import GroundStation

AllSurfaceObjects = Point | GroundStation | RadarStation

__all__ = [
    "AllSurfaceObjects",
    "GroundStation",
    "Point",
    "RadarBand",
    "RadarStation",
    "TerrainMask",
]
