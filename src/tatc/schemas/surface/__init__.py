"""
Object schemas for surface objects.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .platform import SurfacePlatform
from .point import Point
from .radar import RadarBand, RadarStation
from .station import GroundStation
from .terrain import TerrainMask

AllSurfaceObjects = Point | GroundStation | RadarStation | SurfacePlatform

__all__ = [
    "AllSurfaceObjects",
    "GroundStation",
    "Point",
    "RadarBand",
    "RadarStation",
    "SurfacePlatform",
    "TerrainMask",
]
