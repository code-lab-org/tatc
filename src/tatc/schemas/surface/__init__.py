"""
Object schemas for surface objects.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .point import Point
from .station import GroundStation

AllSurfaceObjects = Point | GroundStation

__all__ = [
    "AllSurfaceObjects",
    "GroundStation",
    "Point",
]
