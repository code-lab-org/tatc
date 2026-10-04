"""
Object schemas for instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .simple import Instrument
from .pointed import PointedInstrument
from .conical import ConicalInstrument

AllInstruments = Instrument | PointedInstrument | ConicalInstrument

__all__ = ["AllInstruments", "ConicalInstrument", "Instrument", "PointedInstrument"]
