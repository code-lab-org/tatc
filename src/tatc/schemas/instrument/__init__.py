"""
Object schemas for instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .simple import Instrument
from .pointed import PointedInstrument

AllInstruments = Instrument | PointedInstrument

__all__ = ["AllInstruments", "Instrument", "PointedInstrument"]
