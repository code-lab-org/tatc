"""
Object schemas for instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .nadir import Instrument
from .off_nadir import PointedInstrument

AllInstruments = Instrument | PointedInstrument

__all__ = ["AllInstruments", "Instrument", "PointedInstrument"]
