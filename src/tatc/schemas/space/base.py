"""
Base object schemas for space systems.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from pydantic import BaseModel, Field

from ..instrument import AllInstruments, Instrument
from ..orbit import AllOrbits


class SpaceSystem(BaseModel):
    """
    Base class for space systems.
    """

    name: str = Field(
        ...,
        description="Space system name.",
        examples=["International Space Station"],
    )
    orbit: AllOrbits = Field(..., description="Orbit specification.")
    instruments: list[AllInstruments] = Field(
        [Instrument()], min_length=1, description="List of assigned instruments."
    )
