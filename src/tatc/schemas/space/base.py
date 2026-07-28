"""
Base object schemas for space systems.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from pydantic import BaseModel, Field

from ..instrument import AllInstruments, Instrument


class SpaceSystem(BaseModel):
    """
    Base class for space systems.
    """

    name: str = Field(
        ...,
        description="Space system name.",
        examples=["International Space Station"],
    )
    instruments: list[AllInstruments] = Field(
        default=[Instrument()], min_length=1, description="List of assigned instruments."
    )
