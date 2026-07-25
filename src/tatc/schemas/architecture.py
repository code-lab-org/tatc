"""
Object schemas for mission architectures.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""
from __future__ import annotations

from pydantic import BaseModel, Field

from .space import Satellite
from .surface import GroundStation


class Architecture(BaseModel):
    """
    Mission architecture.
    """

    name: str = Field(..., description="Name of this mission.")
    satellites: list[Satellite] = Field([], description="List of member satellites.")
    stations: list[GroundStation] = Field(
        [], description="List of member ground stations."
    )
