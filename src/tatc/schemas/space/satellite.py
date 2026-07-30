"""
Object schema for satellites.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Literal

from pydantic import Field

from ..orbit import AllOrbits
from .base import SpaceSystem


class Satellite(SpaceSystem):
    """
    Single satellite.
    """

    type: Literal["satellite"] = Field(
        default="satellite", description="Space system type discriminator."
    )
    orbit: AllOrbits = Field(..., description="Orbit specification.")
