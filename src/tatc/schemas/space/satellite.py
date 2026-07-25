"""
Object schema for satellites.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""
from __future__ import annotations

from pydantic import Field
from typing_extensions import Literal

from .base import SpaceSystem


class Satellite(SpaceSystem):
    """
    Single satellite.
    """

    type: Literal["satellite"] = Field(
        "satellite", description="Space system type discriminator."
    )

    def generate_members(self) -> list[Satellite]:
        """
        Generate space system member satellites (returns a list containing this satellite).

        Returns:
            list[Satellite]: the member satellites
        """
        return [self]
