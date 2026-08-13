"""
Base object schemas for constellations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from .base import SpaceSystem
from .satellite import Satellite


class BaseConstellation(SpaceSystem):
    """
    Base class for constellations.
    """

    def generate_members(self) -> list[Satellite]:
        """
        Generates the members of the constellation.

        Returns:
            list[Satellite]: The list of generated members.
        """
        raise NotImplementedError(
            "generate_members() must be implemented in subclasses."
        )
