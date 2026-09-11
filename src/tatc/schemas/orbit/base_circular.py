"""
Base object schemas for circular orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from pydantic import AliasChoices, ConfigDict, Field

from ... import constants
from .base import OrbitBase


class CircularOrbitBase(OrbitBase):
    """
    Base class for circular orbits.
    """

    model_config = ConfigDict(populate_by_name=True)

    mean_altitude: float = Field(
        ...,
        description="Mean altitude (meters).",
        ge=0,
        validation_alias=AliasChoices("mean_altitude", "altitude"),
    )

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis.

        Returns:
            float: the semimajor axis (meters)
        """
        return constants.EARTH_MEAN_RADIUS + self.mean_altitude

    def get_mean_altitude(self) -> float:
        """
        Gets the mean altitude.

        Returns:
            float: the mean altitude (meters)
        """
        return self.mean_altitude

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity, always 0 for a circular orbit.

        Returns:
            float: the eccentricity
        """
        return 0

    def get_perigee_argument(self) -> float:
        """
        Gets the perigee argument, always 0 for a circular orbit (there is
        no perigee to reference).

        Returns:
            float: the perigee argument
        """
        return 0
