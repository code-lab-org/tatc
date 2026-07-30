"""
Object schema for train constellations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import copy
from datetime import timedelta
from typing import Literal

from pydantic import Field

from ...utils.formatting import zero_pad
from ..orbit import AllOrbits
from .base_constellation import BaseConstellation
from .satellite import Satellite


class TrainConstellation(BaseConstellation):
    """
    A constellation that arranges member satellites in sequence.
    """

    type: Literal["train"] = Field(
        default="train", description="Space system type discriminator."
    )
    orbit: AllOrbits = Field(..., description="Lead orbit for this constellation.")
    number_satellites: int = Field(
        1, description="The count of the number of satellites.", ge=1
    )
    interval: timedelta = Field(
        ...,
        description="The local time interval between satellites in a train constellation.",
    )
    repeat_ground_track: bool = Field(
        True,
        description="True, if the train satellites should repeat the same ground track.",
    )

    def get_delta_mean_anomaly(self) -> float:
        """
        Gets the difference in mean anomaly (decimal degrees) for adjacent
        member satellites.

        Returns:
            float: the difference in mean anomaly
        """
        return -360 * self.interval / self.orbit.get_orbit_period()

    def get_delta_raan(self) -> float:
        """
        Gets the difference in right ascension of ascending node (decimal
        degrees) for adjacent member satellites.

        Returns:
            float: the difference in right ascension of ascending node
        """
        if self.repeat_ground_track:
            return 360 * (self.interval / timedelta(days=1))
        return 0

    def generate_members(self) -> list[Satellite]:
        """
        Generate space system member satellites.

        Returns:
            list[Satellite]: the member satellites
        """
        return [
            Satellite(
                name=zero_pad(self.name, self.number_satellites, i + 1),
                orbit=self.orbit.get_derived_orbit(
                    i * self.get_delta_mean_anomaly(), i * self.get_delta_raan()
                ),
                instruments=copy.deepcopy(self.instruments),
            )
            for i in range(self.number_satellites)
        ]
