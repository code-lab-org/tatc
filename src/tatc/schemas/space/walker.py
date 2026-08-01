"""
Object schemas for Walker constellations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import copy
import math
from enum import Enum
from typing import Literal

import numpy as np
from pydantic import Field, model_validator

from tatc.utils.formatting import zero_pad

from ..orbit import AllOrbits
from .base_constellation import BaseConstellation
from .satellite import Satellite


class WalkerConfiguration(str, Enum):
    """
    Enumeration of different Walker constellation configurations.
    """

    DELTA = "delta"
    STAR = "star"


class WalkerConstellation(BaseConstellation):
    """
    A constellation that arranges member satellites following the Walker pattern.
    """

    type: Literal["walker"] = Field(
        default="walker", description="Space system type discriminator."
    )
    configuration: WalkerConfiguration = Field(
        default=WalkerConfiguration.DELTA, description="Walker configuration."
    )
    orbit: AllOrbits = Field(..., description="Lead orbit for this constellation.")
    number_satellites: int = Field(
        default=1, description="Number of satellites in the constellation.", ge=1
    )
    number_planes: int = Field(
        default=1,
        description="The number of equally-spaced planes in a Walker Delta "
        + "constellation. Ranges from 1 to (number of satellites).",
        ge=1,
    )
    relative_spacing: int = Field(
        default=0,
        description="Relative spacing of satellites between planes for a Walker Delta "
        + "constellation. Ranges from 0 for equal mean anomaly to "
        + "(number of planes) - 1. For example, `relative_spacing=1` "
        + "means the mean anomaly is shifted by `360/number_satellites` "
        + "between adjacent planes.",
        ge=0,
    )

    @model_validator(mode="after")
    def number_planes_le_number_satellites(self) -> WalkerConstellation:
        """
        Validates the number of planes given the number of satellites.
        """
        if (
            self.number_planes is not None
            and self.number_satellites is not None
            and self.number_planes > self.number_satellites
        ):
            raise ValueError("number planes exceeds number satellites")
        return self

    @model_validator(mode="after")
    def relative_spacing_lt_number_planes(self) -> WalkerConstellation:
        """
        Validates the relative spacing given the number of planes.
        """
        if (
            self.relative_spacing is not None
            and self.number_planes is not None
            and self.relative_spacing >= self.number_planes
        ):
            raise ValueError("relative spacing exceeds number planes - 1")
        return self

    def get_satellites_per_plane(self) -> int:
        """
        Gets the (max) number of satellites per plane.

        Returns:
            int: number of satellites per plane
        """
        return math.ceil(self.number_satellites / self.number_planes)

    def get_delta_mean_anomaly_within_planes(self) -> float:
        """
        Gets the difference in mean anomaly (decimal degrees) for adjacent
        member satellites within a single plane.

        Returns:
            float: difference in mean anomaly
        """
        return 360 / self.get_satellites_per_plane()

    def get_delta_mean_anomaly_between_planes(self) -> float:
        """
        Gets the difference in mean anomaly (decimal degrees) for adjacent
        member satellites between adjacent planes.

        Returns:
            float: difference in mean anomaly
        """
        return 360 * self.relative_spacing / self.number_satellites

    def get_delta_raan_between_planes(self) -> float:
        """
        Gets the difference in right ascension of ascending node (decimal
        degrees) for adjacent planes of member satellites.

        Returns:
            float: difference in right ascension of ascending node
        """
        if self.configuration == WalkerConfiguration.DELTA:
            return 360 / self.number_planes
        return 180 / self.number_planes

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
                    np.mod(i, self.get_satellites_per_plane())
                    * self.get_delta_mean_anomaly_within_planes()
                    + (i // self.get_satellites_per_plane())
                    * self.get_delta_mean_anomaly_between_planes(),
                    (i // self.get_satellites_per_plane())
                    * self.get_delta_raan_between_planes(),
                ),
                instruments=copy.deepcopy(self.instruments),
            )
            for i in range(self.number_satellites)
        ]
