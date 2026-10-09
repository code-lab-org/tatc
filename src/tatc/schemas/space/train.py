"""
Object schema for train constellations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import copy
from datetime import timedelta
from typing import Literal

import numpy as np
from pydantic import Field

from ... import constants
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

    def _get_secular_rates(self) -> tuple[float, float]:
        """
        Gets the secular rates (degrees/second) at which the lead orbit, as
        propagated (by SGP4, including the Earth's oblateness), advances
        along its orbit and its right ascension of ascending node
        precesses. The rate along the orbit is that of the mean anomaly or,
        for a near-circular orbit (eccentricity below 0.01), whose perigee
        is not meaningful, that of the argument of latitude (the mean
        anomaly plus the argument of perigee).

        Returns:
            tuple[float, float]: the rates along the orbit and of the right
                ascension of ascending node
        """
        satrec = self.orbit.to_gp_orbit().elements[0].to_satrec()
        along = satrec.mdot + (
            satrec.argpdot if self.orbit.get_eccentricity() < 0.01 else 0
        )
        # sgp4 rates are in radians per minute
        return np.degrees(along) / 60, np.degrees(satrec.nodedot) / 60

    def get_delta_mean_anomaly(self) -> float:
        """
        Gets the difference in mean anomaly (decimal degrees) for adjacent
        member satellites: each trails the one ahead of it along the orbit
        by the interval, at the lead orbit's propagated rate.

        Returns:
            float: the difference in mean anomaly
        """
        along, _ = self._get_secular_rates()
        return -along * self.interval.total_seconds()

    def get_delta_raan(self) -> float:
        """
        Gets the difference in right ascension of ascending node (decimal
        degrees) for adjacent member satellites. When repeating the ground
        track, each satellite occupies, in the Earth-fixed frame, the
        position of the one ahead of it one interval earlier: its ascending
        node is shifted by the Earth's rotation (relative to inertial space,
        i.e. the sidereal day) during the interval, less the precession of
        the orbit plane during the interval. (For a sun-synchronous orbit,
        an interval of whole days therefore keeps the satellites in the same
        plane.)

        Returns:
            float: the difference in right ascension of ascending node
        """
        if self.repeat_ground_track:
            _, precession = self._get_secular_rates()
            return (
                360 / constants.EARTH_SIDEREAL_DAY_S - precession
            ) * self.interval.total_seconds()
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
