"""
Object schemas for satellite orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
from pydantic import Field
from typing_extensions import Literal

from ... import constants, utils
from .molniya import MolniyaOrbit


class TundraOrbit(MolniyaOrbit):
    """
    Orbit defined by Tundra parameters.
    """

    type: Literal["tundra"] = Field("tundra", description="Orbit type discriminator.")

    def get_orbit_period(self) -> timedelta:
        """
        Gets the orbit period defined to be approximately 1 sidereal day.

        Returns:
            timedelta: the orbit period
        """
        # TODO this needs to be corrected to account for J2 effects
        return timedelta(seconds=constants.EARTH_SIDEREAL_DAY_S)

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> TundraOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            Tundra Orbit: the derived orbit
        """
        true_anomaly = utils.orbital.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360),
            eccentricity=self.get_eccentricity(),
        )
        raan = np.mod(self.right_ascension_ascending_node + delta_raan, 360)
        return TundraOrbit(
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            right_ascension_ascending_node=raan,
            northern_coverage=self.northern_coverage,
        )
