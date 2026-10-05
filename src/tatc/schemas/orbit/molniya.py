"""
Object schemas for Molniya orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta
from typing import Literal

import numpy as np
from pydantic import Field

from ... import constants, utils
from .base_molniya_tundra import MolniyaTundraOrbitBase


class MolniyaOrbit(MolniyaTundraOrbitBase):
    """
    Orbit defined by Molniya parameters.
    """

    type: Literal["molniya"] = Field(
        default="molniya", description="Orbit type discriminator."
    )

    def get_orbit_period(self) -> timedelta:
        """
        Gets the orbit period, for which the ground track repeats twice daily and corrected for the Earth's oblateness
        (see `_compute_repeat_orbit_period`).

        Returns:
            timedelta: the orbit period
        """
        return self._compute_repeat_orbit_period(2)

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> MolniyaOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            Molniya Orbit: the derived orbit
        """
        true_anomaly = utils.orbital.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360),
            eccentricity=self.get_eccentricity(),
        )
        raan = np.mod(self.right_ascension_ascending_node + delta_raan, 360)
        return MolniyaOrbit(
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            perigee_altitude=self.perigee_altitude,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
            northern_coverage=self.northern_coverage,
        )
