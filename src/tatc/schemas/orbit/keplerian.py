"""
Object schemas for Keplerian orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import Field

from ... import utils
from .base import OrbitBase
from .gp import GeneralPerturbationsOrbit
from .gp_elements import GeneralPerturbationsElements


class KeplerianOrbit(OrbitBase):
    """
    Orbit specification using Keplerian elements for elliptical motion.
    """

    type: Literal["keplerian"] = Field(
        default="keplerian", description="Orbit type discriminator."
    )
    semimajor_axis: float = Field(..., description="Semimajor axis (meters).", gt=0)
    inclination: float = Field(0, description="Inclination (degrees).", ge=0, lt=180)
    right_ascension_ascending_node: float = Field(
        0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )
    eccentricity: float = Field(0, description="Eccentricity.", ge=0, lt=1)
    perigee_argument: float = Field(
        0, description="Perigee argument (degrees).", ge=0, lt=360
    )

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis.

        Returns:
            float: the semimajor axis (meters)
        """
        return self.semimajor_axis

    def get_inclination(self) -> float:
        """
        Gets the inclination.

        Returns:
            float: the inclination (degrees)
        """
        return self.inclination

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node.

        Returns:
            float: the right ascension of ascending node (degrees)
        """
        return self.right_ascension_ascending_node

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity.

        Returns:
            float: the eccentricity
        """
        return self.eccentricity

    def get_perigee_argument(self) -> float:
        """
        Gets the perigee argument.

        Returns:
            float: the perigee argument (degrees)
        """
        return self.perigee_argument

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.orbital.true_anomaly_to_mean_anomaly(
            self.true_anomaly, self.eccentricity
        )

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> KeplerianOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            KeplerianOrbit: the derived orbit
        """
        true_anomaly = utils.orbital.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360),
            eccentricity=self.eccentricity,
        )
        raan = np.mod(self.right_ascension_ascending_node + delta_raan, 360)
        return KeplerianOrbit(
            semimajor_axis=self.semimajor_axis,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
            eccentricity=self.eccentricity,
            perigee_argument=self.perigee_argument,
        )

    def _compute_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Computes a general perturbations orbit representation of this
        orbit.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        return GeneralPerturbationsOrbit(
            elements=[
                GeneralPerturbationsElements(
                    epoch=self.get_epoch(),
                    mean_motion=self.get_mean_motion(),
                    eccentricity=self.get_eccentricity(),
                    inclination=self.get_inclination(),
                    ra_of_asc_node=self.get_right_ascension_ascending_node(),
                    arg_of_pericenter=self.get_perigee_argument(),
                    mean_anomaly=self.get_mean_anomaly(),
                )
            ]
        )
