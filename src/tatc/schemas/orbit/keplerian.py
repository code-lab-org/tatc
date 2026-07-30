"""
Object schemas for Keplerian orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta
from typing import Literal

import numpy as np
from pydantic import Field

from ... import config, utils
from .base import OrbitBase
from .gp import GeneralPerturbationsElements, GeneralPerturbationsOrbit


class KeplerianOrbit(OrbitBase):
    """
    Orbit specification using Keplerian elements for elliptical motion.
    """

    type: Literal["keplerian"] = Field(
        default="keplerian", description="Orbit type discriminator."
    )
    semimajor_axis: float = Field(..., description="Semimajor axis (meters).")
    inclination: float = Field(0, description="Inclination (degrees).", ge=0, lt=180)
    right_ascension_ascending_node: float = Field(
        0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )
    eccentricity: float = Field(0, description="Eccentricity.", ge=0)
    perigee_argument: float = Field(
        0, description="Perigee argument (degrees).", ge=0, lt=360
    )

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.orbital.true_anomaly_to_mean_anomaly(
            self.true_anomaly, self.eccentricity
        )

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion (degrees/second).

        Returns:
            float: the mean motion
        """
        return utils.orbital.semimajor_axis_to_mean_motion(self.semimajor_axis)

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period.

        Returns:
            timedelta: the orbit period
        """
        return timedelta(
            seconds=utils.orbital.semimajor_axis_to_orbit_period(self.semimajor_axis)
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

    def to_gp_orbit(self, lazy_load: bool | None = None) -> GeneralPerturbationsOrbit:
        """
        Converts this orbit to a general perturbations orbit representation.

        Args:
            lazy_load (bool | None): True, if this gp orbit should be lazy-loaded.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        if lazy_load is None:
            lazy_load = config.rc.gp_orbit_lazy_load
        if lazy_load:
            gp_orbit = self.__dict__.get("gp_orbit")
        else:
            gp_orbit = None
        if gp_orbit is None:
            gp_orbit = GeneralPerturbationsOrbit(
                elements=[
                    GeneralPerturbationsElements(
                        epoch=self.epoch,
                        mean_motion=self.get_mean_motion(),
                        eccentricity=self.eccentricity,
                        inclination=self.inclination,
                        ra_of_asc_node=self.right_ascension_ascending_node,
                        arg_of_pericenter=self.perigee_argument,
                        mean_anomaly=self.get_mean_anomaly(),
                    )
                ]
            )
            self.__dict__["gp_orbit"] = gp_orbit  # type: ignore
        return gp_orbit
