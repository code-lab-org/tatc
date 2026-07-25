"""
Object schemas for Molniya orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
from pydantic import Field
from typing_extensions import Literal

from ... import config, constants, utils
from .base import OrbitBase
from .gp import GeneralPerturbationsOrbit
from .keplerian import KeplerianOrbit


class MolniyaOrbit(OrbitBase):
    """
    Orbit defined by Molniya parameters.
    """

    type: Literal["molniya"] = Field("molniya", description="Orbit type discriminator.")
    northern_coverage: bool = Field(
        True, description="True, if the orbit generates northern hemisphere coverage."
    )
    perigee_altitude: float = Field(..., description="Perigee altitude (meters).", ge=0)
    right_ascension_ascending_node: float = Field(
        0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )

    def get_inclination(self) -> float:
        """
        Gets the inclination (degrees) of the frozen orbit.

        Returns:
            float: the inclination
        """
        return constants.EARTH_J2_CRITICAL_INCLINATION
    
    def get_perigee_argument(self) -> float:
        """
        Gets the perigee argument (degrees) of the frozen orbit.

        Returns:
            float: the perigee argument
        """
        return 270 if self.northern_coverage else 90

    def get_orbit_period(self) -> timedelta:
        """
        Gets the orbit period defined to be approximately half a sidereal day.

        Returns:
            timedelta: the orbit period
        """
        # TODO this needs to be corrected to account for J2 effects
        return timedelta(seconds=constants.EARTH_SIDEREAL_DAY_S / 2)

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters).

        Returns:
            float: the semimajor axis
        """
        return np.cbrt(
            (constants.EARTH_MU * self.get_orbit_period().total_seconds()**2) / (4 * np.pi**2)
        )

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity (float between 0 and 1).

        Returns:
            float: the eccentricity
        """
        return 1 - (constants.EARTH_MEAN_RADIUS + self.perigee_altitude) / self.get_semimajor_axis()

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.orbital.true_anomaly_to_mean_anomaly(
            self.true_anomaly, self.get_eccentricity()
        )

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
            right_ascension_ascending_node=raan,
            northern_coverage=self.northern_coverage,
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
            lazy_load = config.rc.orbit_tle_lazy_load
        if lazy_load:
            gp_orbit = self.__dict__.get("gp_orbit")
        else:
            gp_orbit = None
        if gp_orbit is None:
            gp_orbit = KeplerianOrbit(
                semimajor_axis=self.get_semimajor_axis(),
                inclination=self.get_inclination(),
                right_ascension_ascending_node=self.right_ascension_ascending_node,
                true_anomaly=self.true_anomaly,
                epoch=self.epoch,
                eccentricity=self.get_eccentricity(),
                perigee_argument=self.get_perigee_argument(),
            ).to_gp_orbit()
            self.__dict__["gp_orbit"] = gp_orbit
        return gp_orbit
