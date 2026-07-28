"""
Base object schemas for Molniya and Tundra orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
from pydantic import Field

from tatc import config, constants

from ... import utils
from .base import OrbitBase
from .gp import GeneralPerturbationsOrbit
from .keplerian import KeplerianOrbit


class MolniyaTundraOrbitBase(OrbitBase):
    """
    Base class for Molniya and Tundra orbits.
    """

    northern_coverage: bool = Field(
        default=True, description="True, if the orbit generates northern hemisphere coverage."
    )
    perigee_altitude: float = Field(..., description="Perigee altitude (meters).", ge=0)
    right_ascension_ascending_node: float = Field(
        default=0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
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
        Gets the orbit period (seconds).

        Returns:
            float: the orbit period
        """
        raise NotImplementedError("get_orbit_period() must be implemented in subclasses.")

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
            gp_orbit = KeplerianOrbit(
                semimajor_axis=self.get_semimajor_axis(),
                inclination=self.get_inclination(),
                right_ascension_ascending_node=self.right_ascension_ascending_node,
                true_anomaly=self.true_anomaly,
                epoch=self.epoch,
                eccentricity=self.get_eccentricity(),
                perigee_argument=self.get_perigee_argument(),
            ).to_gp_orbit()
            self.__dict__["gp_orbit"] = gp_orbit # type: ignore
        return gp_orbit
