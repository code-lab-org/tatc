"""
Object schemas for circular orbits.

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


class CircularOrbit(OrbitBase):
    """
    Orbit specification using Keplerian elements for elliptical motion -- circular motion case.
    """

    type: Literal["circular"] = Field(
        "circular", description="Orbit type discriminator."
    )
    mean_altitude: float = Field(..., description="Mean altitude (meters).")
    inclination: float = Field(0, description="Inclination (degrees).", ge=0, lt=180)
    right_ascension_ascending_node: float = Field(
        0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )
    
    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis.

        Returns:
            float: the semimajor axis (meters)
        """
        return constants.EARTH_MEAN_RADIUS + self.mean_altitude

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion.

        Returns:
            float: the mean motion (revolutions per day)
        """
        return utils.orbital.semimajor_axis_to_mean_motion(self.get_semimajor_axis())

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period.

        Returns:
            timedelta: the orbit period
        """
        return timedelta(seconds=utils.orbital.semimajor_axis_to_orbit_period(self.get_semimajor_axis()))

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> CircularOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            CircularOrbit: the derived orbit
        """
        true_anomaly = utils.orbital.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360)
        )
        raan = np.mod(self.right_ascension_ascending_node + delta_raan, 360)
        return CircularOrbit(
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
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
                true_anomaly=self.true_anomaly,
                epoch=self.epoch,
                inclination=self.inclination,
                right_ascension_ascending_node=self.right_ascension_ascending_node,
                eccentricity=0,
                perigee_argument=0,
            ).to_gp_orbit()
            self.__dict__["gp_orbit"] = gp_orbit
        return gp_orbit