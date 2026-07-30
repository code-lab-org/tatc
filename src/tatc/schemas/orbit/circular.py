"""
Object schemas for circular orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import Field

from ... import config, utils
from .base_circular import CircularOrbitBase
from .gp import GeneralPerturbationsOrbit
from .keplerian import KeplerianOrbit


class CircularOrbit(CircularOrbitBase):
    """
    Orbit specification using Keplerian elements for elliptical motion -- circular motion case.
    """

    type: Literal["circular"] = Field(
        default="circular", description="Orbit type discriminator."
    )
    inclination: float = Field(default=0, description="Inclination (degrees).", ge=0, lt=180)
    right_ascension_ascending_node: float = Field(
        default=0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )

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
            mean_altitude=self.mean_altitude,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
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
                true_anomaly=self.true_anomaly,
                epoch=self.epoch,
                inclination=self.inclination,
                right_ascension_ascending_node=self.right_ascension_ascending_node,
                eccentricity=0,
                perigee_argument=0,
            ).to_gp_orbit()
            self.__dict__["gp_orbit"] = gp_orbit  # type: ignore
        return gp_orbit
