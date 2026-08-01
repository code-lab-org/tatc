"""
Object schemas for geostationary orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import Field

from ... import constants, utils
from .base_circular import CircularOrbitBase

# mean altitude (meters) yielding an orbit period of exactly one sidereal
# day, via Kepler's third law
_GEOSTATIONARY_MEAN_ALTITUDE = (
    utils.orbital.mean_motion_to_semimajor_axis(360 / constants.EARTH_SIDEREAL_DAY_S)
    - constants.EARTH_MEAN_RADIUS
)


class GeostationaryOrbit(CircularOrbitBase):
    """
    Orbit defined by geostationary parameters: a fixed longitude rather
    than a right ascension of ascending node.
    """

    type: Literal["geostationary"] = Field(
        default="geostationary", description="Orbit type discriminator."
    )
    mean_altitude: float = Field(
        default=_GEOSTATIONARY_MEAN_ALTITUDE,
        description="Mean altitude (meters).",
        ge=0,
    )
    inclination: float = Field(
        default=0, description="Inclination (degrees).", ge=0, lt=180
    )
    longitude: float = Field(
        ...,
        description="Longitude (decimal degrees) of the sub-satellite point "
        + "in the WGS 84 coordinate system.",
        ge=-180,
        le=180,
    )

    def get_inclination(self) -> float:
        """
        Gets the inclination.

        Returns:
            float: the inclination (degrees)
        """
        return self.inclination

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node, derived from this
        orbit's fixed longitude and Earth's rotation angle (Greenwich
        Apparent Sidereal Time) at epoch, since a geostationary satellite
        remains above a fixed Earth-relative longitude rather than a
        fixed inertial-frame RAAN. As with SunSynchronousOrbit's
        equator-crossing-time convention, this is computed independently
        of true_anomaly.

        Returns:
            float: the right ascension of ascending node (degrees)
        """
        epoch_time = constants.timescale.from_datetime(self.epoch)
        earth_rotation_angle = epoch_time.gast * 15  # hours -> degrees
        return (earth_rotation_angle + self.longitude) % 360

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> GeostationaryOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and
        longitude (equivalent, for a geostationary orbit, to a right
        ascension of ascending node perturbation, since both epoch and
        Earth's rotation angle stay fixed).

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            GeostationaryOrbit: the derived orbit
        """
        true_anomaly = utils.orbital.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360)
        )
        longitude = ((self.longitude + delta_raan + 180) % 360) - 180
        return GeostationaryOrbit(
            mean_altitude=self.mean_altitude,
            inclination=self.inclination,
            longitude=longitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
        )
