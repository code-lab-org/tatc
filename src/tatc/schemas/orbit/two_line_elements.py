"""
Object schema for orbits defined by two line elements (TLE).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Literal

from pydantic import Field, model_validator

from .base import OrbitBase
from .gp import GeneralPerturbationsOrbit


class TwoLineElements(OrbitBase):
    """
    Orbit defined by two line elements (TLE).
    """

    type: Literal["tle"] = Field(default="tle", description="Orbit type discriminator.")
    tle: tuple[str, str] = Field(
        ...,
        description="The two TLE lines.",
        examples=[
            (
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            )
        ],
    )

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters).

        Returns:
            float: the semimajor axis
        """
        return self.to_gp_orbit().get_semimajor_axis()

    def get_mean_altitude(self) -> float:
        """
        Gets the mean altitude (meters) above the WGS 84 mean radius.

        Returns:
            float: the mean altitude
        """
        return self.to_gp_orbit().get_mean_altitude()

    def get_altitude(self) -> float:
        """
        Gets the mean altitude (meters) above the WGS 84 mean radius.
        Alias for `get_mean_altitude()`, retained for backward
        compatibility with older scripts.

        Returns:
            float: the mean altitude
        """
        return self.get_mean_altitude()

    def get_inclination(self) -> float:
        """
        Gets the inclination (degrees).

        Returns:
            float: the inclination
        """
        return self.to_gp_orbit().get_inclination()

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity (float between 0 and 1).

        Returns:
            float: the eccentricity
        """
        return self.to_gp_orbit().get_eccentricity()

    def get_epoch(self) -> datetime:
        """
        Gets the epoch of the TLE.

        Returns:
            datetime: the epoch
        """
        return self.to_gp_orbit().get_epoch()

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion (revolutions per day).

        Returns:
            float: the mean motion
        """
        return self.to_gp_orbit().get_mean_motion()

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (degrees).

        Returns:
            float: the mean anomaly
        """
        return self.to_gp_orbit().get_mean_anomaly()

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period.

        Returns:
            timedelta: the orbit period
        """
        return self.to_gp_orbit().get_orbit_period()

    def get_true_anomaly(self) -> float:
        """
        Gets the true anomaly (degrees).

        Returns:
            float: the true anomaly
        """
        return self.to_gp_orbit().get_true_anomaly()

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node (degrees).

        Returns:
            float: the right ascension of ascending node
        """
        return self.to_gp_orbit().get_right_ascension_ascending_node()

    def get_perigee_argument(self) -> float:
        """
        Gets the argument of perigee (degrees).

        Returns:
            float: the argument of perigee
        """
        return self.to_gp_orbit().get_perigee_argument()

    def get_catalog_number(self) -> int:
        """
        Gets the NORAD catalog number.

        Returns:
            int: the NORAD catalog number
        """
        return self.to_gp_orbit().get_catalog_number()

    def get_bstar(self) -> float:
        """
        Gets the starred ballistic coefficient.

        Returns:
            float: the starred ballistic coefficient
        """
        return self.to_gp_orbit().get_bstar()

    def get_mean_motion_dot(self) -> float:
        """
        Gets the first derivative of mean motion (degrees/second^2).

        Returns:
            float: the first derivative of mean motion
        """
        return self.to_gp_orbit().get_mean_motion_dot()

    def get_mean_motion_ddot(self) -> float:
        """
        Gets the second derivative of mean motion (degrees/second^3).

        Returns:
            float: the second derivative of mean motion
        """
        return self.to_gp_orbit().get_mean_motion_ddot()

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> TwoLineElements:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            TwoLineElements: the derived orbit
        """
        derived_element = self.to_gp_orbit().get_derived_orbit(
            delta_mean_anomaly, delta_raan
        ).elements[0]
        return TwoLineElements(tle=derived_element.to_tle())

    def _compute_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Computes a general perturbations orbit representation of this
        orbit by parsing its TLE lines via SGP4.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        return GeneralPerturbationsOrbit.from_tle(list(self.tle))

    @model_validator(mode="after")
    def _validate_tle(self) -> TwoLineElements:
        """
        Validates the TLE lines by parsing them, raising if they cannot
        be converted to a general perturbations orbit. Reuses
        `to_gp_orbit()` (rather than calling `_compute_gp_orbit()`
        directly) so the parsed result is cached, avoiding a redundant
        SGP4 fit the first time this orbit's gp orbit is requested.
        """
        self.to_gp_orbit()
        return self
