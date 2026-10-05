"""
Base object schemas for Molniya and Tundra orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import Field, model_validator
from sgp4.api import WGS72, Satrec

from ... import constants, utils
from .base import OrbitBase


class MolniyaTundraOrbitBase(OrbitBase):
    """
    Base class for Molniya and Tundra orbits.
    """

    northern_coverage: bool = Field(
        default=True,
        description="True, if the orbit generates northern hemisphere coverage.",
    )
    perigee_altitude: float = Field(..., description="Perigee altitude (meters).", ge=0)
    right_ascension_ascending_node: float = Field(
        default=0,
        description="Right ascension of ascending node (degrees).",
        ge=0,
        lt=360,
    )

    def get_inclination(self) -> float:
        """
        Gets the inclination (degrees) of the frozen orbit.

        Returns:
            float: the inclination
        """
        return constants.EARTH_J2_CRITICAL_INCLINATION

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node.

        Returns:
            float: the right ascension of ascending node (degrees)
        """
        return self.right_ascension_ascending_node

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
        raise NotImplementedError(
            "get_orbit_period() must be implemented in subclasses."
        )

    def _compute_j2_corrected_orbit_period(self, target_period_s: float) -> timedelta:
        """
        Corrects a nominal repeat-ground-track target period (e.g. exactly
        half or one full sidereal day, ignoring perturbations) to account
        for Earth's J2 oblateness perturbation to the true rate at which
        mean anomaly advances. Subclasses call this from get_orbit_period()
        with their own target (half a sidereal day for Molniya, one
        sidereal day for Tundra).

        The correction is computed in a single pass (not iteratively):
        the target period first implies a "naive" semimajor axis and
        eccentricity via simple two-body Kepler's third law, which are
        then used to evaluate the J2 mean anomaly rate correction
        (compute_j2_mean_motion_rate). Since that correction is itself
        tiny (on the order of 1e-4 relative to the mean motion), using
        the naive semimajor axis/eccentricity to evaluate it introduces
        only a negligible (second-order, ~1e-8 relative) residual error,
        rather than requiring a full iterative solve.

        The returned period, when passed through get_semimajor_axis()'s
        existing (unmodified) two-body Kepler's third law formula,
        yields a semimajor axis whose true (J2-corrected) mean anomaly
        rate matches the target repeat-ground-track requirement.

        Args:
            target_period_s (float): The nominal, uncorrected target orbit
                period (seconds), ignoring J2 perturbation.

        Returns:
            timedelta: The J2-corrected orbit period.
        """
        naive_semimajor_axis = np.cbrt(
            constants.EARTH_MU * target_period_s**2 / (4 * np.pi**2)
        )
        naive_eccentricity = (
            1
            - (constants.EARTH_MEAN_RADIUS + self.perigee_altitude)
            / naive_semimajor_axis
        )
        mean_motion_correction = utils.orbital.compute_j2_mean_motion_rate(
            naive_semimajor_axis, self.get_inclination(), naive_eccentricity
        )
        target_mean_motion = 360 / target_period_s
        return timedelta(seconds=360 / (target_mean_motion - mean_motion_correction))

    def _compute_repeat_orbit_period(self, revolutions_per_day: int) -> timedelta:
        """
        Computes the orbit period for which the ground track repeats with the
        given number of revolutions per day, as propagated by SGP4. Starting
        from the J2-corrected period for the sidereal day
        (`_compute_j2_corrected_orbit_period`), the period is refined using
        SGP4's own secular rates (including the Earth's oblateness) so that
        the satellite completes the given number of anomalistic revolutions
        (perigee to perigee) while the Earth rotates once relative to the
        orbit's precessing ascending node. Without the node's precession,
        the apogee longitudes of a Molniya orbit drift westward by about 0.1
        degrees per day. The result is cached.

        Args:
            revolutions_per_day (int): The number of revolutions per day.

        Returns:
            timedelta: The orbit period.
        """
        cached = self.__dict__.get("repeat_orbit_period")
        if cached is not None:
            return cached
        period = self._compute_j2_corrected_orbit_period(
            constants.EARTH_SIDEREAL_DAY_S / revolutions_per_day
        ).total_seconds()
        epoch = (self.epoch - datetime(1949, 12, 31, tzinfo=timezone.utc)) / timedelta(
            days=1
        )
        # Earth's rotation rate (radians/minute), as SGP4's rates
        rotation_rate = 2 * np.pi / constants.EARTH_SIDEREAL_DAY_S * 60
        for _ in range(5):
            semimajor_axis = np.cbrt(constants.EARTH_MU * period**2 / (4 * np.pi**2))
            eccentricity = (
                1
                - (constants.EARTH_MEAN_RADIUS + self.perigee_altitude) / semimajor_axis
            )
            satrec = Satrec()
            satrec.sgp4init(
                WGS72,
                "i",
                0,
                epoch,
                0.0,
                0.0,
                0.0,
                eccentricity,
                np.radians(self.get_perigee_argument()),
                np.radians(self.get_inclination()),
                0.0,
                2 * np.pi / (period / 60),
                np.radians(self.right_ascension_ascending_node),
            )
            target = revolutions_per_day * (rotation_rate - satrec.nodedot)
            period *= satrec.mdot / target
        result = timedelta(seconds=period)
        self.__dict__["repeat_orbit_period"] = result  # type: ignore
        return result

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters).

        Returns:
            float: the semimajor axis
        """
        return np.cbrt(
            (constants.EARTH_MU * self.get_orbit_period().total_seconds() ** 2)
            / (4 * np.pi**2)
        )

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity (float between 0 and 1).

        Returns:
            float: the eccentricity
        """
        return (
            1
            - (constants.EARTH_MEAN_RADIUS + self.perigee_altitude)
            / self.get_semimajor_axis()
        )

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.orbital.true_anomaly_to_mean_anomaly(
            self.true_anomaly, self.get_eccentricity()
        )

    @model_validator(mode="after")
    def _validate_eccentricity(self) -> MolniyaTundraOrbitBase:
        """
        Validates that perigee_altitude, combined with this orbit's fixed
        orbit period, yields a physically valid elliptical eccentricity in
        [0, 1). A perigee_altitude above the Kepler-derived ceiling
        implied by the fixed period would otherwise silently produce a
        negative eccentricity. Not implemented on the bare
        MolniyaTundraOrbitBase class (get_orbit_period is abstract there),
        so this is a no-op until a concrete subclass defines the period.
        """
        try:
            eccentricity = self.get_eccentricity()
        except NotImplementedError:
            return self
        if not 0 <= eccentricity < 1:
            raise ValueError(
                f"perigee_altitude={self.perigee_altitude} is invalid for this "
                f"orbit's fixed period: implies eccentricity={eccentricity}, "
                "which is outside the valid elliptical range [0, 1)."
            )
        return self
