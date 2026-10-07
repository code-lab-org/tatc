"""
Base object schemas for Molniya and Tundra orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import Field, model_validator
from sgp4.api import WGS72, Satrec
from skyfield.api import wgs84
from typing_extensions import Self

from ... import constants, utils
from ...utils.cache import get_cached
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
    inclination: float = Field(
        default=constants.EARTH_J2_CRITICAL_INCLINATION,
        description="Inclination (degrees). Defaults to the critical "
        + "inclination (about 63.4 degrees), for which the argument of perigee "
        + "is frozen; at other inclinations, the argument of perigee (and so "
        + "the latitude of apogee) precesses as propagated.",
        ge=0,
        lt=180,
    )
    right_ascension_ascending_node: float = Field(
        default=0,
        description="Right ascension of ascending node (degrees).",
        ge=0,
        lt=360,
    )

    @classmethod
    def from_apogee_longitude(
        cls,
        apogee_longitude: float,
        perigee_altitude: float,
        inclination: float = constants.EARTH_J2_CRITICAL_INCLINATION,
        northern_coverage: bool = True,
        true_anomaly: float = 0,
        epoch: datetime = datetime(2020, 1, 1, tzinfo=timezone.utc),
    ) -> Self:
        """
        Creates an orbit whose apogee is over a longitude, rather than at a
        right ascension of ascending node. The ground track of a Tundra (or
        quasi-zenith) orbit is a figure-8 centered on its apogee longitude
        (and a Molniya orbit's apogees alternate between this longitude and
        the one 180 degrees away). The right ascension of ascending node is
        solved so that the first apogee at or after the epoch, as propagated
        (see `get_apogee_longitude`), is over the longitude.

        Args:
            apogee_longitude (float): Longitude (degrees) of the apogee's
                sub-satellite point in the WGS 84 coordinate system.
            perigee_altitude (float): Perigee altitude (meters).
            inclination (float): Inclination (degrees). Defaults to the
                critical inclination.
            northern_coverage (bool): True, if the orbit generates northern
                hemisphere coverage.
            true_anomaly (float): True anomaly (degrees) at epoch.
            epoch (datetime): Timestamp (epoch) of the initial orbital state.

        Returns:
            Self: the orbit (of the class on which this method is called)
        """
        orbit = None
        raan = 0.0
        # the apogee longitude moves with the right ascension of ascending
        # node, so a few corrections converge
        for _ in range(4):
            orbit = cls(
                perigee_altitude=perigee_altitude,
                inclination=inclination,
                northern_coverage=northern_coverage,
                true_anomaly=true_anomaly,
                epoch=epoch,
                right_ascension_ascending_node=raan,
            )
            error = (apogee_longitude - orbit.get_apogee_longitude() + 180) % 360 - 180
            if abs(error) < 1e-6:
                break
            raan = float((raan + error) % 360)
        return orbit  # type: ignore

    def get_apogee_longitude(self) -> float:
        """
        Gets the longitude (degrees) of the sub-satellite point in the WGS 84
        coordinate system at the first apogee at or after the epoch, as
        propagated by SGP4. For a Tundra (or quasi-zenith) orbit, the ground
        track is a figure-8 centered on this longitude.

        Returns:
            float: the apogee longitude (degrees, -180 to 180)
        """
        gp_orbit = self.to_gp_orbit()
        satrec = gp_orbit.elements[0].to_satrec()
        # minutes from epoch until the mean anomaly reaches 180 degrees
        minutes = ((np.pi - satrec.mo) % (2 * np.pi)) / satrec.mdot
        # refine the time of maximum radius from samples every 30 s
        offsets = minutes + np.arange(-30, 30.5, 0.5)
        track = gp_orbit.get_orbit_track(
            [self.epoch + timedelta(minutes=float(offset)) for offset in offsets],
            try_repeat=False,
        )
        radius = np.linalg.norm(track.position.m, axis=0)
        k = int(np.clip(np.argmax(radius), 1, len(radius) - 2))
        # vertex of the parabola through the samples around the maximum
        curvature = radius[k - 1] - 2 * radius[k] + radius[k + 1]
        shift = (
            0.5 * (radius[k - 1] - radius[k + 1]) / curvature if curvature < 0 else 0
        )
        apogee_time = self.epoch + timedelta(minutes=float(offsets[k] + 0.5 * shift))
        position = wgs84.subpoint_of(
            gp_orbit.get_orbit_track(apogee_time, try_repeat=False)
        )
        return float((position.longitude.degrees + 180) % 360 - 180)

    def get_inclination(self) -> float:
        """
        Gets the inclination (degrees).

        Returns:
            float: the inclination
        """
        return self.inclination

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node.

        Returns:
            float: the right ascension of ascending node (degrees)
        """
        return self.right_ascension_ascending_node

    def get_perigee_argument(self) -> float:
        """
        Gets the perigee argument (degrees), at the southernmost (for
        northern coverage) or northernmost point of the orbit.

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
        apogee, whose right ascension changes with the precession of the
        ascending node and, away from the critical inclination, of the
        argument of perigee. Without the node's precession, the apogee
        longitudes of a Molniya orbit drift westward by about 0.1 degrees
        per day. The result is cached.

        Args:
            revolutions_per_day (int): The number of revolutions per day.

        Returns:
            timedelta: The orbit period.
        """
        # keyed by the inputs, so that a copy with changed fields (e.g. from
        # `model_copy(update=...)`, which copies the cache) is recomputed
        key = (
            revolutions_per_day,
            self.perigee_altitude,
            self.inclination,
            self.northern_coverage,
            self.right_ascension_ascending_node,
            self.epoch,
        )
        return get_cached(
            self,
            "repeat_orbit_period",
            key,
            lambda: self._solve_repeat_orbit_period(revolutions_per_day),
        )

    def _solve_repeat_orbit_period(self, revolutions_per_day: int) -> timedelta:
        """
        Solves for the orbit period for which the ground track repeats with
        the given number of revolutions per day, without caching (see
        `_compute_repeat_orbit_period`).

        Args:
            revolutions_per_day (int): The number of revolutions per day.

        Returns:
            timedelta: The orbit period.
        """
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
            # rate of change of the apogee's right ascension with the
            # argument of perigee (zero at the critical inclination)
            argp = np.radians(self.get_perigee_argument())
            cos_i = np.cos(np.radians(self.get_inclination()))
            denominator = np.cos(argp) ** 2 + (cos_i * np.sin(argp)) ** 2
            apogee_rate = (
                satrec.nodedot + cos_i / denominator * satrec.argpdot
                if denominator > 1e-6
                else satrec.nodedot
            )
            target = revolutions_per_day * (rotation_rate - apogee_rate)
            period *= satrec.mdot / target
        return timedelta(seconds=period)

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
