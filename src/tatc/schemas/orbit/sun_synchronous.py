"""
Object schemas for sunsynchronous orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import date, datetime, time, timedelta, timezone
from typing import Literal

import numpy as np
from pydantic import AliasChoices, Field
from sgp4.api import WGS72, Satrec

from ... import constants, utils
from ...utils.cache import get_cached
from .base_circular import CircularOrbitBase


class SunSynchronousOrbit(CircularOrbitBase):
    """
    Orbit defined by sun synchronous parameters.
    """

    type: Literal["sso"] = Field(default="sso", description="Orbit type discriminator.")
    mean_altitude: float = Field(
        ...,
        description="Mean altitude (meters).",
        ge=0,
        lt=12352000 - constants.EARTH_MEAN_RADIUS,
        validation_alias=AliasChoices("mean_altitude", "altitude"),
    )
    equator_crossing_time: time = Field(
        ..., description="Equator crossing time (local solar time)."
    )
    equator_crossing_ascending: bool = Field(
        default=True,
        description="True, if the equator crossing time is ascending (south-to-north).",
    )

    def get_inclination(self) -> float:
        """
        Gets the inclination (decimal degrees) for which the orbit plane, as
        propagated by SGP4, precesses at the rate of the mean Sun (360 degrees
        per tropical year), so that the local mean solar time of the
        equator crossings stays constant. Starting from the classical J2
        estimate (cos i = -(a / 12,352 km)^(7/2)), the inclination is refined
        with SGP4's secular precession rate of the ascending node (which
        includes the Earth's J2 and J4 oblateness terms) for the orbit's
        general perturbations elements. The result is cached.

        Returns:
            float: the inclination
        """
        semimajor_axis = self.get_semimajor_axis()
        return get_cached(
            self,
            "sso_inclination",
            semimajor_axis,
            lambda: self._compute_inclination(semimajor_axis),
        )

    def _compute_inclination(self, semimajor_axis: float) -> float:
        """
        Computes the sun-synchronous inclination, without caching (see
        `get_inclination`).

        Args:
            semimajor_axis (float): The semimajor axis (meters).

        Returns:
            float: the inclination (degrees)
        """
        inclination = np.arccos(-np.power(semimajor_axis / 12352000, 7 / 2))
        # mean motion (radians/minute), as for the general perturbations elements
        mean_motion = (
            np.radians(utils.orbital.semimajor_axis_to_mean_motion(semimajor_axis)) * 60
        )
        # precession rate (radians/minute) of the mean Sun
        target = 2 * np.pi / constants.TROPICAL_YEAR_S * 60
        epoch = (self.epoch - datetime(1949, 12, 31, tzinfo=timezone.utc)) / timedelta(
            days=1
        )
        for _ in range(5):
            satrec = Satrec()
            satrec.sgp4init(
                WGS72,
                "i",
                0,
                epoch,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                inclination,
                0.0,
                mean_motion,
                0.0,
            )
            # the precession rate is proportional to cos(inclination)
            inclination += (satrec.nodedot - target) / (
                satrec.nodedot * np.tan(inclination)
            )
        return float(np.degrees(inclination))

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node (decimal degrees): the
        right ascension at which the ascending node's local mean solar time,
        at the epoch, equals the equator crossing time (or, for a descending
        equator crossing time, the time 12 hours later). Local mean solar
        time at a longitude is the universal time plus the longitude (15
        degrees per hour), so the node's longitude at the epoch is 15 degrees
        per hour of the difference between its local time and the universal
        time, and its right ascension adds the Greenwich mean sidereal time
        (consistent with the frame of the general perturbations elements).

        Returns:
            float: the right ascension of ascending node
        """
        ect_hours = timedelta(
            hours=self.equator_crossing_time.hour,
            minutes=self.equator_crossing_time.minute,
            seconds=self.equator_crossing_time.second,
            microseconds=self.equator_crossing_time.microsecond,
        ) / timedelta(hours=1)
        # local mean solar time of the ascending node
        node_hours = ect_hours + (0 if self.equator_crossing_ascending else 12)
        epoch_time = constants.timescale.from_datetime(self.epoch)
        universal_hours = np.mod(epoch_time.ut1 + 0.5, 1) * 24
        return float(np.mod(15 * (epoch_time.gmst + node_hours - universal_hours), 360))

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> SunSynchronousOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            SunSynchronousOrbit: the derived orbit
        """
        true_anomaly = utils.orbital.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360)
        )
        # every 15 degrees of raan shift ect by 1 hour
        equator_crossing_time = (
            datetime.combine(date(2000, 1, 1), self.equator_crossing_time)
            + timedelta(hours=delta_raan / 15)
        ).time()
        return SunSynchronousOrbit(
            mean_altitude=self.mean_altitude,
            equator_crossing_time=equator_crossing_time,
            equator_crossing_ascending=self.equator_crossing_ascending,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
        )
