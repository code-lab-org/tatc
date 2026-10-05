"""
Object schemas for sunsynchronous orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import date, datetime, time, timedelta
from typing import Literal

import numpy as np
from pydantic import AliasChoices, Field

from ... import constants, utils
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
        Gets the inclination (decimal degrees).

        Returns:
            float: the inclination
        """
        return np.degrees(
            np.arccos(-np.power(self.get_semimajor_axis() / 12352000, 7 / 2))
        )

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
