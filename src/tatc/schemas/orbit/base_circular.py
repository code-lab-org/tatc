"""
Base object schemas for circular orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta

from pydantic import Field

from ... import constants, utils
from .base import OrbitBase


class CircularOrbitBase(OrbitBase):
    """
    Base class for circular orbits.
    """

    mean_altitude: float = Field(..., description="Mean altitude (meters).")
    
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
