# -*- coding: utf-8 -*-
"""
Object schemas for instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import Optional
from datetime import timedelta

from pydantic import BaseModel, Field
import numpy as np
from skyfield.api import wgs84, EarthSatellite
from skyfield.timelib import Time

from ..constants import de421
from ..utils import compute_min_elevation_angle, field_of_regard_to_swath_width


class Instrument(BaseModel):
    """
    Remote sensing instrument.
    """

    name: str = Field("Default", description="Instrument name.")
    field_of_regard: float = Field(
        180,
        description="Angular field (degrees) of possible observations (with pointing).",
        gt=0,
        le=360,
        examples=[50],
    )
    min_access_time: timedelta = Field(
        timedelta(0),
        description="Minimum access (integration) time to record an observation.",
        examples=[timedelta(seconds=10)],
    )
    req_self_sunlit: Optional[bool] = Field(
        None,
        description="Required instrument sunlit state for valid observation "
        + "(`True`: sunlit, `False`: eclipse, `None`: no requirement).",
    )
    req_target_sunlit: Optional[bool] = Field(
        None,
        description="Required target sunlit state for valid observation "
        + "(`True`: sunlit, `False`: eclipse, `None`: no requirement).",
    )
    access_time_fixed: bool = Field(
        False, description="`True`, if access time is fixed to minimum value."
    )

    def get_swath_width(self, height: float) -> float:
        """
        Gets the instrument swath width projected to the Earth's surface.

        Args:
            height (float): Height (meters) above surface of the observation.

        Returns:
            float: The observation diameter (meters).
        """
        return field_of_regard_to_swath_width(height, self.field_of_regard)

    def get_min_elevation_angle(self, height: float) -> float:
        """
        Get the minimum elevation angle required to observe a point.

        Args:
            height (float): Height (meters) above surface of the observation.

        Returns:
            float: The minimum elevation angle (degrees) for observation.
        """
        return compute_min_elevation_angle(height, self.field_of_regard)

    def is_valid_observation(self, sat: EarthSatellite, time: Time) -> bool:
        """Determines if an instrument can provide a valid observations.

        Args:
            sat (EarthSatellite): Satellite hosting this instrument.
            time (Time): Observation time(s).

        Returns:
            bool: `True` if instrument provides a valid observation, otherwise `False`.
        """
        if isinstance(time.tt, float):
            # scalar
            is_valid = True
            if self.req_self_sunlit is not None:
                # compare requirement to satellite sunlit condition
                is_valid = is_valid and (
                    self.req_self_sunlit == sat.at(time).is_sunlit(de421)
                )
            if self.req_target_sunlit is not None:
                # compute solar altitude angle at sub-satellite point
                solar_alt = (
                    (de421["earth"] + wgs84.subpoint_of(sat.at(time)))
                    .at(time)
                    .observe(de421["sun"])
                    .apparent()
                    .altaz()[0]
                    .degrees
                )
                # compare requirement to sub-satellite point sunlit condition
                is_valid = is_valid and (self.req_target_sunlit == (solar_alt > 0))
            return is_valid
        # vector
        is_valid = np.ones(len(time), dtype=bool)
        if self.req_self_sunlit is not None:
            # compare requirement to satellite sunlit condition
            is_valid = np.logical_and(
                is_valid, self.req_self_sunlit == sat.at(time).is_sunlit(de421)
            )
        if self.req_target_sunlit is not None:
            # compute solar altitude angle at sub-satellite points
            solar_alt = np.array(
                [
                    (de421["earth"] + wgs84.subpoint_of(sat.at(t)))
                    .at(t)
                    .observe(de421["sun"])
                    .apparent()
                    .altaz()[0]
                    .degrees
                    for t in time
                ]
            )
            # compare requirement to sub-satellite point sunlit conditions
            is_valid = np.logical_and(
                is_valid, self.req_target_sunlit == (solar_alt > 0)
            )
        return is_valid
