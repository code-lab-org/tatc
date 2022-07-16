# -*- coding: utf-8 -*-
"""
Object schemas for instruments.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

from enum import Enum
from typing import Optional, Union
from pydantic import BaseModel, Field, confloat
from datetime import timedelta
import pandas as pd
import numpy as np
from skyfield.api import wgs84, EarthSatellite
from skyfield.timelib import Time
from shapely.geometry import Point, Polygon, MultiPolygon

from .. import constants, utils
from ..constants import de421, timescale


class Instrument(BaseModel):
    """
    Representation of a remote sensing instrument.

    :param name: Name of this instrument
    :type  name: str
    :param field_of_regard: Angular field (degrees) of possible observations (with pointing).
    :type field_of_regard: float (0 to 360), default: 180
    :param min_access_time: Minimum access (integration) time to record an observation.
    :type min_access_time: :class:`datetime.timedelta`, default: 0
    :param req_self_sunlit: Option to set the required instrument sunlit state for observation.
    :type req_self_sunlit: bool, optional, default: None
    :param req_target_sunlit: Option to set the required target sunlit state for observation.
    :type req_target_sunlit: bool, optional, default: None
    """

    name: str = Field(..., description="Name of this instrument.")
    field_of_regard: float = Field(
        180,
        description="Angular field (degrees) of possible observations (with pointing).",
        gt=0,
        le=360,
        example=50,
    )
    min_access_time: timedelta = Field(
        timedelta(0),
        description="Minimum access (integration) time to record an observation.",
        example=timedelta(seconds=10),
    )
    req_self_sunlit: Optional[bool] = Field(
        None,
        description="Option to set the required instrument sunlit state for observation.",
    )
    req_target_sunlit: Optional[bool] = Field(
        None,
        description="Option to set the required target sunlit state for observation.",
    )

    def get_swath_width(self, height: float) -> float:
        """
        Gets the instrument swath width projected to the Earth's surface.

        Args:
            height (float): Height (meters) above surface of the observation.

        Returns:
            float: The observation diameter (meters).
        """
        return utils.field_of_regard_to_swath_width(height, self.field_of_regard)

    def get_min_elevation_angle(self, height: float) -> float:
        """
        Get the minimum elevation angle required to observe a point.

        Args:
            height (float): Height (meters) above surface of the observation.

        Returns:
            float: The minimum elevation angle (degrees) for observation.
        """
        return utils.compute_min_elevation_angle(height, self.field_of_regard)

    def is_valid_observation(self, sat: EarthSatellite, time: Time) -> bool:
        """Determines if an instrument can provide a valid observations.

        :param sat: Satellite hosting this instrument (Skyfield).
        :type sat: :obj:`EarthSatellite`
        :param time: Observation time(s) (Skyfield).
        :type time: :obj:`Time`
        :return: True if instrument can provide valid observations, otherwise False
        :rtype: bool
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
        else:
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
