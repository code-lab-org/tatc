# -*- coding: utf-8 -*-
"""
Object schemas for instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import Optional, Union
from datetime import timedelta

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

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

    def is_valid_observation(
        self, orbit_track: Geocentric, target: GeographicPosition = None
    ) -> Union[bool,npt.NDArray]:
        """Determines if an instrument can provide a valid observations.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): orbit track position/velocity from Skyfield
            target (skyfield.toposlib.GeographicPosition): target position from Skyfield

        Returns:
            numpy.typing.NDArray: Array of indicators: `True` if instrument provides a valid observation.
        """
        if target is None:
            # support backwards compatibility
            target = wgs84.subpoint_of(orbit_track)
        is_valid = np.ones(np.size(orbit_track.t), dtype=bool)
        if self.req_self_sunlit is not None:
            # compare requirement to satellite sunlit condition
            is_self_sunlit_valid = orbit_track.is_sunlit(de421) == self.req_self_sunlit
            is_valid = np.logical_and(is_valid, is_self_sunlit_valid)
        if self.req_target_sunlit is not None:
            # compute solar altitude angle at sub-satellite points
            solar_alt = (
                (de421["earth"] + target)
                .at(orbit_track.t)
                .observe(de421["sun"])
                .apparent()
                .altaz()[0]
                .degrees
            )
            # compare requirement to sub-satellite point sunlit conditions
            is_target_sunlit_valid = (solar_alt > 0) == self.req_target_sunlit
            is_valid = np.logical_and(is_valid, is_target_sunlit_valid)
        if np.size(orbit_track.t) > 1:
            return is_valid
        return bool(is_valid[0])


class PointedInstrument(Instrument):
    """
    Remote sensing instrument with off-nadir orientation.
    """

    cross_track_field_of_view: float = Field(
        ...,
        description="Angular field (degrees) of view orthogonal to instrument motion.",
        gt=0,
        le=180,
    )
    along_track_field_of_view: float = Field(
        ...,
        description="Angular field (degrees) of view in direction of instrument motion.",
        gt=0,
        le=180,
    )
    roll_angle: float = Field(
        0,
        description="Left/right look angle (degrees) orthogonal to instrument motion.",
        ge=-180,
        le=180,
    )
    pitch_angle: float = Field(
        0,
        description="Fore/aft look angle (degrees) in direction of instrument motion.",
        ge=-180,
        le=180,
    )
    is_rectangular: bool = Field(
        False, description="True, if this instrument produces a rectangular view."
    )
    cross_track_pixels: int = Field(
        1, description="Number of pixels in cross-track direction.", ge=1
    )
    along_track_pixels: int = Field(
        1, description="Number of pixels in along-track direction.", ge=1
    )
