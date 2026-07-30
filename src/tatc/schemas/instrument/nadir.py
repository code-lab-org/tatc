"""
Object schemas for nadir-pointing instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field
from shapely import MultiPolygon, Polygon
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ...constants import de421
from ...utils.observation import (
    compute_min_elevation_angle,
    field_of_regard_to_swath_width,
)
from ...utils.projection import (
    compute_footprint,
    compute_projected_ray_position,
)


class Instrument(BaseModel):
    """
    Remote sensing instrument.
    """

    name: str = Field(default="Default", description="Instrument name.")
    field_of_regard: float = Field(
        default=180,
        description="Angular field (degrees) of possible observations (with pointing).",
        gt=0,
        le=360,
        examples=[50],
    )
    min_access_time: timedelta = Field(
        default=timedelta(0),
        description="Minimum access (integration) time to record an observation.",
        examples=[timedelta(seconds=10)],
    )
    req_self_sunlit: bool | None = Field(
        default=None,
        description="Required instrument sunlit state for valid observation "
        + "(`True`: sunlit, `False`: eclipse, `None`: no requirement).",
    )
    req_target_sunlit: bool | None = Field(
        default=None,
        description="Required target sunlit state for valid observation "
        + "(`True`: sunlit, `False`: eclipse, `None`: no requirement).",
    )
    access_time_fixed: bool = Field(
        default=False, description="`True`, if access time is fixed to minimum value."
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

    def compute_footprint(
        self,
        orbit_track: Geocentric,
        number_points: int | None = None,
        elevation: float = 0,
    ) -> list[Polygon | MultiPolygon]:
        """
        Compute the instanteous instrument footprint.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            number_points (int | None): The required number of polygon points to generate.
            elevation (float): The elevation (meters) at which project the footprint.

        Returns:
            list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]: The instrument footprint(s).
        """
        return compute_footprint(
            orbit_track=orbit_track,
            cross_track_field_of_view=self.field_of_regard,
            along_track_field_of_view=self.field_of_regard,
            roll_angle=0,
            pitch_angle=0,
            is_rectangular=False,
            number_points=number_points,
            elevation=elevation,
        )

    def compute_footprint_center(
        self,
        orbit_track: Geocentric,
        elevation: float = 0,
    ) -> GeographicPosition:
        """
        Compute the center of an instaneous instrument footprint.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            elevation (float): The elevation (meters) at which project the footprint.

        Returns:
            skyfield.toposlib.GeographicPosition: The instrument footprint center.
        """
        return compute_projected_ray_position(
            orbit_track=orbit_track,
            cross_track_field_of_view=0,
            along_track_field_of_view=0,
            roll_angle=0,
            pitch_angle=0,
            is_rectangular=False,
            angle=0,
            elevation=elevation,
        )

    def is_valid_observation(
        self, orbit_track: Geocentric, target: GeographicPosition | None = None
    ) -> npt.NDArray[np.bool_]:
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
        is_valid = np.ones(np.size(orbit_track.t), dtype=bool)  # type: ignore
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
        return is_valid
