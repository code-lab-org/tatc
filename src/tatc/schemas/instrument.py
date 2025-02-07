# -*- coding: utf-8 -*-
"""
Object schemas for instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import List, Optional, Tuple, Union
from datetime import timedelta

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field
from shapely import Geometry
from shapely.geometry import MultiPoint, Point
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ..constants import de421
from ..utils import (
    compute_footprint,
    compute_min_elevation_angle,
    compute_projected_ray_position,
    field_of_regard_to_swath_width,
)


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

    def compute_footprint(
        self,
        orbit_track: Geocentric,
        number_points: int = None,
        elevation: float = 0,
    ) -> Union[Geometry, List[Geometry]]:
        """
        Compute the instanteous instrument footprint.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            number_points (int): The required number of polygon points to generate.
            elevation (float): The elevation (meters) at which project the footprint.

        Returns:
            Union[shapely.Geometry, List[shapely.Geometry]: The instrument footprint(s).
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
        self, orbit_track: Geocentric, target: GeographicPosition = None
    ) -> Union[bool, npt.NDArray]:
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
    Remote sensing instrument with optional off-nadir orientation.
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
    cross_track_oversampling: float = Field(
        0, description="Fraction of pixel overlap in cross-track diraction.", ge=0, lt=1
    )
    along_track_oversampling: float = Field(
        0, description="Fraction of pixel overlap in along-track diraction.", ge=0, lt=1
    )

    def get_cross_track_instantaneous_field_of_view(self) -> float:
        """
        Gets the instananeous field of view (degrees) for cross-track pixels.

        Returns:
            float: the cross-track instantaneous pixel field of view (degrees)
        """
        return (
            self.cross_track_field_of_view
            / self.cross_track_pixels
            / (1 - self.cross_track_oversampling)
        )

    def get_along_track_instantaneous_field_of_view(self) -> float:
        """
        Gets the instananeous field of view (degrees) for along-track pixels.

        Returns:
            float: the along-track instantaneous pixel field of view (degrees)
        """
        return (
            self.along_track_field_of_view
            / self.along_track_pixels
            / (1 - self.along_track_oversampling)
        )

    def compute_footprint(
        self,
        orbit_track: Geocentric,
        number_points: int = None,
        elevation: float = 0,
    ) -> Union[Geometry, List[Geometry]]:
        """
        Compute the instanteous instrument footprint.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            number_points (int): The required number of polygon points to generate.
            elevation (float): The elevation (meters) at which project the footprint.

        Returns:
            Union[shapely.Geometry, List[shapely.Geometry]: The instrument footprint(s).
        """
        return compute_footprint(
            orbit_track=orbit_track,
            cross_track_field_of_view=self.cross_track_field_of_view,
            along_track_field_of_view=self.along_track_field_of_view,
            roll_angle=self.roll_angle,
            pitch_angle=self.pitch_angle,
            is_rectangular=self.is_rectangular,
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
            roll_angle=self.roll_angle,
            pitch_angle=self.pitch_angle,
            is_rectangular=False,
            angle=0,
            elevation=elevation,
        )

    def compute_projected_pixel_position(
        self,
        orbit_track: Geocentric,
        cross_track_index: int,
        along_track_index: int,
        elevation: float = 0,
    ) -> GeographicPosition:
        """
        Get the location of a projected pixel.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): the satellite orbit track.
            cross_track_index (int): cross-track pixel index (left-to-right).
            along_track_index (int): along-track pixel index (fore-to-aft).
            elevation (float): The elevation (meters) at which project the pixel.

        Returns:
            (skyfield.toposlib.GeographicPosition): the geographic position of the projected pixel
        """
        cone, clock = self.get_pixel_cone_and_clock_angle(
            cross_track_index, along_track_index
        )

        return compute_projected_ray_position(
            orbit_track=orbit_track,
            cross_track_field_of_view=cone,
            along_track_field_of_view=cone,
            roll_angle=self.roll_angle,
            pitch_angle=self.pitch_angle,
            is_rectangular=False,
            angle=clock,
            elevation=elevation,
        )

    def get_pixel_cone_and_clock_angle(
        self, cross_track_index: int, along_track_index: int
    ) -> Tuple[float, float]:
        """
        Gets the cone () and clock angles (degrees) for given pixel.

        Args:
            cross_track_index (int): pixel index in cross-track dimension (left to right).
            along_track_index (int): pixel index in along-track dimension (fore to aft).

        Returns:
            Tuple[float, float]: cone (from nadir-looking) and clock
                (counter-clockwise from right-looking) angles (degrees).
        """
        cross_track_offset = (
            (0.5 + cross_track_index - self.cross_track_pixels / 2)
            * (1 - self.cross_track_oversampling)
            * self.cross_track_field_of_view
            / self.cross_track_pixels
        )
        along_track_offset = (
            (self.along_track_pixels / 2 - 0.5 - along_track_index)
            * (1 - self.along_track_oversampling)
            * self.along_track_field_of_view
            / self.along_track_pixels
        )
        return (
            np.sqrt(cross_track_offset**2 + along_track_offset**2),
            np.arctan2(along_track_offset, cross_track_offset),
        )

    def compute_footprint_pixel_array(
        self,
        orbit_track: Geocentric,
        elevation: float = 0,
    ) -> Union[MultiPoint, List[MultiPoint]]:
        """
        Compute the instanteous footprint pixel array.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            elevation (float): The elevation (meters) at which project the footprint.

        Returns:
            Union[shapely.geometry.MultiPoint, List[shapely.geometry.MultiPoint]]: The instrument pixel array(s).
        """
        points = [
            self.compute_projected_pixel_position(
                orbit_track=orbit_track,
                cross_track_index=i,
                along_track_index=j,
                elevation=elevation,
            )
            for i in range(self.cross_track_pixels)
            for j in range(self.along_track_pixels)
        ]
        if np.size(orbit_track.t) > 1:
            return [
                MultiPoint(
                    [
                        Point(
                            point.longitude.degrees[i],
                            point.latitude.degrees[i],
                            point.elevation.m[i],
                        )
                        for point in points
                    ]
                )
                for i in range(np.size(orbit_track.t))
            ]
        return MultiPoint(
            [
                Point(
                    point.longitude.degrees, point.latitude.degrees, point.elevation.m
                )
                for point in points
            ]
        )
