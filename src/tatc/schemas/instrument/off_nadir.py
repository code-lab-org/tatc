"""
Object schemas for off-nadir pointing instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
from pydantic import Field
from shapely import MultiPolygon, Polygon
from shapely.geometry import MultiPoint, Point
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ...utils.projection import (
    compute_footprint,
    compute_projected_ray_position,
)
from .nadir import Instrument


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
        default=0,
        description="Left/right look angle (degrees) orthogonal to instrument motion.",
        ge=-180,
        le=180,
    )
    pitch_angle: float = Field(
        default=0,
        description="Fore/aft look angle (degrees) in direction of instrument motion.",
        ge=-180,
        le=180,
    )
    is_rectangular: bool = Field(
        default=False,
        description="True, if this instrument produces a rectangular view.",
    )
    cross_track_pixels: int = Field(
        default=1, description="Number of pixels in cross-track direction.", ge=1
    )
    along_track_pixels: int = Field(
        default=1, description="Number of pixels in along-track direction.", ge=1
    )
    cross_track_oversampling: float = Field(
        default=0,
        description="Fraction of pixel overlap in cross-track diraction.",
        ge=0,
        lt=1,
    )
    along_track_oversampling: float = Field(
        default=0,
        description="Fraction of pixel overlap in along-track diraction.",
        ge=0,
        lt=1,
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
    ) -> tuple[float, float]:
        """
        Gets the cone () and clock angles (degrees) for given pixel.

        Args:
            cross_track_index (int): pixel index in cross-track dimension (left to right).
            along_track_index (int): pixel index in along-track dimension (fore to aft).

        Returns:
            tuple[float, float]: cone (from nadir-looking) and clock
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
    ) -> list[MultiPoint]:
        """
        Compute the instanteous footprint pixel array.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            elevation (float): The elevation (meters) at which project the footprint.

        Returns:
            list[shapely.geometry.MultiPoint]: The instrument pixel array(s).
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
            for i in range(np.size(orbit_track.t))  # type: ignore
        ]
