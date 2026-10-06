"""
Object schemas for off-nadir pointing instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from pydantic import Field, field_validator
from shapely import MultiPolygon, Polygon
from shapely.geometry import MultiPoint, Point
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ...utils.orbital import compute_argument_of_latitude
from ...utils.projection import (
    VelocityFrame,
    ViewGeometry,
    compute_footprint,
    compute_projected_ray_position,
    compute_view_angles,
    compute_view_tangents,
)
from .simple import Instrument


class PointedInstrument(Instrument):
    """
    Remote sensing instrument with an optional off-nadir (roll/pitch)
    pointing offset and a rectangular or elliptical pixel array.
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
        description="Left/right look angle (degrees): a rotation of the view about the "
        + "along-track axis, positive to the left of the direction of motion.",
        ge=-180,
        le=180,
    )
    roll_angle_profile: list[tuple[float, float]] | None = Field(
        default=None,
        description="Roll angle (degrees) as a function of the satellite's "
        + "argument of latitude (degrees from the ascending node, 0 to 360, in "
        + "the direction of motion), as (argument of latitude, roll angle) pairs "
        + "interpolated linearly and periodically, for an instrument whose roll "
        + "varies around the orbit (for example, the roll steering of a "
        + "synthetic aperture radar). If given, it overrides `roll_angle`.",
        min_length=1,
    )
    pitch_angle: float = Field(
        default=0,
        description="Fore/aft look angle (degrees): a rotation of the view about the "
        + "rolled cross-track axis (after roll), positive forward.",
        ge=-180,
        le=180,
    )
    view_geometry: ViewGeometry = Field(
        default=ViewGeometry.FRAME,
        description="Geometry in which the fields of view and pixels are defined: "
        + "in a plane perpendicular to the boresight (`frame`, as for a framing "
        + "camera or a pushbroom array), or in cross-track and along-track angles "
        + "(`scan`, as for a cross-track scanner, whose along-track angular extent "
        + "is constant across the scan).",
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
        description="Fraction of pixel overlap in cross-track direction.",
        ge=0,
        lt=1,
    )
    along_track_oversampling: float = Field(
        default=0,
        description="Fraction of pixel overlap in along-track direction.",
        ge=0,
        lt=1,
    )
    velocity_frame: VelocityFrame = Field(
        default=VelocityFrame.EARTH_FIXED,
        description="Reference frame of the velocity vector that defines the "
        + "along-track direction (`earth_fixed`: aligned with the ground track, "
        + "as for a yaw-steered spacecraft; `inertial`: aligned with the orbit "
        + "plane, as for a spacecraft without yaw steering).",
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
            roll_angle=self.get_roll_angle(orbit_track),
            pitch_angle=self.pitch_angle,
            is_rectangular=self.is_rectangular,
            number_points=number_points,
            elevation=elevation,
            velocity_frame=self.velocity_frame,
            nadir_reference=self.nadir_reference,
            view_geometry=self.view_geometry,
        )

    @field_validator("roll_angle_profile")
    @classmethod
    def _validate_roll_angle_profile(
        cls, profile: list[tuple[float, float]] | None
    ) -> list[tuple[float, float]] | None:
        """
        Validates that the profile's arguments of latitude lie in [0, 360)
        and its roll angles in [-180, 180], and sorts it by argument of
        latitude.
        """
        if profile is None:
            return None
        for argument, roll in profile:
            if not 0 <= argument < 360:
                raise ValueError("Arguments of latitude must be in [0, 360).")
            if not -180 <= roll <= 180:
                raise ValueError("Roll angles must be in [-180, 180].")
        return sorted(profile)

    def get_roll_angle(
        self, orbit_track: Geocentric | None = None
    ) -> float | npt.NDArray[np.float64]:
        """
        Gets the roll angle (degrees) at the satellite position(s) of an orbit
        track: interpolated from `roll_angle_profile` at the satellite's
        argument of latitude, if given, or otherwise `roll_angle`.

        Args:
            orbit_track (skyfield.positionlib.Geocentric | None): The
                satellite position/velocity (in an inertial frame).

        Returns:
            float | numpy.typing.NDArray[numpy.float64]: the roll angle(s) (degrees)
        """
        if self.roll_angle_profile is None or orbit_track is None:
            return self.roll_angle
        argument = compute_argument_of_latitude(
            orbit_track.position.m, orbit_track.velocity.m_per_s
        )
        arguments, rolls = zip(*self.roll_angle_profile)
        roll = np.interp(argument, arguments, rolls, period=360)
        return float(roll) if np.ndim(roll) == 0 else roll

    def is_in_field_of_view(
        self,
        orbit_track: Geocentric,
        target: GeographicPosition,
    ) -> npt.NDArray[np.bool_]:
        """
        Determines if a target lies within the instantaneous field of view.
        Does not check whether the target is above the satellite's horizon.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite position/velocity.
            target (skyfield.toposlib.GeographicPosition): The target position.

        Returns:
            numpy.typing.NDArray: Array of indicators: `True` if the target is in the field of view.
        """
        if self.view_geometry == ViewGeometry.SCAN:
            roll, pitch = compute_view_angles(
                orbit_track, target, self.velocity_frame, self.nadir_reference
            )
            # angular offsets from the view center, normalized by the half widths
            cross_offset = (
                (roll - self.get_roll_angle(orbit_track) + 180) % 360 - 180
            ) / (self.cross_track_field_of_view / 2)
            along_offset = (pitch - self.pitch_angle) / (
                self.along_track_field_of_view / 2
            )
        else:
            along, cross = compute_view_tangents(
                orbit_track,
                target,
                self.velocity_frame,
                self.get_roll_angle(orbit_track),
                self.pitch_angle,
                self.nadir_reference,
            )
            # offsets from the view center, normalized by the view half widths
            along_offset = along / np.tan(
                np.radians(self.along_track_field_of_view / 2)
            )
            cross_offset = cross / np.tan(
                np.radians(self.cross_track_field_of_view / 2)
            )
        with np.errstate(invalid="ignore"):
            if self.is_rectangular:
                inside = (np.abs(along_offset) <= 1) & (np.abs(cross_offset) <= 1)
            else:
                inside = along_offset**2 + cross_offset**2 <= 1
        return np.atleast_1d(inside)

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
            roll_angle=self.get_roll_angle(orbit_track),
            pitch_angle=self.pitch_angle,
            is_rectangular=False,
            angle=0,
            elevation=elevation,
            velocity_frame=self.velocity_frame,
            nadir_reference=self.nadir_reference,
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
            cross_track_index (int): cross-track pixel index (right to left:
                index 0 is on the right of the direction of motion).
            along_track_index (int): along-track pixel index (fore-to-aft).
            elevation (float): The elevation (meters) at which project the pixel.

        Returns:
            (skyfield.toposlib.GeographicPosition): the geographic position of the projected pixel
        """
        if self.view_geometry == ViewGeometry.SCAN:
            cross_offset, along_offset = self._get_pixel_offsets(
                cross_track_index, along_track_index
            )
            return compute_projected_ray_position(
                orbit_track=orbit_track,
                cross_track_field_of_view=0,
                along_track_field_of_view=0,
                roll_angle=self.get_roll_angle(orbit_track) + cross_offset,
                pitch_angle=self.pitch_angle + along_offset,
                elevation=elevation,
                velocity_frame=self.velocity_frame,
                nadir_reference=self.nadir_reference,
            )
        cone, clock = self.get_pixel_cone_and_clock_angle(
            cross_track_index, along_track_index
        )

        return compute_projected_ray_position(
            orbit_track=orbit_track,
            # `cone` is the pixel's angular offset from boresight (not a
            # field of view), but compute_projected_ray_position's
            # elliptical ray offsets by half of the field of view it is
            # given; doubling `cone` here cancels that halving so the
            # pixel lands at its true cone-angle offset.
            cross_track_field_of_view=2 * cone,
            along_track_field_of_view=2 * cone,
            roll_angle=self.get_roll_angle(orbit_track),
            pitch_angle=self.pitch_angle,
            is_rectangular=False,
            angle=clock,
            elevation=elevation,
            velocity_frame=self.velocity_frame,
            nadir_reference=self.nadir_reference,
        )

    def _get_pixel_offsets(
        self, cross_track_index: int, along_track_index: int
    ) -> tuple[float, float]:
        """
        Gets a pixel's angular offsets (degrees) from the view center, across
        track (positive left) and along track (positive forward).

        Args:
            cross_track_index (int): pixel index in cross-track dimension
                (right to left: index 0 is on the right of the direction of motion).
            along_track_index (int): pixel index in along-track dimension (fore to aft).

        Returns:
            tuple[float, float]: the cross-track and along-track offsets (degrees).
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
        return cross_track_offset, along_track_offset

    def get_pixel_cone_and_clock_angle(
        self, cross_track_index: int, along_track_index: int
    ) -> tuple[float, float]:
        """
        Gets the cone and clock angles (degrees) for a given pixel: its
        angular offset from the instrument boresight (after roll/pitch),
        expressed in polar form.

        Args:
            cross_track_index (int): pixel index in cross-track dimension
                (right to left: index 0 is on the right of the direction of motion).
            along_track_index (int): pixel index in along-track dimension (fore to aft).

        Returns:
            tuple[float, float]: cone (the pixel's total angular
                displacement from boresight) and clock (about the boresight,
                from the cross-track axis on the left of the direction of
                motion toward the along-track axis forward) angles (degrees).
        """
        cross_track_offset, along_track_offset = self._get_pixel_offsets(
            cross_track_index, along_track_index
        )
        return (
            np.sqrt(cross_track_offset**2 + along_track_offset**2),
            np.degrees(np.arctan2(along_track_offset, cross_track_offset)),
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
