"""
Object schemas for the attitude of space systems: the orientation of a
spacecraft's body relative to its orbit, which orients the views of its
instruments, and schedules of attitudes over time.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime
from typing import Annotated, Literal

import numpy as np
import numpy.typing as npt
from pydantic import AwareDatetime, BaseModel, Field, field_validator
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ..utils.ellipsoid import compute_ellipsoid_intersection, rectangular_to_geodetic
from ..utils.glint import compute_specular_point
from ..utils.orbital import _check_angle_profile, _interpolate_angle_profile
from ..utils.projection import (
    VelocityFrame,
    _compute_line_of_sight,
    _compute_nadir_axes,
    compute_view_angles,
)
from ..utils.time import _index_orbit_track
from .surface.point import Point


class Attitude(BaseModel):
    """
    Attitude of a spacecraft's body relative to its orbit. The body is
    rotated from the orbit's nadir frame (its nadir axis and its along-track
    axis, aligned with the velocity in `velocity_frame`) first by the yaw
    angle about the nadir axis, then by the roll angle about the yawed
    along-track axis, then by the pitch angle about the rolled cross-track
    axis. Each angle is constant or, if its profile is given, varies with the
    satellite's argument of latitude.

    The instruments of a space system point relative to its body: their
    roll and pitch angles, scan azimuths, and pixels are measured from the
    body's axes. For example, a yaw angle of 180 degrees (flying backward)
    swaps the left and right of a cross-track scanner and the fore and aft
    of a conical scan. An instrument's field of regard remains measured from
    nadir.
    """

    type: Literal["attitude"] = Field(
        default="attitude", description="Attitude type discriminator."
    )
    velocity_frame: VelocityFrame = Field(
        default=VelocityFrame.EARTH_FIXED,
        description="Reference frame of the velocity vector that defines the "
        + "along-track axis (`earth_fixed`: aligned with the ground track, as "
        + "for a yaw-steered spacecraft; `inertial`: aligned with the orbit "
        + "plane, as for a spacecraft without yaw steering).",
    )
    yaw_angle: float = Field(
        default=0,
        description="Rotation (degrees) of the body about the nadir axis, "
        + "positive turning the front of the body to the left of the direction "
        + "of motion (180: flying backward).",
        ge=-180,
        le=180,
    )
    roll_angle: float = Field(
        default=0,
        description="Rotation (degrees) of the body about its along-track axis, "
        + "positive tilting its nadir axis to the left of the direction of motion.",
        ge=-180,
        le=180,
    )
    pitch_angle: float = Field(
        default=0,
        description="Rotation (degrees) of the body about its rolled "
        + "cross-track axis, positive tilting its nadir axis forward.",
        ge=-180,
        le=180,
    )
    yaw_angle_profile: list[tuple[float, float]] | None = Field(
        default=None,
        description="Yaw angle (degrees) as a function of the satellite's "
        + "argument of latitude (degrees from the ascending node, 0 to 360), as "
        + "(argument of latitude, yaw angle) pairs interpolated linearly and "
        + "periodically. If given, it overrides `yaw_angle`.",
        min_length=1,
    )
    roll_angle_profile: list[tuple[float, float]] | None = Field(
        default=None,
        description="Roll angle (degrees) as a function of the satellite's "
        + "argument of latitude (see `yaw_angle_profile`), for example the roll "
        + "steering of a synthetic aperture radar. If given, it overrides "
        + "`roll_angle`.",
        min_length=1,
    )
    pitch_angle_profile: list[tuple[float, float]] | None = Field(
        default=None,
        description="Pitch angle (degrees) as a function of the satellite's "
        + "argument of latitude (see `yaw_angle_profile`). If given, it "
        + "overrides `pitch_angle`.",
        min_length=1,
    )

    @field_validator("yaw_angle_profile", "roll_angle_profile", "pitch_angle_profile")
    @classmethod
    def _validate_angle_profile(
        cls, profile: list[tuple[float, float]] | None, info
    ) -> list[tuple[float, float]] | None:
        """
        Validates an angle profile (see `_check_angle_profile`).
        """
        return _check_angle_profile(profile, info.field_name.split("_")[0].title())

    def get_angles(
        self, orbit_track: Geocentric
    ) -> tuple[npt.ArrayLike, npt.ArrayLike, npt.ArrayLike]:
        """
        Gets the yaw, roll, and pitch angles (degrees) of the body at the
        satellite position(s) of an orbit track.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite
                position/velocity (in an inertial frame).

        Returns:
            tuple[numpy.typing.ArrayLike, numpy.typing.ArrayLike, numpy.typing.ArrayLike]:
                the yaw, roll, and pitch angles (degrees)
        """
        return tuple(
            (
                angle
                if profile is None
                else _interpolate_angle_profile(profile, orbit_track)
            )
            for angle, profile in (
                (self.yaw_angle, self.yaw_angle_profile),
                (self.roll_angle, self.roll_angle_profile),
                (self.pitch_angle, self.pitch_angle_profile),
            )
        )

    def get_body_angles(
        self, orbit_track: Geocentric
    ) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
        """
        Gets the yaw, roll, and pitch angles (degrees) of the body and whether
        its along-track axis is aligned with the inertial velocity, at each
        time of an orbit track, as arrays with the shape of its times (as
        used by the projection utilities, see
        `tatc.utils.projection._compute_nadir_axes`).

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite
                position/velocity (in an inertial frame).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
                the yaw, roll, and pitch angles (degrees), and whether the
                velocity frame is inertial
        """
        shape = np.shape(orbit_track.t.tt)
        yaw, roll, pitch = (
            np.broadcast_to(np.asarray(angle, dtype=float), shape)
            for angle in self.get_angles(orbit_track)
        )
        inertial = np.full(shape, self.velocity_frame == VelocityFrame.INERTIAL)
        return yaw, roll, pitch, inertial


class TargetAttitude(BaseModel):
    """
    Attitude of a spacecraft's body that points its nadir axis at ground
    targets, as for an agile satellite imaging sites or a target mode: at
    each time, at the target closest to nadir (with the smallest off-nadir
    angle) among those within `max_off_nadir_angle` of nadir and above the
    horizon, and otherwise along nadir. The body is rotated by the yaw angle
    about the nadir axis and then by the roll and pitch angles that point
    its nadir axis at the target (see `tatc.utils.compute_view_angles`), so
    that its along-track axis stays close to the direction of motion. The
    body turns instantly between targets and nadir (slews are not modeled).

    The fields of regard of its instruments (measured from nadir) should
    include the off-nadir angles at which the body points. As an
    instrument's view stays on a target while the body points at it, the
    epoch of an observation of the target is a time while the body points
    at it (rather than the time a view sweeps over it).
    """

    type: Literal["target"] = Field(
        default="target", description="Attitude type discriminator."
    )
    targets: list[Point] = Field(
        ..., min_length=1, description="The ground targets at which to point."
    )
    max_off_nadir_angle: float = Field(
        default=30,
        description="Maximum angle (degrees) from nadir at which the body "
        + "points at a target.",
        gt=0,
        lt=90,
    )
    velocity_frame: VelocityFrame = Field(
        default=VelocityFrame.EARTH_FIXED,
        description="Reference frame of the velocity vector that defines the "
        + "along-track axis (see `Attitude.velocity_frame`).",
    )
    yaw_angle: float = Field(
        default=0,
        description="Rotation (degrees) of the body about the nadir axis "
        + "before pointing (see `Attitude.yaw_angle`).",
        ge=-180,
        le=180,
    )

    def _get_targets(self) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
        """
        Gets the targets' longitudes, latitudes (degrees), and elevations (meters).
        """
        return tuple(
            np.array([getattr(target, name) for target in self.targets], dtype=float)
            for name in ("longitude", "latitude", "elevation")
        )

    def get_target_index(self, orbit_track: Geocentric) -> npt.NDArray:
        """
        Gets the index of the target at which the body points at each time of
        an orbit track, or -1 if it points along nadir.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite
                position/velocity (in an inertial frame).

        Returns:
            numpy.typing.NDArray: the target indices, with the shape of the
                orbit track's times
        """
        shape = np.shape(orbit_track.t.tt)
        p_m, nadir, _, _ = _compute_nadir_axes(orbit_track, self.velocity_frame)
        p_m, nadir = np.reshape(p_m, (3, -1)), np.reshape(nadir, (3, -1))
        longitude, latitude, elevation = self._get_targets()
        xyz = np.array(wgs84.latlon(latitude, longitude, elevation).itrs_xyz.m)
        up = np.array(
            [
                np.cos(np.radians(latitude)) * np.cos(np.radians(longitude)),
                np.cos(np.radians(latitude)) * np.sin(np.radians(longitude)),
                np.sin(np.radians(latitude)),
            ]
        )
        off_nadir = np.full((len(self.targets), p_m.shape[1]), np.inf)
        for k in range(len(self.targets)):
            los = _compute_line_of_sight(xyz[:, k], p_m)
            distance = np.linalg.norm(los, axis=0)
            angle = np.degrees(
                np.arccos(np.clip(np.sum(los * nadir, axis=0) / distance, -1, 1))
            )
            # above the target's horizon and within the maximum off-nadir angle
            visible = (np.sum(-los * up[:, k : k + 1], axis=0) > 0) & (
                angle <= self.max_off_nadir_angle
            )
            off_nadir[k] = np.where(visible, angle, np.inf)
        index = np.where(
            np.isfinite(off_nadir.min(axis=0)), np.argmin(off_nadir, axis=0), -1
        )
        return np.reshape(index, shape)

    def get_body_angles(
        self, orbit_track: Geocentric
    ) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
        """
        Gets the yaw, roll, and pitch angles (degrees) of the body and whether
        its velocity frame is inertial, at each time of an orbit track (see
        `Attitude.get_body_angles`): the roll and pitch angles that point its
        nadir axis at the target of each time (see `get_target_index`), or
        zero.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite
                position/velocity (in an inertial frame).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
                the yaw, roll, and pitch angles (degrees), and whether the
                velocity frame is inertial
        """
        shape = np.shape(orbit_track.t.tt)
        index = np.reshape(self.get_target_index(orbit_track), -1)
        roll, pitch = np.zeros(index.shape), np.zeros(index.shape)
        pointing = np.flatnonzero(index >= 0)
        if len(pointing) > 0:
            longitude, latitude, elevation = (
                values[index[pointing]] for values in self._get_targets()
            )
            if len(shape) == 0:
                track = orbit_track
                target = wgs84.latlon(latitude[0], longitude[0], elevation[0])
            else:
                track = _index_orbit_track(orbit_track, pointing)
                target = wgs84.latlon(latitude, longitude, elevation)
            # the roll and pitch angles of the target in the yawed body frame
            roll_angle, pitch_angle = compute_view_angles(
                track,
                target,
                self.velocity_frame,
                attitude=Attitude(
                    yaw_angle=self.yaw_angle, velocity_frame=self.velocity_frame
                ),
            )
            roll[pointing] = np.reshape(roll_angle, -1)
            pitch[pointing] = np.reshape(pitch_angle, -1)
        return (
            np.full(shape, float(self.yaw_angle)),
            np.reshape(roll, shape),
            np.reshape(pitch, shape),
            np.full(shape, self.velocity_frame == VelocityFrame.INERTIAL),
        )


class GlintAttitude(BaseModel):
    """
    Attitude of a spacecraft's body that points its nadir axis at the sun
    glint (the specular point, see `tatc.utils.compute_specular_point`), as
    for the glint mode of a greenhouse gas sounder observing the ocean: at
    each time, toward the specular point, or short of it toward nadir by an
    offset angle (in the plane of the nadir axis and the line of sight to
    the specular point), while the Sun is at least `min_solar_elevation`
    above the horizon at the specular point (equal to the satellite's
    elevation angle seen from it) and the body points within
    `max_off_nadir_angle` of nadir, and otherwise along nadir. The body is
    rotated by the yaw angle about the nadir axis and then by the roll and
    pitch angles that point its nadir axis (see
    `tatc.utils.compute_view_angles`), which follows the glint as it moves
    with the satellite and the Sun.

    The fields of regard of its instruments (measured from nadir) should
    include the off-nadir angles at which the body points.
    """

    type: Literal["glint"] = Field(
        default="glint", description="Attitude type discriminator."
    )
    max_off_nadir_angle: float = Field(
        default=60,
        description="Maximum angle (degrees) from nadir at which the body "
        + "points toward the glint.",
        gt=0,
        lt=90,
    )
    offset_angle: float = Field(
        default=0,
        description="Angle (degrees) by which the body points short of the "
        + "specular point, toward nadir, as for an instrument that observes "
        + "near the glint rather than at its center.",
        ge=0,
        lt=90,
    )
    offset_angle_profile: list[tuple[float, float]] | None = Field(
        default=None,
        description="Offset angle (degrees, see `offset_angle`) as a function "
        + "of the specular point's angle from nadir (degrees), as (off-nadir "
        + "angle, offset angle) pairs interpolated linearly (and held at the "
        + "first and last offsets beyond them). If given, it overrides "
        + "`offset_angle`.",
        min_length=1,
    )
    min_solar_elevation: float = Field(
        default=0,
        description="Minimum solar elevation angle (degrees) at the specular "
        + "point at which the body points at it.",
        ge=0,
        lt=90,
    )
    velocity_frame: VelocityFrame = Field(
        default=VelocityFrame.EARTH_FIXED,
        description="Reference frame of the velocity vector that defines the "
        + "along-track axis (see `Attitude.velocity_frame`).",
    )
    yaw_angle: float = Field(
        default=0,
        description="Rotation (degrees) of the body about the nadir axis "
        + "before pointing (see `Attitude.yaw_angle`).",
        ge=-180,
        le=180,
    )

    @field_validator("offset_angle_profile")
    @classmethod
    def _validate_offset_angle_profile(
        cls, profile: list[tuple[float, float]] | None
    ) -> list[tuple[float, float]] | None:
        """
        Sorts an offset angle profile by off-nadir angle and checks its angles.
        """
        if profile is None:
            return None
        for off_nadir, offset in profile:
            if not 0 <= off_nadir < 90:
                raise ValueError("Off-nadir angles must be in [0, 90).")
            if not 0 <= offset < 90:
                raise ValueError("Offset angles must be in [0, 90).")
        return sorted(profile)

    def get_offset_angle(self, off_nadir_angle: npt.ArrayLike) -> npt.NDArray:
        """
        Gets the offset angle (degrees) toward nadir for specular points at
        angles (degrees) from nadir (see `offset_angle_profile`).

        Args:
            off_nadir_angle (numpy.typing.ArrayLike): The specular points'
                angles from nadir (degrees).

        Returns:
            numpy.typing.NDArray: the offset angles (degrees).
        """
        if self.offset_angle_profile is None:
            return np.full(np.shape(off_nadir_angle), float(self.offset_angle))
        angles, offsets = np.array(self.offset_angle_profile, dtype=float).T
        return np.interp(off_nadir_angle, angles, offsets)

    def get_aim_point(self, orbit_track: Geocentric) -> GeographicPosition:
        """
        Gets the point on the WGS 84 ellipsoid at which the body points at
        each time of an orbit track: the specular point, or the point short
        of it toward nadir by the offset angle, where the criteria are met,
        or `NaN` where the body points along nadir.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite
                position/velocity (in an inertial frame).

        Returns:
            skyfield.toposlib.GeographicPosition: the aim points.
        """
        point = compute_specular_point(orbit_track)
        latitude = np.reshape(point.latitude.degrees, -1)
        longitude = np.reshape(point.longitude.degrees, -1)
        p_m, nadir, _, _ = _compute_nadir_axes(orbit_track, self.velocity_frame)
        p_m, nadir = np.reshape(p_m, (3, -1)), np.reshape(nadir, (3, -1))
        los = np.reshape(np.array(point.itrs_xyz.m), (3, -1)) - p_m
        distance = np.linalg.norm(los, axis=0)
        with np.errstate(invalid="ignore"):
            los = los / distance
            off_nadir = np.degrees(
                np.arccos(np.clip(np.sum(los * nadir, axis=0), -1, 1))
            )
            # the elevation of the satellite (and of the Sun) seen from the point
            up = np.array(
                [
                    np.cos(np.radians(latitude)) * np.cos(np.radians(longitude)),
                    np.cos(np.radians(latitude)) * np.sin(np.radians(longitude)),
                    np.sin(np.radians(latitude)),
                ]
            )
            elevation = np.degrees(np.arcsin(np.sum(-los * up, axis=0)))
            # the direction offset toward nadir in the plane of nadir and the line of sight
            angle = np.maximum(off_nadir - self.get_offset_angle(off_nadir), 0)
            across = los - np.sum(los * nadir, axis=0) * nadir
            across = across / np.linalg.norm(across, axis=0)
            direction = np.cos(np.radians(angle)) * nadir + np.sin(
                np.radians(angle)
            ) * np.nan_to_num(across)
            pointing = (
                np.isfinite(latitude)
                & (angle <= self.max_off_nadir_angle)
                & (elevation >= self.min_solar_elevation)
            )
        aim = np.full(p_m.shape, np.nan)
        if np.any(pointing):
            position, found = compute_ellipsoid_intersection(
                p_m[:, pointing], direction[:, pointing]
            )
            aim[:, pointing] = np.where(found, position, np.nan)
        longitude, latitude, _ = rectangular_to_geodetic(aim)
        if np.ndim(orbit_track.t.tt) == 0:
            return wgs84.latlon(latitude[0], longitude[0])
        return wgs84.latlon(latitude, longitude)

    def get_body_angles(
        self, orbit_track: Geocentric
    ) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
        """
        Gets the yaw, roll, and pitch angles (degrees) of the body and whether
        its velocity frame is inertial, at each time of an orbit track (see
        `Attitude.get_body_angles`): the roll and pitch angles that point its
        nadir axis at the aim point (see `get_aim_point`), or zero.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite
                position/velocity (in an inertial frame).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
                the yaw, roll, and pitch angles (degrees), and whether the
                velocity frame is inertial
        """
        shape = np.shape(orbit_track.t.tt)
        point = self.get_aim_point(orbit_track)
        latitude = np.reshape(point.latitude.degrees, -1)
        longitude = np.reshape(point.longitude.degrees, -1)
        roll, pitch = np.zeros(latitude.shape), np.zeros(latitude.shape)
        pointing = np.flatnonzero(np.isfinite(latitude))
        if len(pointing) > 0:
            if len(shape) == 0:
                track = orbit_track
                target = wgs84.latlon(latitude[0], longitude[0])
            else:
                track = _index_orbit_track(orbit_track, pointing)
                target = wgs84.latlon(latitude[pointing], longitude[pointing])
            # the roll and pitch angles of the aim point in the yawed body frame
            roll_angle, pitch_angle = compute_view_angles(
                track,
                target,
                self.velocity_frame,
                attitude=Attitude(
                    yaw_angle=self.yaw_angle, velocity_frame=self.velocity_frame
                ),
            )
            roll[pointing] = np.reshape(roll_angle, -1)
            pitch[pointing] = np.reshape(pitch_angle, -1)
        return (
            np.full(shape, float(self.yaw_angle)),
            np.reshape(roll, shape),
            np.reshape(pitch, shape),
            np.full(shape, self.velocity_frame == VelocityFrame.INERTIAL),
        )


ScheduledAttitudes = Annotated[
    Attitude | TargetAttitude | GlintAttitude, Field(discriminator="type")
]
"""Attitudes that a schedule may switch between."""


class AttitudeSegment(BaseModel):
    """
    An attitude that applies from a start time (until the start of the next
    segment of a schedule).
    """

    start: AwareDatetime = Field(
        ..., description="Time (timezone-aware) from which the attitude applies."
    )
    attitude: ScheduledAttitudes = Field(..., description="The attitude.")


class AttitudeSchedule(BaseModel):
    """
    Schedule of attitudes over time: each segment's attitude applies from its
    start until the start of the next segment, and the first segment's
    attitude also applies before its start. For example, a satellite that
    periodically flips to fly backward (yaw 180 degrees) and forward, or
    that alternates between observing modes.
    """

    type: Literal["schedule"] = Field(
        default="schedule", description="Attitude type discriminator."
    )
    segments: list[AttitudeSegment] = Field(
        ...,
        min_length=1,
        description="The segments, each with a start time and an attitude.",
    )

    @field_validator("segments")
    @classmethod
    def _sort_segments(cls, segments: list[AttitudeSegment]) -> list[AttitudeSegment]:
        """
        Sorts the segments by start time.
        """
        return sorted(segments, key=lambda segment: segment.start)

    @classmethod
    def from_attitudes(
        cls,
        starts: list[datetime],
        attitudes: list[Attitude | TargetAttitude | GlintAttitude],
    ) -> AttitudeSchedule:
        """
        Creates a schedule from the start times and attitudes of its segments.

        Args:
            starts (list[datetime.datetime]): The start times (timezone-aware).
            attitudes (list[Attitude | TargetAttitude | GlintAttitude]): The
                attitude of each segment.

        Returns:
            AttitudeSchedule: the schedule
        """
        if len(starts) != len(attitudes):
            raise ValueError("Each segment needs a start time and an attitude.")
        return cls(
            segments=[
                AttitudeSegment(start=start, attitude=attitude)
                for start, attitude in zip(starts, attitudes)
            ]
        )

    def get_attitude(self, time: datetime) -> Attitude | TargetAttitude | GlintAttitude:
        """
        Gets the attitude that applies at a time.

        Args:
            time (datetime.datetime): The time (timezone-aware).

        Returns:
            Attitude | TargetAttitude | GlintAttitude: the attitude
        """
        index = 0
        for i, segment in enumerate(self.segments):
            if segment.start <= time:
                index = i
        return self.segments[index].attitude

    def get_body_angles(
        self, orbit_track: Geocentric
    ) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
        """
        Gets the yaw, roll, and pitch angles (degrees) of the body and whether
        its velocity frame is inertial, at each time of an orbit track, from
        the attitude of the segment in which each time lies (see
        `Attitude.get_body_angles`).

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The satellite
                position/velocity (in an inertial frame).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
                the yaw, roll, and pitch angles (degrees), and whether the
                velocity frame is inertial
        """
        tt = np.asarray(orbit_track.t.tt, dtype=float)
        shape = np.shape(tt)
        starts = np.atleast_1d(
            orbit_track.t.ts.from_datetimes([s.start for s in self.segments]).tt
        )
        segment = np.maximum(
            np.searchsorted(starts, np.reshape(tt, -1), "right") - 1, 0
        )
        yaw, roll, pitch = (np.zeros(segment.shape) for _ in range(3))
        inertial = np.zeros(segment.shape, dtype=bool)
        for i in np.unique(segment):
            index = np.flatnonzero(segment == i)
            track = (
                orbit_track
                if len(shape) == 0
                else _index_orbit_track(orbit_track, index)
            )
            values = self.segments[i].attitude.get_body_angles(track)
            for array, value in zip((yaw, roll, pitch, inertial), values):
                array[index] = np.reshape(value, -1)
        return tuple(np.reshape(array, shape) for array in (yaw, roll, pitch, inertial))
