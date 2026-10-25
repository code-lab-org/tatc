"""
Object schemas for the attitude of space systems: the orientation of a
spacecraft's body relative to its orbit, which orients the views of its
instruments, and schedules of attitudes over time.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime
from typing import Literal

import numpy as np
import numpy.typing as npt
from pydantic import AwareDatetime, BaseModel, Field, field_validator
from skyfield.positionlib import Geocentric

from ..utils.orbital import compute_argument_of_latitude
from ..utils.projection import VelocityFrame
from ..utils.time import _index_orbit_track


def _check_angle_profile(
    profile: list[tuple[float, float]] | None, kind: str
) -> list[tuple[float, float]] | None:
    """
    Checks that an angle profile's arguments of latitude lie in [0, 360) and
    its angles in [-180, 180], and sorts it by argument of latitude.

    Args:
        profile (list[tuple[float, float]] | None): (argument of latitude,
            angle) pairs (degrees).
        kind (str): The kind of angle (for example, "Roll"), for messages.

    Returns:
        list[tuple[float, float]] | None: the sorted profile
    """
    if profile is None:
        return None
    for argument, angle in profile:
        if not 0 <= argument < 360:
            raise ValueError("Arguments of latitude must be in [0, 360).")
        if not -180 <= angle <= 180:
            raise ValueError(f"{kind} angles must be in [-180, 180].")
    return sorted(profile)


def _interpolate_angle_profile(
    profile: list[tuple[float, float]], orbit_track: Geocentric
) -> float | npt.NDArray[np.float64]:
    """
    Interpolates an angle profile linearly and periodically at the
    satellite's argument(s) of latitude.

    Args:
        profile (list[tuple[float, float]]): (argument of latitude, angle) pairs.
        orbit_track (skyfield.positionlib.Geocentric): The satellite
            position/velocity (in an inertial frame).

    Returns:
        float | numpy.typing.NDArray[numpy.float64]: the angle(s) (degrees)
    """
    argument = compute_argument_of_latitude(
        orbit_track.position.m, orbit_track.velocity.m_per_s
    )
    arguments, angles = zip(*profile)
    angle = np.interp(argument, arguments, angles, period=360)
    return float(angle) if np.ndim(angle) == 0 else angle


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


class AttitudeSegment(BaseModel):
    """
    An attitude that applies from a start time (until the start of the next
    segment of a schedule).
    """

    start: AwareDatetime = Field(
        ..., description="Time (timezone-aware) from which the attitude applies."
    )
    attitude: Attitude = Field(..., description="The attitude.")


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
        cls, starts: list[datetime], attitudes: list[Attitude]
    ) -> AttitudeSchedule:
        """
        Creates a schedule from the start times and attitudes of its segments.

        Args:
            starts (list[datetime.datetime]): The start times (timezone-aware).
            attitudes (list[Attitude]): The attitude of each segment.

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

    def get_attitude(self, time: datetime) -> Attitude:
        """
        Gets the attitude that applies at a time.

        Args:
            time (datetime.datetime): The time (timezone-aware).

        Returns:
            Attitude: the attitude
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
        starts = np.array(
            [orbit_track.t.ts.from_datetime(s.start).tt for s in self.segments]
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
