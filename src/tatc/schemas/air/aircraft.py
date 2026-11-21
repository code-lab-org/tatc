"""
Object schema for aircraft.

@author: Naeun Kim <hkim439@asu.edu>
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, model_validator
from skyfield.positionlib import Geocentric

from ...utils.time import _to_time
from ..instrument import AllInstruments, Instrument, _check_platform_instruments
from .trajectory import Trajectory


def _to_body_angles(
    crab: npt.NDArray, roll: npt.NDArray, pitch: npt.NDArray
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
    """
    Converts the attitude of an aircraft (its heading relative to the track,
    then its pitch, then its roll: the aerospace sequence) to the yaw, roll,
    and pitch angles of a body in TAT-C (yaw about nadir, positive to the
    left, then roll about the along-track axis, positive tilting nadir to the
    left, then pitch about the rolled cross-track axis, positive tilting nadir
    forward; see `tatc.schemas.Attitude`), which orient the same axes.

    Args:
        crab (numpy.typing.NDArray): The headings (degrees clockwise) relative
            to the track.
        roll (numpy.typing.NDArray): The roll angles (degrees, right wing down
            positive).
        pitch (numpy.typing.NDArray): The pitch angles (degrees, nose up
            positive).

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
            the yaw, roll, and pitch angles (degrees) of the body
    """
    d, r, p = np.radians(crab), np.radians(roll), np.radians(pitch)
    # the body's down axis tilts left (cross-track) by the sine of the roll
    # and forward by the pitch; its left axis sets the yaw
    body_roll = np.arcsin(np.clip(np.sin(r) * np.cos(p), -1, 1))
    body_pitch = np.arctan2(np.sin(p), np.cos(r) * np.cos(p))
    body_yaw = np.arctan2(
        np.sin(r) * np.sin(p) * np.cos(d) - np.cos(r) * np.sin(d),
        np.sin(r) * np.sin(p) * np.sin(d) + np.cos(r) * np.cos(d),
    )
    return np.degrees(body_yaw), np.degrees(body_roll), np.degrees(body_pitch)


class _TrajectoryAttitude:
    """
    Attitude of an aircraft's body along its trajectory, from which the views
    of its instruments are rotated (see `tatc.schemas.Attitude`).
    """

    def __init__(self, trajectory: Trajectory):
        self.trajectory = trajectory

    def get_body_angles(
        self, orbit_track: Geocentric
    ) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
        """
        Gets the yaw, roll, and pitch angles (degrees) of the body, relative to
        the track, and whether its along-track axis is aligned with the
        inertial velocity (never), at each time of a track (see
        `tatc.schemas.Attitude.get_body_angles`).

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The aircraft
                position/velocity (see `Trajectory.get_orbit_track_at_time`).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
                the yaw, roll, and pitch angles (degrees), and whether the
                velocity frame is inertial
        """
        yaw, roll, pitch = _to_body_angles(*self.trajectory.get_angles(orbit_track.t))
        return yaw, roll, pitch, np.zeros(np.shape(yaw), dtype=bool)

    def is_observing(self, orbit_track: Geocentric) -> npt.NDArray[np.bool_]:
        """
        Determines if the instruments observe at each time of a track: at the
        times that the trajectory contains (see `Trajectory.contains`), so
        that the footprints swept over a time step are trimmed to the
        trajectory.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The aircraft
                position/velocity.

        Returns:
            numpy.typing.NDArray: Array of indicators: `True` if the instruments observe.
        """
        return self.trajectory.contains(orbit_track.t)


class Aircraft(BaseModel):
    """
    Aircraft hosting one or more instruments along a trajectory.

    Instruments are oriented as on a satellite whose body is the aircraft's:
    a view at zero roll and pitch looks along the aircraft's down axis, a
    positive roll angle looks to its left (toward the left wing), and a
    positive pitch angle looks forward (toward the nose). With level wings and
    nose, along the track, the down axis points to nadir. The aircraft's roll,
    pitch, and heading (see `Trajectory`) rotate the views of a
    `PointedInstrument`; a generic `Instrument` observes a cone about nadir.

    The roll and pitch of an aircraft (right wing down and nose up positive,
    as in the aerospace convention and in flight recordings such as IWG1) and
    those of a satellite's attitude and an instrument's view (see
    `tatc.schemas.Attitude` and `PointedInstrument`) are the same rotations:
    right-handed about the forward and right axes of a frame whose axes point
    forward, right, and down (an aircraft's body axes, or a satellite's
    local-vertical, local-horizontal axes). A positive roll tilts the down
    axis, and so the views, to the left, and a positive pitch tilts them
    forward: an aircraft rolled 10 degrees shifts a nadir view as an
    instrument's `roll_angle` of 10 degrees does. The yaw differs: an
    aircraft's heading is clockwise (positive turning the nose to the right,
    as is standard for aircraft and spacecraft), whereas the yaw angle of a
    satellite's attitude is positive turning its front to the left. The look
    angles of some planning tools are positive to the right (as are the swath
    offset angles of HyPlan's sensors), the opposite of a roll angle (see
    `tatc.preprocess.airborne.instrument_from_sensor`).
    """

    name: str = Field(..., description="Aircraft name.", examples=["NASA P-3B"])
    trajectory: Trajectory = Field(..., description="Trajectory of the aircraft.")
    instruments: list[AllInstruments] = Field(
        default=[Instrument()],
        min_length=1,
        description="List of hosted instruments.",
    )

    @model_validator(mode="after")
    def check_instruments(self) -> Aircraft:
        """
        Checks that each instrument's view is defined on an aircraft: without
        roll or pitch angle profiles (which vary with an orbit's argument of
        latitude) or an inertial velocity frame (aligned with an orbit plane).
        """
        _check_platform_instruments(self.instruments, "an aircraft", "track")
        return self

    def get_instrument(self, index: int = 0) -> Instrument:
        """
        Gets an instrument carried by this aircraft (see
        `Instrument.bind_attitude`), with its views rotated from the axes of
        the aircraft's body (level along the track if the trajectory gives no
        attitude), which observes at the times that the trajectory contains
        (see `Trajectory.contains`).

        Args:
            index (int): The index of the instrument.

        Returns:
            Instrument: the instrument
        """
        return self.instruments[index].bind_attitude(
            _TrajectoryAttitude(self.trajectory)
        )

    def get_track(self, times: datetime | list[datetime] | None = None) -> Geocentric:
        """
        Gets the track of the aircraft (see
        `Trajectory.get_orbit_track_at_time`), from which its instruments
        project their views as from a satellite's orbit track.

        Args:
            times (datetime.datetime | list[datetime.datetime] | None): The
                time(s). Defaults to the times of the trajectory.

        Returns:
            skyfield.positionlib.Geocentric: The position/velocity.
        """
        return self.trajectory.get_orbit_track_at_time(
            _to_time(self.trajectory.time if times is None else times)
        )
