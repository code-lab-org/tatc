"""
Object schema for the trajectories of aircraft.

@author: Naeun Kim <hkim439@asu.edu>
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Annotated

import numpy as np
import numpy.typing as npt
import pandas as pd
from pydantic import AwareDatetime, BaseModel, Field, PrivateAttr, model_validator
from skyfield.constants import DAY_S
from skyfield.positionlib import Geocentric
from skyfield.timelib import Time

from ...utils.ellipsoid import geodetic_to_rectangular, rectangular_to_geodetic
from ...utils.projection import _earth_fixed_to_geocentric
from ...utils.time import _to_time

_REQUIRED_COLUMNS = (
    "time",
    "latitude",
    "longitude",
    "altitude",
    "track",
    "ground_speed",
)
_ATTITUDE_COLUMNS = ("heading", "roll", "pitch")

# speed (meters per second) below which the direction of motion is taken from
# the track rather than from the interpolated velocity (as when stationary)
_MIN_SPEED = 1e-3

Latitude = Annotated[float, Field(ge=-90, le=90)]
Longitude = Annotated[float, Field(ge=-180, le=180)]
Speed = Annotated[float, Field(ge=0)]


def _get_local_axes(
    longitude: npt.NDArray, latitude: npt.NDArray
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Gets the Earth-fixed unit vectors toward the local (geodetic) east and
    north.

    Args:
        longitude (numpy.typing.NDArray): The longitudes (degrees).
        latitude (numpy.typing.NDArray): The geodetic latitudes (degrees).

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the east and north
            unit vectors (shape (3, N))
    """
    lon, lat = np.radians(longitude), np.radians(latitude)
    east = np.array([-np.sin(lon), np.cos(lon), np.zeros_like(lon)])
    north = np.array(
        [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)]
    )
    return east, north


def _interpolate_hermite(
    seconds: npt.NDArray,
    position: npt.NDArray,
    velocity: npt.NDArray,
    times: npt.NDArray,
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Interpolates positions between states along the cubic (Hermite) curves
    through their positions and velocities, and extrapolates them at the
    velocity of the first (or last) state before (or after) the states.

    Args:
        seconds (numpy.typing.NDArray): The times (seconds) of the states,
            in increasing order.
        position (numpy.typing.NDArray): The positions of the states (meters,
            shape (3, N)).
        velocity (numpy.typing.NDArray): The velocities of the states (meters
            per second, shape (3, N)).
        times (numpy.typing.NDArray): The times (seconds) at which to
            interpolate.

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the positions
            (meters) and velocities (meters per second) at the times (shape
            (3, M))
    """
    s, p, v = seconds, position, velocity
    index = np.clip(np.searchsorted(s, times, side="right") - 1, 0, max(len(s) - 2, 0))
    # before the first time (or after the last), continue from its state
    after = times > s[-1]
    index[after] = len(s) - 1
    result = p[:, index] + v[:, index] * (times - s[index])
    rate = v[:, index].copy()
    inside = (times >= s[0]) & ~after & (len(s) > 1)
    if np.any(inside):
        i = index[inside]
        dt = s[i + 1] - s[i]
        u = (times[inside] - s[i]) / dt
        # the cubic Hermite basis functions and their derivatives
        h00, h10 = 2 * u**3 - 3 * u**2 + 1, u**3 - 2 * u**2 + u
        h01, h11 = -2 * u**3 + 3 * u**2, u**3 - u**2
        d00, d10 = 6 * u**2 - 6 * u, 3 * u**2 - 4 * u + 1
        d11 = 3 * u**2 - 2 * u
        result[:, inside] = (
            h00 * p[:, i]
            + h10 * dt * v[:, i]
            + h01 * p[:, i + 1]
            + h11 * dt * v[:, i + 1]
        )
        rate[:, inside] = (
            d00 * (p[:, i] - p[:, i + 1]) / dt + d10 * v[:, i] + d11 * v[:, i + 1]
        )
    return result, rate


class Trajectory(BaseModel):
    """
    Trajectory of an aircraft: its states at a sequence of times (for example,
    from a flight plan, a flight recording, or a simulation). Between the
    times, its position follows the cubic (Hermite) curve through the
    positions and horizontal velocities of the neighboring states, with the
    altitude, roll, pitch, and heading interpolated linearly; before the first
    and after the last time, it continues along the track at the ground speed
    with the attitude held constant. Near a gap between states (as between
    flight legs), the position thus continues along the track.

    The attitude of the aircraft's body follows the aerospace convention:
    from the heading, the body is rotated by the pitch angle (nose up
    positive) and then by the roll angle (right wing down positive). The
    instruments of an aircraft point relative to its body (see
    `tatc.schemas.Aircraft`).
    """

    time: list[AwareDatetime] = Field(
        ...,
        min_length=1,
        description="Times of the states, in increasing order.",
    )
    latitude: list[Latitude] = Field(
        ..., description="Geodetic latitudes (decimal degrees) in WGS 84."
    )
    longitude: list[Longitude] = Field(
        ..., description="Longitudes (decimal degrees) in WGS 84."
    )
    altitude: list[float] = Field(
        ..., description="Altitudes (meters) above the WGS 84 ellipsoid."
    )
    track: list[float] = Field(
        ...,
        description="Directions of motion over the ground (degrees clockwise "
        + "from north).",
    )
    ground_speed: list[Speed] = Field(
        ...,
        description="Speeds (meters per second) over the ground: the "
        + "horizontal speeds of the aircraft (at its altitude, as measured by "
        + "satellite navigation).",
    )
    heading: list[float] | None = Field(
        default=None,
        description="Directions of the aircraft's nose (degrees clockwise from "
        + "north), which differ from the track in a crosswind, or `None` for "
        + "the track.",
    )
    roll: list[float] | None = Field(
        default=None,
        description="Roll angles (degrees, right wing down positive), or `None` "
        + "for level wings.",
    )
    pitch: list[float] | None = Field(
        default=None,
        description="Pitch angles (degrees, nose up positive), or `None` for a "
        + "level nose.",
    )
    # the states as arrays (see `_get_states`)
    _states: dict[str, npt.NDArray] | None = PrivateAttr(default=None)

    @model_validator(mode="after")
    def check_states(self) -> Trajectory:
        """
        Checks that every state has each given quantity and that the times
        increase.
        """
        for name in _REQUIRED_COLUMNS[1:] + _ATTITUDE_COLUMNS:
            values = getattr(self, name)
            if values is not None and len(values) != len(self.time):
                raise ValueError(
                    f"{name} has {len(values)} values for {len(self.time)} times"
                )
        if any(b <= a for a, b in zip(self.time[:-1], self.time[1:])):
            raise ValueError("times must be strictly increasing")
        return self

    @classmethod
    def from_dataframe(cls, data: pd.DataFrame) -> Trajectory:
        """
        Creates a trajectory from a data frame with one row per state and
        columns named as the fields (`time`, `latitude`, `longitude`,
        `altitude`, `track`, and `ground_speed`, and optionally `heading`,
        `roll`, and `pitch`); other columns are ignored. The rows are sorted by
        time, and times without a time zone are taken as UTC.

        Args:
            data (pandas.DataFrame): The states.

        Returns:
            Trajectory: the trajectory
        """
        missing = set(_REQUIRED_COLUMNS).difference(data.columns)
        if missing:
            raise ValueError(f"missing required columns: {sorted(missing)}")
        columns = list(_REQUIRED_COLUMNS) + [
            name for name in _ATTITUDE_COLUMNS if name in data.columns
        ]
        states = data[columns].copy()
        states["time"] = pd.to_datetime(states["time"], utc=True)
        states = states.sort_values("time", kind="stable")
        values = {
            name: states[name].to_numpy(dtype=float).tolist() for name in columns[1:]
        }
        return cls(time=states["time"].dt.to_pydatetime().tolist(), **values)

    def has_attitude(self) -> bool:
        """
        Checks whether the trajectory gives the aircraft's attitude: a heading,
        roll, or pitch angle.

        Returns:
            bool: `True` if any attitude angle is given
        """
        return any(getattr(self, name) is not None for name in _ATTITUDE_COLUMNS)

    def _get_states(self) -> dict[str, npt.NDArray]:
        """
        Gets the states as arrays: the seconds (TT) since the first time and
        the Earth-fixed positions (meters) and horizontal velocities (meters
        per second), and the altitudes, unit vectors along the track, and
        attitude angles (degrees, the heading relative to the track).
        """
        if self._states is None:
            t = _to_time(self.time)
            longitude = np.asarray(self.longitude, dtype=float)
            latitude = np.asarray(self.latitude, dtype=float)
            track = np.radians(self.track)
            east, north = _get_local_axes(longitude, latitude)
            direction = np.sin(track) * east + np.cos(track) * north
            # the heading relative to the track, continuous in time
            crab = (
                np.zeros(len(self.time))
                if self.heading is None
                else np.degrees(
                    np.unwrap(
                        np.radians(
                            np.mod(np.subtract(self.heading, self.track) + 180, 360)
                            - 180
                        )
                    )
                )
            )
            self._states = {
                "seconds": (t.whole - t.whole[0]) * DAY_S
                + (t.tt_fraction - t.tt_fraction[0]) * DAY_S,
                "position": geodetic_to_rectangular(longitude, latitude, self.altitude),
                "velocity": direction * np.asarray(self.ground_speed, dtype=float),
                "direction": direction,
                "altitude": np.asarray(self.altitude, dtype=float),
                "crab": crab,
                "roll": (
                    np.zeros(len(self.time))
                    if self.roll is None
                    else np.asarray(self.roll, dtype=float)
                ),
                "pitch": (
                    np.zeros(len(self.time))
                    if self.pitch is None
                    else np.asarray(self.pitch, dtype=float)
                ),
            }
        return self._states

    def _get_seconds(self, t: Time) -> npt.NDArray:
        """
        Gets the seconds (TT) of times since the first time of the trajectory.

        Args:
            t (skyfield.timelib.Time): The time(s).

        Returns:
            numpy.typing.NDArray: the seconds (with the shape of the times)
        """
        t0 = _to_time(self.time[0])
        return (np.asarray(t.whole) - t0.whole) * DAY_S + (
            np.asarray(t.tt_fraction) - t0.tt_fraction
        ) * DAY_S

    def get_angles(self, t: Time) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray]:
        """
        Gets the attitude angles (degrees) at times, interpolated linearly
        (and held constant before the first and after the last time).

        Args:
            t (skyfield.timelib.Time): The time(s).

        Returns:
            tuple[numpy.typing.NDArray, numpy.typing.NDArray, numpy.typing.NDArray]:
                the heading relative to the track (degrees clockwise), the
                roll angle, and the pitch angle, each with the shape of the
                times
        """
        states = self._get_states()
        seconds = self._get_seconds(t)
        return tuple(
            np.interp(seconds, states["seconds"], states[name])
            for name in ("crab", "roll", "pitch")
        )

    def get_orbit_track_at_time(self, t: Time) -> Geocentric:
        """
        Gets the position and velocity of the aircraft at times, in the form of
        a satellite's orbit track (named as for an orbit, see
        `tatc.schemas.GeneralPerturbationsOrbit.get_orbit_track_at_time`, so
        that the ground track methods accept a trajectory in place of an
        orbit). The velocity is relative to the Earth-fixed frame, so that it
        orients the instruments' views along the track.

        Args:
            t (skyfield.timelib.Time): The time(s), a scalar or a vector.

        Returns:
            skyfield.positionlib.Geocentric: The position/velocity (with shape
                (3,) for a scalar time or (3, N) for N times).
        """
        states = self._get_states()
        seconds = np.atleast_1d(self._get_seconds(t))
        s = states["seconds"]
        position, velocity = _interpolate_hermite(
            s, states["position"], states["velocity"], seconds
        )
        # set the altitude, interpolated linearly
        longitude, latitude, _ = rectangular_to_geodetic(position)
        altitude = np.interp(seconds, s, states["altitude"])
        position = geodetic_to_rectangular(longitude, latitude, altitude)
        # where (nearly) stationary, move slowly along the interpolated track
        slow = np.linalg.norm(velocity, axis=0) < _MIN_SPEED
        if np.any(slow):
            direction = np.array(
                [np.interp(seconds[slow], s, states["direction"][k]) for k in range(3)]
            )
            velocity[:, slow] = (
                _MIN_SPEED * direction / np.linalg.norm(direction, axis=0)
            )
        if np.ndim(t.tt) == 0:
            position, velocity = position[:, 0], velocity[:, 0]
        return _earth_fixed_to_geocentric(position, velocity, t)
