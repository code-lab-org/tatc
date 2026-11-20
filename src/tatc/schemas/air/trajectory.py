"""
Object schema for the trajectories of aircraft.

@author: Naeun Kim <hkim439@asu.edu>
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Sequence
from datetime import datetime, timedelta
from typing import Annotated

import numpy as np
import numpy.typing as npt
import pandas as pd
from pydantic import (
    AwareDatetime,
    BaseModel,
    ConfigDict,
    Field,
    PrivateAttr,
    model_validator,
)
from pyproj import Geod
from shapely.geometry import LineString
from skyfield.constants import DAY_S
from skyfield.positionlib import Geocentric
from skyfield.timelib import Time

from ...constants import EARTH_ECCENTRICITY, EARTH_EQUATORIAL_RADIUS
from ...utils.ellipsoid import geodetic_to_rectangular, rectangular_to_geodetic
from ...utils.projection import _earth_fixed_to_geocentric
from ...utils.time import _to_time

_REQUIRED_COLUMNS = ("time", "latitude", "longitude", "altitude")
_MOTION_COLUMNS = ("track", "ground_speed")
_ATTITUDE_COLUMNS = ("heading", "roll", "pitch")

# speed (meters per second) below which the direction of motion is taken from
# the track rather than from the interpolated velocity (as when stationary)
_MIN_SPEED = 1e-3

# tolerance (seconds) of the times at the ends of the parts of a trajectory,
# which absorbs their rounding
_TIME_TOLERANCE = 1e-6

_GEOD = Geod(ellps="WGS84")

Finite = Annotated[float, Field(allow_inf_nan=False)]
Latitude = Annotated[float, Field(ge=-90, le=90, allow_inf_nan=False)]
Longitude = Annotated[float, Field(ge=-180, le=180, allow_inf_nan=False)]
Speed = Annotated[float, Field(ge=0, allow_inf_nan=False)]


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


def _get_radius_of_curvature(
    latitude: npt.NDArray, azimuth: npt.NDArray
) -> npt.NDArray:
    """
    Gets the radius of curvature (meters) of the WGS 84 ellipsoid in an
    azimuth at a latitude (by Euler's theorem, from its meridian and prime
    vertical radii of curvature), by which a horizontal distance at an
    altitude scales to the distance over the ellipsoid below.

    Args:
        latitude (numpy.typing.NDArray): The geodetic latitudes (degrees).
        azimuth (numpy.typing.NDArray): The azimuths (degrees clockwise from
            north).

    Returns:
        numpy.typing.NDArray: the radii of curvature (meters)
    """
    lat, azi = np.radians(latitude), np.radians(azimuth)
    e2 = EARTH_ECCENTRICITY**2
    w = 1 - e2 * np.sin(lat) ** 2
    meridian = EARTH_EQUATORIAL_RADIUS * (1 - e2) / w**1.5
    prime = EARTH_EQUATORIAL_RADIUS / np.sqrt(w)
    return 1 / (np.cos(azi) ** 2 / meridian + np.sin(azi) ** 2 / prime)


class Trajectory(BaseModel):
    """
    Trajectory of an aircraft: its states at a sequence of times (for example,
    from a flight plan, a flight recording, or a simulation, read from a
    table with `from_dataframe`, or flown through waypoints at a constant
    speed with `from_waypoints`). Between the
    times, its position follows the cubic (Hermite) curve through the
    positions and horizontal velocities of the neighboring states, with the
    altitude, roll, pitch, and heading interpolated linearly; before the first
    and after the last time, it continues along the track at the ground speed
    with the attitude held constant. The trajectory gives the aircraft's
    state only from its first to its last time and outside the gaps between
    states longer than `max_gap` (see `contains`), as between flight legs or
    over missing data: the analyses omit other times.

    The track and ground speed, if not given, are those of the horizontal
    velocity estimated from the positions of the neighboring states (by
    finite differences within each part of the trajectory between gaps).

    The attitude of the aircraft's body follows the aerospace convention:
    from the heading, the body is rotated by the pitch angle (nose up
    positive) and then by the roll angle (right wing down positive). The
    instruments of an aircraft point relative to its body (see
    `tatc.schemas.Aircraft`, which compares these angles with those of a
    satellite's attitude).
    """

    model_config = ConfigDict(validate_assignment=True)

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
    altitude: list[Finite] = Field(
        ...,
        description="Altitudes (meters) above the WGS 84 ellipsoid: heights "
        + "above mean sea level plus the geoid height (see "
        + "`tatc.utils.compute_geoid_height`).",
    )
    track: list[Finite] | None = Field(
        default=None,
        description="Directions of motion over the ground (degrees clockwise "
        + "from north), or `None` to estimate them from the positions.",
    )
    ground_speed: list[Speed] | None = Field(
        default=None,
        description="Speeds (meters per second) over the ground: the "
        + "horizontal speeds of the aircraft (at its altitude, as measured by "
        + "satellite navigation), or `None` to estimate them from the "
        + "positions.",
    )
    heading: list[Finite] | None = Field(
        default=None,
        description="Directions of the aircraft's nose (degrees clockwise from "
        + "north), which differ from the track in a crosswind, or `None` for "
        + "the track.",
    )
    roll: list[Finite] | None = Field(
        default=None,
        description="Roll angles (degrees, right wing down positive), or `None` "
        + "for level wings.",
    )
    pitch: list[Finite] | None = Field(
        default=None,
        description="Pitch angles (degrees, nose up positive), or `None` for a "
        + "level nose.",
    )
    max_gap: timedelta | None = Field(
        default=timedelta(seconds=60),
        description="Longest time between consecutive states over which the "
        + "trajectory is interpolated: a longer gap (as between flight legs, or "
        + "over missing data) separates the parts of the trajectory, and the "
        + "times within it are not part of the trajectory. `None` for no gaps.",
    )
    # the states as arrays (see `_get_states`)
    _states: dict[str, npt.NDArray] | None = PrivateAttr(default=None)

    @model_validator(mode="after")
    def check_states(self) -> Trajectory:
        """
        Checks that every state has each given quantity, that the times
        increase, that the maximum gap is positive, and that a state whose
        track or ground speed is estimated has a neighbor within the maximum
        gap. Clears the states computed from the fields (see `_get_states`),
        which are computed again after a field is assigned.
        """
        for name in _REQUIRED_COLUMNS[1:] + _MOTION_COLUMNS + _ATTITUDE_COLUMNS:
            values = getattr(self, name)
            if values is not None and len(values) != len(self.time):
                raise ValueError(
                    f"{name} has {len(values)} values for {len(self.time)} times"
                )
        if any(b <= a for a, b in zip(self.time[:-1], self.time[1:])):
            raise ValueError("times must be strictly increasing")
        if self.max_gap is not None and self.max_gap <= timedelta(0):
            raise ValueError(f"max_gap must be positive, not {self.max_gap}")
        if self.track is None or self.ground_speed is None:
            for part in self._get_parts():
                if len(part) == 1:
                    raise ValueError(
                        f"the track and ground speed at {self.time[part[0]]} "
                        "cannot be estimated without a neighboring state within "
                        "max_gap: give track and ground_speed"
                    )
        self._states = None
        return self

    @classmethod
    def from_dataframe(
        cls,
        data: pd.DataFrame,
        max_gap: timedelta | None = timedelta(seconds=60),
    ) -> Trajectory:
        """
        Creates a trajectory from a data frame with one row per state and
        columns named as the fields (`time`, `latitude`, `longitude`, and
        `altitude`, and optionally `track`, `ground_speed`, `heading`, `roll`,
        and `pitch`); other columns are ignored. The rows are sorted by time,
        and times without a time zone are taken as UTC. An optional column
        without any value is ignored, and the rows missing a value of the
        other columns (as in a recording with missing data) are omitted.

        Args:
            data (pandas.DataFrame): The states.
            max_gap (datetime.timedelta | None): The longest time between
                consecutive states over which the trajectory is interpolated
                (see `Trajectory.max_gap`).

        Returns:
            Trajectory: the trajectory
        """
        missing = set(_REQUIRED_COLUMNS).difference(data.columns)
        if missing:
            raise ValueError(f"missing required columns: {sorted(missing)}")
        columns = list(_REQUIRED_COLUMNS) + [
            name
            for name in _MOTION_COLUMNS + _ATTITUDE_COLUMNS
            if name in data.columns and data[name].notna().any()
        ]
        states = data[columns].copy()
        states["time"] = pd.to_datetime(states["time"], utc=True)
        states = states.dropna().sort_values("time", kind="stable")
        values = {
            name: states[name].to_numpy(dtype=float).tolist() for name in columns[1:]
        }
        return cls(
            time=states["time"].dt.to_pydatetime().tolist(), max_gap=max_gap, **values
        )

    @classmethod
    def from_waypoints(
        cls,
        waypoints: Sequence[tuple[float, float]] | LineString,
        start: datetime,
        speed: float,
        altitude: float | Sequence[float],
        time_step: timedelta = timedelta(seconds=1),
    ) -> Trajectory:
        """
        Creates the trajectory of a flight through waypoints, as planned for
        an underflight or a survey: from the first waypoint at the `start`,
        along the geodesic (on the WGS 84 ellipsoid) to each next waypoint at
        a constant ground speed, turning at each waypoint without slowing. Its
        states are sampled at each time step from the start, and at the last
        waypoint, with the track along the geodesic and the ground speed. The
        trajectory has no gaps (`max_gap` is `None`) and no attitude (level
        wings and nose, heading along the track): for turns, climbs, or an
        attitude, see `tatc.preprocess.sample_flight_plan` and
        `tatc.preprocess.estimate_attitude`.

        Args:
            waypoints (Sequence[tuple[float, float]] | shapely.geometry.LineString):
                The waypoints, as (longitude, latitude) pairs (decimal degrees)
                or the vertices of a line, at least two distinct.
            start (datetime.datetime): The time (timezone-aware) at the first
                waypoint.
            speed (float): The ground speed (meters per second): the
                horizontal speed at the aircraft's altitude (see
                `Trajectory.ground_speed`).
            altitude (float | Sequence[float]): The altitude (meters above the
                WGS 84 ellipsoid), or that at each waypoint, interpolated
                linearly over each leg.
            time_step (datetime.timedelta): The time between states.

        Returns:
            Trajectory: the trajectory
        """
        if start.utcoffset() is None:
            raise ValueError("start must be a timezone-aware datetime")
        if not speed > 0:
            raise ValueError(f"speed must be positive, not {speed}")
        if not time_step > timedelta(0):
            raise ValueError(f"time_step must be positive, not {time_step}")
        points = np.asarray(
            waypoints.coords if isinstance(waypoints, LineString) else waypoints,
            dtype=float,
        )[:, :2]
        altitudes = np.broadcast_to(np.asarray(altitude, dtype=float), len(points))
        # omit repeated waypoints, which do not move the aircraft
        distinct = np.r_[True, np.any(points[1:] != points[:-1], axis=1)]
        points, altitudes = points[distinct], altitudes[distinct]
        if len(points) < 2:
            raise ValueError("waypoints must include at least two distinct points")
        longitude, latitude = points[:, 0], points[:, 1]
        azimuth, _, length = _GEOD.inv(
            longitude[:-1], latitude[:-1], longitude[1:], latitude[1:]
        )
        # the time to fly each leg at the ground speed at its altitude: the
        # distance over the ellipsoid below scales with its radius of
        # curvature in the leg's azimuth
        radius = _get_radius_of_curvature(latitude[:-1], azimuth)
        mean_altitude = (altitudes[:-1] + altitudes[1:]) / 2
        leg_seconds = length * (radius + mean_altitude) / radius / speed
        leg_start = np.concatenate([[0], np.cumsum(leg_seconds)])
        seconds = np.arange(0, leg_start[-1], time_step.total_seconds())
        if leg_start[-1] - seconds[-1] > _TIME_TOLERANCE:
            seconds = np.append(seconds, leg_start[-1])
        leg = np.clip(
            np.searchsorted(leg_start, seconds, side="right") - 1, 0, len(length) - 1
        )
        fraction = np.clip((seconds - leg_start[leg]) / leg_seconds[leg], 0, 1)
        state_longitude, state_latitude, back_azimuth = _GEOD.fwd(
            longitude[leg], latitude[leg], azimuth[leg], fraction * length[leg]
        )
        return cls(
            time=[start + timedelta(seconds=float(value)) for value in seconds],
            latitude=np.asarray(state_latitude).tolist(),
            longitude=np.asarray(state_longitude).tolist(),
            altitude=(
                altitudes[leg] + fraction * (altitudes[leg + 1] - altitudes[leg])
            ).tolist(),
            track=((np.asarray(back_azimuth) + 180) % 360).tolist(),
            ground_speed=[float(speed)] * len(seconds),
            max_gap=None,
        )

    def has_attitude(self) -> bool:
        """
        Checks whether the trajectory gives the aircraft's attitude: a heading,
        roll, or pitch angle.

        Returns:
            bool: `True` if any attitude angle is given
        """
        return any(getattr(self, name) is not None for name in _ATTITUDE_COLUMNS)

    def _get_parts(self) -> list[npt.NDArray]:
        """
        Gets the indices of the states of each part of the trajectory between
        gaps longer than `max_gap`.

        Returns:
            list[numpy.typing.NDArray]: the indices of each part
        """
        index = np.arange(len(self.time))
        if self.max_gap is None:
            return [index]
        gaps = [
            i + 1
            for i, (a, b) in enumerate(zip(self.time[:-1], self.time[1:]))
            if b - a > self.max_gap
        ]
        return np.split(index, gaps)

    def contains(
        self, times: datetime | list[datetime] | Time
    ) -> npt.NDArray[np.bool_]:
        """
        Determines whether the trajectory gives the aircraft's state at times:
        from its first to its last time, and not within a gap between states
        longer than `max_gap`.

        Args:
            times (datetime.datetime | list[datetime.datetime] | skyfield.timelib.Time):
                The time(s).

        Returns:
            numpy.typing.NDArray: Array of indicators (with the shape of the
                times): `True` if the trajectory contains the time.
        """
        t = times if isinstance(times, Time) else _to_time(times)
        seconds = self._get_seconds(t)
        states = self._get_states()
        first, last = states["part_first"], states["part_last"]
        part = np.searchsorted(first, seconds, side="right") - 1
        return (part >= 0) & (seconds <= last[np.maximum(part, 0)])

    def _get_states(self) -> dict[str, npt.NDArray]:
        """
        Gets the states as arrays: the seconds (TT) since the first time and
        the Earth-fixed positions (meters) and horizontal velocities (meters
        per second), and the altitudes, unit vectors along the track, and
        attitude angles (degrees, the heading relative to the track), and the
        first and last seconds of each part of the trajectory between gaps
        (within a small tolerance that absorbs the rounding of times).
        """
        if self._states is None:
            t = _to_time(self.time)
            seconds = np.atleast_1d(
                (t.whole - t.whole[0]) * DAY_S
                + (t.tt_fraction - t.tt_fraction[0]) * DAY_S
            )
            longitude = np.asarray(self.longitude, dtype=float)
            latitude = np.asarray(self.latitude, dtype=float)
            position = geodetic_to_rectangular(longitude, latitude, self.altitude)
            east, north = _get_local_axes(longitude, latitude)
            if self.track is None or self.ground_speed is None:
                # the horizontal velocity, by finite differences of the
                # positions within each part of the trajectory
                velocity = np.zeros_like(position)
                for part in self._get_parts():
                    velocity[:, part] = np.gradient(
                        position[:, part], seconds[part], axis=1
                    )
                v_east = np.sum(velocity * east, axis=0)
                v_north = np.sum(velocity * north, axis=0)
            track = (
                np.arctan2(v_east, v_north)
                if self.track is None
                else np.radians(self.track)
            )
            speed = (
                np.hypot(v_east, v_north)
                if self.ground_speed is None
                else np.asarray(self.ground_speed, dtype=float)
            )
            direction = np.sin(track) * east + np.cos(track) * north
            # the heading relative to the track, continuous in time
            crab = (
                np.zeros(len(self.time))
                if self.heading is None
                else np.degrees(
                    np.unwrap(
                        np.radians(
                            np.mod(
                                np.subtract(self.heading, np.degrees(track)) + 180, 360
                            )
                            - 180
                        )
                    )
                )
            )
            self._states = {
                "seconds": seconds,
                "part_first": np.array([seconds[p[0]] for p in self._get_parts()])
                - _TIME_TOLERANCE,
                "part_last": np.array([seconds[p[-1]] for p in self._get_parts()])
                + _TIME_TOLERANCE,
                "position": position,
                "velocity": direction * speed,
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
