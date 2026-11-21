"""
Preprocessing utilities that derive the inputs of airborne analyses (the
states of a `tatc.schemas.Trajectory` and the `tatc.schemas.PointedInstrument`
of a sensor) from flight plans, flight recordings, and sensor models, such as
those of HyPlan (https://github.com/ryanpavlick/hyplan, installed with
`pip install tatc[airborne]`). These utilities read HyPlan's outputs (data
frames and sensor objects) without importing it.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from pyproj import Geod

from ..schemas import PointedInstrument
from ..utils.geoid import compute_geoid_height

_GEOD = Geod(ellps="WGS84")
_METERS_PER_FOOT = 0.3048
_METERS_PER_SECOND_PER_KNOT = 1852 / 3600
_GRAVITY = 9.80665  # standard gravity (meters per second squared)
_MIN_SPEED = 1e-6  # meters per second, below which the pitch is level

_STATE_COLUMNS = ["time", "latitude", "longitude", "altitude", "track", "ground_speed"]


def sample_flight_plan(
    plan: pd.DataFrame,
    takeoff_time: datetime,
    time_step: timedelta = timedelta(seconds=1),
    geoid_height: float | None = None,
) -> pd.DataFrame:
    """
    Samples the states of an aircraft at regular times along a flight plan
    computed by HyPlan (`hyplan.compute_flight_plan`: one row per segment,
    with its path as a line string, its duration in minutes, and its
    altitudes in feet above mean sea level). Each segment is flown at a
    constant speed along its path (including the turns of connecting
    segments), with the altitude changing linearly. The track is the
    direction of the path at each state, and the heading the track plus the
    crab angle that HyPlan computes for a flight line in a wind (or the
    track elsewhere). The roll and pitch are not given (see
    `estimate_attitude`).

    Args:
        plan (pandas.DataFrame): The flight plan, with columns `geometry`,
            `time_to_segment`, `start_altitude`, and `end_altitude`, and
            optionally `segment_type`, `segment_name`, and `crab_angle_deg`.
            Segments without altitudes are skipped (their time elapses).
        takeoff_time (datetime.datetime): The time at which the plan starts
            (the takeoff time given to `hyplan.compute_flight_plan`).
        time_step (datetime.timedelta): The time between states.
        geoid_height (float | None): The height (meters) of the geoid above
            the WGS 84 ellipsoid, added to the altitudes above mean sea level
            to give altitudes above the ellipsoid: a constant (for example,
            about -30 meters near Phoenix, Arizona), or `None` for the
            EGM2008 geoid at each state (see
            `tatc.utils.compute_geoid_height`).

    Returns:
        pandas.DataFrame: The states, with columns `time`, `latitude`,
            `longitude`, `altitude` (meters above the WGS 84 ellipsoid),
            `track`, `ground_speed` (meters per second), `heading`,
            `segment_type`, and `segment_name`, read by
            `tatc.schemas.Trajectory.from_dataframe`.
    """
    step = time_step.total_seconds()
    if step <= 0:
        raise ValueError(f"time_step must be positive, not {time_step}")
    duration = pd.to_numeric(plan["time_to_segment"]).fillna(0.0).to_numpy() * 60.0
    ends = np.cumsum(duration)
    starts = ends - duration
    seconds = np.arange(0.0, ends[-1] + 1e-9, step) if len(ends) > 0 else []
    # the segment of each time (the later one at a boundary, except the end)
    segment = np.clip(np.searchsorted(ends, seconds, side="right"), 0, len(plan) - 1)
    parts = []
    for i, row in enumerate(plan.itertuples(index=False)):
        t = seconds[segment == i]
        if t.size == 0 or duration[i] <= 0:
            continue
        if pd.isna(row.start_altitude) or pd.isna(row.end_altitude):
            continue
        fraction = (t - starts[i]) / duration[i]
        lon, lat, track, length = _sample_path(row.geometry, fraction)
        altitude = _METERS_PER_FOOT * (
            row.start_altitude + fraction * (row.end_altitude - row.start_altitude)
        ) + _get_geoid_height(lon, lat, geoid_height)
        crab = getattr(row, "crab_angle_deg", np.nan)
        heading = track + (0.0 if pd.isna(crab) else float(crab))
        parts.append(
            pd.DataFrame(
                {
                    "time": pd.Timestamp(takeoff_time) + pd.to_timedelta(t, unit="s"),
                    "latitude": lat,
                    "longitude": lon,
                    "altitude": altitude,
                    "track": track,
                    "ground_speed": length / duration[i],
                    "heading": np.mod(heading, 360.0),
                    "segment_type": getattr(row, "segment_type", None),
                    "segment_name": getattr(row, "segment_name", None),
                }
            )
        )
    if len(parts) == 0:
        return pd.DataFrame(
            columns=_STATE_COLUMNS + ["heading", "segment_type", "segment_name"]
        )
    states = pd.concat(parts, ignore_index=True)
    states["time"] = pd.to_datetime(states["time"], utc=True)
    return states


def _get_geoid_height(
    longitude: npt.ArrayLike, latitude: npt.ArrayLike, geoid_height: float | None
) -> float | npt.NDArray:
    """
    Gets the height (meters) of the geoid above the WGS 84 ellipsoid: a
    constant, or, if `None`, that of the EGM2008 geoid at positions (where
    given; see `tatc.utils.compute_geoid_height`).

    Args:
        longitude (numpy.typing.ArrayLike): The longitudes (degrees).
        latitude (numpy.typing.ArrayLike): The latitudes (degrees).
        geoid_height (float | None): The constant height, or `None`.

    Returns:
        float | numpy.typing.NDArray: the geoid height(s)
    """
    if geoid_height is not None:
        return geoid_height
    longitude = np.asarray(longitude, dtype=float)
    latitude = np.asarray(latitude, dtype=float)
    known = np.isfinite(longitude) & np.isfinite(latitude)
    result = np.full(longitude.shape, np.nan)
    result[known] = compute_geoid_height(longitude[known], latitude[known])
    return result


def _sample_path(
    path: Any, fraction: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Samples a path (a line string of longitudes and latitudes) at fractions
    of its length along the geodesics between its vertices.

    Args:
        path (shapely.LineString): The path.
        fraction (numpy.ndarray): The fractions (0 to 1) of its length.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, float]: the
            longitudes, latitudes, and directions (degrees clockwise from
            north) at the fractions, and the length (meters) of the path
    """
    lon, lat = np.asarray(path.coords)[:, :2].T
    azimuth, _, length = _GEOD.inv(lon[:-1], lat[:-1], lon[1:], lat[1:])
    keep = np.asarray(length) > 0
    if not np.any(keep):
        # a path without length (as a hold): stationary at its start
        n = len(fraction)
        return np.full(n, lon[0]), np.full(n, lat[0]), np.zeros(n), 0.0
    azimuth, length = np.asarray(azimuth)[keep], np.asarray(length)[keep]
    start_lon, start_lat = lon[:-1][keep], lat[:-1][keep]
    cumulative = np.concatenate([[0.0], np.cumsum(length)])
    along = fraction * cumulative[-1]
    k = np.clip(
        np.searchsorted(cumulative, along, side="right") - 1, 0, len(length) - 1
    )
    sample_lon, sample_lat, back = _GEOD.fwd(
        start_lon[k], start_lat[k], azimuth[k], along - cumulative[k]
    )
    return (
        np.asarray(sample_lon),
        np.asarray(sample_lat),
        np.mod(np.asarray(back) + 180.0, 360.0),
        float(cumulative[-1]),
    )


def estimate_attitude(
    states: pd.DataFrame,
    smoothing: timedelta | None = timedelta(seconds=15),
    max_gap: timedelta | None = timedelta(seconds=60),
) -> pd.DataFrame:
    """
    Estimates the roll and pitch of an aircraft from its states (for
    example, of a flight plan, see `sample_flight_plan`, or of a recording
    without attitude): the roll of a coordinated turn in still air,
    `atan(V dψ/dt / g)` for a ground speed `V`, a rate of change of the
    track `dψ/dt`, and the standard gravity `g`, and the pitch of the flight
    path, `atan2(dh/dt, V)` for a rate of climb `dh/dt` (without an angle of
    attack). The rates are computed within each part of the states between
    gaps longer than `max_gap` (a state without a neighbor within `max_gap`
    is level), and smoothed by a centered moving average of duration
    `smoothing` applied twice (a triangular kernel), which reduces the
    fluctuations of sampled tracks but blurs short transitions. In a wind,
    the roll of a turn depends on the rate of change of the heading at the
    airspeed instead, which this does not model. The pitch is that of the
    flight path, so that a climb or descent much faster than the ground
    speed (as of a hovering aircraft) pitches the nose toward vertical.

    Args:
        states (pandas.DataFrame): The states, with columns `time`,
            `altitude`, `track`, and `ground_speed`, without missing values
            (rows missing a value can be omitted, which leaves a gap).
        smoothing (datetime.timedelta | None): The duration of the moving
            average, or `None` for no smoothing.
        max_gap (datetime.timedelta | None): The longest time between
            consecutive states over which rates are computed (see
            `tatc.schemas.Trajectory.max_gap`), or `None` for no gaps.

    Returns:
        pandas.DataFrame: A copy of the states, with the columns `roll`
            (degrees, right wing down positive) and `pitch` (degrees, nose
            up positive).
    """
    missing = {"time", "altitude", "track", "ground_speed"}.difference(states.columns)
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    incomplete = [
        name
        for name in ("time", "altitude", "track", "ground_speed")
        if states[name].isna().any()
    ]
    if incomplete:
        raise ValueError(f"missing values in columns: {incomplete}")
    result = states.copy()
    time = pd.to_datetime(result["time"], utc=True)
    if not time.is_monotonic_increasing or time.duplicated().any():
        raise ValueError("times must be strictly increasing")
    seconds = (time - time.iloc[0]).dt.total_seconds().to_numpy()
    turn_rate, climb_rate = np.zeros(len(result)), np.zeros(len(result))
    breaks = (
        np.flatnonzero(np.diff(seconds) > max_gap.total_seconds()) + 1
        if max_gap is not None
        else []
    )
    for part in np.split(np.arange(len(result)), breaks):
        if len(part) < 2:
            continue
        t = seconds[part]
        track = np.unwrap(np.radians(result["track"].to_numpy(dtype=float)[part]))
        altitude = result["altitude"].to_numpy(dtype=float)[part]
        turn_rate[part] = _smooth(np.gradient(track, t), time.iloc[part], smoothing)
        climb_rate[part] = _smooth(np.gradient(altitude, t), time.iloc[part], smoothing)
    speed = result["ground_speed"].to_numpy(dtype=float)
    result["roll"] = np.degrees(np.arctan(speed * turn_rate / _GRAVITY))
    result["pitch"] = np.degrees(np.arctan2(climb_rate, np.maximum(speed, _MIN_SPEED)))
    return result


def _smooth(
    values: np.ndarray, time: pd.Series, smoothing: timedelta | None
) -> np.ndarray:
    """
    Smooths values by a centered moving average over a duration, applied
    twice (a triangular kernel).

    Args:
        values (numpy.ndarray): The values.
        time (pandas.Series): Their times.
        smoothing (datetime.timedelta | None): The duration, or `None`.

    Returns:
        numpy.ndarray: the smoothed values
    """
    if smoothing is None or smoothing <= timedelta(0):
        return values
    series = pd.Series(values, index=pd.DatetimeIndex(time))
    for _ in range(2):
        series = series.rolling(smoothing, center=True, min_periods=1).mean()
    return series.to_numpy()


def states_from_iwg1(
    data: pd.DataFrame, geoid_height: float | None = None
) -> pd.DataFrame:
    """
    Gets the states of an aircraft from a flight recording in the IWG1
    format (NASA's interchange format of aircraft navigation and state data),
    as read by HyPlan (`hyplan.aircraft.iwg1.load_iwg1`): its GPS position
    and altitude, ground speed and track, and measured heading, roll, and
    pitch. Rows missing a value are omitted, and an attitude angle without
    any value is left out.

    The altitude above the WGS 84 ellipsoid is that recorded (IWG1's "WGS84
    Altitude", in a column `altitude_wgs84_ft`, which HyPlan does not read:
    add it, in feet, to use it), or otherwise the GPS altitude above mean
    sea level plus the geoid height. The recorded altitude avoids the
    difference between the geoid model of the receiver, which gave the
    altitude above mean sea level, and that added back.

    Args:
        data (pandas.DataFrame): The recording, with columns `timestamp`
            (UTC), `latitude`, `longitude`, `altitude_gps_ft` (feet above
            mean sea level) or `altitude_wgs84_ft` (feet above the WGS 84
            ellipsoid), or both, `groundspeed` (knots), and `track`, and
            optionally `true_heading`, `roll_deg`, and `pitch_deg`.
        geoid_height (float | None): The height (meters) of the geoid above
            the WGS 84 ellipsoid, added to the altitudes above mean sea level
            (where no altitude above the ellipsoid is recorded) to give
            altitudes above the ellipsoid: a constant (for example, about -30
            meters near Phoenix, Arizona), or `None` for the EGM2008 geoid at
            each state (see `tatc.utils.compute_geoid_height`).

    Returns:
        pandas.DataFrame: The states, with columns `time`, `latitude`,
            `longitude`, `altitude` (meters above the WGS 84 ellipsoid),
            `track`, `ground_speed` (meters per second), and those of the
            attitude angles given (`heading`, `roll`, and `pitch`), read by
            `tatc.schemas.Trajectory.from_dataframe`.
    """
    missing = {
        "timestamp",
        "latitude",
        "longitude",
        "groundspeed",
        "track",
    }.difference(data.columns)
    if not {"altitude_gps_ft", "altitude_wgs84_ft"} & set(data.columns):
        missing.add("altitude_gps_ft")
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    states = pd.DataFrame(
        {
            "time": pd.to_datetime(data["timestamp"], utc=True),
            "latitude": data["latitude"],
            "longitude": data["longitude"],
            "altitude": _get_iwg1_altitude(data, geoid_height),
            "track": data["track"],
            "ground_speed": data["groundspeed"] * _METERS_PER_SECOND_PER_KNOT,
        }
    )
    for name, column in (
        ("heading", "true_heading"),
        ("roll", "roll_deg"),
        ("pitch", "pitch_deg"),
    ):
        if column in data.columns and data[column].notna().any():
            states[name] = data[column]
    return states.dropna().reset_index(drop=True)


def _get_iwg1_altitude(data: pd.DataFrame, geoid_height: float | None) -> pd.Series:
    """
    Gets the altitudes (meters) above the WGS 84 ellipsoid of a flight
    recording (see `states_from_iwg1`): those recorded above the ellipsoid
    (`altitude_wgs84_ft`), where given, or otherwise those above mean sea
    level (`altitude_gps_ft`) plus the geoid height.

    Args:
        data (pandas.DataFrame): The recording.
        geoid_height (float | None): The height (meters) of the geoid above
            the ellipsoid, or `None` for that of the EGM2008 geoid.

    Returns:
        pandas.Series: the altitudes (meters), `NaN` where none is recorded
    """
    altitude = pd.Series(np.nan, index=data.index)
    if "altitude_gps_ft" in data.columns:
        altitude = data["altitude_gps_ft"].astype(float) * _METERS_PER_FOOT + (
            _get_geoid_height(data["longitude"], data["latitude"], geoid_height)
        )
    if "altitude_wgs84_ft" in data.columns:
        recorded = data["altitude_wgs84_ft"].astype(float) * _METERS_PER_FOOT
        altitude = recorded.fillna(altitude)
    return altitude


def instrument_from_sensor(
    sensor: Any,
    name: str | None = None,
    along_track_field_of_view: float | None = None,
) -> PointedInstrument:
    """
    Gets the instrument of a sensor model of HyPlan, as carried by an
    aircraft (see `tatc.schemas.Aircraft`):

    - A scanning sensor, whose swath lies between the cross-track angles
      of its `swath_offset_angles` (positive to the right; as a line
      scanner, scanning lidar, or side-looking radar), is a swept
      `PointedInstrument` that scans across the track (see
      `PointedInstrument.is_swept`): its cross-track field of view spans the
      swath and its roll angle (positive to the left, the opposite sign)
      centers it.
    - A frame camera (with fields of view `fov_x` across and `fov_y` along
      the track, tilted by `tilt_angle` toward `tilt_direction`) is a
      rectangular `PointedInstrument` whose roll and pitch angles point its
      center as HyPlan projects it (a tilt direction of 0 degrees forward and
      of 90 degrees to the left).

    Args:
        sensor (Any): The sensor (for example, `hyplan.create_sensor("HyTES")`).
        name (str | None): The name of the instrument. Defaults to the
            sensor's name.
        along_track_field_of_view (float | None): The along-track field of
            view (degrees) of a scanning sensor, which does not change its
            swept footprints. Defaults to the sensor's instantaneous field of
            view (`ifov`), the divergence of a lidar's beam, or the beam width
            of a radar's antenna (its wavelength over its length, or over
            twice its single-look azimuth resolution).

    Returns:
        PointedInstrument: the instrument
    """
    name = getattr(sensor, "name", "Default") if name is None else name
    if callable(getattr(sensor, "swath_offset_angles", None)):
        port, starboard = (float(angle) for angle in sensor.swath_offset_angles())
        if along_track_field_of_view is None:
            along_track_field_of_view = _get_along_track_field_of_view(sensor)
        return PointedInstrument(
            name=name,
            cross_track_field_of_view=starboard - port,
            along_track_field_of_view=along_track_field_of_view,
            roll_angle=-(port + starboard) / 2 + 0.0,  # not -0.0
            view_geometry="scan",
            is_swept=True,
        )
    if all(hasattr(sensor, key) for key in ("fov_x", "fov_y", "tilt_angle")):
        tilt = np.radians(float(sensor.tilt_angle))
        direction = np.radians(float(getattr(sensor, "tilt_direction", 0.0)))
        # the view center's components to the left, forward, and down, of
        # which the roll and then the pitch give `cos(pitch) sin(roll)`,
        # `sin(pitch)`, and `cos(pitch) cos(roll)`
        left = np.sin(tilt) * np.sin(direction)
        forward = np.sin(tilt) * np.cos(direction)
        down = np.cos(tilt)
        return PointedInstrument(
            name=name,
            cross_track_field_of_view=float(sensor.fov_x),
            along_track_field_of_view=float(sensor.fov_y),
            roll_angle=float(np.degrees(np.arctan2(left, down))),
            pitch_angle=float(np.degrees(np.arcsin(np.clip(forward, -1, 1)))),
            is_rectangular=True,
        )
    raise TypeError(
        f"sensor {name!r} ({type(sensor).__name__}) is not supported: give a "
        "scanning sensor (with swath_offset_angles) or a frame camera (with "
        "fov_x, fov_y, and tilt_angle)"
    )


def _get_along_track_field_of_view(sensor: Any) -> float:
    """
    Gets the along-track field of view (degrees) of a scanning sensor (see
    `instrument_from_sensor`).

    Args:
        sensor (Any): The sensor.

    Returns:
        float: the along-track field of view
    """
    if getattr(sensor, "ifov", None) is not None:
        return float(sensor.ifov)
    if getattr(sensor, "beam_divergence", None) is not None:
        return float(sensor.beam_divergence.m_as("degree"))
    lens = getattr(sensor, "lens", None)
    if getattr(lens, "divergence_mrad", None) is not None:
        return float(np.degrees(lens.divergence_mrad / 1000))
    # a radar's beam width: its wavelength over its antenna's length, which
    # is twice the single-look azimuth resolution of a stripmap radar
    length = getattr(sensor, "antenna_length", None)
    if length is None and getattr(sensor, "azimuth_resolution", None) is not None:
        length = 2 * sensor.azimuth_resolution
    if length is not None and getattr(sensor, "wavelength", None) is not None:
        return float(np.degrees((sensor.wavelength / length).m_as("dimensionless")))
    raise ValueError(
        f"the along-track field of view of sensor {getattr(sensor, 'name', '')!r} "
        "is unknown: give along_track_field_of_view"
    )
