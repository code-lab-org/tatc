"""
Methods to perform radio occultation (RO) coverage analysis.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""
from __future__ import annotations

from datetime import datetime, timedelta
from itertools import chain

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import MultiPoint, Point
from skyfield.api import Distance, wgs84
from skyfield.positionlib import Geocentric
from skyfield.searchlib import find_discrete

from ..constants import timescale
from ..schemas import Satellite


def _tangent_point_geometry(
    tx_pv: Geocentric,
    rx_pv: Geocentric,
    rx_v_u: list[float],
    rx_n_u: list[float],
    rx_b_u: list[float],
    compute_velocity: bool = False,
):
    """
    Computes tangent point position (and, optionally, velocity) and
    receiver-frame pitch/yaw angles of the transmitter, as seen from the
    receiver, at one or more times.

    Tangent point velocity is not needed for geodetic position or azimuth
    computations (Skyfield ignores it there), so it is skipped by default;
    pass `compute_velocity=True` to compute it anyway.
    """
    # relative position, velocity of transmitter from receiver
    # x_(rx,tx) = x_tx - x_rx; v_(rx,tx) = v_tx - v_rx
    rx_tx_pv = tx_pv - rx_pv
    # tangent point position (m)
    # x_tp = x_tx - x_(rx,tx) . [ x_tx . x_(rx,tx) ] / || x_(rx,tx) ||
    tp_p = tx_pv.position.m - np.einsum(
        "ij,j->ij",
        rx_tx_pv.position.m,
        np.divide(
            np.einsum("ij,ij->j", tx_pv.position.m, rx_tx_pv.position.m),
            np.einsum("ij,ij->j", rx_tx_pv.position.m, rx_tx_pv.position.m),
        ),
    )
    if compute_velocity:
        # tangent point velocity (m/s) - derived using chain rule
        # v_tp = v_tx - v_(rx,tx) . [ x_tx . x_(rx,tx) ] / || x_(rx,tx) ||
        #        - x_(rx,tx) . [
        #           [ v_tx . x_(rx,tx) ] + [ x_tx . v_(rx,tx) ] ] / || x_(rx,tx) || ]
        #           - 2 * [ v_(rx,tx) . x_(rx,tx) ] * [ x_tx . x_(rx,tx) ] / || x_(rx,tx) ||^2
        #        ]
        tp_v = (
            tx_pv.velocity.m_per_s
            - np.einsum(
                "ij,j->ij",
                rx_tx_pv.velocity.m_per_s,
                np.divide(
                    np.einsum("ij,ij->j", tx_pv.position.m, rx_tx_pv.position.m),
                    np.einsum("ij,ij->j", rx_tx_pv.position.m, rx_tx_pv.position.m),
                ),
            )
            - np.einsum(
                "ij,j->ij",
                rx_tx_pv.position.m,
                (
                    np.divide(
                        (
                            np.einsum(
                                "ij,ij->j", tx_pv.velocity.m_per_s, rx_tx_pv.position.m
                            )
                            + np.einsum(
                                "ij,ij->j", tx_pv.position.m, rx_tx_pv.velocity.m_per_s
                            )
                        ),
                        np.einsum(
                            "ij,ij->j", rx_tx_pv.position.m, rx_tx_pv.position.m
                        ),
                    )
                    - 2
                    * np.divide(
                        np.multiply(
                            np.einsum(
                                "ij,ij->j",
                                rx_tx_pv.velocity.m_per_s,
                                rx_tx_pv.position.m,
                            ),
                            np.einsum(
                                "ij,ij->j", tx_pv.position.m, rx_tx_pv.position.m
                            ),
                        ),
                        np.power(
                            np.einsum(
                                "ij,ij->j", rx_tx_pv.position.m, rx_tx_pv.position.m
                            ),
                            2,
                        ),
                    )
                ),
            )
        )
    else:
        tp_v = None
    # intersecting (-1) or parallel (+1) view of tangent point
    tp_sign = np.sign(
        np.einsum("ij,ij->j", tp_p - tx_pv.position.m, tp_p - rx_pv.position.m)
    )
    # relative transmitter position from receiver in plane normal to receiver orbit
    rx_tx_p_rx_n_plane = rx_tx_pv.position.m - np.einsum(
        "ij,j->ij", rx_n_u, np.einsum("ij,ij->j", rx_n_u, rx_tx_pv.position.m)
    )
    # relative transmitter position from receiver in plane binormal to receiver orbit
    rx_tx_p_rx_t_plane = rx_tx_pv.position.m - np.einsum(
        "ij,j->ij", rx_b_u, np.einsum("ij,ij->j", rx_b_u, rx_tx_pv.position.m)
    )
    # transmitter pitch angle in receiver body-fixed frame
    rx_tx_pitch = np.degrees(
        np.arctan2(
            np.einsum("ij,ij->j", rx_tx_p_rx_n_plane, rx_b_u),
            np.einsum("ij,ij->j", rx_tx_p_rx_n_plane, rx_v_u),
        )
    )
    # transmitter yaw angle in receiver body-fixed frame
    rx_tx_yaw = np.degrees(
        np.arctan2(
            np.einsum("ij,ij->j", rx_tx_p_rx_t_plane, rx_n_u),
            np.einsum("ij,ij->j", rx_tx_p_rx_t_plane, rx_v_u),
        )
    )
    return tp_p, tp_v, tp_sign, rx_tx_pitch, rx_tx_yaw


def _receiver_frame_vectors(rx_pv: Geocentric):
    """
    Computes the receiver body-fixed (VNB) frame unit vectors.
    """
    # unit vector tangent to receiver orbit plane (VNB x-axis)
    rx_v_u = np.divide(
        rx_pv.velocity.m_per_s, np.linalg.norm(rx_pv.velocity.m_per_s, axis=0)
    )
    # unit vector normal to receiver orbit plane (VNB y-axis)
    rx_n_u = np.cross(rx_pv.position.m, rx_pv.velocity.m_per_s, 0, 0, -1).T
    rx_n_u = np.divide(rx_n_u, np.linalg.norm(rx_n_u, axis=0))
    # unit vector orthogonal to receiver orbit plane (VNB z-axis)
    rx_b_u = np.divide(rx_pv.position.m, np.linalg.norm(rx_pv.position.m, axis=0))
    return rx_v_u, rx_n_u, rx_b_u


def _make_ro_validity_function(
    transmitter: Satellite, receiver: Satellite, max_yaw: float, step_days: float
):
    """
    Builds a Skyfield-compatible discrete function of time returning whether the
    tangent point intersects the Earth and the transmitter yaw is within bounds,
    for use with `skyfield.searchlib.find_discrete`.
    """

    def f(t):
        # get_orbit_track_at_time always does true (directly propagated) inertial
        # propagation, which is required here since receiver and transmitter
        # positions are compared directly in the inertial frame.
        rx_pv = receiver.orbit.to_gp_orbit().get_orbit_track_at_time(t)
        tx_pv = transmitter.orbit.to_gp_orbit().get_orbit_track_at_time(t)
        rx_v_u, rx_n_u, rx_b_u = _receiver_frame_vectors(rx_pv)
        _, _, tp_sign, _, rx_tx_yaw = _tangent_point_geometry(
            tx_pv, rx_pv, rx_v_u, rx_n_u, rx_b_u
        )
        # valid if tangent point intersects and yaw angle below maximum
        valid = np.logical_and(
            tp_sign < 0, np.abs(rx_tx_yaw) % (180 - max_yaw) < max_yaw
        )
        return valid.astype(int)

    f.step_days = step_days
    return f


def _tangent_point_tx_azimuth(
    transmitter: Satellite,
    times: list[datetime],
    t,
    tp_p: np.ndarray,
) -> np.ndarray:
    """
    Computes the transmitter azimuth (deg, clockwise from North) as viewed from
    each point of a tangent point track, vectorized per distinct TLE element used
    across the track (almost always a single element, given how short RO arcs are).

    Only the tangent point position (not velocity) is needed: Skyfield's
    geodetic and azimuth computations do not use it.
    """
    orbit = transmitter.orbit.to_gp_orbit()
    element_indices = np.asarray(orbit.get_closest_element_index(times))
    azimuth = np.empty(len(times))
    for element_index in np.unique(element_indices):
        mask = element_indices == element_index
        sat = orbit.elements[element_index].to_skyfield()
        tpp_geo = wgs84.geographic_position_of(
            Geocentric(Distance(m=tp_p[:, mask]).au, None, t[mask])
        )
        azimuth[mask] = (sat - tpp_geo).at(t[mask]).altaz()[1].degrees
    return azimuth


def _sample_ro_arc(
    transmitter: Satellite,
    receiver: Satellite,
    arc_start: datetime,
    arc_end: datetime,
    time_step: timedelta,
    range_elevation: tuple[float],
) -> list[dict]:
    """
    Samples tangent point observations across a single valid RO arc, splitting it
    into one or more observations if the tangent point elevation leaves the
    specified range.
    """
    # sample the arc at (at most) the specified time step, including both endpoints
    steps = max(int(np.ceil((arc_end - arc_start) / time_step)), 1)
    times = [arc_start + i * (arc_end - arc_start) / steps for i in range(steps + 1)]
    t = timescale.from_datetimes(times)

    rx_pv = receiver.orbit.to_gp_orbit().get_orbit_track_at_time(t)
    rx_v_u, rx_n_u, rx_b_u = _receiver_frame_vectors(rx_pv)
    tx_pv = transmitter.orbit.to_gp_orbit().get_orbit_track_at_time(t)
    tp_p, _, _, rx_tx_pitch, rx_tx_yaw = _tangent_point_geometry(
        tx_pv, rx_pv, rx_v_u, rx_n_u, rx_b_u
    )

    # tangent point geodetic position, computed once for the whole arc
    tpp_geo = wgs84.geographic_position_of(Geocentric(Distance(m=tp_p).au, None, t))
    longitude = tpp_geo.longitude.degrees
    latitude = tpp_geo.latitude.degrees
    elevation = tpp_geo.elevation.m
    # azimuth of transmitter from geodetic tangent point (clockwise from North)
    tp_tx_azimuth = _tangent_point_tx_azimuth(transmitter, times, t, tp_p)
    # tangent point height within elevation range
    in_range = np.logical_and(
        elevation > range_elevation[0], elevation < range_elevation[1]
    )

    # occultation observations
    occ_obs = []
    # occultation arc
    occ_arc = None
    for j in range(len(times)):
        if in_range[j]:
            if occ_arc is None:
                # start of new RO observation
                occ_arc = {
                    "tx": transmitter.name,
                    "is_rising": rx_tx_pitch[j] > -90,
                    "points": [],
                }
            occ_arc["points"].append(
                {
                    "time": times[j],
                    "longitude": longitude[j],
                    "latitude": latitude[j],
                    "elevation": elevation[j],
                    "rx_tx_pitch": rx_tx_pitch[j],
                    "rx_tx_yaw": rx_tx_yaw[j],
                    "tp_tx_azimuth": tp_tx_azimuth[j],
                }
            )
            if j + 1 >= len(times):
                # end of RO observation due to arc boundary
                occ_obs.append(occ_arc)
                occ_arc = None
        elif occ_arc is not None:
            # end of RO observation due to elevation constraints
            occ_obs.append(occ_arc)
            occ_arc = None
    return occ_obs


def _collect_ro_series(
    transmitter: Satellite,
    receiver: Satellite,
    start: datetime,
    end: datetime,
    time_step: timedelta,
    max_yaw: float,
    range_elevation: tuple[float],
    min_profile_duration: timedelta,
) -> list[dict]:
    # discrete function of time: 1 if the tangent point intersects and the
    # transmitter yaw angle is below maximum, 0 otherwise
    # scan at half the shortest profile duration we must not skip, decoupled
    # from time_step so long mission durations don't blow up the coarse scan
    is_valid = _make_ro_validity_function(
        transmitter, receiver, max_yaw, (min_profile_duration / 2) / timedelta(days=1)
    )
    # find the precise times at which validity changes
    transition_times, transition_values = find_discrete(
        timescale.from_datetime(start), timescale.from_datetime(end), is_valid
    )
    initial_valid = bool(is_valid(timescale.from_datetimes([start]))[0])
    final_valid = bool(transition_values[-1]) if len(transition_values) else initial_valid
    # boundary times/values delimiting alternating valid/invalid segments
    boundary_times = [start] + list(transition_times.utc_datetime()) + [end]
    boundary_values = (
        [initial_valid] + [bool(value) for value in transition_values] + [final_valid]
    )
    # keep only the segments where validity holds
    arcs = [
        (boundary_times[i], boundary_times[i + 1])
        for i in range(len(boundary_times) - 1)
        if boundary_values[i] and boundary_times[i + 1] > boundary_times[i]
    ]
    return list(
        chain.from_iterable(
            _sample_ro_arc(
                transmitter, receiver, arc_start, arc_end, time_step, range_elevation
            )
            for arc_start, arc_end in arcs
        )
    )


def _interpolate_ro_point(points: list[dict], sample_elevation: float) -> dict:
    """
    Interpolates RO observation attributes at the specified tangent point elevation.

    Args:
        points (list[dict]): the RO observation points (ordered by time).
        sample_elevation (float): the tangent point elevation (m) at which to interpolate.

    Returns:
        dict: interpolated longitude (deg), latitude (deg), elevation (m), rx_tx_pitch (deg),
            rx_tx_yaw (deg), tp_tx_azimuth (deg), and time.
    """
    elevations = np.array([point["elevation"] for point in points])
    diffs = elevations - sample_elevation
    # bracketing indices where the tangent point elevation crosses the sample elevation
    crossings = np.nonzero(np.diff(np.sign(diffs)))[0]
    if len(crossings) > 0:
        i = crossings[0]
        p0, p1 = points[i], points[i + 1]
        denom = diffs[i] - diffs[i + 1]
        frac = diffs[i] / denom if denom != 0 else 0.0
    else:
        # sample elevation is outside the observed range: clamp to the nearest endpoint
        i = 0 if abs(diffs[0]) <= abs(diffs[-1]) else len(points) - 1
        p0 = p1 = points[i]
        frac = 0.0

    def lerp(a, b):
        return a + frac * (b - a)

    def lerp_angle(a, b, low=-180.0):
        # interpolate along the shortest angular path, then wrap to [low, low + 360)
        diff = ((b - a + 180) % 360) - 180
        return (a + frac * diff - low) % 360 + low

    return {
        "longitude": lerp_angle(p0["longitude"], p1["longitude"]),
        "latitude": lerp(p0["latitude"], p1["latitude"]),
        "elevation": lerp(p0["elevation"], p1["elevation"]),
        "rx_tx_pitch": lerp_angle(p0["rx_tx_pitch"], p1["rx_tx_pitch"]),
        "rx_tx_yaw": lerp_angle(p0["rx_tx_yaw"], p1["rx_tx_yaw"]),
        "tp_tx_azimuth": lerp_angle(p0["tp_tx_azimuth"], p1["tp_tx_azimuth"], low=0.0),
        "time": p0["time"] + frac * (p1["time"] - p0["time"]),
    }


def _get_empty_ro_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for ro results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "receiver": pd.Series([], dtype="str"),
        "transmitter": pd.Series([], dtype="str"),
        "is_rising": pd.Series([], dtype="float"),
        "geometry": pd.Series([], dtype="object"),
        "position": pd.Series([], dtype="object"),
        "rx_tx_pitch": pd.Series([], dtype="float"),
        "rx_tx_yaw": pd.Series([], dtype="float"),
        "tp_tx_azimuth": pd.Series([], dtype="float"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "time": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_ro_observations(
    receiver: Satellite,
    transmitters: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
    time_step: timedelta = timedelta(seconds=10),
    sample_elevation: float = -80e3,
    max_yaw: float = 65,
    range_elevation: tuple[float] = (-200e3, 60e3),
    min_profile_duration: timedelta = timedelta(seconds=30),
) -> gpd.GeoDataFrame:
    """
    Collects Radio Occultation (RO) observations.

    Args:
        receiver (Satellite): the satellite with a RO receiver.
        transmitters (Satellite | list[Satellite]]): the satellite(s) with a RO transmitter.
        start (datetime.datetime): the start of the analysis period.
        end (datetime.datetime): the end of the analysis period.
        time_step (datetime.timedelta): the time step used to sample tangent point
            tracks within each observation period, once its bounds are found.
        sample_elevation: (float): the elevation (m) at which to interpolate observation attributes.
        max_yaw (float): the maximum transmitter yaw angle (from receiver body-fixed frame) for a valid obsevation.
        range_elevation: (tuple[float]): the lower and upper bound on tangent point elevation (m) for a valid observation.
        min_profile_duration (datetime.timedelta): the shortest RO observation period
            guaranteed to be detected. Sets the coarse scan resolution used to search
            for observation periods (via `skyfield.searchlib.find_discrete`),
            independent of `time_step` and of the overall analysis duration. Set this
            no larger than the shortest profile you expect; a smaller value costs more
            computation but guards against silently skipping brief observation periods.
    """
    # generate observations
    obs = list(
        chain.from_iterable(
            _collect_ro_series(
                transmitter,
                receiver,
                start,
                end,
                time_step,
                max_yaw,
                range_elevation,
                min_profile_duration,
            )
            for transmitter in (
                transmitters if isinstance(transmitters, list) else [transmitters]
            )
        )
    )
    if len(obs) == 0:
        return _get_empty_ro_frame()
    # format results
    return gpd.GeoDataFrame(
        [
            {
                "receiver": receiver.name,
                "transmitter": o["tx"],
                "is_rising": o["is_rising"],
                "geometry": MultiPoint(
                    [
                        [point["longitude"], point["latitude"], point["elevation"]]
                        for point in o["points"]
                    ]
                ),
                "position": Point(
                    sample["longitude"], sample["latitude"], sample["elevation"]
                ),
                "rx_tx_pitch": sample["rx_tx_pitch"],
                "rx_tx_yaw": sample["rx_tx_yaw"],
                "tp_tx_azimuth": sample["tp_tx_azimuth"],
                "start": o["points"][0]["time"],
                "end": o["points"][-1]["time"],
                "time": sample["time"],
            }
            for o in obs
            for sample in [_interpolate_ro_point(o["points"], sample_elevation)]
        ],
        crs="EPSG:4326",
    ).sort_values("time", ignore_index=True)
