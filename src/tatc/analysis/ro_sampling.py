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
from skyfield.timelib import Time

from ..constants import timescale
from ..schemas import Satellite
from ..utils.ellipsoid import _ellipsoidal_tangent_point, _itrs_rotation
from ..utils.orbital import compute_vnb_frame
from ..utils.propagation import _index_orbit_track
from .check import _check_satellite, _check_satellites


def _tangent_point_geometry(
    tx_pv: Geocentric,
    rx_pv: Geocentric,
    rx_v_u: np.ndarray,
    rx_n_u: np.ndarray,
    rx_b_u: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes tangent point position and receiver-frame pitch/yaw angles of
    the transmitter, as seen from the receiver, at one or more times.

    The tangent point is the point on the receiver-transmitter line with
    minimum WGS 84 geodetic altitude (see
    `tatc.utils.compute_tangent_point`), the
    convention used to geolocate operationally processed RO profiles. It
    differs from the line's closest approach to the Earth's center by up to
    about 20 km horizontally at middle latitudes (the two coincide at the
    equator and poles), but by only meters in altitude.
    """
    # relative position of transmitter from receiver: x_(rx,tx) = x_tx - x_rx
    rx_tx_pv = tx_pv - rx_pv
    rx_tx_p_m = np.array(rx_tx_pv.position.m)
    rx_p_m = np.array(rx_pv.position.m)
    tx_p_m = np.array(tx_pv.position.m)
    # tangent point position (m)
    tp_p = _ellipsoidal_tangent_point(rx_p_m, rx_tx_p_m, rx_pv.t)
    # intersecting (-1) or parallel (+1) view of tangent point
    tp_sign = np.sign(np.einsum("ij,ij->j", tp_p - tx_p_m, tp_p - rx_p_m))
    # relative transmitter position from receiver in plane normal to receiver orbit
    rx_tx_p_rx_n_plane = rx_tx_pv.position.m - np.einsum(
        "ij,j->ij", rx_n_u, np.einsum("ij,ij->j", rx_n_u, rx_tx_p_m)
    )
    # relative transmitter position from receiver in plane binormal to receiver orbit
    rx_tx_p_rx_t_plane = rx_tx_pv.position.m - np.einsum(
        "ij,j->ij", rx_b_u, np.einsum("ij,ij->j", rx_b_u, rx_tx_p_m)
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
    return tp_p, tp_sign, rx_tx_pitch, rx_tx_yaw


def _get_ro_validity(
    transmitters: list[Satellite],
    receiver: Satellite,
    t: Time,
    parts: list[tuple[int, np.ndarray | None]],
    max_yaw: float,
) -> list[np.ndarray]:
    """
    Evaluates whether each of several transmitters' tangent points intersect
    the Earth with the transmitter yaw within bounds (1 if so, 0 otherwise),
    each at its own part of Skyfield times `t`: (transmitter index, index of
    its times in `t`, or None for all of `t`). The receiver is propagated
    once for all of `t`, and the parts share the Earth orientation
    quantities cached on `t` (see `tatc.utils.propagation._index_time`).
    """
    # receiver and transmitter positions are compared directly in the
    # inertial frame, in which repeat tracks are expressed at the true time
    rx_track = receiver.orbit.to_gp_orbit().get_orbit_track_at_time(t)
    values = []
    for k, index in parts:
        rx_pv = rx_track if index is None else _index_orbit_track(rx_track, index)
        tx_pv = transmitters[k].orbit.to_gp_orbit().get_orbit_track_at_time(rx_pv.t)
        rx_v_u, rx_n_u, rx_b_u = compute_vnb_frame(rx_pv)
        _, tp_sign, _, rx_tx_yaw = _tangent_point_geometry(
            tx_pv, rx_pv, rx_v_u, rx_n_u, rx_b_u
        )
        # valid if tangent point intersects and yaw angle below maximum
        valid = np.logical_and(
            tp_sign < 0, np.abs(rx_tx_yaw) % (180 - max_yaw) < max_yaw
        )
        values.append(valid.astype(int))
    return values


def _find_ro_arcs(
    transmitters: list[Satellite],
    receiver: Satellite,
    start: datetime,
    end: datetime,
    max_yaw: float,
    step_days: float,
    epsilon: float = 1e-3 / 86400,
    num: int = 12,
) -> list[tuple[int, datetime, datetime]]:
    """
    Finds the arcs (periods) when each transmitter's RO observation is valid
    (see `_get_ro_validity`) as by Skyfield's `find_discrete` for each
    transmitter (sampled at most `step_days` apart, and refined by dividing
    the brackets of each change of validity into `num` samples until they
    are at most `epsilon` days long), but for all transmitters together:
    they share the initial samples, and the refined samples of all
    transmitters are evaluated together at each step.

    Returns:
        list[tuple[int, datetime, datetime]]: each arc's transmitter index,
            start, and end, by transmitter and in time order.
    """
    if len(transmitters) == 0:
        return []
    t_start, t_end = timescale.from_datetime(start), timescale.from_datetime(end)
    jd0, jd1 = t_start.tt, t_end.tt
    if jd0 >= jd1:
        raise ValueError(
            f"your start_time {t_start} is later than your end_time {t_end}"
        )
    everything = [(k, None) for k in range(len(transmitters))]
    # validity at the start of the first segment, from the start time itself
    initial = _get_ro_validity(
        transmitters, receiver, timescale.from_datetimes([start]), everything, max_yaw
    )
    # initial samples (with at least the end points), shared by all transmitters
    jd = np.linspace(jd0, jd1, int((jd1 - jd0) / step_days) + 2)
    samples = dict.fromkeys(range(len(transmitters)), jd)
    values = dict(
        zip(
            samples,
            _get_ro_validity(
                transmitters, receiver, timescale.tt_jd(jd), everything, max_yaw
            ),
        )
    )
    end_mask = np.linspace(0.0, 1.0, num)
    start_mask = end_mask[::-1]
    transitions = {}
    while len(samples) > 0:
        for k in list(samples):
            jd, y = samples[k], values[k]
            indices = np.flatnonzero(np.diff(y))
            if len(indices) == 0:
                # no change of validity
                transitions[k] = (jd[indices], y[indices])
                del samples[k]
                continue
            starts, ends = jd[indices], jd[indices + 1]
            # brackets narrow at the same rate (from equal initial intervals),
            # so only the first is tested
            if ends[0] - starts[0] <= epsilon:
                # keep only the last of changes less than epsilon apart
                mask = np.concatenate((np.diff(ends) > 3.0 * epsilon, [True]))
                transitions[k] = (ends[mask], y[indices + 1][mask])
                del samples[k]
                continue
            samples[k] = (
                np.multiply.outer(starts, start_mask).flatten()
                + np.multiply.outer(ends, end_mask).flatten()
            )
        if len(samples) > 0:
            # evaluate the refined samples of all transmitters together
            bounds = np.cumsum([0] + [len(jd) for jd in samples.values()])
            values = dict(
                zip(
                    samples,
                    _get_ro_validity(
                        transmitters,
                        receiver,
                        timescale.tt_jd(np.concatenate(list(samples.values()))),
                        [
                            (k, np.arange(bounds[i], bounds[i + 1]))
                            for i, k in enumerate(samples)
                        ],
                        max_yaw,
                    ),
                )
            )
    # times of all changes of validity, converted together
    order = sorted(transitions)
    counts = np.cumsum([0] + [len(transitions[k][0]) for k in order])
    utc = np.atleast_1d(
        timescale.tt_jd(
            np.concatenate([transitions[k][0] for k in order] + [np.array([])])
        ).utc_datetime()
    )
    arcs = []
    for i, k in enumerate(order):
        # boundary times delimiting N+1 alternating valid/invalid segments (N
        # = number of transitions); each segment's validity is the value that
        # becomes active at its start (the initial validity for the first
        # segment, else the corresponding transition value) -- the final
        # boundary time (`end`) is a pure endpoint with no segment-start
        # value of its own
        boundary_times = [start] + list(utc[counts[i] : counts[i + 1]]) + [end]
        boundary_values = [bool(initial[k][0])] + [
            bool(value) for value in transitions[k][1]
        ]
        # keep only the segments where validity holds
        arcs.extend(
            (k, boundary_times[j], boundary_times[j + 1])
            for j in range(len(boundary_times) - 1)
            if boundary_values[j] and boundary_times[j + 1] > boundary_times[j]
        )
    return arcs


def _tangent_point_tx_azimuth(
    tp_p: np.ndarray,
    tx_p: np.ndarray,
    t: Time,
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> np.ndarray:
    """
    Computes the transmitter azimuth (deg, clockwise from North) as viewed from
    each point of a tangent point track, in the local horizontal (east, north)
    frame of the geodetic tangent point (as Skyfield's
    `(satellite - geographic_position).at(t).altaz()` would).

    Works directly from the tangent point and transmitter positions (m, shape
    (3, N), GCRS) already computed for the track, rotated to the Earth-fixed
    frame, rather than re-propagating the transmitter: the rotation reuses the
    Earth orientation quantities cached on `t` by that propagation, which
    Skyfield would otherwise recompute for every new (sliced) time object.
    """
    tp_tx = np.einsum("ij...,j...->i...", _itrs_rotation(t), tx_p - tp_p)
    lat, lon = np.radians(latitude), np.radians(longitude)
    east = -np.sin(lon) * tp_tx[0] + np.cos(lon) * tp_tx[1]
    north = (
        -np.sin(lat) * np.cos(lon) * tp_tx[0]
        - np.sin(lat) * np.sin(lon) * tp_tx[1]
        + np.cos(lat) * tp_tx[2]
    )
    return np.degrees(np.arctan2(east, north)) % 360


def _sample_ro_arcs(
    transmitters: list[Satellite],
    receiver: Satellite,
    arcs: list[tuple[int, datetime, datetime]],
    time_step: timedelta,
    range_elevation: tuple[float, float],
) -> list[dict]:
    """
    Samples tangent point observations across valid RO arcs (see
    `_find_ro_arcs`: each with its transmitter index, start, and end),
    splitting each into one or more observations if the tangent point
    elevation leaves the specified range. The samples of all arcs are
    computed together.
    """
    if len(arcs) == 0:
        return []
    # sample each arc at (at most) the specified time step, including both endpoints
    arc_times = []
    for _, arc_start, arc_end in arcs:
        steps = max(int(np.ceil((arc_end - arc_start) / time_step)), 1)
        arc_times.append(
            [arc_start + i * (arc_end - arc_start) / steps for i in range(steps + 1)]
        )
    bounds = np.cumsum([0] + [len(times) for times in arc_times])
    t = timescale.from_datetimes(list(chain.from_iterable(arc_times)))

    rx_track = receiver.orbit.to_gp_orbit().get_orbit_track_at_time(t)
    tp_p = np.empty((3, bounds[-1]))
    tx_p = np.empty((3, bounds[-1]))
    rx_tx_pitch = np.empty(bounds[-1])
    rx_tx_yaw = np.empty(bounds[-1])
    for k in sorted(set(k for k, _, _ in arcs)):
        # the samples of the arcs of each transmitter
        index = np.concatenate(
            [
                np.arange(bounds[i], bounds[i + 1])
                for i, (arc_k, _, _) in enumerate(arcs)
                if arc_k == k
            ]
        )
        rx_pv = _index_orbit_track(rx_track, index)
        rx_v_u, rx_n_u, rx_b_u = compute_vnb_frame(rx_pv)
        tx_pv = transmitters[k].orbit.to_gp_orbit().get_orbit_track_at_time(rx_pv.t)
        tp_p[:, index], _, rx_tx_pitch[index], rx_tx_yaw[index] = (
            _tangent_point_geometry(tx_pv, rx_pv, rx_v_u, rx_n_u, rx_b_u)
        )
        tx_p[:, index] = tx_pv.position.m

    # tangent point geodetic position, computed once for all arcs
    tpp_geo = wgs84.geographic_position_of(Geocentric(Distance(m=tp_p).au, None, t))
    longitude = np.array(tpp_geo.longitude.degrees)
    latitude = np.array(tpp_geo.latitude.degrees)
    elevation = np.array(tpp_geo.elevation.m)
    # azimuth of transmitter from geodetic tangent point (clockwise from North)
    tp_tx_azimuth = _tangent_point_tx_azimuth(tp_p, tx_p, t, latitude, longitude)
    # tangent point height within elevation range
    in_range = np.logical_and(
        elevation > range_elevation[0], elevation < range_elevation[1]
    )

    # occultation observations
    occ_obs = []
    for (k, _, _), times, offset in zip(arcs, arc_times, bounds):
        # occultation arc
        occ_arc = None
        for j, time in enumerate(times, start=offset):
            if in_range[j]:
                if occ_arc is None:
                    # start of new RO observation. rx_tx_pitch is the
                    # transmitter's pitch angle relative to the receiver,
                    # where -90 deg points at the geocenter (never actually
                    # reached by a real RO profile, since the signal must
                    # pass through the atmosphere); pitch above -90 deg means
                    # the transmitter is "ahead" of the receiver (a
                    # rising/emersion occultation), below -90 deg means
                    # "behind" (a setting/immersion one). This is an
                    # approximation of the more direct (but more expensive)
                    # definition -- the sign of the tangent point's own
                    # elevation rate -- using this arc's first sample only.
                    occ_arc = {
                        "tx": transmitters[k].name,
                        "is_rising": rx_tx_pitch[j] > -90,
                        "points": [],
                    }
                occ_arc["points"].append(
                    {
                        "time": time,
                        "longitude": longitude[j],
                        "latitude": latitude[j],
                        "elevation": elevation[j],
                        "rx_tx_pitch": rx_tx_pitch[j],
                        "rx_tx_yaw": rx_tx_yaw[j],
                        "tp_tx_azimuth": tp_tx_azimuth[j],
                    }
                )
                if j + 1 >= offset + len(times):
                    # end of RO observation due to arc boundary
                    occ_obs.append(occ_arc)
                    occ_arc = None
            elif occ_arc is not None:
                # end of RO observation due to elevation constraints
                occ_obs.append(occ_arc)
                occ_arc = None
    return occ_obs


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
    range_elevation: tuple[float, float] = (-200e3, 60e3),
    min_profile_duration: timedelta = timedelta(seconds=30),
) -> gpd.GeoDataFrame:
    """
    Collects Radio Occultation (RO) observations.

    The tangent point (the geometric basis for every reported position and
    elevation) is the point on the straight receiver-transmitter line with
    minimum WGS 84 geodetic altitude, the convention used to geolocate
    operationally processed RO profiles. Atmospheric refraction is not
    modeled, so tangent point elevations are those of the straight line,
    well below those of the refracted signal in the lower atmosphere (hence
    the negative default `sample_elevation`). See `_tangent_point_geometry`
    for details.

    Args:
        receiver (Satellite): the satellite with a RO receiver.
        transmitters (Satellite | list[Satellite]]): the satellite(s) with a RO transmitter.
        start (datetime.datetime): the start of the analysis period.
        end (datetime.datetime): the end of the analysis period.
        time_step (datetime.timedelta): the time step used to sample tangent point
            tracks within each observation period, once its bounds are found.
        sample_elevation: (float): the elevation (m) at which to interpolate observation attributes.
        max_yaw (float): the maximum transmitter yaw angle (from receiver body-fixed frame)
            for a valid obsevation.
        range_elevation: (tuple[float, float]): the lower and upper bound on tangent
            point elevation (m) for a valid observation.
        min_profile_duration (datetime.timedelta): the shortest RO observation period
            guaranteed to be detected. Sets the coarse scan resolution used to search
            for observation periods (via `skyfield.searchlib.find_discrete`),
            independent of `time_step` and of the overall analysis duration. Set this
            no larger than the shortest profile you expect; a smaller value costs more
            computation but guards against silently skipping brief observation periods.
    """
    _check_satellite(receiver, "receiver")
    transmitters = _check_satellites(transmitters, "transmitters")
    # find the valid observation periods of all transmitters together,
    # scanning at half the shortest profile duration we must not skip
    # (decoupled from time_step so long mission durations don't blow up the
    # coarse scan), and sample them together
    obs = _sample_ro_arcs(
        transmitters,
        receiver,
        _find_ro_arcs(
            transmitters,
            receiver,
            start,
            end,
            max_yaw,
            (min_profile_duration / 2) / timedelta(days=1),
        ),
        time_step,
        range_elevation,
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
