"""
Methods to perform space-based (satellite-to-satellite) coverage analysis:
the periods when one satellite can observe another, as by a sensor or a
radio frequency (RF) link, subject to constraints on their line of sight.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from skyfield.timelib import Time

from ..constants import EARTH_EQUATORIAL_RADIUS, EARTH_POLAR_RADIUS
from ..schemas import Satellite
from ..utils.computation import TimeRequest, _run
from ..utils.earth_orientation import _interpolate_nutation
from ..utils.ellipsoid import (
    _find_tangent_distance,
    _itrs_rotation,
    rectangular_to_geodetic,
)
from ..utils.geometry import hash_geometry
from ..utils.time import _index_time, _to_time_from_offsets
from .check import _check_satellites, _check_time_window
from .sampling import _assemble_periods, _find_crossings


def _compute_min_altitude(from_p: np.ndarray, to_p: np.ndarray) -> np.ndarray:
    """
    Computes the minimum WGS 84 geodetic altitude of each line-of-sight
    segment between two positions: the altitude of its tangent point (see
    `tatc.utils.ellipsoid.compute_tangent_point`) if it lies between them,
    or else of the nearer position.

    Args:
        from_p (numpy.ndarray): The positions at one end (meters, Earth-fixed,
            shape (3, N)).
        to_p (numpy.ndarray): The positions at the other end (meters,
            Earth-fixed, shape (3, N)).

    Returns:
        numpy.ndarray: the minimum altitudes (meters)
    """
    d = to_p - from_p
    s = np.clip(_find_tangent_distance(from_p, d), 0, 1)
    return rectangular_to_geodetic(from_p + s * d)[2]


def _get_space_residual(
    from_p: np.ndarray,
    to_p: np.ndarray,
    min_range: float | None = None,
    max_range: float | None = None,
    min_grazing_altitude: float | None = 0,
) -> np.ndarray:
    """
    Computes a residual function of the line of sight from one satellite to
    another that is not positive when the second is observable from the
    first: the largest violation (meters) of each constraint, which are
    - a slant range of at least `min_range` and at most `max_range`, and
    - a minimum WGS 84 geodetic altitude of the line-of-sight segment (see
      `_compute_min_altitude`) of at least `min_grazing_altitude`, so that
      it is not occluded by the Earth (or by its atmosphere below that
      altitude).

    The minimum altitude is bounded by the closest approach of the segment
    to the Earth's center, less the equatorial or polar radius: it is
    computed only where these bounds do not decide the constraint, and the
    bound that decides it is used otherwise (with the correct sign, but
    within about 21 km of the minimum altitude).

    Args:
        from_p (numpy.ndarray): The observing satellite positions (meters,
            Earth-fixed, shape (3, N)).
        to_p (numpy.ndarray): The observed satellite positions (meters,
            Earth-fixed, shape (3, N)).
        min_range (float | None): The minimum slant range (meters), if any.
        max_range (float | None): The maximum slant range (meters), if any.
        min_grazing_altitude (float | None): The minimum altitude (meters)
            of the line of sight above the WGS 84 ellipsoid, or None to
            ignore Earth occlusion.

    Returns:
        numpy.ndarray: the residuals (meters, or -1 without constraints)
    """
    terms = []
    d = to_p - from_p
    if min_range is not None or max_range is not None:
        slant_range = np.linalg.norm(d, axis=0)
        if min_range is not None:
            terms.append(min_range - slant_range)
        if max_range is not None:
            terms.append(slant_range - max_range)
    if min_grazing_altitude is not None:
        # closest approach of the segment to the Earth's center
        s = np.clip(
            -np.einsum("ij,ij->j", from_p, d) / np.einsum("ij,ij->j", d, d), 0, 1
        )
        distance = np.linalg.norm(from_p + s * d, axis=0)
        # bounds of the minimum altitude, which decide the constraint unless
        # the lower bound is below it and the upper bound above it
        lower = distance - EARTH_EQUATORIAL_RADIUS
        upper = distance - EARTH_POLAR_RADIUS
        altitude = np.where(lower > min_grazing_altitude, lower, upper)
        undecided = (lower <= min_grazing_altitude) & (upper >= min_grazing_altitude)
        if np.any(undecided):
            altitude[undecided] = _compute_min_altitude(
                from_p[:, undecided], to_p[:, undecided]
            )
        terms.append(min_grazing_altitude - altitude)
    if len(terms) == 0:
        return np.full(from_p.shape[1], -1.0)
    return np.max(terms, axis=0)


def _get_positions(satellite: Satellite, t: Time, rotation: np.ndarray) -> np.ndarray:
    """
    Gets the Earth-fixed positions (meters, shape (3, N)) of a satellite at
    Skyfield times `t`, given their GCRS -> ITRS rotation matrices (see
    `tatc.utils.ellipsoid._itrs_rotation`).
    """
    track = satellite.orbit.to_gp_orbit().get_orbit_track_at_time(t)
    return np.einsum("ij...,j...->i...", rotation, track.position.m)


def _get_pair_positions(
    satellites: list[Satellite],
    from_index: np.ndarray,
    to_index: np.ndarray,
    t: Time,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the Earth-fixed positions (meters, shape (3, N)) of the observing
    and observed satellites (indices `from_index` and `to_index` of
    `satellites`) at each of Skyfield times `t`, propagating each satellite
    once for all of its times; they share the Earth orientation quantities
    cached on `t` (see `tatc.utils.time._index_time`).
    """
    size = len(t.tt)
    # interpolate the costly nutation angles (see _interpolate_nutation)
    _interpolate_nutation(t)
    rotation = _itrs_rotation(t)
    from_p = np.empty((3, size))
    to_p = np.empty((3, size))
    for u in np.unique(np.concatenate((from_index, to_index))):
        index = np.flatnonzero((from_index == u) | (to_index == u))
        position = _get_positions(
            satellites[u],
            t if len(index) == size else _index_time(t, index),
            rotation[:, :, index],
        )
        is_from = from_index[index] == u
        is_to = to_index[index] == u
        from_p[:, index[is_from]] = position[:, is_from]
        to_p[:, index[is_to]] = position[:, is_to]
    return from_p, to_p


def _get_empty_space_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for space-based coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "target_hash": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        "from_satellite": pd.Series([], dtype="str"),
        "to_satellite": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_space_observations(
    from_satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
    to_satellites: Satellite | list[Satellite] | None = None,
    min_range: float | None = None,
    max_range: float | None = None,
    min_grazing_altitude: float | None = 0,
    min_duration: timedelta = timedelta(seconds=60),
) -> gpd.GeoDataFrame:
    """
    Collects the periods when satellites can observe other satellites, as by
    a sensor or a radio frequency (RF) link, subject to constraints on the
    line of sight between them: its slant range and its occlusion by the
    Earth (the minimum WGS 84 geodetic altitude of the line-of-sight
    segment). Satellite instruments are not used.

    Each ordered pair of an observing ("from") and an observed ("to")
    satellite is analyzed, except a satellite with itself, so that by
    default (with the same satellites observing and observed) each period
    is reported for both directions.

    Each observation records as its geometry the lines of sight from the
    observing to the observed satellite at the start and end of the period
    (a `MultiLineString` of the satellites' longitude, latitude, and
    altitude), identified by the hash of the geometry (`target_hash`, see
    `tatc.utils.hash_geometry`). The lines join only the satellite positions,
    so they do not follow the line of sight across the Earth's surface and
    are not divided at the antimeridian.

    Args:
        from_satellites (Satellite | list[Satellite]): The observing satellite(s).
        start (datetime.datetime): The start of the analysis period.
        end (datetime.datetime): The end of the analysis period.
        to_satellites (Satellite | list[Satellite] | None): The observed
            satellite(s), or None for the observing satellites.
        min_range (float | None): The minimum slant range (meters), if any.
        max_range (float | None): The maximum slant range (meters), if any.
        min_grazing_altitude (float | None): The minimum altitude (meters)
            of the line of sight above the WGS 84 ellipsoid (for example, to
            avoid the atmosphere), or None to ignore Earth occlusion.
        min_duration (datetime.timedelta): The shortest observation period
            guaranteed to be detected (sets the resolution of the search for
            observation periods): a smaller value costs more computation but
            guards against skipping brief periods (or gaps between periods).

    Returns:
        geopandas.GeoDataFrame: The observation periods, sorted by start.
    """
    from_satellites = _check_satellites(from_satellites, "from_satellites")
    to_satellites = (
        from_satellites
        if to_satellites is None
        else _check_satellites(to_satellites, "to_satellites")
    )
    _check_time_window(start, end)
    if min_range is not None and max_range is not None and min_range > max_range:
        raise ValueError(
            f"min_range ({min_range}) is greater than max_range ({max_range})"
        )
    if min_duration <= timedelta(0):
        raise ValueError(f"min_duration ({min_duration}) must be positive")
    # each satellite is propagated once, whether observing, observed, or both
    satellites: list[Satellite] = []
    unique: dict[int, int] = {}
    for satellite in from_satellites + to_satellites:
        if id(satellite) not in unique:
            unique[id(satellite)] = len(satellites)
            satellites.append(satellite)
    pairs = [
        (unique[id(from_satellite)], unique[id(to_satellite)])
        for from_satellite in from_satellites
        for to_satellite in to_satellites
        if from_satellite != to_satellite
    ]
    if len(pairs) == 0 or end == start:
        return _get_empty_space_frame()
    pair_from = np.array([u for u, _ in pairs])
    pair_to = np.array([v for _, v in pairs])

    def residual(from_p: np.ndarray, to_p: np.ndarray) -> np.ndarray:
        return _get_space_residual(
            from_p, to_p, min_range, max_range, min_grazing_altitude
        )

    # residual of each pair at samples (with at least the end points) shared
    # by all pairs, at which each satellite is propagated once
    duration = (end - start).total_seconds()
    samples = np.linspace(
        0, duration, int(duration / (min_duration / 2).total_seconds()) + 2
    )
    t = _to_time_from_offsets(start, samples)
    _interpolate_nutation(t)
    rotation = _itrs_rotation(t)
    positions = [_get_positions(satellite, t, rotation) for satellite in satellites]
    values = []
    for u in np.unique(pair_from):
        # the pairs observed from each satellite together
        observed = pair_to[pair_from == u]
        to_p = np.stack([positions[v] for v in observed], axis=1)
        from_p = np.broadcast_to(positions[u][:, np.newaxis, :], to_p.shape)
        values.extend(
            np.reshape(
                residual(from_p.reshape(3, -1), to_p.reshape(3, -1)),
                (len(observed), len(samples)),
            )
        )
    # brackets of each change of sign between samples, refined together
    brackets = [
        (p, j)
        for p, value in enumerate(values)
        for j in np.flatnonzero((value[:-1] <= 0) != (value[1:] <= 0))
    ]
    crossings = {}
    if len(brackets) > 0:
        bracket_pair = np.array([p for p, _ in brackets])

        def evaluate(seconds: np.ndarray, index: np.ndarray) -> TimeRequest:
            t = _to_time_from_offsets(start, seconds)
            yield t
            return residual(
                *_get_pair_positions(
                    satellites,
                    pair_from[bracket_pair[index]],
                    pair_to[bracket_pair[index]],
                    t,
                )
            )

        crossing, _ = _run(
            _find_crossings(
                evaluate,
                np.array([samples[j] for _, j in brackets]),
                np.array([samples[j + 1] for _, j in brackets]),
            )
        )
        crossings = dict(zip(brackets, crossing))
    periods = [
        (p, left, right)
        for p, parts in enumerate(
            _assemble_periods([samples] * len(pairs), values, crossings)
        )
        for left, right in parts
        if right > left
    ]
    if len(periods) == 0:
        return _get_empty_space_frame()
    # positions of each pair at the start and end of each period, computed together
    period_pair = np.array([p for p, _, _ in periods])
    seconds = np.array([[left, right] for _, left, right in periods]).flatten()
    from_p, to_p = _get_pair_positions(
        satellites,
        np.repeat(pair_from[period_pair], 2),
        np.repeat(pair_to[period_pair], 2),
        _to_time_from_offsets(start, seconds),
    )
    from_geo = np.array(rectangular_to_geodetic(from_p)).T
    to_geo = np.array(rectangular_to_geodetic(to_p)).T
    # lines of sight at the start and end of each period
    geometry = shapely.multilinestrings(
        shapely.linestrings(np.stack((from_geo, to_geo), axis=1)),
        indices=np.repeat(np.arange(len(periods)), 2),
    )
    reference = pd.Timestamp(start).tz_convert("UTC")
    return (
        gpd.GeoDataFrame(
            {
                "target_hash": [hash_geometry(line) for line in geometry],
                "geometry": geometry,
                "from_satellite": [satellites[pair_from[p]].name for p in period_pair],
                "to_satellite": [satellites[pair_to[p]].name for p in period_pair],
                "start": reference + pd.to_timedelta(seconds[0::2], unit="s"),
                "end": reference + pd.to_timedelta(seconds[1::2], unit="s"),
            },
            crs="EPSG:4326",
        )
        .sort_values("start", kind="stable")
        .reset_index(drop=True)
    )
