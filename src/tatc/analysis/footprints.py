"""
Methods to collect, dissolve, and sample the footprints of any observing
platforms together: satellites, aircraft, fixed surface platforms, and
radar stations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import replace
from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely import geometry as geo

from .._warnings import warn
from ..schemas import AllPlatforms, Point
from ..utils.geometry import _get_point_coordinates, _hash_geometries
from ._tracks import Mask, _collect_tracks, _FootprintRequest
from .airborne_track import _collect_aircraft
from .check import (
    _PLATFORM_KINDS,
    AIRCRAFT,
    RADAR_STATION,
    SATELLITE,
    SURFACE_PLATFORM,
    _dissolve_footprints,
    _get_empty_observation_frame,
    _get_instruments,
    _get_platform_columns,
    _get_platforms,
    _group_platforms,
    _is_single,
    _PlatformKind,
)
from .ground_track import _collect_satellite
from .radar_track import _collect_radar_station, _get_empty_radar_track
from .surface_track import _collect_surface_platform

_COLLECTORS = {
    SATELLITE: (_collect_satellite,),
    AIRCRAFT: (_collect_aircraft,),
    SURFACE_PLATFORM: (_collect_surface_platform,),
    RADAR_STATION: (_collect_radar_station, _get_empty_radar_track),
}
"""
The collector of the footprints of a platform of each kind (see
`tatc.analysis._tracks._collect_tracks`), with the empty result of its kind
if not the default.
"""


def collect_footprints(
    platforms: AllPlatforms | list[AllPlatforms],
    times: list[datetime] | None = None,
    instrument_index: int | None = 0,
    elevation: float | Mapping[str, float | None] | None = None,
    mask: Mask | None = None,
    platform_altaz: bool = False,
    solar_altaz: bool = False,
    time_step: timedelta | None = None,
    pixels: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the footprints of instruments on any observing platforms
    together: satellites (see `tatc.analysis.ground_track.collect_ground_track`),
    aircraft (see `tatc.analysis.airborne_track.collect_airborne_ground_track`),
    fixed surface platforms (see
    `tatc.analysis.surface_track.collect_surface_track`), and radar stations
    (see `tatc.analysis.radar_track.collect_radar_track`), for example to
    find their coincident observations (see
    `tatc.analysis.coincidence.find_coincident_ground_track`) or their
    combined coverage (see `dissolve_footprints`).

    Each footprint names its platform in the column of its kind
    (`satellite`, `aircraft`, `surface_platform`, or `radar_station`), and
    has the other columns of its kind's results (missing for the
    footprints of other kinds). A footprint without a time (of a surface
    platform or a radar station without `times`, whose footprint does not
    change) is static.

    Args:
        platforms (AllPlatforms | list[AllPlatforms]): The observing
                platform(s), of any kinds (see `tatc.schemas.AllPlatforms`).
        times (list[datetime.datetime] | None): The times to sample:
                required for satellites; for an aircraft, those that its
                trajectory contains (by default, the times of its
                trajectory); for a surface platform or a radar station,
                `None` for a static footprint.
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` for every instrument of each
                platform (a radar station has none).
        elevation (float | Mapping[str, float | None] | None): The
                elevation (meters) above the WGS 84 ellipsoid at which to
                compute footprints, or that of each kind of platform, by its
                column (for example, `{"radar_station": 3000}`, the others
                by default): by default, the ellipsoid (zero) for satellites
                and aircraft, and the ground elevation of each surface
                platform. Required for radar stations (see
                `collect_radar_track`). Footprints of several kinds at
                different elevations warn: they are compared in longitude
                and latitude (as by `find_coincident_ground_track`), and a
                view off nadir is offset between elevations.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.
        platform_altaz (bool): `True` to include the altitude/azimuth angles
                of each platform seen from the footprint center
                (`platform_alt` and `platform_az`; `NaN` for radar
                stations, whose footprints are at many elevation angles).
        solar_altaz (bool): `True` to include the solar altitude/azimuth
                angles (`solar_alt` and `solar_az`) and the glint angle
                (`glint_angle`) at the footprint center: requires `times` for
                surface platforms.
        time_step (datetime.timedelta | None): The time step that each time
                represents, over which the footprint of a swept instrument
                on a satellite or an aircraft sweeps (see
                `collect_ground_track`).
        pixels (bool): `True` to collect the centers of the pixels of each
                instrument (a rectangular `PointedInstrument`, see
                `tatc.analysis.ground_track.collect_ground_pixels`) rather
                than its footprint: one row per pixel, with its validity and
                angles seen from it. Radar stations have no pixels.

    Returns:
        geopandas.GeoDataFrame: The data frame of footprints (or pixels):
            for a single platform and instrument (an integer
            `instrument_index`), those of its kind's analysis; otherwise,
            those of each kind of platform, concatenated and sorted by time
            (footprints without a time last, in the order of their
            platforms).
    """
    single = _is_single(platforms, instrument_index=instrument_index)
    groups = _group_platforms(platforms)
    elevations = _get_kind_elevations(elevation)
    if groups[SATELLITE] and times is None:
        raise ValueError("times are required to collect the footprints of satellites")
    if groups[RADAR_STATION] and pixels:
        raise ValueError("radar stations have no pixels")
    if groups[RADAR_STATION] and elevations[RADAR_STATION] is None:
        raise ValueError(
            "elevation is required to collect the footprints of radar stations "
            "(see collect_radar_track)"
        )
    request = _FootprintRequest.create(
        times,
        instrument_index,
        None,
        mask,
        platform_altaz,
        solar_altaz,
        time_step,
        pixels,
    )
    frames = [
        _collect_tracks(
            members,
            kind,
            replace(request, elevation=elevations[kind]),
            *_COLLECTORS[kind],
        )
        for kind, members in groups.items()
        if members
    ]
    if len(frames) == 0:
        return gpd.GeoDataFrame(geometry=gpd.GeoSeries([], crs="EPSG:4326"))
    if single or len(frames) == 1:
        return frames[0]
    _warn_kind_elevations(frames)
    filled = [frame for frame in frames if not frame.empty]
    combined = pd.concat(filled if filled else frames, ignore_index=True)
    if "time" in combined.columns:
        combined = combined.sort_values("time", kind="stable", na_position="last")
    return gpd.GeoDataFrame(combined, crs="EPSG:4326").reset_index(drop=True)


def dissolve_footprints(
    footprints: gpd.GeoDataFrame, by: str | list[str] | None = None
) -> gpd.GeoDataFrame:
    """
    Dissolve the footprints of valid observations (those whose `valid_obs`
    is `True`, if the column is given) of any platforms (see
    `collect_footprints`) into one geometry (for example, the area observed
    over the times) or, with `by`, one geometry for each group (for example,
    `by="aircraft"` for the area each aircraft observes). The `compute_*`
    track functions (such as `compute_surface_track`) dissolve their tracks
    this way.

    Args:
        footprints (geopandas.GeoDataFrame): The footprints.
        by (str | list[str] | None): The column(s) by which to group the
                footprints (missing values form their own group), or `None`
                for one geometry.

    Returns:
        geopandas.GeoDataFrame: The dissolved footprints: one row (with
            column `geometry`), or one row for each group (with the columns
            of `by`), or none if no observation is valid.
    """
    return _dissolve_footprints(footprints, by)


def collect_footprint_observations(
    points: (
        Point | geo.Point | list[Point | geo.Point] | gpd.GeoDataFrame | gpd.GeoSeries
    ),
    footprints: gpd.GeoDataFrame,
    max_gap: timedelta | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the observations of points by the footprints of any platforms
    (see `collect_footprints`): the periods when each point lies within the
    valid footprints of a platform's instrument (`valid_obs`, if given) at
    consecutive times of its track. Unlike
    `tatc.analysis.point_sampling.collect_observations`, which finds the
    periods of a satellite's observations of points from its orbit, the
    periods are sampled at the times of the footprints: each spans the
    times of a run of consecutive footprints of an instrument that contain
    the point, from the first (`start`) to the last (`end`), with `epoch` at
    its midpoint, so that a point within a single footprint is observed at
    an instant. A footprint should tile the time step of its time (for
    example, with an along-track field of view tailored to the time step,
    see `tatc.utils.compute_along_track_field_of_view`) for the runs to be
    unbroken.

    Footprints are consecutive if no other footprint of the track lies
    between them and their times are at most `max_gap` apart, so that a run
    does not span the times without footprints (as those a mask culled, or
    a gap in an aircraft's trajectory). By default, `max_gap` is 1.5 times
    the shortest time between the footprints of each track (its time step).

    A point is compared with the footprints in longitude and latitude: a
    footprint projected to another elevation than the point's (its
    z-coordinate, or zero) is offset from it along the line of sight, which
    warns (see `collect_footprints` for the elevation of footprints).

    The observations have the columns of `collect_observations` (without
    angles), with each platform named in the column of its kind, so that
    `aggregate_observations`, `reduce_observations`, and
    `find_coincident_observations` compute the access, revisit, and
    coincidences of the observations of any platforms (also together with
    those of `collect_observations`, for the same points).

    Args:
        points (Point | shapely.geometry.Point | list[Point | shapely.geometry.Point] | geopandas.GeoDataFrame | geopandas.GeoSeries):
                The point(s) of interest (see `collect_observations`).
        footprints (geopandas.GeoDataFrame): The footprints, each with a
                time: of satellites, aircraft, or surface platforms with
                times.
        max_gap (datetime.timedelta | None): The longest time between
                consecutive footprints of a run, or `None` for 1.5 times the
                time step of each track.

    Returns:
        geopandas.GeoDataFrame: The data frame of observations, with columns
            `target_hash`, `geometry` (the point), the platform columns of
            the footprints, `instrument`, `start`, `epoch`, and `end`,
            sorted by start time.
    """
    if isinstance(points, gpd.GeoDataFrame):
        points = points.geometry
    if isinstance(points, gpd.GeoSeries):
        points = list(points)
    points = list(points) if isinstance(points, (list, tuple)) else [points]
    for point in points:
        if not isinstance(point, (Point, geo.Point)):
            raise TypeError(
                f"point must be a Point or shapely Point, not a {type(point).__name__}"
            )
    platform_columns = _get_platform_columns(footprints)
    if len(points) == 0 or footprints.empty:
        return _get_empty_observation_frame(platform_columns)
    if "time" not in footprints.columns or footprints["time"].isna().any():
        raise ValueError(
            "each footprint must have a time: collect the footprints of surface "
            "platforms and radar stations with times"
        )
    frame = footprints.reset_index(drop=True)
    platforms = _get_platforms(frame)
    times = pd.DatetimeIndex(frame["time"]).tz_convert("UTC").as_unit("ns")
    # the instrument of each footprint (none for a radar station)
    instruments = _get_instruments(frame)
    # the track of each platform's instrument, and the position of each
    # footprint in it (in time order)
    track, _ = pd.factorize(
        pd.Series(platforms) + "\x00" + pd.Series(instruments).astype(str).to_numpy()
    )
    order = np.lexsort((times.asi8, track))
    position = np.empty(len(frame), dtype=int)
    position[order] = np.arange(len(frame))
    limit = _get_max_gaps(track[order], times.asi8[order], max_gap)
    # the valid footprints that contain each point
    valid = (
        frame["valid_obs"].eq(True).to_numpy()
        if "valid_obs" in frame.columns
        else np.ones(len(frame), dtype=bool)
    )
    rows = np.flatnonzero(valid)
    coordinates = [_get_point_coordinates(point) for point in points]
    geometries = [geo.Point(x, y, z) for x, y, z in coordinates]
    tree = shapely.STRtree(shapely.force_2d(frame.geometry.to_numpy()[rows]))
    point, row = tree.query(
        shapely.force_2d(np.array(geometries, dtype=object)), predicate="intersects"
    )
    row = rows[row]
    if len(point) == 0:
        return _get_empty_observation_frame(platform_columns)
    _warn_elevations(
        np.array([z for _, _, z in coordinates])[point],
        _get_footprint_elevations(frame.geometry.to_numpy()[row]),
    )
    # the runs of consecutive footprints of each track that contain each point
    pairs = np.lexsort((position[row], point))
    point, row = point[pairs], row[pairs]
    new = np.ones(len(row), dtype=bool)
    new[1:] = (
        (point[1:] != point[:-1])
        | (track[row[1:]] != track[row[:-1]])
        | (position[row[1:]] != position[row[:-1]] + 1)
        | (times.asi8[row[1:]] - times.asi8[row[:-1]] > limit[track[row[1:]]])
    )
    first = np.flatnonzero(new)
    last = np.append(first[1:], len(row)) - 1
    start = times[row[first]]
    end = times[row[last]]
    observed = point[first]
    hashes = _hash_geometries([geometries[p] for p in np.unique(observed)])
    target_hash = dict(zip(np.unique(observed), hashes))
    observations = gpd.GeoDataFrame(
        {
            "target_hash": [target_hash[p] for p in observed],
            "geometry": [geometries[p] for p in observed],
            **{
                column: frame[column].to_numpy()[row[first]]
                for column in platform_columns
            },
            "instrument": instruments[row[first]],
            "start": start,
            "epoch": start + (end - start) / 2,
            "end": end,
        },
        crs="EPSG:4326",
    )
    return observations.sort_values("start", kind="stable").reset_index(drop=True)


_MAX_GAP_STEPS = 1.5
"""
The default longest time between consecutive footprints of a run, in time
steps of their track (see `collect_footprint_observations`).
"""

_ELEVATION_TOLERANCE = 1.0
"""
The difference (meters) between the elevations of a point and of a footprint
that contains it above which `collect_footprint_observations` warns.
"""


def _get_max_gaps(
    track: np.ndarray, nanoseconds: np.ndarray, max_gap: timedelta | None
) -> np.ndarray:
    """
    Gets the longest time (nanoseconds) between consecutive footprints of
    each track: `max_gap`, or 1.5 times the shortest positive time between
    the footprints of the track (infinite without one).

    Args:
        track (numpy.ndarray): The track code (from zero) of each footprint,
            in order of track and time.
        nanoseconds (numpy.ndarray): The time (nanoseconds) of each
            footprint, in the same order.
        max_gap (datetime.timedelta | None): The longest time, or `None` for
            that of each track's time step.

    Returns:
        numpy.ndarray: the longest time (nanoseconds) of each track code
    """
    count = int(track.max()) + 1 if len(track) > 0 else 0
    if max_gap is not None:
        if max_gap < timedelta(0):
            raise ValueError(f"max_gap must be non-negative, not {max_gap}")
        return np.full(count, pd.Timedelta(max_gap).value, dtype=float)
    step = np.full(count, np.inf)
    gaps = np.diff(nanoseconds).astype(float)
    within = (track[1:] == track[:-1]) & (gaps > 0)
    np.minimum.at(step, track[1:][within], gaps[within])
    return _MAX_GAP_STEPS * step


def _get_footprint_elevations(geometries: np.ndarray) -> np.ndarray:
    """
    Gets the elevation of each footprint: the z-coordinate of its first
    vertex (footprints are projected to a constant elevation), or `NaN` if
    it is empty or has none.

    Args:
        geometries (numpy.ndarray): The footprints.

    Returns:
        numpy.ndarray: the elevation (meters) of each footprint
    """
    elevations = np.full(len(geometries), np.nan)
    if len(geometries) == 0:
        return elevations
    coordinates, index = shapely.get_coordinates(
        geometries, include_z=True, return_index=True
    )
    rows, first = np.unique(index, return_index=True)
    elevations[rows] = coordinates[first, 2]
    return elevations


def _warn_elevations(points: np.ndarray, footprints: np.ndarray) -> None:
    """
    Warns if points lie within footprints projected to other elevations
    (see `collect_footprint_observations`).

    Args:
        points (numpy.ndarray): The elevation (meters) of each point.
        footprints (numpy.ndarray): The elevation (meters) of the footprint
            that contains each point (`NaN` if unknown).
    """
    differs = np.abs(points - footprints) > _ELEVATION_TOLERANCE
    if not np.any(differs):
        return
    i = np.flatnonzero(differs)[0]
    warn(
        f"points lie within footprints at other elevations (such as a point at "
        f"{points[i]:.0f} m within a footprint at {footprints[i]:.0f} m), which "
        "are offset from them along the line of sight: give the points the "
        "elevation of the footprints (as a z-coordinate), or collect the "
        "footprints at the elevation of the points (see collect_footprints)"
    )


def _get_kind_elevations(
    elevation: float | Mapping[str, float | None] | None,
) -> dict[_PlatformKind, float | None]:
    """
    Gets the elevation at which to compute the footprints of each kind of
    platform (see `collect_footprints`): the elevation, or that given for
    the kind's column, or `None` for the default of the kind.

    Args:
        elevation (float | Mapping[str, float | None] | None): The elevation
            (meters), or that of each kind of platform by its column.

    Returns:
        dict[_PlatformKind, float | None]: the elevation of each kind
    """
    if not isinstance(elevation, Mapping):
        return {kind: elevation for kind in _PLATFORM_KINDS}
    columns = [kind.column for kind in _PLATFORM_KINDS]
    for column in elevation:
        if column not in columns:
            raise ValueError(
                f"elevation must be given for platform columns ({', '.join(columns)}), "
                f"not {column!r}"
            )
    return {kind: elevation.get(kind.column) for kind in _PLATFORM_KINDS}


def _warn_kind_elevations(frames: list[gpd.GeoDataFrame]) -> None:
    """
    Warns if the footprints of different kinds of platforms are at different
    elevations (see `collect_footprints`).

    Args:
        frames (list[geopandas.GeoDataFrame]): The footprints of each kind.
    """
    ranges = {}
    for frame in frames:
        elevations = _get_footprint_elevations(frame.geometry.to_numpy())
        if frame.empty or np.all(np.isnan(elevations)):
            continue
        ranges[_get_platform_columns(frame)[0]] = (
            np.nanmin(elevations),
            np.nanmax(elevations),
        )
    if len(ranges) < 2:
        return
    lows, highs = zip(*ranges.values())
    if max(highs) - min(lows) <= _ELEVATION_TOLERANCE:
        return
    described = ", ".join(
        (
            f"{column} at {low:.0f} m"
            if high - low <= _ELEVATION_TOLERANCE
            else f"{column} at {low:.0f} to {high:.0f} m"
        )
        for column, (low, high) in ranges.items()
    )
    warn(
        f"footprints of different kinds of platforms are at different elevations "
        f"({described}), which are offset from each other along their lines of "
        "sight: give elevation (or one for each kind of platform) to compare "
        "them at the same elevation"
    )
