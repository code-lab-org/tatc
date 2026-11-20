"""
Methods to collect, dissolve, and sample the footprints of any observing
platforms together: satellites, aircraft, fixed surface platforms, and
radar stations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely import geometry as geo

from ..schemas import AllPlatforms, Point
from ..utils.geometry import _get_point_coordinates, _hash_geometries
from ._tracks import Mask
from .airborne_track import _collect_airborne
from .check import (
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
)
from .ground_track import _collect_satellites
from .radar_track import collect_radar_track
from .surface_track import _collect_surface


def collect_footprints(
    platforms: AllPlatforms | list[AllPlatforms],
    times: list[datetime] | None = None,
    instrument_index: int | None = 0,
    elevation: float | None = None,
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
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid at which to compute footprints: by default, the
                ellipsoid (zero) for satellites and aircraft, and the ground
                elevation of each surface platform. Required for radar
                stations, for which it is in the vertical reference of
                their elevations (see `collect_radar_track`).
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
    if groups[SATELLITE] and times is None:
        raise ValueError("times are required to collect the footprints of satellites")
    if groups[RADAR_STATION] and pixels:
        raise ValueError("radar stations have no pixels")
    if groups[RADAR_STATION] and elevation is None:
        raise ValueError(
            "elevation is required to collect the footprints of radar stations "
            "(see collect_radar_track)"
        )
    frames = []
    for kind, members in groups.items():
        if not members:
            continue
        if kind is RADAR_STATION:
            frames.append(
                collect_radar_track(
                    members, elevation, mask, times=times  # type: ignore
                )
            )
        elif kind is SURFACE_PLATFORM:
            # by default, at the ground elevation of each platform
            frames.append(
                _collect_surface(
                    members,
                    times,
                    instrument_index,
                    elevation,
                    mask,
                    platform_altaz,
                    solar_altaz,
                    pixels,
                )
            )
        else:
            collect = _collect_satellites if kind is SATELLITE else _collect_airborne
            frames.append(
                collect(
                    members,
                    times,  # type: ignore
                    instrument_index,
                    0 if elevation is None else elevation,
                    mask,
                    platform_altaz,
                    solar_altaz,
                    time_step,
                    pixels,
                )
            )
    if len(frames) == 0:
        return gpd.GeoDataFrame(geometry=gpd.GeoSeries([], crs="EPSG:4326"))
    if single or len(frames) == 1:
        return frames[0]
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
    # the runs of consecutive footprints of each track that contain each point
    pairs = np.lexsort((position[row], point))
    point, row = point[pairs], row[pairs]
    new = np.ones(len(row), dtype=bool)
    new[1:] = (
        (point[1:] != point[:-1])
        | (track[row[1:]] != track[row[:-1]])
        | (position[row[1:]] != position[row[:-1]] + 1)
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
