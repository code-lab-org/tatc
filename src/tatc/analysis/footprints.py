"""
Methods to collect, dissolve, and sample the footprints of any observing
platforms together: satellites, aircraft, fixed surface platforms, and
radar stations, including their observations of points and regions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely import geometry as geo
from skyfield.api import wgs84

from .._warnings import warn
from ..schemas import AllPlatforms, Point
from ..utils.geometry import (
    _get_point_coordinates,
    _hash_geometries,
    project_polygon_to_elevation,
    split_polygon,
)
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
    _get_platform_key,
    _get_platforms,
    _group_platforms,
    _is_single,
    _PlatformKind,
)
from .ground_track import _collect_satellite
from .radar_track import _collect_radar_station, _get_empty_radar_track
from .region_sampling import _get_polygons, _get_region_elevation
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
                by default): by default, that of each kind (see
                `tatc.analysis.check._PlatformKind.default_elevation`): the
                ellipsoid (zero) for satellites and aircraft, and the ground
                elevation of each surface platform. Required for radar
                stations (see `collect_radar_track`). Footprints of several kinds at
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
    missing = _get_kind_without_elevation(groups, elevations)
    if missing is not None:
        raise ValueError(
            f"elevation is required to collect the footprints of {missing.label}s: "
            f"give elevation, or one for the {missing.column} column"
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
    targets: (
        Point
        | geo.Point
        | geo.Polygon
        | geo.MultiPolygon
        | list[Point | geo.Point | geo.Polygon | geo.MultiPolygon]
        | gpd.GeoDataFrame
        | gpd.GeoSeries
    ),
    footprints: gpd.GeoDataFrame,
    max_gap: timedelta | None = None,
    platforms: AllPlatforms | list[AllPlatforms] | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the observations of points and regions by the footprints of any
    platforms (see `collect_footprints`): the periods when a target lies
    within (a point) or overlaps (a region) the valid footprints of a
    platform's instrument (`valid_obs`, if given) at consecutive times of its
    track. Unlike `tatc.analysis.point_sampling.collect_observations` and
    `tatc.analysis.region_sampling.collect_region_observations`, which find
    the periods of a satellite's observations from its orbit, the periods
    are sampled at the times of the footprints, for satellites, aircraft,
    fixed surface platforms, and radar stations alike: each spans the times
    of a run of consecutive footprints of an instrument that contain the
    point (or overlap the region in a positive area), from the first
    (`start`) to the last (`end`), with `epoch` at its midpoint, so that a
    target within a single footprint is observed at an instant. A footprint
    should tile the time step of its time (for example, with an
    along-track field of view tailored to the time step, see
    `tatc.utils.compute_along_track_field_of_view`) for the runs to be
    unbroken.

    Footprints are consecutive if no other footprint of the track lies
    between them and their times are at most `max_gap` apart, so that a run
    does not span the times without footprints (as those a mask culled, or
    a gap in an aircraft's trajectory). By default, `max_gap` is 1.5 times
    the shortest time between the footprints of each track (its time step).

    Each point is a TAT-C `Point` or a shapely `Point` (at the elevation of
    its z-coordinate, or zero), and each region a shapely `Polygon` or
    `MultiPolygon` (at the elevation of the mean of its z-coordinates, or
    zero), split along the anti-meridian and poles (see
    `tatc.utils.geometry.split_polygon`), as for `collect_observations` and
    `collect_region_observations`. The geometry of a point's observation is
    the point, and that of a region's observation its swath: the part of the
    region that the footprints of the run cover (at the region's elevation,
    if it has z-coordinates). A target is compared with the footprints in
    longitude and latitude: a footprint projected to another elevation than
    the target's is offset from it along the line of sight, which warns
    (see `collect_footprints` for the elevation of footprints).

    The observations have the columns of `collect_observations` (and of
    `collect_region_observations`), with each platform named in the column
    of its kind and the `target_hash` of each target as they give it, so
    that `aggregate_observations`, `reduce_observations`, and
    `find_coincident_observations` compute the access, revisit, and
    coincidences of the observations of any platforms (also together with
    those of satellites from `collect_observations` and
    `collect_region_observations`, for the same targets). With the
    `platforms` of the footprints, the observations of points also have
    their platform's altitude and azimuth angles seen from the point at the
    epoch (`platform_alt` and `platform_az`, degrees, as for
    `collect_observations`; `NaN` for regions).

    Args:
        targets (Point | shapely.geometry.Point | shapely.geometry.Polygon | shapely.geometry.MultiPolygon | list[Point | shapely.geometry.Point | shapely.geometry.Polygon | shapely.geometry.MultiPolygon] | geopandas.GeoDataFrame | geopandas.GeoSeries):
                The point(s) and region(s) of interest (the geometries of a
                data frame or series).
        footprints (geopandas.GeoDataFrame): The footprints, each with a
                time: of satellites, aircraft, or surface platforms with
                times.
        max_gap (datetime.timedelta | None): The longest time between
                consecutive footprints of a run, or `None` for 1.5 times the
                time step of each track.
        platforms (AllPlatforms | list[AllPlatforms] | None): The platforms
                of the footprints (and possibly others), to compute the
                angles of the observations of points, or `None` to omit them.

    Returns:
        geopandas.GeoDataFrame: The data frame of observations, with columns
            `target_hash`, `geometry` (the point, or the region's swath),
            the platform columns of the footprints, `instrument`, `start`,
            `epoch`, and `end` (and, with `platforms`, `platform_alt` and
            `platform_az`), sorted by start time.
    """
    targets = _get_targets(targets)
    lookup = None if platforms is None else _get_platform_lookup(platforms)
    platform_columns = _get_platform_columns(footprints)
    empty = _get_empty_footprint_observations(platform_columns, lookup is not None)
    if len(targets) == 0 or footprints.empty:
        return empty
    if "time" not in footprints.columns or footprints["time"].isna().any():
        raise ValueError(
            "each footprint must have a time: collect the footprints of surface "
            "platforms and radar stations with times"
        )
    frame = footprints.reset_index(drop=True)
    keys = _get_platforms(frame)
    if lookup is not None:
        missing = sorted(set(keys) - set(lookup))
        if missing:
            raise ValueError(
                "platforms must include the platform of each footprint, such as "
                f"{missing[0]!r}"
            )
    times = pd.DatetimeIndex(frame["time"]).tz_convert("UTC").as_unit("ns")
    # the instrument of each footprint (none for a radar station)
    instruments = _get_instruments(frame)
    # the track of each platform's instrument, and the position of each
    # footprint in it (in time order)
    track, _ = pd.factorize(
        pd.Series(keys) + "\x00" + pd.Series(instruments).astype(str).to_numpy()
    )
    order = np.lexsort((times.asi8, track))
    position = np.empty(len(frame), dtype=int)
    position[order] = np.arange(len(frame))
    limit = _get_max_gaps(track[order], times.asi8[order], max_gap)
    # the valid footprints that contain (or overlap) each target
    valid = (
        frame["valid_obs"].eq(True).to_numpy()
        if "valid_obs" in frame.columns
        else np.ones(len(frame), dtype=bool)
    )
    rows = np.flatnonzero(valid)
    shapes = shapely.force_2d(frame.geometry.to_numpy())
    target, row = shapely.STRtree(shapes[rows]).query(
        targets.shapes, predicate="intersects"
    )
    row = rows[row]
    # the part of each region within each footprint, of positive area
    pieces = np.full(len(target), None, dtype=object)
    regional = ~targets.is_point[target]
    pieces[regional] = shapely.intersection(
        targets.shapes[target[regional]], shapes[row[regional]]
    )
    keep = ~regional
    keep[regional] = shapely.area(pieces[regional]) > 0
    target, row, pieces = target[keep], row[keep], pieces[keep]
    if len(target) == 0:
        return empty
    _warn_elevations(
        targets.elevations[target],
        _get_footprint_elevations(frame.geometry.to_numpy()[row]),
    )
    # the runs of consecutive footprints of each track that contain (or
    # overlap) each target
    pairs = np.lexsort((position[row], target))
    target, row, pieces = target[pairs], row[pairs], pieces[pairs]
    new = np.ones(len(row), dtype=bool)
    new[1:] = (
        (target[1:] != target[:-1])
        | (track[row[1:]] != track[row[:-1]])
        | (position[row[1:]] != position[row[:-1]] + 1)
        | (times.asi8[row[1:]] - times.asi8[row[:-1]] > limit[track[row[1:]]])
    )
    first = np.flatnonzero(new)
    last = np.append(first[1:], len(row)) - 1
    start = times[row[first]]
    end = times[row[last]]
    epoch = start + (end - start) / 2
    observed = target[first]
    observations = gpd.GeoDataFrame(
        {
            "target_hash": targets.hashes[observed],
            "geometry": [
                (
                    targets.geometries[t]
                    if targets.is_point[t]
                    else _get_swath(pieces[a : b + 1], targets.geometries[t])
                )
                for t, a, b in zip(observed, first, last)
            ],
            **{
                column: frame[column].to_numpy()[row[first]]
                for column in platform_columns
            },
            "instrument": instruments[row[first]],
            "start": start,
            "epoch": epoch,
            "end": end,
        },
        crs="EPSG:4326",
    )
    if lookup is not None:
        _add_point_angles(
            observations, targets, observed, keys[row[first]], epoch, lookup
        )
    return observations.sort_values("start", kind="stable").reset_index(drop=True)


@dataclass(frozen=True)
class _Targets:
    """
    The points and regions observed by footprints (see
    `collect_footprint_observations`).
    """

    geometries: np.ndarray
    """The geometry of each target recorded with its observations."""
    shapes: np.ndarray
    """The two-dimensional shape of each target compared with footprints."""
    is_point: np.ndarray
    """Whether each target is a point (otherwise a region)."""
    elevations: np.ndarray
    """The elevation (meters) of each target."""
    hashes: np.ndarray
    """The hash of each target (see `tatc.utils.geometry.hash_geometry`)."""

    def __len__(self) -> int:
        return len(self.geometries)


def _get_targets(targets: object) -> _Targets:
    """
    Gets the points and regions observed by footprints, as
    `collect_observations` and `collect_region_observations` take them.

    Args:
        targets (object): A point or region, a list of them, or a data frame
            or series of their geometries.

    Returns:
        _Targets: the targets
    """
    if isinstance(targets, gpd.GeoDataFrame):
        targets = targets.geometry
    if isinstance(targets, gpd.GeoSeries):
        targets = list(targets)
    targets = list(targets) if isinstance(targets, (list, tuple)) else [targets]
    geometries, shapes, elevations = [], [], []
    for value in targets:
        if isinstance(value, (Point, geo.Point)):
            x, y, z = _get_point_coordinates(value)
            geometries.append(geo.Point(x, y, z))
            shapes.append(geo.Point(x, y))
            elevations.append(z)
        elif isinstance(value, (geo.Polygon, geo.MultiPolygon)):
            geometries.append(value)
            shapes.append(split_polygon(value))
            elevations.append(_get_region_elevation(value))
        else:
            raise TypeError(
                "target must be a Point, shapely Point, Polygon, or MultiPolygon, "
                f"not a {type(value).__name__}"
            )
    return _Targets(
        geometries=_to_object_array(geometries),
        shapes=_to_object_array(shapes),
        is_point=np.array([isinstance(g, geo.Point) for g in geometries], dtype=bool),
        elevations=np.array(elevations, dtype=float),
        hashes=_to_object_array(_hash_geometries(geometries) if geometries else []),
    )


def _to_object_array(values: list) -> np.ndarray:
    """
    Gets a one-dimensional array of objects (such as geometries), which
    numpy would otherwise unpack if they are sequences.

    Args:
        values (list): The objects.

    Returns:
        numpy.ndarray: the array of the objects
    """
    array = np.empty(len(values), dtype=object)
    array[:] = values
    return array


def _get_swath(
    pieces: np.ndarray, region: geo.Polygon | geo.MultiPolygon
) -> geo.Polygon | geo.MultiPolygon:
    """
    Gets the swath of a run of footprints over a region: the union of the
    parts of the region within each footprint, at the region's elevation if
    it has z-coordinates (as of `collect_region_observations`).

    Args:
        pieces (numpy.ndarray): The (two-dimensional) part of the region
            within each footprint of the run.
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
            The region.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: the swath
    """
    swath = _get_polygons(shapely.union_all(pieces))
    if region.has_z:
        return project_polygon_to_elevation(swath, _get_region_elevation(region))
    return swath


def _get_platform_lookup(
    platforms: AllPlatforms | list[AllPlatforms],
) -> dict[str, AllPlatforms]:
    """
    Gets each platform by its key in results (see
    `tatc.analysis.check._get_platforms`).

    Args:
        platforms (AllPlatforms | list[AllPlatforms]): The platforms.

    Returns:
        dict[str, AllPlatforms]: the platform of each key
    """
    return {
        _get_platform_key(kind.column, platform.name): platform
        for kind, members in _group_platforms(platforms).items()
        for platform in members
    }


_ANGLE_COLUMNS = ("platform_alt", "platform_az")
"""
The columns of the angles of observations of points by footprints (see
`collect_footprint_observations`), as of `collect_observations`.
"""


def _get_empty_footprint_observations(
    platform_columns: list[str], angles: bool
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame of the observations of targets by footprints
    (see `collect_footprint_observations`).

    Args:
        platform_columns (list[str]): The columns that name the platforms.
        angles (bool): `True` to include the platform's angles.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    frame = _get_empty_observation_frame(platform_columns)
    if angles:
        for column in _ANGLE_COLUMNS:
            frame[column] = pd.Series([], dtype="float")
    return frame


def _add_point_angles(
    observations: gpd.GeoDataFrame,
    targets: _Targets,
    observed: np.ndarray,
    keys: np.ndarray,
    epoch: pd.DatetimeIndex,
    lookup: Mapping[str, AllPlatforms],
) -> None:
    """
    Adds the altitude and azimuth angles (degrees) of the platform of each
    observation of a point, seen from the point at its epoch (`NaN` for the
    observations of regions).

    Args:
        observations (geopandas.GeoDataFrame): The observations, modified in
            place.
        targets (_Targets): The targets.
        observed (numpy.ndarray): The target of each observation.
        keys (numpy.ndarray): The platform (key) of each observation.
        epoch (pandas.DatetimeIndex): The epoch of each observation.
        lookup (Mapping[str, AllPlatforms]): The platform of each key.
    """
    angles = {column: np.full(len(observations), np.nan) for column in _ANGLE_COLUMNS}
    points = targets.is_point[observed]
    for key in pd.unique(keys[points]):
        index = np.flatnonzero(points & (keys == key))
        geometries = targets.geometries[observed[index]]
        platform_track = lookup[key].get_track(list(epoch[index].to_pydatetime()))
        topos = wgs84.latlon(
            shapely.get_y(geometries),
            shapely.get_x(geometries),
            targets.elevations[observed[index]],
        )
        altitude, azimuth, _ = (platform_track - topos.at(platform_track.t)).altaz()
        angles["platform_alt"][index] = np.atleast_1d(altitude.degrees)
        angles["platform_az"][index] = np.atleast_1d(azimuth.degrees)
    for column, values in angles.items():
        observations[column] = values


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


def _warn_elevations(targets: np.ndarray, footprints: np.ndarray) -> None:
    """
    Warns if targets (points or regions) lie within footprints projected to
    other elevations (see `collect_footprint_observations`).

    Args:
        targets (numpy.ndarray): The elevation (meters) of each target.
        footprints (numpy.ndarray): The elevation (meters) of the footprint
            that contains (or overlaps) each target (`NaN` if unknown).
    """
    differs = np.abs(targets - footprints) > _ELEVATION_TOLERANCE
    if not np.any(differs):
        return
    i = np.flatnonzero(differs)[0]
    warn(
        f"targets lie within footprints at other elevations (such as a target at "
        f"{targets[i]:.0f} m within a footprint at {footprints[i]:.0f} m), which "
        "are offset from them along the line of sight: give the targets the "
        "elevation of the footprints (as a z-coordinate), or collect the "
        "footprints at the elevation of the targets (see collect_footprints)"
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


def _get_kind_without_elevation(
    groups: Mapping[_PlatformKind, list],
    elevations: Mapping[_PlatformKind, float | None],
) -> _PlatformKind | None:
    """
    Gets the first kind of the platforms given without an elevation at
    which to compute their footprints and without a default (see
    `tatc.analysis.check._PlatformKind.default_elevation`), as radar
    stations, if any.

    Args:
        groups (Mapping[_PlatformKind, list]): The platforms of each kind
            (see `tatc.analysis.check._group_platforms`).
        elevations (Mapping[_PlatformKind, float | None]): The elevation of
            each kind (see `_get_kind_elevations`).

    Returns:
        _PlatformKind | None: the kind, or `None` if every kind has one
    """
    return next(
        (
            kind
            for kind, members in groups.items()
            if members and elevations[kind] is None and kind.default_elevation is None
        ),
        None,
    )


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
