"""
Methods to perform coverage analysis of regions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely import geometry as geo
from skyfield.api import wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric

from ..constants import (
    EARTH_ECCENTRICITY,
    EARTH_EQUATORIAL_RADIUS,
    EARTH_FLATTENING,
    EARTH_MU,
    EARTH_POLAR_RADIUS,
    EARTH_ROTATION_RATE,
    timescale,
)
from ..schemas import (
    ConicalInstrument,
    GeneralPerturbationsOrbit,
    Instrument,
    PointedInstrument,
    Satellite,
)
from ..utils.ellipsoid import (
    _get_geodetic_coordinates,
    _get_surface_directions,
    _get_surface_positions,
)
from ..utils.geometry import (
    _get_angular_distance_to_arcs,
    _get_boundary_arcs,
    _get_nearest_arc_points,
    hash_geometry,
    project_polygon_to_elevation,
    split_polygon,
)
from ..utils.projection import NadirReference, VelocityFrame, _compute_view_frame
from ..utils.propagation import _to_time_from_offsets
from .check import _check_satellite, _check_satellites
from .sampling import _refine_access_periods


def _get_visible_polygon_interval_series(
    geometry: geo.Polygon | geo.MultiPolygon,
    satellite: Satellite,
    field_of_regard: float,
    start: datetime,
    end: datetime,
    elevation: float = 0,
    margin: float = 0.1,
    coarse_step: float = 10,
) -> pd.Series:
    """
    Get the series of periods when an instrument's field of regard (a cone
    about nadir) may observe any part of a region: a conservative superset
    of the observation periods, for culling, from which no observation is
    missed however brief.

    The region is observable when the angular distance `d(t)` from the
    satellite's geocentric direction to the region (zero within it) is at
    most the Earth central angle `lambda` from nadir to the edge of the field
    of regard (or the horizon), so the periods are those when

        g(t) = d(t) - lambda - margin

    is not positive. For a single point, this is equivalent to a minimum
    elevation angle condition (see
    `tatc.analysis.point_sampling.compute_access_periods`). The
    central angle `lambda` is evaluated at the orbit's apoapsis (its largest
    over the orbit) on a sphere of the region's smallest geocentric radius
    (larger than on the ellipsoid), and widened by the largest separation of
    the geodetic and geocentric nadir directions as seen on the ground, so
    that it is conservative over the whole orbit.

    `d(t)` changes no faster than the angular rate `omega` of the satellite's
    direction in the Earth-fixed frame (bounded by the orbital angular rate
    at periapsis plus the Earth's rotation), so no period lies within an
    interval `[t_a, t_b]` in which `g(t_a) + g(t_b) > omega * (t_b - t_a)`, and
    the whole interval lies within a period if `g(t_a) + g(t_b) < -omega *
    (t_b - t_a)`. Starting from samples spaced by `coarse_step` degrees of the
    satellite's motion, intervals that are neither are halved until they
    are, or until they are shorter than a millisecond (kept as part of a
    period).

    The region's boundary is treated as great circle arcs between points at
    most a degree apart along its edges (see
    `tatc.utils.geometry._get_boundary_arcs`). The region is split along the
    anti-meridian and poles (see `tatc.utils.geometry.split_polygon`), so
    that its parts can be tested in longitude and latitude.

    Args:
        geometry (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region (with EPSG:4326 coordinates) to observe.
        satellite (Satellite): Satellite doing the observation.
        field_of_regard (float): The instrument's field of regard (degrees).
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        elevation (float): Elevation (meters) of the region above the WGS 84 ellipsoid.
        margin (float): Additional central angle (degrees) to widen the field of regard.
        coarse_step (float): Angle (degrees) of the satellite's motion between initial samples.

    Returns:
        pandas.Series: Series of observation intervals.
    """
    geometry = split_polygon(geometry)
    shapely.prepare(geometry)
    arcs = _get_boundary_arcs(geometry, elevation)
    orbit = satellite.orbit.to_gp_orbit()
    shapes = [
        (element.get_semimajor_axis(), element.eccentricity)
        for element in orbit.elements
    ]
    # bounds on the orbit radius, with a margin for the short-period
    # perturbations of the osculating radius about the mean elements, and on
    # the angular rate of the satellite's direction in the Earth-fixed frame
    # (radians/day), with a 2 percent margin, over all elements
    apoapsis = max(a * (1 + e) for a, e in shapes) + 25e3
    rate = 86400 * (
        1.02
        * max(
            np.sqrt(EARTH_MU * (1 + e) / (a * (1 - e))) / (a * (1 - e))
            for a, e in shapes
        )
        + EARTH_ROTATION_RATE
    )
    # earth central angle from nadir to the edge of the field of regard, or
    # to the horizon, at apoapsis, on a sphere of the smallest (geocentric)
    # radius of the region
    _, min_latitude, _, max_latitude = geometry.bounds
    latitude = np.radians(max(abs(min_latitude), abs(max_latitude)))
    radius = elevation + np.sqrt(
        (
            (EARTH_EQUATORIAL_RADIUS**2 * np.cos(latitude)) ** 2
            + (EARTH_POLAR_RADIUS**2 * np.sin(latitude)) ** 2
        )
        / (
            (EARTH_EQUATORIAL_RADIUS * np.cos(latitude)) ** 2
            + (EARTH_POLAR_RADIUS * np.sin(latitude)) ** 2
        )
    )
    half_angle = np.radians(min(field_of_regard / 2, 90))
    ratio = apoapsis / radius * np.sin(half_angle)
    central_angle = (
        np.arcsin(ratio) - half_angle if ratio < 1 else np.arccos(radius / apoapsis)
    )
    # separation on the ground of the geodetic and geocentric nadir
    # directions (which differ by up to about the flattening, in radians)
    nadir_offset = (apoapsis - radius) / radius * EARTH_FLATTENING
    threshold = central_angle + nadir_offset + np.radians(margin)
    e2 = EARTH_ECCENTRICITY**2

    def excess(x: np.ndarray) -> np.ndarray:
        position = np.array(
            orbit.get_orbit_track_at_time(timescale.tt_jd(x)).frame_xyz(itrs).m
        ).reshape(3, -1)
        u = position / np.linalg.norm(position, axis=0)
        distance = _get_angular_distance_to_arcs(u, arcs)
        # geodetic latitude of the surface point in the satellite's direction
        longitude = np.degrees(np.arctan2(u[1], u[0]))
        latitude = np.degrees(np.arctan2(u[2], (1 - e2) * np.hypot(u[0], u[1])))
        distance[shapely.contains_xy(geometry, longitude, latitude)] = 0
        return distance - threshold

    t_start = timescale.from_datetime(start).tt
    t_end = max(timescale.from_datetime(end).tt, t_start)
    count = int(np.ceil((t_end - t_start) * rate / np.radians(coarse_step))) + 1
    x = np.linspace(t_start, t_end, max(count, 2))
    f = excess(x)
    lower, upper, f_lower, f_upper = x[:-1], x[1:], f[:-1], f[1:]
    tolerance = 1e-3 / 86400
    windows = []
    while len(lower) > 0:
        span = rate * (upper - lower)
        inside = f_lower + f_upper < -span
        unresolved = ~inside & (f_lower + f_upper <= span)
        # keep intervals within a period, and those too short to resolve
        small = unresolved & (upper - lower <= tolerance)
        windows.extend(zip(lower[inside | small], upper[inside | small]))
        halve = unresolved & ~small
        lower, upper = lower[halve], upper[halve]
        f_lower, f_upper = f_lower[halve], f_upper[halve]
        if len(lower) > 0:
            middle = (lower + upper) / 2
            f_middle = excess(middle)
            lower, upper = np.concatenate((lower, middle)), np.concatenate(
                (middle, upper)
            )
            f_lower = np.concatenate((f_lower, f_middle))
            f_upper = np.concatenate((f_middle, f_upper))
    # merge adjoining intervals into periods
    periods = []
    for left, right in sorted(windows):
        if len(periods) > 0 and left <= periods[-1][1]:
            periods[-1][1] = max(periods[-1][1], right)
        else:
            periods.append([left, right])
    if len(periods) == 0:
        return pd.Series([], dtype="interval")
    bounds = timescale.tt_jd(np.array(periods).flatten()).utc_datetime()
    first = pd.Timestamp(start.astimezone(tz=timezone.utc))
    last = pd.Timestamp(end.astimezone(tz=timezone.utc))
    # widen the periods by a millisecond to cover the precision of the
    # (Julian date) times at which they were found
    epsilon = pd.Timedelta(milliseconds=1)
    return pd.Series(
        [
            pd.Interval(
                left=max(pd.Timestamp(left) - epsilon, first),
                right=min(pd.Timestamp(right) + epsilon, last),
                closed="both",
            )
            for left, right in zip(bounds[::2], bounds[1::2])
        ],
        dtype="interval",
    )


def _get_region_elevation(region: geo.Polygon | geo.MultiPolygon) -> float:
    """
    Gets the elevation of a region: the mean of its z coordinates (meters
    above the WGS 84 ellipsoid), or zero if it has none.

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The region.

    Returns:
        float: the elevation (meters)
    """
    if not region.has_z:
        return 0.0
    return float(np.mean(shapely.get_coordinates(region, include_z=True)[:, 2]))


def compute_region_access_periods(
    region: geo.Polygon | geo.MultiPolygon,
    satellite: Satellite,
    start: datetime,
    end: datetime,
    field_of_regard: float = 180,
    elevation: float | None = None,
    margin: float = 0.1,
) -> pd.Series:
    """
    Compute a conservative superset of the periods when an instrument's
    field of regard (a cone about nadir, or the horizon) may observe any
    part of a region, from which no observation is missed however brief:
    for example, to cull the times at which to propagate an orbit or compute
    footprints. The field of regard is evaluated conservatively (at the
    orbit's apoapsis, on a sphere of the region's smallest radius, and
    widened by a margin), so the periods may begin somewhat before and end
    somewhat after the region is observed. For the periods when the region
    is observed, see `collect_region_observations`.

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region (longitude and latitude in degrees).
        satellite (Satellite): The satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        field_of_regard (float): The instrument's field of regard (degrees;
                180, the default, for the horizon).
        elevation (float | None): The region's elevation (meters) above the
                WGS 84 ellipsoid; by default, the mean of its z coordinates,
                if any, or otherwise zero.
        margin (float): Additional central angle (degrees) to widen the field of regard.

    Returns:
        pandas.Series: the access periods (`pandas.Interval` of UTC
            timestamps), in time order.
    """
    _check_satellite(satellite)
    if not isinstance(region, (geo.Polygon, geo.MultiPolygon)):
        raise TypeError(
            f"region must be a Polygon or MultiPolygon, not a {type(region).__name__}"
        )
    if elevation is None:
        elevation = _get_region_elevation(region)
    return _get_visible_polygon_interval_series(
        region, satellite, field_of_regard, start, end, elevation, margin
    )


def _get_region_view(
    region: geo.Polygon | geo.MultiPolygon,
    arcs: tuple[np.ndarray, np.ndarray],
    orbit_track: Geocentric,
    nadir_reference: NadirReference,
    elevation: float = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Gets the view from a satellite of a region's point nearest to nadir:
    the nadir point itself, if the region contains it, or otherwise the
    point of its boundary at the least angular distance from the nadir
    point (see `tatc.utils.geometry._get_nearest_arc_points`).

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region, split along the anti-meridian and poles (see
                `tatc.utils.geometry.split_polygon`).
        arcs (tuple[numpy.ndarray, numpy.ndarray]): The region's boundary
                arcs (see `tatc.utils.geometry._get_boundary_arcs`).
        orbit_track (skyfield.positionlib.Geocentric): The satellite orbit track.
        nadir_reference (NadirReference): The definition of the nadir direction.
        elevation (float): The region's elevation (meters) above the WGS 84 ellipsoid.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:
            the point's angle from nadir (degrees), the satellite's elevation
            angle seen from the point (degrees), and the point's longitude
            and latitude (degrees)
    """
    position, nadir, _, _ = _compute_view_frame(
        orbit_track, 0, 0, VelocityFrame.EARTH_FIXED, nadir_reference
    )
    position, nadir = np.reshape(position, (3, -1)), np.reshape(nadir, (3, -1))
    # direction of the nadir point
    if nadir_reference == NadirReference.GEODETIC:
        subpoint = wgs84.geographic_position_of(orbit_track)
        origin = _get_surface_directions(
            np.reshape(subpoint.longitude.degrees, -1),
            np.reshape(subpoint.latitude.degrees, -1),
            elevation,
        )
    else:
        origin = position / np.linalg.norm(position, axis=0)
    longitude, latitude = _get_geodetic_coordinates(origin)
    inside = shapely.contains_xy(region, longitude, latitude)
    _, nearest = _get_nearest_arc_points(origin, arcs)
    nearest[:, inside] = origin[:, inside]
    longitude, latitude = _get_geodetic_coordinates(nearest)
    # line of sight from the satellite to the point
    los = _get_surface_positions(nearest, elevation) - position
    distance = np.linalg.norm(los, axis=0)
    angle = np.degrees(
        np.arccos(np.clip(np.sum(los * nadir, axis=0) / distance, -1, 1))
    )
    angle[inside] = 0
    # satellite elevation angle above the point's (geodetic) horizon
    lon, lat = np.radians(longitude), np.radians(latitude)
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    sat_elevation = np.degrees(
        np.arcsin(np.clip(-np.sum(los * up, axis=0) / distance, -1, 1))
    )
    return angle, sat_elevation, longitude, latitude


def _get_empty_region_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for region observations.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "target_hash": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _find_footprint_periods(
    region: geo.Polygon | geo.MultiPolygon,
    orbit: GeneralPerturbationsOrbit,
    instrument: PointedInstrument | ConicalInstrument,
    windows: list[pd.Interval],
    elevation: float = 0,
    time_step: timedelta = timedelta(seconds=10),
    tolerance: timedelta = timedelta(milliseconds=1),
) -> list[pd.Interval]:
    """
    Finds the periods within windows (which must contain them) when an
    instrument's footprint intersects a region. Footprints are sampled at
    most `time_step` apart within each window. Where neither of two
    consecutive footprints intersects the region but the convex hull of both
    (which contains the area the footprint sweeps between them) does, and
    their distances to the region do not exceed the distance a footprint
    moves between them, the interval is divided until a
    footprint intersects the region or the interval is shorter than
    `tolerance`, so that brief observations between samples are found. The
    start and end of each period are refined by bisection to within
    `tolerance`.

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region, split along the anti-meridian and poles (see
                `tatc.utils.geometry.split_polygon`).
        orbit (GeneralPerturbationsOrbit): The orbit.
        instrument (PointedInstrument | ConicalInstrument): The instrument.
        windows (list[pandas.Interval]): The windows to search.
        elevation (float): The elevation (meters) of the region above the
                WGS 84 ellipsoid, at which to project footprints.
        time_step (datetime.timedelta): The maximum time between samples.
        tolerance (datetime.timedelta): The precision of the period bounds.

    Returns:
        list[pandas.Interval]: the periods
    """
    if len(windows) == 0:
        return []
    reference = windows[0].left
    resolution = tolerance.total_seconds()

    def footprints(seconds: np.ndarray) -> np.ndarray:
        if len(seconds) == 0:
            return np.array([], dtype=object)
        orbit_track = orbit.get_orbit_track_at_time(
            _to_time_from_offsets(reference, seconds)
        )
        return np.array(
            instrument.compute_footprint(orbit_track, elevation=elevation),
            dtype=object,
        )

    def observes(views: np.ndarray) -> np.ndarray:
        return shapely.intersects(views, region)

    # bound on the angular rate (radians/second) at which footprints move
    # over the ground: the orbital angular rate at periapsis plus the
    # Earth's rotation, with a 50 percent margin for off-nadir views
    rate = 1.5 * (
        max(
            np.sqrt(EARTH_MU * (1 + e.eccentricity) / e.get_semimajor_axis() ** 3)
            / (1 - e.eccentricity) ** 1.5
            for e in orbit.elements
        )
        + EARTH_ROTATION_RATE
    )
    # orthographic projection (of the unit sphere) centered on the region,
    # under which distances on the hemisphere facing it do not exceed those
    # on the sphere (radians), at any latitude and across the anti-meridian;
    # boundaries are first divided into segments of at most 1 degree, as
    # they are straight in longitude and latitude
    center = region.representative_point()
    sin_lat0, cos_lat0 = np.sin(np.radians(center.y)), np.cos(np.radians(center.y))

    def project(coords: np.ndarray) -> np.ndarray:
        longitude = np.radians(coords[:, 0] - center.x)
        latitude = np.radians(coords[:, 1])
        return np.column_stack(
            [
                np.cos(latitude) * np.sin(longitude),
                cos_lat0 * np.sin(latitude)
                - sin_lat0 * np.cos(latitude) * np.cos(longitude),
            ]
        )

    def facing(geometries: np.ndarray) -> np.ndarray:
        # whether every vertex of each geometry is on the facing hemisphere
        coords, index = shapely.get_coordinates(geometries, return_index=True)
        cosine = sin_lat0 * np.sin(np.radians(coords[:, 1])) + cos_lat0 * np.cos(
            np.radians(coords[:, 1])
        ) * np.cos(np.radians(coords[:, 0] - center.x))
        result = np.ones(len(geometries), dtype=bool)
        np.logical_and.at(result, index, cosine > 0)
        return result

    projected_region = shapely.transform(shapely.segmentize(region, 1), project)

    def distance(views: np.ndarray) -> np.ndarray:
        segmented = shapely.segmentize(views, 1)
        return np.where(
            facing(segmented),
            np.nan_to_num(
                shapely.distance(
                    shapely.transform(segmented, project), projected_region
                )
            ),
            0,
        )

    def may_observe(
        views_a: np.ndarray, views_b: np.ndarray, seconds: np.ndarray
    ) -> np.ndarray:
        # the area swept between two footprints lies within their convex hull
        possible = shapely.intersects(
            shapely.convex_hull(shapely.union(views_a, views_b)), region
        )
        # and within the distance a footprint moves of either footprint
        return possible & (distance(views_a) + distance(views_b) <= rate * seconds)

    # samples within each window
    bounds = [
        ((w.left - reference).total_seconds(), (w.right - reference).total_seconds())
        for w in windows
    ]
    samples = [
        np.linspace(
            lo, hi, max(2, int(np.ceil((hi - lo) / time_step.total_seconds())) + 1)
        )
        for lo, hi in bounds
    ]
    counts = np.cumsum([len(s) for s in samples])[:-1]
    views = np.split(footprints(np.concatenate(samples)), counts)
    seen = [observes(v) for v in views]
    # search intervals where the footprint may have swept over the region
    # between samples that do not observe it
    times = [list(zip(s, p)) for s, p in zip(samples, seen)]
    search = [
        (k, s[i], s[i + 1], v[i], v[i + 1])
        for k, (s, v, p) in enumerate(zip(samples, views, seen))
        for i in np.flatnonzero(
            ~p[:-1] & ~p[1:] & may_observe(v[:-1], v[1:], np.diff(s))
        )
    ]
    while len(search) > 0:
        middle = np.array([(lo + hi) / 2 for _, lo, hi, _, _ in search])
        views_middle = footprints(middle)
        seen_middle = observes(views_middle)
        divided = []
        for (k, lo, hi, view_lo, view_hi), t, view, found in zip(
            search, middle, views_middle, seen_middle
        ):
            if found:
                times[k].append((t, True))
            elif hi - lo > 2 * resolution:
                divided += [(k, lo, t, view_lo, view), (k, t, hi, view, view_hi)]
        search = [
            interval
            for interval, keep in zip(
                divided,
                (
                    may_observe(
                        np.array([i[3] for i in divided], dtype=object),
                        np.array([i[4] for i in divided], dtype=object),
                        np.array([i[2] - i[1] for i in divided]),
                    )
                    if divided
                    else []
                ),
            )
            if keep
        ]
    # brackets of each change of state between (sorted) samples
    times = [sorted(t) for t in times]
    brackets = [
        (k, i)
        for k, t in enumerate(times)
        for i in range(len(t) - 1)
        if t[i][1] != t[i + 1][1]
    ]
    lower = np.array([times[k][i][0] for k, i in brackets])
    upper = np.array([times[k][i + 1][0] for k, i in brackets])
    state = np.array([times[k][i][1] for k, i in brackets], dtype=bool)
    while len(lower) > 0 and np.any(upper - lower > resolution):
        middle = (lower + upper) / 2
        same = observes(footprints(middle)) == state
        lower, upper = np.where(same, middle, lower), np.where(same, upper, middle)
    crossings = dict(zip(brackets, (lower + upper) / 2))
    periods = []
    for k, t in enumerate(times):
        left = t[0][0] if t[0][1] else None
        for i in range(len(t) - 1):
            if t[i][1] == t[i + 1][1]:
                continue
            if t[i + 1][1]:
                left = crossings[(k, i)]
            else:
                periods.append((left, crossings[(k, i)]))
                left = None
        if left is not None:
            periods.append((left, t[-1][0]))
    return [
        pd.Interval(
            left=reference + pd.Timedelta(seconds=float(left)),
            right=reference + pd.Timedelta(seconds=float(right)),
        )
        for left, right in periods
    ]


def _get_swaths(
    region: geo.Polygon | geo.MultiPolygon,
    orbit: GeneralPerturbationsOrbit,
    instrument: Instrument,
    periods: list[pd.Interval],
    elevation: float = 0,
    time_step: timedelta = timedelta(seconds=10),
    min_time_step: timedelta = timedelta(milliseconds=100),
) -> list[geo.Polygon | geo.MultiPolygon]:
    """
    Gets the part of a region swept by an instrument's footprint (see
    `compute_footprint`) during each of a set of periods: the union of
    footprints sampled at most `time_step` apart, where the interval between
    two consecutive footprints that do not intersect is divided until they
    do (or it is shorter than `min_time_step`), clipped to the region.
    Except for a `ConicalInstrument` (whose footprint is an arc), footprints
    are convex, so the convex hull of two consecutive footprints (a single
    polygon spanning less than 180 degrees of longitude) fills the area
    swept between them. The footprints of all periods are computed together.

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region, split along the anti-meridian and poles (see
                `tatc.utils.geometry.split_polygon`).
        orbit (GeneralPerturbationsOrbit): The orbit.
        instrument (Instrument): The instrument.
        periods (list[pandas.Interval]): The periods.
        elevation (float): The elevation (meters) of the region above the
                WGS 84 ellipsoid, at which to project footprints.
        time_step (datetime.timedelta): The maximum time between samples.
        min_time_step (datetime.timedelta): The minimum time between samples.

    Returns:
        list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]: the
            swath of each period
    """
    if len(periods) == 0:
        return []
    reference = periods[0].left

    def footprints(seconds: np.ndarray) -> np.ndarray:
        orbit_track = orbit.get_orbit_track_at_time(
            _to_time_from_offsets(reference, seconds)
        )
        return np.array(
            instrument.compute_footprint(orbit_track, elevation=elevation),
            dtype=object,
        )

    # samples (seconds from the reference) within each period
    seconds = []
    for period in periods:
        duration = (period.right - period.left).total_seconds()
        seconds.append(
            (period.left - reference).total_seconds()
            + np.linspace(
                0,
                duration,
                max(2, int(np.ceil(duration / time_step.total_seconds())) + 1),
            )
        )
    views = np.split(
        footprints(np.concatenate(seconds)), np.cumsum([len(s) for s in seconds])[:-1]
    )
    # divide intervals between consecutive footprints that do not intersect
    while True:
        gaps = [
            np.flatnonzero(
                ~shapely.intersects(v[:-1], v[1:])
                & (np.diff(s) > 2 * min_time_step.total_seconds())
            )
            for s, v in zip(seconds, views)
        ]
        counts = [len(g) for g in gaps]
        if sum(counts) == 0:
            break
        middles = [(s[g] + s[g + 1]) / 2 for s, g in zip(seconds, gaps)]
        views_middle = np.split(
            footprints(np.concatenate(middles)), np.cumsum(counts)[:-1]
        )
        seconds = [np.insert(s, g + 1, m) for s, g, m in zip(seconds, gaps, middles)]
        views = [np.insert(v, g + 1, w) for v, g, w in zip(views, gaps, views_middle)]
    swaths = []
    for period_views in views:
        parts = [period_views]
        if not isinstance(instrument, ConicalInstrument):
            # fill the area swept between consecutive (convex) footprints
            pairs = shapely.union(period_views[:-1], period_views[1:])
            bounds = shapely.bounds(pairs)
            single = (shapely.get_type_id(pairs) == 3) & (
                bounds[:, 2] - bounds[:, 0] < 180
            )
            parts.append(shapely.convex_hull(pairs[single]))
        swath = shapely.intersection(shapely.union_all(np.concatenate(parts)), region)
        # keep the polygonal parts of the swath, without z coordinates
        polygons = [
            shapely.force_2d(part)
            for part in shapely.get_parts(swath)
            if isinstance(part, geo.Polygon) and not part.is_empty
        ]
        swaths.append(polygons[0] if len(polygons) == 1 else geo.MultiPolygon(polygons))
    return swaths


def collect_region_observations(
    region: geo.Polygon | geo.MultiPolygon,
    satellite: Satellite,
    start: datetime,
    end: datetime,
    instrument_index: int = 0,
) -> gpd.GeoDataFrame:
    """
    Collect single satellite observations of a region of interest: a
    shapely `Polygon` or `MultiPolygon` in longitude and latitude (degrees),
    at the elevation of the mean of its z coordinates (meters above the
    WGS 84 ellipsoid), if any, or otherwise zero. The region is split along
    the anti-meridian and poles (see `tatc.utils.geometry.split_polygon`).
    For a point, see `tatc.analysis.point_sampling.collect_observations`.

    Each observation spans a period when the instrument can observe any part
    of the region, from its start to its end (refined to a millisecond):

    - for a `PointedInstrument` or a `ConicalInstrument`, while its
      instantaneous footprint (see `compute_footprint`) intersects the
      region, searched within the periods when its field of regard (which
      must contain its footprint) may observe the region (see
      `compute_region_access_periods`);
    - otherwise, while any part of the region lies within the field of
      regard (a cone about nadir) and above the satellite's horizon.

    An observation is recorded if it lasts at least the instrument's minimum
    access time and the instrument's requirements (such as illumination) are
    met at its midpoint, for a point of the region within the instrument's
    view. The instrument's fixed access time (`access_time_fixed`) does not
    apply to regions.

    The geometry of each observation is its swath: the part of the region
    swept by the instrument's footprint (see `compute_footprint`) during the
    observation. The region is identified by its `target_hash` (see
    `tatc.utils.geometry.hash_geometry`), with which `aggregate_observations`
    and `reduce_observations` group observations of the same region.

    If the orbit is propagated with a repeat cycle (see
    `GeneralPerturbationsOrbit.repeat_cycle`), it is modeled as maintained on
    its repeat ground track before its first and after its last element's
    epoch (see `GeneralPerturbationsOrbit.get_orbit_track_at_time`).

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region of interest.
        satellite (Satellite): The observing satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        instrument_index (int): The index of the observing instrument in satellite.

    Returns:
        geopandas.GeoDataFrame: The data frame of observations, with the
            region's `target_hash`, the swath (`geometry`), the `satellite`
            and `instrument` names, and the `start` and `end` of each
            observation.
    """
    _check_satellite(satellite)
    if not isinstance(region, (geo.Polygon, geo.MultiPolygon)):
        raise TypeError(
            "region must be a Polygon or MultiPolygon, not a "
            f"{type(region).__name__} (see collect_observations for a point)"
        )
    elevation = _get_region_elevation(region)
    target_hash = hash_geometry(region)
    # split along the anti-meridian and poles (which discards z coordinates),
    # restoring the region's elevation, if any
    geometry = split_polygon(region)
    if region.has_z:
        geometry = project_polygon_to_elevation(geometry, elevation)
    instrument = satellite.instruments[instrument_index]
    orbit = satellite.orbit.to_gp_orbit()
    shapely.prepare(geometry)
    windows = list(
        _get_visible_polygon_interval_series(
            geometry, satellite, instrument.field_of_regard, start, end, elevation
        )
    )
    if isinstance(instrument, (PointedInstrument, ConicalInstrument)):
        periods = _find_footprint_periods(
            geometry, orbit, instrument, windows, elevation
        )
    else:
        # refine the windows to the field of regard: when the angle from
        # nadir of the region's point nearest to nadir is at most half the
        # field of regard, and the satellite is above that point's horizon
        arcs = _get_boundary_arcs(geometry, elevation)
        half_angle = instrument.field_of_regard / 2

        def residual(orbit_track: Geocentric) -> np.ndarray:
            angle, sat_elevation, _, _ = _get_region_view(
                geometry, arcs, orbit_track, instrument.nadir_reference, elevation
            )
            return np.maximum(angle - half_angle, -sat_elevation)

        periods = _refine_access_periods(
            residual, orbit, windows, max_step=timedelta(seconds=10)
        )
    periods = [
        period
        for period in periods
        if period.right - period.left >= instrument.min_access_time
    ]
    if len(periods) == 0:
        return _get_empty_region_frame()
    # instrument validity (illumination) is only checked at each period's
    # midpoint, for a point of the region within the instrument's view, as an
    # approximation of the whole period (for all periods at once)
    orbit_track = orbit.get_orbit_track([period.mid for period in periods])
    views = np.array(
        instrument.compute_footprint(orbit_track, elevation=elevation), dtype=object
    )
    observed = shapely.intersection(views, geometry)
    points = shapely.point_on_surface(
        np.where(shapely.is_empty(observed), geometry, observed)
    )
    target = wgs84.latlon(shapely.get_y(points), shapely.get_x(points), elevation)
    valid = np.atleast_1d(instrument.is_valid_observation(orbit_track, target))
    periods = [period for period, is_valid in zip(periods, valid) if is_valid]
    if len(periods) == 0:
        return _get_empty_region_frame()
    swaths = _get_swaths(geometry, orbit, instrument, periods, elevation)
    return gpd.GeoDataFrame(
        [
            {
                "target_hash": target_hash,
                "geometry": (
                    project_polygon_to_elevation(swath, elevation)
                    if region.has_z
                    else swath
                ),
                "satellite": satellite.name,
                "instrument": instrument.name,
                "start": period.left,
                "end": period.right,
            }
            for period, swath in zip(periods, swaths)
        ],
        crs="EPSG:4326",
    )


def collect_multi_region_observations(
    region: geo.Polygon | geo.MultiPolygon,
    satellites: Satellite | list[Satellite],
    start: datetime,
    end: datetime,
) -> gpd.GeoDataFrame:
    """
    Collect multiple satellite observations of a region of interest: calls
    `collect_region_observations` for every instrument on every satellite in
    `satellites`, and concatenates the results into one data frame.

    Args:
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon):
                The region of interest (see `collect_region_observations`).
        satellites (Satellite | list[Satellite]): The observing satellite(s),
                each contributing an observation per instrument it carries.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.

    Returns:
        geopandas.GeoDataFrame: The data frame with all recorded observations.
    """
    gdfs = [
        collect_region_observations(region, satellite, start, end, instrument_index)
        for satellite in _check_satellites(satellites)
        for instrument_index in range(len(satellite.instruments))
    ]
    if len(gdfs) == 0:
        # an empty `satellites` list leaves nothing to concatenate
        return _get_empty_region_frame()
    # concatenate into one data frame, sort by start time, and re-index
    return pd.concat(gdfs).sort_values("start").reset_index(drop=True)
