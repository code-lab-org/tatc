"""
Ground-based radar utility functions.

Heights in this module may use any vertical reference (e.g. mean sea
level or the WGS 84 ellipsoid), as long as station, target, and terrain
heights all use the same one: the geometry depends only on height
differences.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from numba import njit
from pyproj import Transformer
from shapely import make_valid
from shapely.geometry import GeometryCollection, MultiPolygon, Point, Polygon
from shapely.ops import transform

from .. import constants
from .geometry import geodesic_destination, project_polygon_to_elevation, split_polygon


@njit
def compute_radar_beam_height(
    slant_range: float, elevation_angle: float, station_height: float = 0
) -> float:
    """
    Fast computation of the radar beam height at a
    specified slant range along a beam at a specified elevation angle,
    using the standard-atmosphere "4/3 Earth radius" refraction
    approximation (a spherical Earth with an effective radius scaled by
    `constants.EFFECTIVE_EARTH_RADIUS_FACTOR`).

    Args:
        slant_range (float): Slant range (meters) along the radar beam.
        elevation_angle (float): Radar beam elevation angle (degrees)
            above local horizontal.
        station_height (float): Radar antenna height (meters).

    Returns:
        float: The beam height (meters), in the same vertical reference as
        `station_height`.
    """
    effective_radius = (
        constants.EFFECTIVE_EARTH_RADIUS_FACTOR * constants.EARTH_MEAN_RADIUS
    )
    theta = np.radians(elevation_angle)
    return (
        np.sqrt(
            slant_range**2
            + effective_radius**2
            + 2 * slant_range * effective_radius * np.sin(theta)
        )
        - effective_radius
        + station_height
    )


@njit
def compute_radar_ground_range(
    slant_range: float, elevation_angle: float, station_height: float = 0
) -> float:
    """
    Fast computation of the ground (surface) range to the point reached by
    a specified slant range along a radar beam at a specified elevation
    angle, using the standard-atmosphere "4/3 Earth radius" refraction
    approximation.

    Args:
        slant_range (float): Slant range (meters) along the radar beam.
        elevation_angle (float): Radar beam elevation angle (degrees)
            above local horizontal.
        station_height (float): Radar antenna height (meters).

    Returns:
        float: The ground range (meters), measured along the Earth's
        surface from the station, to the point reached by `slant_range`.
    """
    effective_radius = (
        constants.EFFECTIVE_EARTH_RADIUS_FACTOR * constants.EARTH_MEAN_RADIUS
    )
    theta = np.radians(elevation_angle)
    height = compute_radar_beam_height(slant_range, elevation_angle, station_height)
    sin_arg = slant_range * np.cos(theta) / (effective_radius + height - station_height)
    # guard against floating-point overshoot past 1 (mathematically sin_arg <= 1)
    sin_arg = min(sin_arg, 1.0)
    return effective_radius * np.arcsin(sin_arg)


@njit
def compute_radar_slant_range(
    elevation_angle: float, target_height: float, station_height: float = 0
) -> float:
    """
    Fast computation of the slant range at which a radar beam at a
    specified elevation angle reaches a specified target height, using the
    standard-atmosphere "4/3 Earth radius" refraction approximation. This
    is the inverse of `compute_radar_beam_height`.

    Args:
        elevation_angle (float): Radar beam elevation angle (degrees)
            above local horizontal.
        target_height (float): Target height (meters), in the same vertical
            reference as `station_height`.
        station_height (float): Radar antenna height (meters).

    Returns:
        float: The slant range (meters) at which the beam reaches
        `target_height`, or `numpy.nan` if there is no real, non-negative
        solution. A beam departing at or above local horizontal
        (`elevation_angle >= 0`) only climbs with increasing range under
        this model, so `target_height < station_height` always yields
        `numpy.nan`; `target_height == station_height` yields `0`. A beam
        departing below local horizontal (`elevation_angle < 0`) first
        descends and then climbs, so it can cross a height twice; this
        returns the farther (climbing) crossing.
    """
    effective_radius = (
        constants.EFFECTIVE_EARTH_RADIUS_FACTOR * constants.EARTH_MEAN_RADIUS
    )
    theta = np.radians(elevation_angle)
    discriminant = (target_height - station_height + effective_radius) ** 2 - (
        effective_radius * np.cos(theta)
    ) ** 2
    if discriminant < 0:
        return np.nan
    slant_range = -effective_radius * np.sin(theta) + np.sqrt(discriminant)
    if slant_range < -1e-6:
        # a mathematically real but non-physical (negative range) root,
        # which can arise for target_height < station_height
        return np.nan
    # clamp away floating-point noise (e.g. near elevation_angle = 90) that
    # can otherwise leave an exact-zero solution slightly negative
    return max(slant_range, 0.0)


def compute_radar_ground_range_bounds(
    min_elevation_angle: float,
    max_elevation_angle: float,
    max_range: float,
    target_elevation: float,
    station_elevation: float = 0,
) -> tuple[float, float] | None:
    """
    Computes the ground-range annulus (meters) within which a ground-based
    radar can detect a target at a specified elevation, using the
    standard-atmosphere "4/3 Earth radius" refraction approximation.

    The inner bound is the ground range at which the highest scanned
    elevation angle first reaches `target_elevation` (closer in, no
    scanned angle has climbed that high -- an overhead "cone of
    silence"). The outer bound is the ground range at which the lowest
    scanned elevation angle reaches `target_elevation`, capped by
    `max_range`: if the lowest tilt reaches `target_elevation` within
    `max_range`, the beam has overshot the target height before running
    out of range.

    A target at or below `station_elevation` is observable only if the
    lowest elevation angle is below local horizontal: such a beam first
    descends and then climbs, so it passes below the target between its
    descending and climbing crossings of `target_elevation`, which become
    the inner and outer bounds (higher beams, assumed at or above local
    horizontal, always pass above the target). Otherwise, every beam
    stays above the antenna and the target is not observable.

    Args:
        min_elevation_angle (float): Lowest scanned elevation angle (degrees).
        max_elevation_angle (float): Highest scanned elevation angle (degrees).
        max_range (float): Maximum unambiguous slant range (meters).
        target_elevation (float): Target height (meters), in the same
            vertical reference as `station_elevation`.
        station_elevation (float): Radar antenna height (meters).

    Returns:
        tuple[float, float] | None: The `(inner_ground_range,
        outer_ground_range)` bounds (meters), or `None` if there is no
        ground range at which the target is observable (including the
        degenerate case `min_elevation_angle == max_elevation_angle`).
    """
    if target_elevation <= station_elevation:
        if min_elevation_angle >= 0:
            return None
        outer_slant_range = compute_radar_slant_range(
            min_elevation_angle, target_elevation, station_elevation
        )
        if np.isnan(outer_slant_range):
            # the lowest beam never descends as low as the target
            return None
        # the descending and climbing crossings are symmetric about the
        # beam's lowest point, at a slant range of -R sin(theta)
        effective_radius = (
            constants.EFFECTIVE_EARTH_RADIUS_FACTOR * constants.EARTH_MEAN_RADIUS
        )
        inner_slant_range = max(
            -2 * effective_radius * np.sin(np.radians(min_elevation_angle))
            - outer_slant_range,
            0.0,
        )
        if inner_slant_range >= max_range:
            return None
        inner_ground_range = compute_radar_ground_range(
            inner_slant_range, min_elevation_angle, station_elevation
        )
        outer_ground_range = compute_radar_ground_range(
            min(outer_slant_range, max_range), min_elevation_angle, station_elevation
        )
        if inner_ground_range >= outer_ground_range:
            return None
        return (inner_ground_range, outer_ground_range)
    inner_slant_range = compute_radar_slant_range(
        max_elevation_angle, target_elevation, station_elevation
    )
    if np.isnan(inner_slant_range):
        # not expected for target_elevation > station_elevation; defensive fallback
        inner_ground_range = 0.0
    else:
        inner_ground_range = compute_radar_ground_range(
            inner_slant_range, max_elevation_angle, station_elevation
        )
    outer_slant_range = compute_radar_slant_range(
        min_elevation_angle, target_elevation, station_elevation
    )
    if np.isnan(outer_slant_range):
        # not expected for target_elevation > station_elevation; defensive fallback
        return None
    outer_slant_range = min(outer_slant_range, max_range)
    outer_ground_range = compute_radar_ground_range(
        outer_slant_range, min_elevation_angle, station_elevation
    )
    if inner_ground_range >= outer_ground_range:
        return None
    return (inner_ground_range, outer_ground_range)


@njit
def compute_terrain_elevation_angle(
    ground_distance: float, terrain_elevation: float, station_elevation: float = 0
) -> float:
    """
    Fast computation of the curvature-corrected elevation angle, as seen
    from a ground-based station, to a terrain point at a specified ground
    distance and elevation. Uses the standard-atmosphere "4/3 Earth
    radius" refraction approximation (the same effective-Earth-radius
    model used elsewhere in this module) to account for both the Earth's
    curvature and atmospheric refraction bending the line of sight: this
    is the small-angle "curvature drop" correction standard to radio/radar
    line-of-sight terrain analysis, distinct from (and simpler than) the
    exact beam-height geometry used by `compute_radar_beam_height`, since
    here the target (a terrain point) is given directly by its ground
    distance and elevation, not by a slant range and elevation angle.

    Args:
        ground_distance (float): Ground (surface) distance (meters) to
            the terrain point.
        terrain_elevation (float): Elevation (meters) of the terrain point,
            in the same vertical reference as `station_elevation`.
        station_elevation (float): Elevation (meters) of the station
            (antenna).

    Returns:
        float: The elevation angle (degrees) to the terrain point, as seen
        from the station; positive above local horizontal, negative below.
    """
    effective_radius = (
        constants.EFFECTIVE_EARTH_RADIUS_FACTOR * constants.EARTH_MEAN_RADIUS
    )
    curvature_drop = ground_distance**2 / (2 * effective_radius)
    return np.degrees(
        np.arctan2(
            (terrain_elevation - station_elevation) - curvature_drop, ground_distance
        )
    )


def compute_radar_footprint(
    longitude: float,
    latitude: float,
    inner_ground_range: float,
    outer_ground_range: float,
    elevation: float = 0,
) -> Polygon | MultiPolygon:
    """
    Builds a ground-based radar coverage footprint (a disk, or an annulus
    if `inner_ground_range` is positive) centered on the specified
    longitude/latitude, by reprojecting to a distance-preserving CRS,
    buffering by the ground ranges, and reprojecting back.

    Args:
        longitude (float): Longitude (degrees) of the radar station.
        latitude (float): Latitude (degrees) of the radar station.
        inner_ground_range (float): The inner ground range (meters) of the
            footprint annulus; `0` for a full disk.
        outer_ground_range (float): The outer ground range (meters) of the
            footprint.
        elevation (float): The elevation (meters) at which to project the footprint.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The radar footprint.
    """
    distance_crs = f"+proj=eqc +lat_ts={latitude} +datum=WGS84 +units=m"
    to_crs = Transformer.from_crs("EPSG:4326", distance_crs, always_xy=True)
    from_crs = Transformer.from_crs(distance_crs, "EPSG:4326", always_xy=True)
    center = transform(to_crs.transform, Point(longitude, latitude))  # type: ignore
    outer = center.buffer(outer_ground_range)
    footprint = (
        outer.difference(center.buffer(inner_ground_range))
        if inner_ground_range > 0
        else outer
    )
    return project_polygon_to_elevation(
        split_polygon(transform(from_crs.transform, footprint)),  # type: ignore
        elevation,
    )


def _profile_polygon(
    longitude: float,
    latitude: float,
    azimuths: list[float],
    ground_ranges: Iterable[float],
) -> Polygon | MultiPolygon:
    """
    Builds a polygon through the points at sampled ground ranges along
    sampled azimuths from a center point, repairing any self-intersections
    that can arise from zero-width (e.g. fully blocked) azimuth samples.
    """
    polygon = Polygon(
        [
            geodesic_destination(longitude, latitude, azimuth, ground_range)
            for azimuth, ground_range in zip(azimuths, ground_ranges)
        ]
    )
    if polygon.is_valid:
        return polygon
    polygon = make_valid(polygon)
    if isinstance(polygon, GeometryCollection):
        polygons = [g for g in polygon.geoms if isinstance(g, Polygon)] + [
            p for g in polygon.geoms if isinstance(g, MultiPolygon) for p in g.geoms
        ]
        polygon = polygons[0] if len(polygons) == 1 else MultiPolygon(polygons)
    return polygon


def compute_radar_footprint_profile(
    longitude: float,
    latitude: float,
    azimuths: Iterable[float],
    outer_ground_ranges: Iterable[float],
    inner_ground_range: float | Iterable[float],
    elevation: float = 0,
) -> Polygon | MultiPolygon:
    """
    Builds an azimuthally irregular ground-based radar coverage footprint
    from a sampled outer-boundary profile (e.g. reflecting per-azimuth
    terrain blockage), optionally subtracting an inner hole: either a
    uniform circle (e.g. the overhead "cone of silence" for a target above
    the antenna, governed by the antenna's maximum scan elevation angle
    rather than terrain) or a sampled inner-boundary profile (e.g. for a
    target below the antenna, where the inner bound depends on the lowest
    usable elevation angle and therefore on terrain).

    Args:
        longitude (float): Longitude (degrees) of the radar station.
        latitude (float): Latitude (degrees) of the radar station.
        azimuths (Iterable[float]): Azimuth samples (degrees, clockwise
            from north), in increasing order, spanning one full revolution.
        outer_ground_ranges (Iterable[float]): The outer ground range
            (meters) of the footprint at each corresponding azimuth.
        inner_ground_range (float | Iterable[float]): The inner ground
            range (meters) of the footprint annulus, either uniform (`0`
            for no hole) or at each corresponding azimuth.
        elevation (float): The elevation (meters) at which to project the footprint.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The radar footprint.
    """
    azimuths = list(azimuths)
    footprint = _profile_polygon(longitude, latitude, azimuths, outer_ground_ranges)
    if not isinstance(inner_ground_range, (int, float)):
        inner = _profile_polygon(longitude, latitude, azimuths, inner_ground_range)
        footprint = footprint.difference(inner)
    elif inner_ground_range > 0:
        distance_crs = f"+proj=eqc +lat_ts={latitude} +datum=WGS84 +units=m"
        to_crs = Transformer.from_crs("EPSG:4326", distance_crs, always_xy=True)
        from_crs = Transformer.from_crs(distance_crs, "EPSG:4326", always_xy=True)
        center = transform(to_crs.transform, Point(longitude, latitude))  # type: ignore
        inner_circle = transform(
            from_crs.transform, center.buffer(inner_ground_range)  # type: ignore
        )
        footprint = footprint.difference(inner_circle)
    return project_polygon_to_elevation(split_polygon(footprint), elevation)
