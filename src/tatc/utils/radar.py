"""
Ground-based radar utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
from numba import njit

from .. import constants


@njit
def compute_radar_beam_height(
    slant_range: float, elevation_angle: float, station_height: float = 0
) -> float:
    """
    Fast computation of the radar beam height above the WGS 84 datum at a
    specified slant range along a beam at a specified elevation angle,
    using the standard-atmosphere "4/3 Earth radius" refraction
    approximation (a spherical Earth with an effective radius scaled by
    `constants.EFFECTIVE_EARTH_RADIUS_FACTOR`).

    Args:
        slant_range (float): Slant range (meters) along the radar beam.
        elevation_angle (float): Radar beam elevation angle (degrees)
            above local horizontal.
        station_height (float): Radar antenna height (meters) above the
            WGS 84 datum.

    Returns:
        float: The beam height (meters) above the WGS 84 datum.
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
        station_height (float): Radar antenna height (meters) above the
            WGS 84 datum.

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
        target_height (float): Target height (meters) above the WGS 84 datum.
        station_height (float): Radar antenna height (meters) above the
            WGS 84 datum.

    Returns:
        float: The slant range (meters) at which the beam reaches
        `target_height`, or `numpy.nan` if there is no real, non-negative
        solution. A beam departing at or above local horizontal
        (`elevation_angle >= 0`) only climbs with increasing range under
        this model, so `target_height < station_height` always yields
        `numpy.nan`; `target_height == station_height` yields `0`.
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

    A target at or below `station_elevation` is treated as a full disk out
    to `max_range` (zero inner bound), since the elevation-angle geometry
    used otherwise is not meaningful for a target at/below the antenna.

    Args:
        min_elevation_angle (float): Lowest scanned elevation angle (degrees).
        max_elevation_angle (float): Highest scanned elevation angle (degrees).
        max_range (float): Maximum unambiguous slant range (meters).
        target_elevation (float): Target height (meters) above the WGS 84 datum.
        station_elevation (float): Radar antenna height (meters) above the
            WGS 84 datum.

    Returns:
        tuple[float, float] | None: The `(inner_ground_range,
        outer_ground_range)` bounds (meters), or `None` if there is no
        ground range at which the target is observable (including the
        degenerate case `min_elevation_angle == max_elevation_angle`).
    """
    if target_elevation <= station_elevation:
        return (
            0.0,
            compute_radar_ground_range(
                max_range, min_elevation_angle, station_elevation
            ),
        )
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
