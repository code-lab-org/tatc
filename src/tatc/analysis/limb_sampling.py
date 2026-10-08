"""
Methods to perform limb sounding coverage analysis.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta
from enum import Enum

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Point
from skyfield.api import Distance, wgs84
from skyfield.positionlib import Geocentric

from ..constants import EARTH_MEAN_RADIUS
from ..schemas import Satellite
from ..utils.ellipsoid import (
    _ellipsoidal_tangent_distance,
    _itrs_rotation,
    rectangular_to_geodetic,
)
from ..utils.orbital import compute_vnb_frame
from .check import _check_satellite


class ScanDirection(str, Enum):
    """
    Enumeration of the two broad classes of limb sounder vertical scan:
    sweeping from the lowest to the highest requested tangent point
    elevation (`UPWARD`), or the reverse (`DOWNWARD`).
    """

    UPWARD = "upward"
    DOWNWARD = "downward"


def _default_scan_direction(scan_azimuth: float) -> ScanDirection:
    """
    Chooses a default scan direction from the sensor's viewing azimuth: a
    forward-looking sensor (`scan_azimuth` closer to 0 deg than to 180 deg)
    defaults to an upward scan, and a rearward-looking one (closer to 180
    deg than to 0 deg) defaults to a downward scan. Exactly sideways
    (90 or 270 deg, equidistant from both) defaults to upward. Compared as
    exact angular distances (not trigonometrically), so this tie is exact
    rather than subject to floating-point rounding noise.
    """
    normalized = scan_azimuth % 360.0
    distance_to_forward = min(normalized, 360.0 - normalized)
    distance_to_rearward = abs(normalized - 180.0)
    return (
        ScanDirection.UPWARD
        if distance_to_forward <= distance_to_rearward
        else ScanDirection.DOWNWARD
    )


def _limb_tangent_point(
    sat_pv: Geocentric,
    scan_azimuth: float,
    scan_elevations: list[float],
    tolerance: float = 1e-3,
    max_iterations: int = 10,
    scans: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes the limb sounder's tangent point position and per-sample
    viewing-domain validity, given the satellite's position/velocity at
    each sample (one `scan_elevations` entry per sample/column of `sat_pv`).

    The tangent point is defined with respect to the WGS 84 reference
    ellipsoid -- the point along the look direction `d` with minimum
    geodetic altitude (see `_ellipsoidal_tangent_distance`) -- and the
    viewing (depression) angle is solved iteratively so that this point's
    geodetic altitude matches the requested elevation (to within
    `tolerance` meters). The samples of each scan (given by `scans`, the
    integer index of each sample's scan; by default, all samples are one
    scan) are refined until all of them converge, independently of other
    scans.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: tangent point positions (m,
            shape (3, N), GCRS) and a boolean mask (shape (N,)) that is
            `False` wherever the requested elevation is above the
            satellite's own altitude at that sample (no valid viewing
            angle exists).
    """
    v_u, n_u, b_u = compute_vnb_frame(sat_pv)

    # satellite position and geocentric radius at each sample
    sat_p = np.array(sat_pv.position.m)
    r_sat = np.linalg.norm(sat_p, axis=0)
    target_elevations = np.asarray(scan_elevations, dtype=float)

    # GCRS -> ITRS rotation at each sample, to express rays relative to the
    # Earth-fixed WGS 84 ellipsoid (a single matrix for a scalar time)
    rotation = _itrs_rotation(sat_pv.t)
    sat_p_itrs = np.einsum("ij...,j...->i...", rotation, sat_p)

    # a valid viewing angle exists only for elevations up to the
    # satellite's own geodetic altitude
    in_domain = target_elevations <= rectangular_to_geodetic(sat_p_itrs)[2] + tolerance

    # look direction: rotate from the velocity direction toward the
    # orbit-normal by scan_azimuth (staying in the local horizontal
    # plane), then tilt down toward the quasi-nadir direction by el
    az = np.radians(scan_azimuth)
    horizontal = np.cos(az) * v_u + np.sin(az) * n_u

    # solve for the viewing angle el: a spherical-Earth relationship
    # cos(el) = (R + h) / r_sat, with an effective reference radius R per
    # sample, iteratively corrected by the achieved geodetic altitude error
    # (starting from the mean Earth radius; converges in a few iterations)
    reference_radius = np.full(target_elevations.shape, EARTH_MEAN_RADIUS)
    scans = np.zeros(target_elevations.shape, dtype=int) if scans is None else scans
    for _ in range(max_iterations):
        cos_el = (reference_radius + target_elevations) / r_sat
        el = np.arccos(np.clip(cos_el, -1, 1))
        d = np.cos(el) * horizontal - np.sin(el) * b_u
        d_itrs = np.einsum("ij...,j...->i...", rotation, d)
        s = _ellipsoidal_tangent_distance(sat_p_itrs, d_itrs, target_elevations)
        error = rectangular_to_geodetic(sat_p_itrs + s * d_itrs)[2] - target_elevations
        # scans with a sample (in the viewing domain) not yet converged
        unconverged = np.zeros(np.max(scans, initial=-1) + 1, dtype=bool)
        np.logical_or.at(
            unconverged, scans[in_domain], np.abs(error[in_domain]) > tolerance
        )
        if not np.any(unconverged):
            break
        # refine only those scans (the others are unchanged, so their next
        # look directions and tangent points are the same)
        reference_radius = np.where(
            unconverged[scans], reference_radius - error, reference_radius
        )

    # tangent point: the same distance s along the (inertial) ray
    tp_p = sat_p + s * d
    return tp_p, in_domain


def _constant_rate_scan_fractions(values: np.ndarray) -> np.ndarray:
    """
    Computes the elapsed-time fraction (0 to 1) at which each value in a
    sequence is reached by a process that moves through the sequence at a
    constant rate -- i.e. time proportional to cumulative absolute change
    between consecutive values. Used to time a vertical scan's samples as
    a constant-angular-rate scan mirror (the simplest and most common real
    scan mechanization) would reach them.
    """
    if len(values) <= 1:
        return np.zeros(len(values))
    cumulative = np.concatenate(([0.0], np.cumsum(np.abs(np.diff(values)))))
    total = cumulative[-1]
    if total == 0:
        # degenerate case (e.g. every requested value is identical):
        # no meaningful rate to infer, fall back to even time spacing
        return np.linspace(0, 1, len(values))
    return cumulative / total


def _sample_limb_scan(
    satellite: Satellite,
    start: datetime,
    scan_azimuth: float,
    scan_elevations: list[float],
    scan_duration: timedelta,
) -> list[dict]:
    """
    Samples tangent points across a single vertical scan: one sample per
    requested elevation, timed across `scan_duration` as a constant
    angular-rate scan mirror would reach them (see
    `_constant_rate_scan_fractions`), and computed using that sample's own
    satellite position (so a long scan, with significant along-track
    motion, stays accurate throughout). See `_sample_limb_scans`.
    """
    return _sample_limb_scans(
        satellite, [start], scan_azimuth, scan_elevations, scan_duration
    )[0]


def _sample_limb_scans(
    satellite: Satellite,
    starts: list[datetime],
    scan_azimuth: float,
    scan_elevations: list[float],
    scan_duration: timedelta,
) -> list[list[dict]]:
    """
    Samples tangent points across vertical scans starting at each of
    `starts` (see `_sample_limb_scan`), computing the samples of all scans
    together.
    """
    if len(starts) == 0:
        return []
    orbit = satellite.orbit.to_gp_orbit()
    # reference satellite radius at each scan's start, used only to convert
    # target elevations to reference viewing angles for timing purposes
    # (altitude changes negligibly over one scan's duration; the tangent
    # points actually reported are still computed exactly, per sample,
    # below)
    r_sat0 = [
        np.linalg.norm(position)
        for position in np.reshape(
            orbit.get_orbit_track(list(starts)).position.m, (3, -1)
        ).T
    ]
    target_elevations = np.asarray(scan_elevations, dtype=float)
    num_samples = len(scan_elevations)
    sample_times = []
    for start, radius in zip(starts, r_sat0):
        angle0 = np.arccos(
            np.clip((EARTH_MEAN_RADIUS + target_elevations) / radius, -1, 1)
        )
        fractions = _constant_rate_scan_fractions(angle0)
        sample_times.extend(start + f * scan_duration for f in fractions)

    sat_pv = orbit.get_orbit_track(sample_times)
    tp_p, in_domain = _limb_tangent_point(
        sat_pv,
        scan_azimuth,
        np.tile(target_elevations, len(starts)),
        scans=np.repeat(np.arange(len(starts)), num_samples),
    )

    tpp_geo = wgs84.geographic_position_of(
        Geocentric(Distance(m=tp_p).au, None, sat_pv.t)
    )
    longitude = np.array(tpp_geo.longitude.degrees)
    latitude = np.array(tpp_geo.latitude.degrees)
    elevation = np.array(tpp_geo.elevation.m)

    return [
        [
            {
                "time": sample_times[i],
                "longitude": longitude[i],
                "latitude": latitude[i],
                "elevation": elevation[i],
            }
            for i in range(k * num_samples, (k + 1) * num_samples)
            if in_domain[i]
        ]
        for k in range(len(starts))
    ]


def _interpolate_limb_point(points: list[dict], sample_elevation: float) -> dict:
    """
    Interpolates limb scan attributes at the specified tangent point elevation.

    Args:
        points (list[dict]): the limb scan points (ordered by time).
        sample_elevation (float): the tangent point elevation (m) at which to interpolate.

    Returns:
        dict: interpolated longitude (deg), latitude (deg), elevation (m), and time.
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
        "time": p0["time"] + frac * (p1["time"] - p0["time"]),
    }


def _get_empty_limb_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for limb sounding results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "satellite": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        "position": pd.Series([], dtype="object"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "time": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_limb_observations(
    satellite: Satellite,
    times: list[datetime],
    scan_azimuth: float,
    scan_elevations: list[float],
    scan_duration: timedelta,
    scan_direction: ScanDirection | None = None,
    sample_elevation: float | None = None,
) -> gpd.GeoDataFrame:
    """
    Collects limb sounding observations.

    The tangent point (the geometric basis for every reported position and
    elevation) is the point along the sensor's look direction with minimum
    WGS 84 geodetic altitude, and the viewing angle is solved so that this
    altitude matches each requested elevation. See `_limb_tangent_point`
    for details.

    Args:
        satellite (Satellite): the satellite carrying the limb-sounding instrument.
        times (list[datetime.datetime]): the start time of each vertical scan.
        scan_azimuth (float): the sensor's fixed viewing azimuth (deg),
            relative to the satellite's velocity direction, measured in
            the local horizontal plane (0 = forward along velocity,
            90 = cross-track).
        scan_elevations (list[float]): the target tangent point elevations
            (m) that make up the vertical scan. Order does not matter here
            -- `scan_direction` determines the sweep order -- so these can
            be listed in any order (e.g. a natural ascending pressure- or
            altitude-based grid). Samples are timed across `scan_duration`
            as a constant angular-rate scan mirror would reach them (time
            proportional to each elevation's change in viewing angle, not
            to its position in the list), so unevenly-spaced elevations
            are not visited at evenly-spaced times. Science applications
            that specify levels by atmospheric pressure rather than
            altitude can convert with `tatc.utils.pressure_to_altitude`.
        scan_duration (datetime.timedelta): the total time to sweep
            through `scan_elevations`, starting at each requested time.
        scan_direction (ScanDirection | None): whether the scan sweeps
            from the lowest to the highest requested elevation (`UPWARD`,
            e.g. MLS) or the reverse (`DOWNWARD`) -- the two broad classes
            of real limb sounder vertical scans. Instruments that alternate
            between them (e.g. SABER) can be modeled with one call per
            direction. Defaults to
            `UPWARD` for a forward-looking sensor (`scan_azimuth` closer
            to 0 deg than 180 deg) and `DOWNWARD` for a rearward-looking
            one, when `None` (see `_default_scan_direction`).
        sample_elevation (float | None): the tangent point elevation (m)
            at which to interpolate a single representative point/time for
            each scan. Defaults to the midpoint of `scan_elevations`
            (`(min + max) / 2`) when `None`.
    """
    _check_satellite(satellite)
    if scan_direction is None:
        scan_direction = _default_scan_direction(scan_azimuth)
    scan_elevations = sorted(
        scan_elevations, reverse=(scan_direction == ScanDirection.DOWNWARD)
    )
    if sample_elevation is None:
        sample_elevation = (min(scan_elevations) + max(scan_elevations)) / 2

    # sample all scans together
    scans = [
        points
        for points in _sample_limb_scans(
            satellite, times, scan_azimuth, scan_elevations, scan_duration
        )
        if len(points) > 0
    ]
    if len(scans) == 0:
        return _get_empty_limb_frame()
    # the sampled points of each scan, built at once
    geometry = shapely.multipoints(
        [
            [point["longitude"], point["latitude"], point["elevation"]]
            for scan in scans
            for point in scan
        ],
        indices=np.repeat(np.arange(len(scans)), [len(scan) for scan in scans]),
    )
    # format results
    return gpd.GeoDataFrame(
        [
            {
                "satellite": satellite.name,
                "geometry": points,
                "position": Point(
                    sample["longitude"], sample["latitude"], sample["elevation"]
                ),
                "start": scan[0]["time"],
                "end": scan[-1]["time"],
                "time": sample["time"],
            }
            for scan, points in zip(scans, geometry)
            for sample in [_interpolate_limb_point(scan, sample_elevation)]
        ],
        crs="EPSG:4326",
    ).sort_values("time", ignore_index=True)
