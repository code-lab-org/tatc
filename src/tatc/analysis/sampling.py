"""
Methods shared by the sampling analyses: observation data frames and the
refinement of access periods (see `point_sampling`, `region_sampling`, and
`space_sampling`), and the interpolation of tangent point profiles (see `ro_sampling` and
`limb_sampling`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable
from datetime import timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import geometry as geo
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric

from ..constants import de421
from ..schemas import GeneralPerturbationsOrbit, Instrument, Satellite
from ..utils.computation import TimeRequest
from ..utils.time import _to_time_from_offsets


def _get_empty_coverage_frame(omit_solar: bool) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for coverage analysis results.

    Args:
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "target_hash": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "sat_alt": pd.Series(dtype="float"),
        "sat_az": pd.Series(dtype="float"),
    }
    if not omit_solar:
        columns = {
            **columns,
            "sat_sunlit": pd.Series(dtype="bool"),
            "solar_alt": pd.Series(dtype="float"),
            "solar_az": pd.Series(dtype="float"),
            "solar_time": pd.Series(dtype="float"),
        }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _build_observation_frame(
    observations: list[tuple[pd.Interval, pd.Timestamp, tuple[float, float, float]]],
    target_hash: str,
    geometry: geo.base.BaseGeometry,
    satellite: Satellite,
    instrument: Instrument,
    omit_solar: bool,
    orbit_track: Geocentric | None = None,
) -> gpd.GeoDataFrame:
    """
    Builds the data frame of observations of a point or region, with the
    satellite (and solar) angles of each observed point at its epoch.

    Args:
        observations (list[tuple[pandas.Interval, pandas.Timestamp, tuple[float, float, float]]]):
                Each observation's period, epoch, and observed point's longitude
                (degrees), latitude (degrees), and elevation (meters).
        target_hash (str): The hash of the target recorded with each
                observation (see `tatc.utils.geometry.hash_geometry`).
        geometry (shapely.geometry.base.BaseGeometry): The geometry recorded with
                each observation.
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        omit_solar (bool): `True`, to omit solar angles to improve performance.
        orbit_track (skyfield.positionlib.Geocentric | None): The satellite's
                orbit track at each observation's epoch, if already computed.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations.
    """
    if len(observations) == 0:
        return _get_empty_coverage_frame(omit_solar)
    gdf = gpd.GeoDataFrame(
        [
            {
                "target_hash": target_hash,
                "geometry": geometry,
                "satellite": satellite.name,
                "instrument": instrument.name,
                "start": (
                    period.left
                    if not instrument.access_time_fixed
                    else epoch - instrument.min_access_time / 2
                ),
                "end": (
                    period.right
                    if not instrument.access_time_fixed
                    else epoch + instrument.min_access_time / 2
                ),
                "epoch": epoch,
            }
            for period, epoch, _ in observations
        ],
        crs="EPSG:4326",
    )
    # observed point of each observation
    longitude, latitude, elevation = np.array(
        [target for _, _, target in observations]
    ).T
    topos = wgs84.latlon(latitude, longitude, elevation)
    if orbit_track is None:
        orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(gdf.epoch.tolist())
    ts = orbit_track.t
    # append satellite altitude/azimuth columns
    sat_altaz = (orbit_track - topos.at(ts)).altaz()
    gdf["sat_alt"] = sat_altaz[0].degrees  # type: ignore
    gdf["sat_az"] = sat_altaz[1].degrees  # type: ignore
    if not omit_solar:
        # append satellite sunlit column
        gdf["sat_sunlit"] = orbit_track.is_sunlit(de421)
        # append solar altitude/azimuth columns
        sun_altaz = (
            (de421["earth"] + topos).at(ts).observe(de421["sun"]).apparent().altaz()
        )
        gdf["solar_alt"] = sun_altaz[0].degrees
        gdf["solar_az"] = sun_altaz[1].degrees
        # append local solar time column
        gdf["solar_time"] = (de421["earth"] + topos).at(ts).observe(
            de421["sun"]
        ).apparent().hadec()[0].hours + 12
    return gdf


def _find_crossings(
    residual: Callable[[np.ndarray, np.ndarray], TimeRequest],
    lower: np.ndarray,
    upper: np.ndarray,
) -> TimeRequest:
    """
    Find a zero of a residual function within each of a set of intervals by
    the Illinois variant of the regula falsi method, vectorized across the
    intervals, as a computation (see `tatc.utils.computation.TimeRequest`).

    Args:
        residual (Callable[[numpy.ndarray, numpy.ndarray], TimeRequest]):
                The computation of the residual function, evaluated at an
                array of times (seconds) within the intervals with the given
                indices.
        lower (numpy.ndarray): The lower ends of the intervals (seconds).
        upper (numpy.ndarray): The upper ends of the intervals (seconds).

    Returns:
        TimeRequest: the computation of the zero in each interval (or, if
            the residual does not change sign, the end with the smaller
            residual) and whether the interval brackets a zero
            (tuple[numpy.ndarray, numpy.ndarray]).
    """
    lower, upper = np.array(lower, dtype=float), np.array(upper, dtype=float)
    index = np.arange(len(lower))
    f_lower = yield from residual(lower, index)
    f_upper = yield from residual(upper, index)
    # without a sign change, use the end with the smaller residual
    crossing = np.where(np.abs(f_lower) <= np.abs(f_upper), lower, upper)
    bracketed = np.sign(f_lower) * np.sign(f_upper) < 0
    # side of the bracket replaced in the previous iteration (-1: lower, 1: upper)
    side = np.zeros(len(lower))
    for _ in range(50):
        active = bracketed & (upper - lower > 1e-3)
        if not np.any(active):
            break
        with np.errstate(divide="ignore", invalid="ignore"):
            x = np.where(
                active, (lower * f_upper - upper * f_lower) / (f_upper - f_lower), lower
            )
        f_x = np.zeros(len(lower))
        f_x[active] = yield from residual(x[active], index[active])
        replace_lower = active & (np.sign(f_x) == np.sign(f_lower))
        replace_upper = active & ~replace_lower
        # Illinois modification: halve the residual of an end retained twice
        f_upper = np.where(replace_lower & (side == -1), f_upper / 2, f_upper)
        f_lower = np.where(replace_upper & (side == 1), f_lower / 2, f_lower)
        lower = np.where(replace_lower, x, lower)
        f_lower = np.where(replace_lower, f_x, f_lower)
        upper = np.where(replace_upper, x, upper)
        f_upper = np.where(replace_upper, f_x, f_upper)
        side = np.where(replace_lower, -1, np.where(replace_upper, 1, side))
        crossing = np.where(active, x, crossing)
    return crossing, bracketed


def _refine_access_periods(
    residual: Callable[[Geocentric], np.ndarray],
    orbit: GeneralPerturbationsOrbit,
    periods: list[pd.Interval],
    max_step: timedelta | None = None,
) -> TimeRequest:
    """
    Refine visible periods to the times when a residual function of the
    orbit track is not positive: for example, a target's angle from nadir
    less half an instrument's field of regard. The visible periods, from a
    conservative condition, bracket these times. Each period is sampled at
    21 times (or more, if needed to sample at least every `max_step`), and
    each change of sign of the residual between samples is refined, so that
    a period may be divided into several parts (for example, as a
    satellite passes over separate parts of a region). Periods in which the
    residual is positive at every sample are removed; period ends at which
    it is not positive (for example, at the ends of the analysis period) are
    kept. As a computation (see `tatc.utils.computation.TimeRequest`).

    Args:
        residual (Callable[[skyfield.positionlib.Geocentric], numpy.ndarray]):
                The residual function of an orbit track (at one or more times).
        orbit (GeneralPerturbationsOrbit): The orbit.
        periods (list[pandas.Interval]): The visible periods.
        max_step (datetime.timedelta | None): The maximum time between samples.

    Returns:
        TimeRequest: the computation of the refined periods
            (list[pandas.Interval]).
    """
    if len(periods) == 0:
        return periods
    reference = periods[0].left

    def evaluate(seconds: np.ndarray, _index: np.ndarray) -> TimeRequest:
        t = _to_time_from_offsets(reference, seconds)
        yield t
        return np.reshape(residual(orbit.get_orbit_track_at_time(t)), -1)

    lower = np.array([(period.left - reference).total_seconds() for period in periods])
    upper = np.array([(period.right - reference).total_seconds() for period in periods])
    counts = np.full(len(periods), 21)
    if max_step is not None:
        counts = np.maximum(
            counts, np.ceil((upper - lower) / max_step.total_seconds()).astype(int) + 1
        )
    samples = [
        np.linspace(lo, hi, count) for lo, hi, count in zip(lower, upper, counts)
    ]
    values = np.split(
        (yield from evaluate(np.concatenate(samples), np.array([]))),
        np.cumsum(counts)[:-1],
    )
    # brackets of each change of sign between samples
    brackets = [
        (i, j)
        for i, value in enumerate(values)
        for j in np.flatnonzero((value[:-1] <= 0) != (value[1:] <= 0))
    ]
    crossings = {}
    if len(brackets) > 0:
        crossing, _ = yield from _find_crossings(
            evaluate,
            np.array([samples[i][j] for i, j in brackets]),
            np.array([samples[i][j + 1] for i, j in brackets]),
        )
        crossings = dict(zip(brackets, crossing))
    return [
        pd.Interval(
            left=reference + pd.Timedelta(seconds=float(left)),
            right=reference + pd.Timedelta(seconds=float(right)),
        )
        for parts in _assemble_periods(samples, values, crossings)
        for left, right in parts
    ]


def _assemble_periods(
    samples: list[np.ndarray],
    values: list[np.ndarray],
    crossings: dict[tuple[int, int], float],
) -> list[list[tuple[float, float]]]:
    """
    Assembles the periods when a residual function is not positive from its
    values at several sets of samples (for example, of the periods or pairs
    of satellites analyzed) and the refined crossings at each of its changes
    of sign between samples. A period starts at the first sample (or ends at
    the last sample) of a set if the residual is not positive there.

    Args:
        samples (list[numpy.ndarray]): The times (seconds) of each set of samples.
        values (list[numpy.ndarray]): The residual at each set of samples.
        crossings (dict[tuple[int, int], float]): The time (seconds) of the
            crossing between samples `j` and `j + 1` of set `i`, by `(i, j)`,
            for each change of sign.

    Returns:
        list[list[tuple[float, float]]]: the start and end (seconds) of the
            periods of each set of samples
    """
    periods = []
    for i, (sample, value) in enumerate(zip(samples, values)):
        parts = []
        inside = value <= 0
        left = sample[0] if inside[0] else None
        for j in np.flatnonzero(inside[:-1] != inside[1:]):
            if inside[j + 1]:
                left = crossings[(i, j)]
            else:
                parts.append((left, crossings[(i, j)]))
                left = None
        if left is not None:
            parts.append((left, sample[-1]))
        periods.append(parts)
    return periods


def _interpolate_profile_point(
    points: list[dict],
    sample_elevation: float,
    angles: dict[str, float] | None = None,
) -> dict:
    """
    Interpolates the attributes of a tangent point profile (as of a radio
    occultation or a limb scan) at the specified tangent point elevation:
    linearly between the points that bracket the first crossing of the
    elevation, or at the endpoint nearest to it if it is not crossed.

    Args:
        points (list[dict]): the profile points (ordered by time), with
            `longitude` (deg), `latitude` (deg), `elevation` (m), and `time`,
            and any `angles`.
        sample_elevation (float): the tangent point elevation (m) at which
            to interpolate.
        angles (dict[str, float] | None): the other attributes to
            interpolate, as angles (deg) along the shortest angular path,
            each with the lower end of the range to which it is wrapped
            (e.g. -180 or 0).

    Returns:
        dict: interpolated longitude (deg), latitude (deg), elevation (m),
            any `angles` (deg), and time.
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
        **{
            name: lerp_angle(p0[name], p1[name], low=low)
            for name, low in (angles or {}).items()
        },
        "time": p0["time"] + frac * (p1["time"] - p0["time"]),
    }
