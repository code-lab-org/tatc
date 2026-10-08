"""
Methods shared by the point and region sampling analyses (see
`point_sampling` and `region_sampling`): observation data frames and the
refinement of access periods.

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
from ..utils.propagation import TimeRequest, _run, _to_time_from_offsets, _value


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
    residual: Callable[[np.ndarray, np.ndarray], np.ndarray],
    lower: np.ndarray,
    upper: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Find a zero of a residual function within each of a set of intervals by
    the Illinois variant of the regula falsi method, vectorized across the
    intervals.

    Args:
        residual (Callable[[numpy.ndarray, numpy.ndarray], numpy.ndarray]):
                The residual function, evaluated at an array of times
                (seconds) within the intervals with the given indices.
        lower (numpy.ndarray): The lower ends of the intervals (seconds).
        upper (numpy.ndarray): The upper ends of the intervals (seconds).

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: The zero in each interval (or,
            if the residual does not change sign, the end with the smaller
            residual) and whether the interval brackets a zero.
    """
    return _run(
        _find_crossings_steps(
            lambda seconds, index: _value(residual(seconds, index)), lower, upper
        )
    )


def _find_crossings_steps(
    residual: Callable[[np.ndarray, np.ndarray], TimeRequest],
    lower: np.ndarray,
    upper: np.ndarray,
) -> TimeRequest:
    """
    Find a zero of a residual function within each of a set of intervals, as
    a computation (see `_find_crossings` and
    `tatc.utils.propagation.TimeRequest`) of a residual function that is
    itself a computation.
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
) -> list[pd.Interval]:
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
    kept.

    Args:
        residual (Callable[[skyfield.positionlib.Geocentric], numpy.ndarray]):
                The residual function of an orbit track (at one or more times).
        orbit (GeneralPerturbationsOrbit): The orbit.
        periods (list[pandas.Interval]): The visible periods.
        max_step (datetime.timedelta | None): The maximum time between samples.

    Returns:
        list[pandas.Interval]: The refined periods.
    """
    return _run(_refine_access_periods_steps(residual, orbit, periods, max_step))


def _refine_access_periods_steps(
    residual: Callable[[Geocentric], np.ndarray],
    orbit: GeneralPerturbationsOrbit,
    periods: list[pd.Interval],
    max_step: timedelta | None = None,
) -> TimeRequest:
    """
    Refine visible periods to the times when a residual function of the
    orbit track is not positive, as a computation (see
    `_refine_access_periods` and `tatc.utils.propagation.TimeRequest`).
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
        crossing, _ = yield from _find_crossings_steps(
            evaluate,
            np.array([samples[i][j] for i, j in brackets]),
            np.array([samples[i][j + 1] for i, j in brackets]),
        )
        crossings = dict(zip(brackets, crossing))
    refined = []
    for i, (sample, value) in enumerate(zip(samples, values)):
        inside = value <= 0
        left = sample[0] if inside[0] else None
        for j in range(len(sample) - 1):
            if inside[j] == inside[j + 1]:
                continue
            if inside[j + 1]:
                left = crossings[(i, j)]
            else:
                refined.append((left, crossings[(i, j)]))
                left = None
        if left is not None:
            refined.append((left, sample[-1]))
    return [
        pd.Interval(
            left=reference + pd.Timedelta(seconds=float(left)),
            right=reference + pd.Timedelta(seconds=float(right)),
        )
        for left, right in refined
    ]
