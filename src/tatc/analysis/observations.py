"""
Methods shared by the point and region coverage analyses: observation data
frames, the refinement of access periods, and the aggregation, reduction,
and gridding of observations.

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

from ..constants import de421, timescale
from ..schemas import GeneralPerturbationsOrbit, Instrument, Satellite


def _get_empty_coverage_frame(omit_solar: bool) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for coverage analysis results.

    Args:
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
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
    point_id: int,
    geometry: geo.base.BaseGeometry,
    satellite: Satellite,
    instrument: Instrument,
    omit_solar: bool,
) -> gpd.GeoDataFrame:
    """
    Builds the data frame of observations of a point or region, with the
    satellite (and solar) angles of each observed point at its epoch.

    Args:
        observations (list[tuple[pandas.Interval, pandas.Timestamp, tuple[float, float, float]]]):
                Each observation's period, epoch, and observed point's longitude
                (degrees), latitude (degrees), and elevation (meters).
        point_id (int): The identifier recorded with each observation.
        geometry (shapely.geometry.base.BaseGeometry): The geometry recorded with
                each observation.
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        omit_solar (bool): `True`, to omit solar angles to improve performance.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations.
    """
    if len(observations) == 0:
        return _get_empty_coverage_frame(omit_solar)
    gdf = gpd.GeoDataFrame(
        [
            {
                "point_id": point_id,
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
    ts = timescale.from_datetimes(gdf.epoch)
    orbit_track = satellite.orbit.to_gp_orbit().get_orbit_track(gdf.epoch.tolist())
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
    lower, upper = np.array(lower, dtype=float), np.array(upper, dtype=float)
    index = np.arange(len(lower))
    f_lower, f_upper = residual(lower, index), residual(upper, index)
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
        f_x[active] = residual(x[active], index[active])
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
    if len(periods) == 0:
        return periods
    reference = periods[0].left

    def evaluate(seconds: np.ndarray, _index: np.ndarray) -> np.ndarray:
        return np.reshape(
            residual(
                orbit.get_orbit_track(
                    [reference + pd.Timedelta(seconds=float(x)) for x in seconds]
                )
            ),
            -1,
        )

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
        evaluate(np.concatenate(samples), np.array([])), np.cumsum(counts)[:-1]
    )
    # brackets of each change of sign between samples
    brackets = [
        (i, j)
        for i, value in enumerate(values)
        for j in np.flatnonzero((value[:-1] <= 0) != (value[1:] <= 0))
    ]
    crossings = {}
    if len(brackets) > 0:
        crossing, _ = _find_crossings(
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


def _get_target_keys(gdf: gpd.GeoDataFrame) -> list[pd.Series]:
    """
    Gets the keys that identify the target (point or region) of each
    observation: for a region, its `target_hash` (see
    `tatc.utils.geometry.hash_geometry`), as the geometry of each observation
    is the part of the region observed; for a point, its `point_id` and its
    geometry (as well-known binary), so that points with distinct geometries
    are kept apart even if they share a `point_id` (for example, shapely
    points with the default identifier of 0).

    Args:
        gdf (geopandas.GeoDataFrame): The observations.

    Returns:
        list[pandas.Series]: the keys
    """
    if "target_hash" in gdf.columns:
        return [gdf["target_hash"]]
    return [gdf["point_id"], gdf.geometry.to_wkb().rename("geometry_key")]


def _get_target_columns(gdf: gpd.GeoDataFrame) -> dict[str, pd.Series]:
    """
    Gets empty columns that identify the target (point or region) of each
    observation, as in `gdf` (see `_get_target_keys`).

    Args:
        gdf (geopandas.GeoDataFrame): The observations.

    Returns:
        dict[str, pandas.Series]: the empty columns
    """
    if "target_hash" in gdf.columns:
        return {"target_hash": pd.Series([], dtype="str")}
    return {"point_id": pd.Series([], dtype="int")}


def _get_empty_aggregate_frame(
    target_columns: dict[str, pd.Series] | None = None,
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for aggregated coverage analysis results.

    Args:
        target_columns (dict[str, pandas.Series]): The columns that identify
                the target (see `_get_target_columns`), by default `point_id`.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        **(target_columns or {"point_id": pd.Series([], dtype="int")}),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def aggregate_observations(observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Aggregate constellation observations. Interleaves observations by multiple
    satellites to compute aggregate performance metrics including access
    (observation duration) and revisit (duration between observations).
    Overlapping (including fully nested) observations of the same target
    (the same `point_id` and geometry of a point, or the same `target_hash`
    of a region, see `_get_target_keys`), possibly from different
    satellites/instruments, are merged into a single continuous coverage
    period; `satellite`/`instrument` record every contributor to that
    period, comma-separated, and the geometry is the union of theirs (for a
    region, the part of it observed during the period). `epoch` is reassigned to
    the midpoint of the merged period's `start`/`end` (a representative
    instant), not the mean of the constituent observations' own epochs.
    Per-observation columns that lose their meaning once merged across
    satellites and over a potentially much longer period -- e.g. `sat_alt`,
    `sat_az`, `sat_sunlit`, `solar_alt`, `solar_az`, `solar_time` -- are
    intentionally dropped, even if present on `observations`.

    Args:
        observations (geopandas.GeoDataFrame): The collected observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with aggregated observations.
    """
    if observations.empty:
        return _get_empty_aggregate_frame(_get_target_columns(observations))
    gdfs = []
    # split into constituent data frames for each target
    for _, gdf in observations.groupby(_get_target_keys(observations)):
        # sort the values by start datetime
        gdf = gdf.sort_values("start")
        # assign the observation group number based on overlapping start/end times
        gdf["obs"] = (gdf["start"] > gdf["end"].shift().cummax()).cumsum()
        # perform the aggregation to group overlapping observations
        gdf = gdf.dissolve(
            "obs",
            aggfunc={
                **{key: "first" for key in _get_target_columns(gdf)},
                "satellite": ", ".join,  # type: ignore
                "instrument": ", ".join,  # type: ignore
                "start": "min",
                "end": "max",
            },
        )
        # reassign epoch to the midpoint of the merged period, as a single
        # representative instant, rather than the mean of the constituent
        # observations' own (pre-merge) epochs
        gdf["epoch"] = gdf["start"] + (gdf["end"] - gdf["start"]) / 2
        # compute access and revisit metrics
        gdf["access"] = gdf["end"] - gdf["start"]
        gdf["revisit"] = gdf["start"] - gdf["end"].shift()
        # append to the list of data frames
        gdfs.append(gdf)
    # return a concatenated data frame and re-index
    return pd.concat(gdfs).reset_index(drop=True)


def _get_empty_reduce_frame(
    target_columns: dict[str, pd.Series] | None = None,
) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for reduced coverage analysis results.

    Args:
        target_columns (dict[str, pandas.Series]): The columns that identify
                the target (see `_get_target_columns`), by default `point_id`.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        **(target_columns or {"point_id": pd.Series([], dtype="int")}),
        "geometry": pd.Series([], dtype="object"),
        "access": pd.Series([], dtype="timedelta64[ns]"),
        "revisit": pd.Series([], dtype="timedelta64[ns]"),
        "samples": pd.Series([], dtype="int"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def reduce_observations(aggregated_observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reduce constellation observations: for each unique target (the
    `point_id` and geometry of a point, or the `target_hash` of a region, see
    `_get_target_keys`) in `aggregated_observations`, computes the mean
    access period, the mean revisit period, and the total number of samples
    (aggregated periods) over the analysis period. For a region, the
    geometry is the union of the parts of it observed. The first sample's revisit is undefined (no
    prior observation to measure a gap from) and is excluded from the mean
    rather than counted as zero, which would otherwise bias the mean
    downward; a point with only one sample accordingly has an undefined
    (NaT) mean revisit.

    Args:
        aggregated_observations (geopandas.GeoDataFrame): The aggregated observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced observations.
    """
    if aggregated_observations.empty:
        return _get_empty_reduce_frame(_get_target_columns(aggregated_observations))
    # operate on a copy of the data frame
    gdf = aggregated_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"].dt.total_seconds()
    gdf["revisit"] = gdf["revisit"].dt.total_seconds()
    # assign each record to one observation
    gdf["samples"] = 1
    # perform the aggregation operation for each target
    gdf = (
        gdf.dissolve(
            _get_target_keys(gdf),
            aggfunc={
                "access": "mean",
                "revisit": "mean",
                "samples": "sum",
            },
        )
        .reset_index()
        .drop(columns="geometry_key", errors="ignore")
    )
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = pd.to_timedelta(gdf["access"], unit="s")
    gdf["revisit"] = pd.to_timedelta(gdf["revisit"], unit="s")
    return gdf


def grid_observations(
    reduced_observations: gpd.GeoDataFrame, cells: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Grid reduced observations to cells: for every cell, sums the number of
    samples across every point it contains, and combines those points'
    access/revisit statistics into a single representative value per cell.
    Both access (a per-event duration) and revisit (a time-between-events
    duration, i.e. the reciprocal of a sampling rate) use a sample-weighted
    mean -- arithmetic for access, harmonic for revisit, since revisit
    needs to be averaged as a rate to stay a representative statistic.

    Args:
        reduced_observations (geopandas.GeoDataFrame): The reduced observations.
        cells (geopandas.GeoDataFrame): The cell specification.

    Returns:
        geopandas.GeoDataFrame: The data frame with gridded observations.
    """
    if reduced_observations.empty:
        gdf = cells.copy()
        gdf["samples"] = 0
        gdf["access"] = None
        gdf["revisit"] = None
        return gdf
    # operate on a copy of the data frame
    gdf = reduced_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"].dt.total_seconds()
    gdf["revisit"] = gdf["revisit"].dt.total_seconds()
    # pre-transform so the means below reduce to plain sums: groupby().agg()
    # with a dict of {column: function} only ever hands a custom callable
    # its own column's Series, never a sibling column like "samples" needed
    # to compute a weighted statistic within the callable. The weighted
    # harmonic mean of revisit is sum(samples) / sum(samples/revisit).
    gdf["access_x_samples"] = gdf["access"] * gdf["samples"]
    gdf["samples_over_revisit"] = gdf["samples"] / gdf["revisit"]
    gdf = (
        cells.sjoin(gdf, how="inner", predicate="contains")
        .dissolve(
            by="cell_id",
            aggfunc={
                "samples": "sum",
                "access_x_samples": "sum",
                "samples_over_revisit": "sum",
            },
        )
        .reset_index()
    )
    # finish the aggregation: sample-weighted arithmetic mean for access,
    # sample-weighted harmonic mean for revisit
    gdf["access"] = gdf["access_x_samples"] / gdf["samples"]
    gdf["revisit"] = gdf["samples"] / gdf["samples_over_revisit"]
    gdf = gdf.drop(columns=["access_x_samples", "samples_over_revisit"])
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = pd.to_timedelta(gdf["access"], unit="s")
    gdf["revisit"] = pd.to_timedelta(gdf["revisit"], unit="s")
    return gdf
