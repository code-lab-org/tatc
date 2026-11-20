"""
Methods to compute coverage metrics from observations of points and
regions: the aggregation, reduction, and gridding of observations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from ..utils.geometry import _hash_geometries
from .check import (
    _aggregate_names,
    _check_time_window,
    _get_empty_observation_frame,
    _get_platform_columns,
)


def _get_target_keys(gdf: gpd.GeoDataFrame) -> list[pd.Series]:
    """
    Gets the keys that identify the target (point or region) of each
    observation: its `target_hash`, the hash of the point or region (see
    `tatc.utils.geometry.hash_geometry`).

    Args:
        gdf (geopandas.GeoDataFrame): The observations.

    Returns:
        list[pandas.Series]: the keys
    """
    return [gdf["target_hash"]]


def _group_overlapping(
    target: np.ndarray, start: np.ndarray, end: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Groups the overlapping (including touching or nested) observation periods
    of each target: sorts the observations by target and start (stably), and
    starts a group with each target and with each observation that starts
    after every earlier observation of the target has ended.

    Args:
        target (numpy.ndarray): The target code of each observation.
        start (numpy.ndarray): The start of each observation (comparable
            values, such as naive `numpy.datetime64`).
        end (numpy.ndarray): The end of each observation.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: the order of the observations
            (by target and start), and the group (from 0, in that order) of
            each observation in that order
    """
    order = np.lexsort((start, target))
    target, start, end = target[order], start[order], end[order]
    first = np.concatenate([[True], target[1:] != target[:-1]])
    reach = pd.Series(end).groupby(target).cummax().to_numpy()
    group = np.cumsum(first | (start > np.concatenate([reach[:1], reach[:-1]]))) - 1
    return order, group


def _get_naive_utc(times: pd.Series) -> np.ndarray:
    """
    Gets times (timezone-aware) as naive UTC `numpy.datetime64` values.

    Args:
        times (pandas.Series): The times.

    Returns:
        numpy.ndarray: the naive UTC times
    """
    return pd.to_datetime(times, utc=True).dt.tz_localize(None).to_numpy()


def _get_targets(
    observations: pd.DataFrame, targets: gpd.GeoDataFrame | gpd.GeoSeries | None
) -> gpd.GeoDataFrame:
    """
    Gets the targets of a summary of observations: those given (identified
    by the hash of their geometry, see `tatc.utils.geometry.hash_geometry`,
    as observations identify them, with the other columns of a data frame of
    targets, such as a `point_id` or a weight), or those observed.

    Args:
        observations (pandas.DataFrame): The observations, with the
            `target_hash` and geometry of each.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed, or `None` for those
            observed.

    Returns:
        geopandas.GeoDataFrame: the `target_hash` and geometry of each target
    """
    if targets is None:
        first = ~observations["target_hash"].duplicated().to_numpy()
        result = observations.loc[first, ["target_hash", "geometry"]]
        result = gpd.GeoDataFrame(result.reset_index(drop=True), geometry="geometry")
        return result if result.crs is not None else result.set_crs("EPSG:4326")
    geometry = targets if isinstance(targets, gpd.GeoSeries) else targets.geometry
    attributes = (
        {}
        if isinstance(targets, gpd.GeoSeries)
        else {
            column: targets[column].to_numpy()
            for column in targets.columns
            if column not in ("target_hash", targets.geometry.name)
        }
    )
    result = gpd.GeoDataFrame(
        {"target_hash": _hash_geometries(geometry)},
        geometry=geometry.to_numpy(),
        crs=geometry.crs if geometry.crs is not None else "EPSG:4326",
    )
    # the hash and geometry, then the attributes
    return result.assign(**attributes)


def _include_targets(
    reduced: gpd.GeoDataFrame,
    targets: gpd.GeoDataFrame | gpd.GeoSeries,
    counts: list[str],
) -> gpd.GeoDataFrame:
    """
    Gives the reduced observations of every target, in the order of the
    targets: those of an observed target (with its observed geometry), or,
    for a target not observed, its geometry, no samples (counts of zero), and
    missing statistics.

    Args:
        reduced (geopandas.GeoDataFrame): The reduced observations of the
            observed targets.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries): The targets.
        counts (list[str]): The columns of counts, zero for a target not
            observed.

    Returns:
        geopandas.GeoDataFrame: the reduced observations of every target
    """
    result = _get_targets(reduced, targets)
    observed = reduced.drop_duplicates("target_hash").set_index("target_hash")
    index = observed.index.get_indexer(result["target_hash"])
    found = index >= 0
    geometry = result.geometry.to_numpy().copy()
    geometry[found] = observed.geometry.to_numpy()[index[found]]
    result = result.set_geometry(gpd.GeoSeries(geometry, crs=result.crs))
    for column in observed.columns:
        if column == observed.geometry.name:
            continue
        values = observed[column].reindex(result["target_hash"]).to_numpy()
        if column in counts:
            values = np.where(found, values, 0).astype(int)
        result[column] = values
    # the target's attributes after its hash and geometry, then the statistics
    first = [
        column
        for column in reduced.columns
        if column in ("target_hash", reduced.geometry.name)
    ]
    attributes = [column for column in result.columns if column not in reduced.columns]
    statistics = [column for column in reduced.columns if column not in first]
    return result[[*first, *attributes, *statistics]]


def _sum_by_cell(
    points: gpd.GeoDataFrame, cells: gpd.GeoDataFrame, columns: list[str]
) -> pd.DataFrame:
    """
    Sums columns of points over the cells that contain them, for every cell
    (zero for a cell that contains none).

    Args:
        points (geopandas.GeoDataFrame): The points.
        cells (geopandas.GeoDataFrame): The cells, identified by `cell_id`.
        columns (list[str]): The columns to sum.

    Returns:
        pandas.DataFrame: the sums of each cell, in the order of the cells
    """
    joined = cells[["cell_id", cells.geometry.name]].sjoin(
        points[[*columns, points.geometry.name]], how="inner", predicate="contains"
    )
    sums = joined.groupby("cell_id")[columns].sum()
    return sums.reindex(cells["cell_id"].to_numpy()).fillna(0).astype(float)


def aggregate_observations(observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Aggregate constellation observations. Interleaves observations by multiple
    platforms to compute aggregate performance metrics including access
    (observation duration) and revisit (duration between observations).
    Overlapping (including fully nested) observations of the same target
    (the same `target_hash`, see `_get_target_keys`), possibly from different
    platforms/instruments, are merged into a single continuous coverage
    period; the platform columns (`satellite`, or, as for the observations
    of `tatc.analysis.footprints.collect_footprint_observations`, `aircraft`,
    `surface_platform`, or `radar_station`) and `instrument` record every
    contributor to that period, comma-separated (each platform column only
    those of its kind, missing if there are none), and the geometry is the
    union of theirs (for a
    region, the part of it observed during the period). `epoch` is reassigned to
    the midpoint of the merged period's `start`/`end` (a representative
    instant), not the mean of the constituent observations' own epochs.
    Per-observation columns that lose their meaning once merged across
    satellites and over a potentially much longer period -- e.g.
    `platform_alt`, `platform_az`, `sat_sunlit`, `solar_alt`, `solar_az`, `solar_time`,
    `glint_angle` -- are
    intentionally dropped, even if present on `observations`.

    Args:
        observations (geopandas.GeoDataFrame): The collected observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with aggregated observations.
    """
    platform_columns = _get_platform_columns(observations)
    if observations.empty:
        return _get_empty_observation_frame(platform_columns)
    # targets (in sorted order of their keys), with each target's observations
    # sorted by start (stably), for all targets together, grouped by their
    # overlapping periods
    target, _ = pd.MultiIndex.from_arrays(_get_target_keys(observations)).factorize(
        sort=True
    )
    order, group = _group_overlapping(
        target,
        _get_naive_utc(observations["start"]),
        _get_naive_utc(observations["end"]),
    )
    gdf = observations.iloc[order].reset_index(drop=True)
    target = target[order]
    grouped = gdf.groupby(group, sort=False)
    # the geometry of each group: the union of its observations' (for a region,
    # the part of it observed during the period)
    geometry = grouped.geometry.first().to_numpy()
    counts = np.bincount(group)
    members = np.split(gdf.geometry.to_numpy(), np.cumsum(counts)[:-1])
    for g in np.flatnonzero(counts > 1):
        geometry[g] = shapely.union_all(members[g])
    aggregated = gpd.GeoDataFrame(
        {
            "target_hash": grouped["target_hash"].first(),
            **{
                column: _aggregate_names(grouped, column).astype(gdf[column].dtype)
                for column in [*platform_columns, "instrument"]
            },
            "start": grouped["start"].min(),
            "end": grouped["end"].max(),
        },
        geometry=gpd.GeoSeries(
            geometry, index=grouped["start"].min().index, crs=observations.crs
        ),
    )
    aggregated = aggregated[
        ["geometry", "target_hash", *platform_columns, "instrument", "start", "end"]
    ]
    # reassign epoch to the midpoint of the merged period, as a single
    # representative instant, rather than the mean of the constituent
    # observations' own (pre-merge) epochs
    aggregated["epoch"] = (
        aggregated["start"] + (aggregated["end"] - aggregated["start"]) / 2
    )
    # compute access and revisit metrics (the revisit of each target's first
    # period is undefined)
    aggregated["access"] = aggregated["end"] - aggregated["start"]
    group_target = target[np.concatenate([[True], group[1:] != group[:-1]])]
    revisit = aggregated["start"] - aggregated["end"].shift()
    revisit[np.concatenate([[True], group_target[1:] != group_target[:-1]])] = pd.NaT
    aggregated["revisit"] = revisit
    return aggregated.reset_index(drop=True)


def _get_empty_reduce_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for reduced coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "target_hash": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        "access": pd.Series([], dtype="timedelta64[ns]"),
        "revisit": pd.Series([], dtype="timedelta64[ns]"),
        "samples": pd.Series([], dtype="int"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def reduce_observations(
    aggregated_observations: gpd.GeoDataFrame,
    targets: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Reduce constellation observations: for each unique target (the same
    `target_hash`, see `_get_target_keys`) in `aggregated_observations`,
    computes the mean
    access period, the mean revisit period, and the total number of samples
    (aggregated periods) over the analysis period. For a region, the
    geometry is the union of the parts of it observed. The first sample's revisit is undefined (no
    prior observation to measure a gap from) and is excluded from the mean
    rather than counted as zero, which would otherwise bias the mean
    downward; a point with only one sample accordingly has an undefined
    (NaT) mean revisit. The revisit is the time from the end of one
    aggregated period to the start of the next; for the times between
    samples (such as their epochs), their quantiles, and the longest gap,
    see `tatc.analysis.sufficiency.compute_revisit_statistics`.

    With `targets` (such as the points given to `collect_observations`), the
    result has a row for every target, in their order: a target not observed
    has no samples (zero) and undefined (NaT) access and revisit, so that
    unobserved targets appear in maps and statistics rather than being
    omitted. Targets given as a data frame keep their other columns (such as
    a `point_id` or a weight, after the geometry), except those named as the
    statistics.

    Args:
        aggregated_observations (geopandas.GeoDataFrame): The aggregated observations.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed (identified by the hash of
            their geometry, as the observations' `target_hash`), with their
            other columns, or `None` for those observed.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced observations.
    """
    if targets is not None:
        return _include_targets(
            reduce_observations(aggregated_observations), targets, ["samples"]
        )
    if aggregated_observations.empty:
        return _get_empty_reduce_frame()
    # operate on a copy of the data frame
    gdf = aggregated_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"].dt.total_seconds()
    gdf["revisit"] = gdf["revisit"].dt.total_seconds()
    # assign each record to one observation
    gdf["samples"] = 1
    # perform the aggregation operation for each target
    gdf = gdf.dissolve(
        _get_target_keys(gdf),
        aggfunc={
            "access": "mean",
            "revisit": "mean",
            "samples": "sum",
        },
    ).reset_index()
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
    needs to be averaged as a rate to stay a representative statistic. The
    revisit is averaged over the points with a revisit (at least two
    samples), and is undefined (NaT) for a cell with none. Every cell is
    gridded: a cell that contains no observed point has no samples (zero)
    and undefined (NaT) access and revisit.

    Args:
        reduced_observations (geopandas.GeoDataFrame): The reduced observations.
        cells (geopandas.GeoDataFrame): The cell specification, with a
            `cell_id` column.

    Returns:
        geopandas.GeoDataFrame: The cells, with the `samples`, `access`, and
            `revisit` of each.
    """
    # operate on a copy of the data frame
    gdf = reduced_observations.copy()
    # convert access and revisit to numeric values before aggregation
    access = pd.to_timedelta(gdf["access"]).dt.total_seconds()
    revisit = pd.to_timedelta(gdf["revisit"]).dt.total_seconds()
    # pre-transform so the means below reduce to plain sums. The weighted
    # harmonic mean of revisit is sum(samples) / sum(samples/revisit), over
    # the points with a revisit.
    gdf["access_x_samples"] = access * gdf["samples"]
    gdf["revisit_samples"] = gdf["samples"].where(revisit.notna(), 0)
    gdf["samples_over_revisit"] = (gdf["samples"] / revisit).fillna(0)
    sums = _sum_by_cell(
        gdf,
        cells,
        ["samples", "access_x_samples", "revisit_samples", "samples_over_revisit"],
    )
    result = cells.copy()
    samples = sums["samples"].to_numpy()
    result["samples"] = samples.astype(int)
    with np.errstate(divide="ignore", invalid="ignore"):
        # finish the aggregation: sample-weighted arithmetic mean for access,
        # sample-weighted harmonic mean for revisit
        mean_access = np.where(
            samples > 0, sums["access_x_samples"].to_numpy() / samples, np.nan
        )
        rate = sums["samples_over_revisit"].to_numpy()
        mean_revisit = np.where(
            rate > 0, sums["revisit_samples"].to_numpy() / rate, np.nan
        )
    result["access"] = pd.to_timedelta(mean_access, unit="s").as_unit("ns")
    result["revisit"] = pd.to_timedelta(mean_revisit, unit="s").as_unit("ns")
    return result


def count_observations(
    observations: pd.DataFrame,
    start: datetime,
    end: datetime,
    step: timedelta,
    by: str | list[str] | None = None,
) -> pd.DataFrame:
    """
    Count the observations in progress over time: for each step of a period
    (from `start`, every `step`, to `end`), the number of observations whose
    period (from `start` to `end`, inclusive) overlaps it, as the number of
    satellites that observe at once, the points that a constellation
    observes in each minute, or the links available to a satellite (see
    `collect_space_observations`). A short step counts the observations in
    progress at each time. Observations of any kind with `start` and `end`
    columns may be counted.

    Args:
        observations (pandas.DataFrame): The observations, with their
                `start` and `end`.
        start (datetime.datetime): The start of the period.
        end (datetime.datetime): The end of the period.
        step (datetime.timedelta): The duration of each step.
        by (str | list[str] | None): The column(s) whose values group the
                counts, or `None` for one count.

    Returns:
        pandas.DataFrame: The counts, indexed by the start of each step
            (`start`): one column (`observations`), or one for each group
            (named by its value, or its values for several columns).
    """
    _check_time_window(start, end)
    if step <= timedelta(0):
        raise ValueError(f"step ({step}) must be positive")
    count = int(np.ceil((end - start) / step))
    index = pd.DatetimeIndex(
        [pd.Timestamp(start) + k * pd.Timedelta(step) for k in range(count)],
        name="start",
    ).tz_convert("UTC")
    columns = [] if by is None else [by] if isinstance(by, str) else list(by)
    # the steps that each observation overlaps (from the step of its start to
    # that of its end), within the period
    origin = pd.Timestamp(start).tz_convert("UTC")
    first = np.floor(
        (pd.to_datetime(observations["start"], utc=True) - origin) / pd.Timedelta(step)
    ).to_numpy()
    last = np.floor(
        (pd.to_datetime(observations["end"], utc=True) - origin) / pd.Timedelta(step)
    ).to_numpy()
    inside = (last >= 0) & (first < count)
    first = np.clip(first[inside], 0, count - 1).astype(int)
    last = np.clip(last[inside], 0, count - 1).astype(int)
    if columns:
        groups = pd.MultiIndex.from_frame(observations.loc[inside, columns])
        codes, labels = groups.factorize(sort=True)
        labels = labels.get_level_values(0) if len(columns) == 1 else pd.Index(labels)
    else:
        codes = np.zeros(len(first), dtype=int)
        labels = pd.Index(["observations"])
    # count the starts and ends of the overlapped steps, and accumulate them
    changes = np.zeros((count + 1, len(labels)), dtype=int)
    np.add.at(changes, (first, codes), 1)
    np.add.at(changes, (last + 1, codes), -1)
    counts = np.cumsum(changes[:-1], axis=0)
    return pd.DataFrame(counts, index=index, columns=labels)
