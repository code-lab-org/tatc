"""
Methods to compute coverage metrics from observations of points and
regions: the aggregation, reduction, and gridding of observations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd


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


def _get_empty_aggregate_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for aggregated coverage analysis results.

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
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def aggregate_observations(observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Aggregate constellation observations. Interleaves observations by multiple
    satellites to compute aggregate performance metrics including access
    (observation duration) and revisit (duration between observations).
    Overlapping (including fully nested) observations of the same target
    (the same `target_hash`, see `_get_target_keys`), possibly from different
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
        return _get_empty_aggregate_frame()
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
                "target_hash": "first",
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


def reduce_observations(aggregated_observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
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
    (NaT) mean revisit.

    Args:
        aggregated_observations (geopandas.GeoDataFrame): The aggregated observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced observations.
    """
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
