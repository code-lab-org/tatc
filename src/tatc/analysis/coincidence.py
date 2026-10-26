"""
Methods to find coincident observations by two or more observing
platforms: the periods when platforms observe the same point and the
areas that their ground tracks observe in common.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta
from itertools import combinations

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely

from .coverage_metrics import _get_target_keys


def _check_max_time_difference(max_time_difference: timedelta) -> pd.Timedelta:
    """
    Checks the maximum time difference of coincident observations.

    Args:
        max_time_difference (datetime.timedelta): The maximum time difference.

    Returns:
        pandas.Timedelta: the maximum time difference
    """
    value = pd.Timedelta(max_time_difference)
    if value < pd.Timedelta(0):
        raise ValueError(
            f"max_time_difference must be non-negative, not {max_time_difference}"
        )
    return value


def _get_empty_coincidence_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for coincident observations.

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


def _find_coincident_periods(
    target: np.ndarray,
    platform: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    min_platforms: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Finds the maximal periods when at least a minimum number of platforms
    observe each target at once, sweeping over the start and end of each
    platform's (closed) observation periods, for all targets together.

    Args:
        target (numpy.ndarray): The target code of each observation period.
        platform (numpy.ndarray): The platform code of each observation period.
        start (numpy.ndarray): The start of each observation period
                (`numpy.datetime64`).
        end (numpy.ndarray): The end of each observation period
                (`numpy.datetime64`).
        min_platforms (int): The minimum number of platforms.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: the target code,
            start, and end of each coincident period, in order of target and
            time
    """
    # merge the overlapping periods of each platform at each target so that
    # it counts once
    order = np.lexsort((start, platform, target))
    target, platform = target[order], platform[order]
    start, end = start[order], end[order]
    first = np.r_[True, (target[1:] != target[:-1]) | (platform[1:] != platform[:-1])]
    reach = pd.Series(end).groupby(np.cumsum(first)).cummax().to_numpy()
    merged = np.cumsum(first | (start > np.r_[reach[:1], reach[:-1]]))
    merged_target = target[np.r_[True, merged[1:] != merged[:-1]]]
    merged_start = pd.Series(start).groupby(merged).min().to_numpy()
    merged_end = pd.Series(end).groupby(merged).max().to_numpy()
    # sweep the events of each target in time order, with starts before ends
    # at the same time so that periods which touch count as coincident at
    # that instant
    event_target = np.concatenate((merged_target, merged_target))
    event_time = np.concatenate((merged_start, merged_end))
    change = np.repeat([1, -1], len(merged_target))
    order = np.lexsort((-change, event_time, event_target))
    event_target, event_time, change = (
        event_target[order],
        event_time[order],
        change[order],
    )
    active = (
        pd.Series(change).groupby(event_target).cumsum().to_numpy() >= min_platforms
    )
    was_active = (
        np.r_[False, active[:-1]] & np.r_[False, event_target[1:] == event_target[:-1]]
    )
    return (
        event_target[active & ~was_active],
        event_time[active & ~was_active],
        event_time[~active & was_active],
    )


def find_coincident_observations(
    observations: gpd.GeoDataFrame,
    min_platforms: int = 2,
    max_time_difference: timedelta = timedelta(0),
    instantaneous: bool = False,
) -> gpd.GeoDataFrame:
    """
    Find coincident observations of points: the periods when at least
    `min_platforms` satellites observe the same target (the same
    `target_hash`, see `_get_target_keys`) at once, as collected by
    `tatc.analysis.point_sampling.collect_observations`. Observations by
    different instruments on the same satellite are not coincident with
    each other.

    Each observation period is widened by half of `max_time_difference`
    before and after, so that two observations coincide if the gap between
    their periods is at most `max_time_difference` (for instantaneous
    observations, if their epochs differ by at most `max_time_difference`).
    Each coincident period is maximal: it spans the time that at least
    `min_platforms` widened periods overlap, even as the satellites
    observing change; `satellite`/`instrument` record every contributor to
    that period, comma-separated and paired, in order of their start time.
    `epoch` is the midpoint of the coincident period.

    With `instantaneous`, each observation is reduced to its `epoch` (for a
    scanning instrument, such as a `ConicalInstrument` or a
    `PointedInstrument`, the time it samples the point, while its period is
    the time the point is within its field of regard), so that observations
    coincide if their epochs differ by at most `max_time_difference`. Each
    coincident observation then spans the epochs of the observations that
    contribute to it, from the first to the last (for two observations,
    `end` less `start` is their time difference).

    The result has the columns of `aggregate_observations` inputs, so
    `aggregate_observations` and `reduce_observations` compute the access
    and revisit of coincident observations.

    Args:
        observations (geopandas.GeoDataFrame): The collected observations of
                points.
        min_platforms (int): The minimum number of satellites (at least 2)
                observing at once.
        max_time_difference (datetime.timedelta): The maximum time
                difference between coincident observations.
        instantaneous (bool): True, to compare the observations' epochs
                rather than their periods.

    Returns:
        geopandas.GeoDataFrame: The data frame of coincident observations,
            sorted by start time.
    """
    if min_platforms < 2:
        raise ValueError(f"min_platforms must be at least 2, not {min_platforms}")
    half_difference = _check_max_time_difference(max_time_difference) / 2
    if observations.empty:
        return _get_empty_coincidence_frame()
    if not (observations.geom_type == "Point").all():
        raise ValueError(
            "observations must be of points (see find_coincident_ground_track "
            "for the footprints of ground tracks)"
        )
    gdf = observations.reset_index(drop=True)
    # targets (in sorted order of their keys) and platforms of the observations
    target, _ = pd.MultiIndex.from_arrays(_get_target_keys(gdf)).factorize(sort=True)
    platform, _ = pd.factorize(gdf["satellite"])

    def naive_utc(times: pd.Series) -> np.ndarray:
        return times.dt.tz_convert("UTC").dt.tz_localize(None).to_numpy()

    # observation times: each period, or each epoch if instantaneous
    observed_start = gdf["epoch"] if instantaneous else gdf["start"]
    observed_end = gdf["epoch"] if instantaneous else gdf["end"]
    start = naive_utc(observed_start - half_difference)
    end = naive_utc(observed_end + half_difference)
    period_target, period_start, period_end = _find_coincident_periods(
        target, platform, start, end, min_platforms
    )
    if len(period_target) == 0:
        return _get_empty_coincidence_frame()
    # each observation overlapping a period contributes to it, in order of
    # its start time
    candidates = pd.DataFrame(
        {
            "target": target,
            "widened_start": start,
            "widened_end": end,
            "start": naive_utc(observed_start),
            "satellite": gdf["satellite"].to_numpy(),
            "instrument": gdf["instrument"].to_numpy(),
        }
    )
    contributors = pd.DataFrame(
        {
            "period": np.arange(len(period_target)),
            "target": period_target,
            "period_start": period_start,
            "period_end": period_end,
        }
    ).merge(candidates, on="target")
    contributors = contributors[
        (contributors["widened_start"] <= contributors["period_end"])
        & (contributors["widened_end"] >= contributors["period_start"])
    ]
    if instantaneous:
        # each coincident observation spans the epochs of its contributors
        period_start = contributors.groupby("period")["start"].min().to_numpy()
        period_end = contributors.groupby("period")["start"].max().to_numpy()
    contributors = (
        contributors.sort_values(["period", "start"], kind="stable")
        .drop_duplicates(["period", "satellite", "instrument"])
        .groupby("period")
        .agg(satellite=("satellite", ", ".join), instrument=("instrument", ", ".join))
    )
    # each target's hash and geometry, from its first observation
    first = gdf.groupby(target).head(1).set_index(pd.Index(target).unique())
    start = pd.Series(period_start).dt.tz_localize("UTC")
    end = pd.Series(period_end).dt.tz_localize("UTC")
    return (
        gpd.GeoDataFrame(
            {
                "target_hash": first["target_hash"].loc[period_target].to_numpy(),
                "geometry": first.geometry.loc[period_target].to_numpy(),
                "satellite": contributors["satellite"].to_numpy(),
                "instrument": contributors["instrument"].to_numpy(),
                "start": start,
                "epoch": start + (end - start) / 2,
                "end": end,
            },
            crs=observations.crs,
        )
        .sort_values("start", kind="stable")
        .reset_index(drop=True)
    )


def find_coincident_ground_track(
    ground_track: gpd.GeoDataFrame,
    max_time_difference: timedelta = timedelta(0),
) -> gpd.GeoDataFrame:
    """
    Find coincident observations of ground tracks: the pairs of footprints
    of different satellites, as collected by
    `tatc.analysis.ground_track.collect_ground_track`, whose times differ
    by at most `max_time_difference` and which overlap (share a positive
    area, in longitude/latitude coordinates). Footprints of different
    instruments on the same satellite are not coincident with each other.

    Each pair records the columns of both footprints, suffixed `_1` and
    `_2` (for example, `time_1`, `satellite_1`, and `valid_obs_1`) and
    ordered so that `time_1` is no later than `time_2`, the time difference
    `time_difference` (`time_2` less `time_1`), and the overlap of the
    footprints as its geometry. Dissolving the result (for example, after
    selecting the pairs of valid observations) gives the area observed in
    common.

    Footprints are first matched by overlap, then by time, so the cost
    grows with the number of times each area is observed.

    Args:
        ground_track (geopandas.GeoDataFrame): The collected ground track.
        max_time_difference (datetime.timedelta): The maximum time
                difference between coincident footprints.

    Returns:
        geopandas.GeoDataFrame: The data frame of coincident footprints,
            sorted by `time_1` and `time_2`.
    """
    max_time_difference = _check_max_time_difference(max_time_difference)
    track = ground_track.sort_values("time", kind="stable").reset_index(drop=True)
    geometries = track.geometry.to_numpy()
    satellites = track["satellite"].to_numpy()
    times = pd.DatetimeIndex(track["time"])
    # find the overlapping footprints of each pair of satellites
    first, second = [np.array([], dtype=int)], [np.array([], dtype=int)]
    for satellite_1, satellite_2 in combinations(pd.unique(satellites), 2):
        index_1 = np.flatnonzero(satellites == satellite_1)
        index_2 = np.flatnonzero(satellites == satellite_2)
        tree = shapely.STRtree(geometries[index_2])
        query_1, query_2 = tree.query(geometries[index_1], predicate="intersects")
        first.append(index_1[query_1])
        second.append(index_2[query_2])
    first, second = np.concatenate(first), np.concatenate(second)
    # order each pair by time (rows are sorted by time)
    first, second = np.minimum(first, second), np.maximum(first, second)
    time_difference = times[second] - times[first]
    selected = np.asarray(time_difference <= max_time_difference)
    first, second = first[selected], second[selected]
    overlap = shapely.intersection(geometries[first], geometries[second])
    # exclude footprints that only touch
    selected = shapely.area(overlap) > 0
    first, second = first[selected], second[selected]
    columns = [column for column in track.columns if column != track.geometry.name]
    pairs = gpd.GeoDataFrame(
        {
            **{
                f"{column}_{suffix}": track[column].iloc[index].reset_index(drop=True)
                for suffix, index in ((1, first), (2, second))
                for column in columns
            },
            "time_difference": (times[second] - times[first]).to_series(
                index=pd.RangeIndex(len(first))
            ),
        },
        geometry=gpd.GeoSeries(overlap[selected], crs=track.crs),
    )
    return pairs.sort_values(["time_1", "time_2"], kind="stable").reset_index(drop=True)
