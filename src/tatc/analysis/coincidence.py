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
    periods: pd.DataFrame, min_platforms: int
) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """
    Finds the maximal periods when at least a minimum number of platforms
    observe at once, sweeping over the start and end of each platform's
    (closed) observation periods.

    Args:
        periods (pandas.DataFrame): The observation periods, with columns
                `platform`, `start`, and `end`.
        min_platforms (int): The minimum number of platforms.

    Returns:
        list[tuple[pandas.Timestamp, pandas.Timestamp]]: the start and end
            of each coincident period, in time order
    """
    starts, ends = [], []
    for _, platform in periods.groupby("platform"):
        # merge the platform's overlapping periods so that it counts once
        platform = platform.sort_values("start")
        group = (platform["start"] > platform["end"].shift().cummax()).cumsum()
        merged = platform.groupby(group).agg(start=("start", "min"), end=("end", "max"))
        starts.append(merged["start"])
        ends.append(merged["end"])
    times = pd.concat(starts + ends, ignore_index=True)
    events = pd.DataFrame(
        {"time": times, "change": np.repeat([1, -1], len(times) // 2)}
    )
    # sweep the events in time order, with starts before ends at the same
    # time so that periods which touch count as coincident at that instant
    events = events.sort_values(["time", "change"], ascending=[True, False])
    active = events["change"].cumsum().to_numpy() >= min_platforms
    was_active = np.concatenate(([False], active[:-1]))
    return list(
        zip(events["time"][active & ~was_active], events["time"][~active & was_active])
    )


def find_coincident_observations(
    observations: gpd.GeoDataFrame,
    min_platforms: int = 2,
    max_time_difference: timedelta = timedelta(0),
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
    records = []
    for _, gdf in observations.groupby(_get_target_keys(observations)):
        periods = pd.DataFrame(
            {
                "platform": gdf["satellite"],
                "start": gdf["start"] - half_difference,
                "end": gdf["end"] + half_difference,
            }
        )
        for start, end in _find_coincident_periods(periods, min_platforms):
            # each observation overlapping the period contributes to it
            contributors = gdf[
                (periods["start"] <= end) & (periods["end"] >= start)
            ].sort_values("start", kind="stable")
            pairs = list(
                dict.fromkeys(
                    zip(contributors["satellite"], contributors["instrument"])
                )
            )
            records.append(
                {
                    "target_hash": gdf["target_hash"].iloc[0],
                    "geometry": gdf.geometry.iloc[0],
                    "satellite": ", ".join(satellite for satellite, _ in pairs),
                    "instrument": ", ".join(instrument for _, instrument in pairs),
                    "start": start,
                    "epoch": start + (end - start) / 2,
                    "end": end,
                }
            )
    if len(records) == 0:
        return _get_empty_coincidence_frame()
    return (
        gpd.GeoDataFrame(records, crs=observations.crs)
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
