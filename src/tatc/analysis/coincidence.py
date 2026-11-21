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

from .check import (
    _PLATFORM_COLUMNS,
    _PLATFORM_KINDS,
    _aggregate_names,
    _get_empty_observation_frame,
    _get_instruments,
    _get_platform_columns,
    _get_platform_kinds,
    _get_platforms,
)
from .coverage_metrics import _get_naive_utc, _get_target_keys, _group_overlapping


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
    order, merged = _group_overlapping(
        target * (int(platform.max()) + 1) + platform, start, end
    )
    merged_target = pd.Series(target[order]).groupby(merged).first().to_numpy()
    merged_start = pd.Series(start[order]).groupby(merged).min().to_numpy()
    merged_end = pd.Series(end[order]).groupby(merged).max().to_numpy()
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
    `min_platforms` platforms observe the same target (the same
    `target_hash`, see `_get_target_keys`) at once, as collected by
    `tatc.analysis.point_sampling.collect_observations` for satellites or
    `tatc.analysis.footprints.collect_footprint_observations` for the
    footprints of any platform (concatenated, for example, to find when an
    aircraft and a satellite observe the same point). Each observation's
    platform is named by its `satellite`, `aircraft`, `surface_platform`, or
    `radar_station` column. Observations by different instruments on the
    same platform are not coincident with each other.

    Each observation period is widened by half of `max_time_difference`
    before and after, so that two observations coincide if the gap between
    their periods is at most `max_time_difference` (for instantaneous
    observations, if their epochs differ by at most `max_time_difference`).
    Each coincident period is maximal: it spans the time that at least
    `min_platforms` widened periods overlap, even as the platforms
    observing change; the platform columns (such as `satellite`) and
    `instrument` record every contributor to that period, comma-separated,
    in order of their kind (in the order of the columns) and then of their
    start time: each platform column names the contributors of its kind
    (missing if there are none), so that the instruments pair with the
    names of the platform columns in turn (radar stations, last, have none).
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
        min_platforms (int): The minimum number of platforms (at least 2)
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
    platform_columns = _get_platform_columns(observations)
    if observations.empty:
        return _get_empty_observation_frame(platform_columns)
    if not (observations.geom_type == "Point").all():
        raise ValueError(
            "observations must be of points (see find_coincident_ground_track "
            "for the footprints of ground tracks)"
        )
    gdf = observations.reset_index(drop=True)
    # targets (in sorted order of their keys) and platforms of the observations
    target, _ = pd.MultiIndex.from_arrays(_get_target_keys(gdf)).factorize(sort=True)
    platform, _ = pd.factorize(_get_platforms(gdf))
    # observation times: each period, or each epoch if instantaneous
    observed_start = gdf["epoch"] if instantaneous else gdf["start"]
    observed_end = gdf["epoch"] if instantaneous else gdf["end"]
    start = _get_naive_utc(observed_start - half_difference)
    end = _get_naive_utc(observed_end + half_difference)
    period_target, period_start, period_end = _find_coincident_periods(
        target, platform, start, end, min_platforms
    )
    if len(period_target) == 0:
        return _get_empty_observation_frame(platform_columns)
    # each observation overlapping a period contributes to it, in order of
    # its start time
    candidates = pd.DataFrame(
        {
            "target": target,
            "widened_start": start,
            "widened_end": end,
            "start": _get_naive_utc(observed_start),
            "rank": pd.Index(_PLATFORM_COLUMNS).get_indexer(_get_platform_kinds(gdf)),
            "platform": platform,
            **{column: gdf[column].to_numpy() for column in platform_columns},
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
        contributors.sort_values(["period", "rank", "start"], kind="stable")
        .drop_duplicates(["period", "platform", "instrument"])
        .groupby("period")
    )
    contributors = pd.DataFrame(
        {
            column: _aggregate_names(contributors, column)
            for column in [*platform_columns, "instrument"]
        }
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
                **{
                    column: contributors[column].to_numpy()
                    for column in [*platform_columns, "instrument"]
                },
                "start": start,
                "epoch": start + (end - start) / 2,
                "end": end,
            },
            crs=observations.crs,
        )
        .sort_values("start", kind="stable")
        .reset_index(drop=True)
    )


def _get_overlaps(
    geometries_1: np.ndarray, geometries_2: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the overlap of each pair of footprints, and whether it has a
    positive area (footprints that only touch do not overlap).

    Args:
        geometries_1 (numpy.ndarray): The first footprint of each pair.
        geometries_2 (numpy.ndarray): The second footprint of each pair.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: the overlap of each pair (empty
            if they do not intersect), and whether it has a positive area
    """
    overlaps = np.full(len(geometries_1), None, dtype=object)
    meet = shapely.intersects(geometries_1, geometries_2)
    overlaps[meet] = shapely.intersection(geometries_1[meet], geometries_2[meet])
    selected = np.zeros(len(geometries_1), dtype=bool)
    selected[meet] = shapely.area(overlaps[meet]) > 0
    return overlaps, selected


def _keep_nearest_copies(
    other: np.ndarray,
    copy: np.ndarray,
    valid: np.ndarray,
    instruments: np.ndarray,
    times: pd.DatetimeIndex,
) -> np.ndarray:
    """
    Selects, of the overlapping pairs of a footprint of a platform (`other`)
    and a timed footprint of a fixed platform (`copy`, see
    `_PlatformKind.fixed`), the copy of each instrument nearest in time to
    the other footprint: a fixed platform's footprint is the same at every
    time, so its copies at other times would pair with the other footprint
    again (see `find_coincident_ground_track`). A valid copy (see
    `valid_obs`) is preferred to an invalid one, and of two as near, the
    earlier. Pairs with a static footprint are kept.

    Args:
        other (numpy.ndarray): The row of the other footprint of each pair.
        copy (numpy.ndarray): The row of the fixed platform's footprint of
            each pair.
        valid (numpy.ndarray): Whether the footprint of each row is valid.
        instruments (numpy.ndarray): The instrument code of each row.
        times (pandas.DatetimeIndex): The time of each row (`NaT` if static).

    Returns:
        numpy.ndarray: whether to keep each pair
    """
    keep = np.ones(len(other), dtype=bool)
    static = np.asarray(times.isna())
    index = np.flatnonzero(~static[copy] & ~static[other])
    if len(index) == 0:
        return keep
    copy, other = copy[index], other[index]
    copy_times = times.asi8[copy]
    distance = np.abs(copy_times - times.asi8[other])
    order = np.lexsort((copy_times, distance, ~valid[copy], instruments[copy], other))
    group_other, group_instrument = other[order], instruments[copy][order]
    nearest = np.ones(len(order), dtype=bool)
    nearest[1:] = (group_other[1:] != group_other[:-1]) | (
        group_instrument[1:] != group_instrument[:-1]
    )
    keep[index[order[~nearest]]] = False
    return keep


def _query_fixed_pairs(
    others: np.ndarray,
    copies: np.ndarray,
    geometries: np.ndarray,
    instruments: np.ndarray,
    valid: np.ndarray,
    times: pd.DatetimeIndex,
    max_time_difference: pd.Timedelta,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Finds the candidate pairs of the footprints of a platform (`others`) and
    of a fixed platform (`copies`, see `_PlatformKind.fixed`): those whose
    bounds overlap at times within the maximum time difference (or static),
    as for other platforms (see `find_coincident_ground_track`), but of the
    timed copies of the fixed platform's footprint only the nearest in time
    of each instrument, geometry, and validity (which `_keep_nearest_copies`
    then reduces, of those that overlap, to the nearest of each instrument).
    The copies of each instrument, geometry, and validity share one entry of
    a spatial index, so that the number of candidates grows with the number
    of footprints rather than with the product of their numbers.

    Args:
        others (numpy.ndarray): The rows of the footprints of the platform.
        copies (numpy.ndarray): The rows of the footprints of the fixed
            platform.
        geometries (numpy.ndarray): The geometry of each row.
        instruments (numpy.ndarray): The instrument code of each row.
        valid (numpy.ndarray): Whether the footprint of each row is valid.
        times (pandas.DatetimeIndex): The time of each row (`NaT` if static).
        max_time_difference (pandas.Timedelta): The maximum time difference.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: the rows of the other footprint
            and of the copy of each candidate pair
    """
    static = np.asarray(times.isna())
    nanoseconds = times.as_unit("ns").asi8
    pair_others, pair_copies = [np.array([], dtype=int)], [np.array([], dtype=int)]

    def query(query_rows: np.ndarray, tree_rows: np.ndarray) -> None:
        # every pair whose bounds overlap
        if len(query_rows) > 0 and len(tree_rows) > 0:
            found, tree = shapely.STRtree(geometries[tree_rows]).query(
                geometries[query_rows]
            )
            pair_others.append(query_rows[found])
            pair_copies.append(tree_rows[tree])

    # a static footprint is coincident at any time: a static copy with every
    # footprint of the other platform, and a static footprint of the other
    # platform with every timed copy
    timed_copies = copies[~static[copies]]
    timed_others = others[~static[others]]
    query(others, copies[static[copies]])
    query(others[static[others]], timed_copies)
    if len(timed_copies) > 0 and len(timed_others) > 0:
        # the groups of copies of the same instrument, geometry, and
        # validity, each in order of time
        shape, _ = pd.factorize(shapely.to_wkb(geometries[timed_copies]))
        group, _ = pd.factorize(
            pd.MultiIndex.from_arrays(
                [instruments[timed_copies], shape, valid[timed_copies]]
            )
        )
        groups = int(group.max()) + 1
        order = np.lexsort((nanoseconds[timed_copies], group))
        sorted_rows = timed_copies[order]
        bounds = np.searchsorted(group[order], np.arange(groups + 1))
        # the groups whose bounds overlap each timed footprint of the other
        # platform, and the nearest copy of each (the earlier of two as near)
        found, found_group = shapely.STRtree(
            geometries[sorted_rows[bounds[:-1]]]
        ).query(geometries[timed_others])
        other = timed_others[found]
        nearest = np.empty(len(other), dtype=int)
        by_group = np.argsort(found_group, kind="stable")
        splits = np.searchsorted(found_group[by_group], np.arange(groups + 1))
        for g in np.flatnonzero(np.diff(splits) > 0):
            index = by_group[splits[g] : splits[g + 1]]
            members = sorted_rows[bounds[g] : bounds[g + 1]]
            member_times = nanoseconds[members]
            other_times = nanoseconds[other[index]]
            position = np.searchsorted(member_times, other_times)
            before = np.clip(position - 1, 0, len(members) - 1)
            after = np.clip(position, 0, len(members) - 1)
            later = np.abs(member_times[after] - other_times) < np.abs(
                member_times[before] - other_times
            )
            nearest[index] = members[np.where(later, after, before)]
        within = (
            np.abs(nanoseconds[nearest] - nanoseconds[other])
            <= max_time_difference.value
        )
        pair_others.append(other[within])
        pair_copies.append(nearest[within])
    return np.concatenate(pair_others), np.concatenate(pair_copies)


def _get_copy_side(
    code_1: int,
    code_2: int,
    fixed: bool,
    fixed_2: bool,
    ranks: np.ndarray,
    names: pd.Index,
) -> int:
    """
    Gets the platform whose timed footprints are the copies paired with the
    footprints of the other (see `_keep_nearest_copies`): the fixed one or,
    if both are fixed, the later in the order of kinds (of `between`, or
    else of `_PLATFORM_COLUMNS`) and then of names, whatever the order of
    the rows.

    Args:
        code_1 (int): The code of the first platform.
        code_2 (int): The code of the second platform.
        fixed (bool): `True` if the first platform is fixed.
        fixed_2 (bool): `True` if the second platform is fixed.
        ranks (numpy.ndarray): The rank of the kind of each platform code.
        names (pandas.Index): The key of each platform code.

    Returns:
        int: the code of the platform of the copies
    """
    if fixed and fixed_2:
        return max(code_1, code_2, key=lambda code: (ranks[code], names[code]))
    return code_1 if fixed else code_2


def _check_between(between: tuple[str, str] | None) -> tuple[str, str] | None:
    """
    Checks the kinds of platforms between which to find coincident footprints.

    Args:
        between (tuple[str, str] | None): The platform columns of the two
            kinds, or `None` for any.

    Returns:
        tuple[str, str] | None: the platform columns
    """
    if between is None:
        return None
    if isinstance(between, str) or len(between) != 2:
        raise ValueError(f"between must be a pair of platform columns, not {between!r}")
    for column in between:
        if column not in _PLATFORM_COLUMNS:
            raise ValueError(
                f"between must name platform columns ({', '.join(_PLATFORM_COLUMNS)}), "
                f"not {column!r}"
            )
    return (between[0], between[1])


def find_coincident_ground_track(
    ground_track: gpd.GeoDataFrame,
    max_time_difference: timedelta = timedelta(0),
    between: tuple[str, str] | None = None,
) -> gpd.GeoDataFrame:
    """
    Find coincident observations of ground tracks: the pairs of footprints
    of different platforms, as collected by
    `tatc.analysis.footprints.collect_footprints` for any platforms together,
    or by `tatc.analysis.ground_track.collect_ground_track` for satellites,
    `tatc.analysis.airborne_track.collect_airborne_ground_track` for
    aircraft, `tatc.analysis.surface_track.collect_surface_track` for fixed
    surface platforms, and `tatc.analysis.radar_track.collect_radar_track`
    for radar stations, concatenated (for example, to find the underflights
    of a satellite by an aircraft, or the overpasses of a tower's view or a
    radar's coverage), whose times differ by at most `max_time_difference`
    and which overlap (share a positive area, in longitude/latitude
    coordinates). Each footprint's platform is named by its `satellite`,
    `aircraft`, `surface_platform`, or `radar_station` column. Footprints of
    different instruments on the same platform are not coincident with each
    other. A footprint without a time (as of a fixed platform or radar
    station collected without times, whose footprint does not change) is
    static: it is coincident with the overlapping footprints of other
    platforms at any time.

    The footprint of a fixed surface platform or radar station collected at
    times is the same at every time, so a footprint of another platform
    pairs only with its overlapping copy nearest in time (of each
    instrument; a valid observation in preference to an invalid one, and the
    earlier of two as near), rather than with every copy within
    `max_time_difference`. For two fixed platforms, each footprint of the
    one of the first kind (of `between`, or else in the order of the
    platform columns), or of two of a kind the first by name, pairs with the
    nearest copy of the other.

    Each pair records the columns of both footprints, suffixed `_1` and
    `_2` (for example, `time_1`, `satellite_1`, and `valid_obs_1`) and
    ordered so that `time_1` is no later than `time_2` (a static footprint
    second), the time difference `time_difference` (`time_2` less `time_1`,
    or `NaT` with a static footprint), whether both are valid observations
    (`valid_obs`, if the footprints have one), and the overlap of the
    footprints as its geometry. Dissolving the result (see
    `tatc.analysis.footprints.dissolve_footprints`, which selects the pairs
    of valid observations) gives the area observed in common.

    With `between`, only the pairs of a footprint of one kind of platform
    and one of another (or of two platforms of one kind) are found, for
    example `("satellite", "aircraft")` for the underflights of satellites
    by aircraft, without the pairs of two satellites. The footprint of the
    first kind is then first (suffixed `_1`, for example `satellite_1` and
    `aircraft_2`), whatever its time, so that `time_difference` is negative
    if it is later.

    Footprints whose bounds overlap are first matched by time, so that
    only those within `max_time_difference` are tested for overlap. The
    timed copies of a fixed platform's footprint are indexed once (for each
    instrument and geometry), so that the time to pair them grows with the
    number of footprints rather than with the product of the numbers of
    footprints of the two platforms.

    Args:
        ground_track (geopandas.GeoDataFrame): The collected ground track.
        max_time_difference (datetime.timedelta): The maximum time
                difference between coincident footprints.
        between (tuple[str, str] | None): The platform columns (`satellite`,
                `aircraft`, `surface_platform`, or `radar_station`) of the
                two kinds of platforms between which to find coincident
                footprints, or `None` for any platforms.

    Returns:
        geopandas.GeoDataFrame: The data frame of coincident footprints,
            sorted by `time_1` and `time_2` (static footprints last).
    """
    max_time_difference = _check_max_time_difference(max_time_difference)
    between = _check_between(between)
    track = ground_track
    if "time" not in track.columns:
        # every footprint is static
        track = track.assign(
            time=pd.Series(pd.NaT, index=track.index, dtype="datetime64[ns, UTC]")
        )
    # in order of time, with static footprints (without a time) last
    track = track.sort_values("time", kind="stable", na_position="last")
    track = track.reset_index(drop=True)
    geometries = track.geometry.to_numpy()
    kinds = _get_platform_kinds(track)
    fixed = np.isin(kinds, [kind.column for kind in _PLATFORM_KINDS if kind.fixed])
    instruments, _ = pd.factorize(pd.Series(_get_instruments(track)).astype(str))
    platforms, names = pd.factorize(_get_platforms(track))
    times = pd.DatetimeIndex(track["time"])
    valid = (
        track["valid_obs"].eq(True).to_numpy()
        if "valid_obs" in track.columns
        else np.ones(len(track), dtype=bool)
    )
    # the footprints of each platform (in time order), and their spatial
    # index (built once)
    order = np.argsort(platforms, kind="stable")
    rows = np.split(order, np.cumsum(np.bincount(platforms, minlength=len(names)))[:-1])
    # the order of the kinds of platforms: that of `between`, or else of the
    # platform columns
    ranks = pd.Index(list(dict.fromkeys(between or _PLATFORM_COLUMNS))).get_indexer(
        [kinds[members[0]] if len(members) > 0 else None for members in rows]
    )
    trees: dict[int, shapely.STRtree] = {}

    def order_pairs(
        row_1: np.ndarray, row_2: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        # each pair in order of time (rows are sorted by time, static ones
        # last) or, with `between`, of kind
        first, second = np.minimum(row_1, row_2), np.maximum(row_1, row_2)
        if between is None:
            return first, second
        swap = kinds[first] != between[0]
        return np.where(swap, second, first), np.where(swap, first, second)

    # find the pairs of each pair of platforms (of the kinds between which to
    # find them): footprints whose bounds overlap, at times within the
    # maximum time difference (or static), that overlap in a positive area
    # (not only touching)
    first, second = [np.array([], dtype=int)], [np.array([], dtype=int)]
    overlaps = [np.array([], dtype=object)]
    for code_1, code_2 in combinations(range(len(rows)), 2):
        pair = (kinds[rows[code_1][0]], kinds[rows[code_2][0]])
        if between is not None and pair not in (between, between[::-1]):
            continue
        is_fixed = fixed[rows[code_1][0]] or fixed[rows[code_2][0]]
        if is_fixed:
            # the copies of a fixed platform's footprint nearest in time
            side = _get_copy_side(
                code_1,
                code_2,
                bool(fixed[rows[code_1][0]]),
                bool(fixed[rows[code_2][0]]),
                ranks,
                names,
            )
            row_1, row_2 = _query_fixed_pairs(
                rows[code_2 if side == code_1 else code_1],
                rows[side],
                geometries,
                instruments,
                valid,
                times,
                max_time_difference,
            )
        else:
            if code_2 not in trees:
                trees[code_2] = shapely.STRtree(geometries[rows[code_2]])
            query_1, query_2 = trees[code_2].query(geometries[rows[code_1]])
            row_1, row_2 = rows[code_1][query_1], rows[code_2][query_2]
            time_difference = times[row_2] - times[row_1]
            # a static footprint is coincident at any time
            selected = np.asarray(
                abs(time_difference) <= max_time_difference
            ) | np.asarray(time_difference.isna())
            row_1, row_2 = row_1[selected], row_2[selected]
        # each pair in order (see `order_pairs`), whose overlap is that of
        # its first and second footprints
        pair_1, pair_2 = order_pairs(row_1, row_2)
        overlap, selected = _get_overlaps(geometries[pair_1], geometries[pair_2])
        if is_fixed:
            # of the overlapping copies, the nearest (row_2 are the copies)
            nearest = _keep_nearest_copies(
                row_1[selected], row_2[selected], valid, instruments, times
            )
            selected[np.flatnonzero(selected)[~nearest]] = False
        first.append(pair_1[selected])
        second.append(pair_2[selected])
        overlaps.append(overlap[selected])
    first, second = np.concatenate(first), np.concatenate(second)
    overlap = np.concatenate(overlaps)
    columns = [column for column in track.columns if column != track.geometry.name]
    validity = (
        {"valid_obs": valid[first] & valid[second]}
        if "valid_obs" in track.columns
        else {}
    )
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
            **validity,
        },
        geometry=gpd.GeoSeries(overlap, crs=track.crs),
    )
    return pairs.sort_values(["time_1", "time_2"], kind="stable").reset_index(drop=True)
