"""
Methods to measure how quickly observations respond to events: the time from
an event's onset (such as a fire's ignition or a storm's formation) to the
first observation of its place (its first look), and to the delivery of that
observation's data.

Events may be given (`compute_event_responses`: each with its place, onset,
and optional expiry), or taken to begin at any time with equal chance
(`compute_response_statistics`: the distribution of the responses to events
whose onsets are uniform over a period, computed exactly from the times of the
observations of each target).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Sequence
from datetime import datetime, timedelta

import geopandas as gpd
import numpy as np
import pandas as pd

from ..utils.geometry import _hash_geometries
from .check import _get_observation_times
from .coverage_metrics import _get_targets

RESPONSE_COLUMNS = ("first_look", "response", "looks", "delivered", "delivery")
"""tuple[str, ...]: The columns that `compute_event_responses` adds to events
(`delivered` and `delivery` only for observations with downlinks)."""


def _get_times(values: pd.Series) -> pd.Series:
    """
    Gets times as timezone-aware (UTC) times, of times without time zones in
    UTC.
    """
    return pd.to_datetime(values, utc=True)


def compute_event_responses(
    observations: pd.DataFrame,
    events: gpd.GeoDataFrame,
    *,
    time_column: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Computes the response of observations to events: for each event (a place,
    such as a point or a region, observed as targets are, and its `onset`),
    the time of the first observation of its place at or after its onset
    (`first_look`), the time from its onset to then (`response`), and the
    number of its observations (`looks`) from its onset to its `expires`
    time, if given (after which its observations do not count, as for an
    event that has ended). With the downlinks of the observations (their
    `downlinked` times, see `compute_latencies`, `compute_data_flow`, and
    `evaluate_architecture`), it also gives the earliest delivery of the data
    of those observations (`delivered`, which may come from a later look than
    the first) and the time from the onset to then (`delivery`).

    Args:
        observations (pandas.DataFrame): The observations, with the
            `target_hash` and time (`time_column`) of each, and optionally
            their `downlinked` times.
        events (geopandas.GeoDataFrame): The events: the place of each (its
            geometry, identified by its hash as observations identify their
            targets, so that observations of the same places are those of the
            events), its `onset`, optionally its `expires` time (`NaT` for
            none), and other columns (such as an identifier or priority),
            which the result keeps. Times without time zones are UTC.
        time_column (str | None): The column of the times of the
            observations, or `None` (default) for the first of `epoch`,
            `time`, or `closest_time`, or else the midpoint of `start` and
            `end`.

    Returns:
        geopandas.GeoDataFrame: the events, with the columns of their
            responses (see `RESPONSE_COLUMNS`): `NaT` (and no `looks`) for
            events without observations

    Raises:
        ValueError: If the events have no `onset` column, or an onset is
            missing.
    """
    if "onset" not in events.columns:
        raise ValueError("events require an onset column")
    onset = _get_times(events["onset"])
    if onset.isna().any():
        raise ValueError("every event requires an onset")
    expires = (
        _get_times(events["expires"])
        if "expires" in events.columns
        else pd.Series(pd.NaT, index=events.index, dtype="datetime64[ns, UTC]")
    )
    delivered = "downlinked" in observations.columns
    looks = pd.DataFrame(
        {
            "target_hash": observations["target_hash"].to_numpy(),
            "time": _get_observation_times(observations, time_column).array,
            **(
                {"downlinked": _get_times(observations["downlinked"]).array}
                if delivered
                else {}
            ),
        }
    )
    places = pd.DataFrame(
        {
            "event": np.arange(len(events)),
            "target_hash": _hash_geometries(events.geometry),
            "onset": onset.array,
            "expires": expires.array,
        }
    )
    matched = places.merge(looks, on="target_hash")
    matched = matched[
        (matched["time"] >= matched["onset"])
        & (matched["expires"].isna() | (matched["time"] <= matched["expires"]))
    ]
    grouped = matched.groupby("event")
    first = grouped["time"].min().reindex(places["event"])
    counts = grouped.size().reindex(places["event"], fill_value=0)
    onsets = places["onset"].array
    columns = {
        "first_look": first.array,
        "response": first.array - onsets,
        "looks": counts.to_numpy().astype(int),
    }
    if delivered:
        earliest = grouped["downlinked"].min().reindex(places["event"])
        columns["delivered"] = earliest.array
        columns["delivery"] = earliest.array - onsets
    result = events.copy()
    for name, values in columns.items():
        result[name] = pd.Series(values, index=events.index)
    return result


def _get_waits(
    look: np.ndarray, answer: np.ndarray, length: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the waits for events whose onsets are uniform over a period (from 0
    to its length, seconds) of a target's looks (their times, seconds, in
    order) and the answers to them (the time each answers an event that
    begins before it: its own time, or the delivery of its data, infinite if
    never): an onset after one look and before the next (or after the start
    of the period, before the first) is answered by the earliest answer of
    the next look and those after it, so that its wait falls uniformly
    within an interval.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: the low and high ends of the
            interval of the waits of each interval of onsets that is answered
            (the length of each is that of its onsets)
    """
    # the earliest answer of each look and those after it
    earliest = np.minimum.accumulate(answer[::-1])[::-1]
    previous = np.concatenate([[0.0], look[:-1]])
    low_onset = np.clip(previous, 0, length)
    high_onset = np.clip(look, 0, length)
    kept = (high_onset > low_onset) & np.isfinite(earliest)
    return (earliest - high_onset)[kept], (earliest - low_onset)[kept]


def _get_wait_statistics(
    low: np.ndarray,
    high: np.ndarray,
    length: float,
    quantiles: Sequence[float],
    deadline: float | None,
) -> dict[str, float]:
    """
    Gets the statistics of the waits of a mixture of uniform intervals of
    waits (see `_get_waits`), each of density `1 / length` (an onset per
    second of the period): the fraction of onsets answered, the mean,
    quantiles, and longest of their waits, and the fraction of onsets
    answered within a deadline.
    """
    answered = float(np.sum(high - low)) / length
    if answered <= 0:
        return {
            "response": np.nan,
            **{f"q{q}": np.nan for q in quantiles},
            "max": np.nan,
            "answered": 0.0,
            "within": 0.0,
        }
    mean = float(np.sum((high**2 - low**2) / 2)) / length / answered
    # the distribution function of the waits: piecewise linear, whose slope
    # is the number of intervals that contain the wait (over the length)
    points = np.concatenate([low, high])
    steps = np.concatenate([np.ones(len(low)), -np.ones(len(high))])
    order = np.argsort(points, kind="stable")
    points, steps = points[order], steps[order]
    slopes = np.cumsum(steps)[:-1] / length
    values = np.concatenate([[0.0], np.cumsum(slopes * np.diff(points))])

    def cumulative(wait: float) -> float:
        return float(np.interp(wait, points, values))

    statistics = {"response": mean, "max": float(points[-1]), "answered": answered}
    for q in quantiles:
        target = q * answered
        j = int(np.searchsorted(values, target, side="left"))
        if j == 0:
            statistics[f"q{q}"] = float(points[0])
        else:
            slope = slopes[j - 1]
            statistics[f"q{q}"] = float(
                points[j - 1] + ((target - values[j - 1]) / slope if slope > 0 else 0.0)
            )
    statistics["within"] = (
        np.nan if deadline is None else cumulative(deadline) if deadline >= 0 else 0.0
    )
    return statistics


def compute_response_statistics(
    observations: pd.DataFrame,
    start: datetime,
    end: datetime,
    *,
    quantiles: Sequence[float] = (0.5, 0.9),
    deadline: timedelta | None = None,
    delivered: bool = False,
    targets: gpd.GeoDataFrame | gpd.GeoSeries | None = None,
    time_column: str | None = None,
) -> gpd.GeoDataFrame:
    """
    Computes statistics of the responses of each target's observations to
    events whose onsets are equally likely at any time of a period: the time
    from an event's onset to the next observation of its place (or, with
    `delivered`, to the earliest delivery of the data of an observation at or
    after its onset, by their `downlinked` times). Its distribution follows
    exactly from the times of the observations: an onset in a gap between
    observations waits for the next, so that the waits of each gap are
    uniform, and their mixture is weighted by the lengths of the gaps. Long
    gaps thus count more than in a mean revisit (the mean wait for a first
    look is the sum of the squares of the gaps over twice the period).

    Onsets after a target's last observation (or, with `delivered`, the last
    one delivered) are not answered: `responded` is the fraction of onsets
    answered, and the mean (`response`), quantiles (`response p50`, ...),
    and longest (`max_response`) waits are of those answered. Observations
    after the end of the period answer onsets before it, if given.

    Args:
        observations (pandas.DataFrame): The observations, with the
            `target_hash` and time (`time_column`) of each, and, with
            `delivered`, their `downlinked` times (see `compute_latencies`).
        start (datetime.datetime): The start of the period of onsets.
        end (datetime.datetime): The end of the period of onsets.
        quantiles (Sequence[float]): The quantiles of the waits.
        deadline (datetime.timedelta | None): A deadline, for the fraction of
            all onsets answered within it (`within_deadline`), or `None`.
        delivered (bool): `True` for the waits for the delivery of data,
            rather than for the observations.
        targets (geopandas.GeoDataFrame | geopandas.GeoSeries | None): The
            targets, including those not observed (whose onsets are not
            answered), with the other columns of a data frame of targets,
            or `None` for those observed.
        time_column (str | None): The column of the times of the
            observations, or `None` (default) for the first of `epoch`,
            `time`, or `closest_time`, or else the midpoint of `start` and
            `end`.

    Returns:
        geopandas.GeoDataFrame: a row per target, with `responded`,
            `response`, its quantiles, `max_response`, and, with a deadline,
            `within_deadline`

    Raises:
        ValueError: If the period is empty, or `delivered` observations have
            no `downlinked` column.
    """
    length = (end - start).total_seconds()
    if length <= 0:
        raise ValueError("the period must end after it starts")
    if delivered and "downlinked" not in observations.columns:
        raise ValueError("delivered responses require downlinked times")
    result = _get_targets(observations, targets)
    times = _get_observation_times(observations, time_column)
    seconds = np.asarray((times - pd.Timestamp(start)).dt.total_seconds(), dtype=float)
    if delivered:
        downlinked = _get_times(observations["downlinked"])
        answers = np.asarray(
            (downlinked - pd.Timestamp(start)).dt.total_seconds(), dtype=float
        )
        answers = np.where(np.isnan(answers), np.inf, answers)
    else:
        answers = seconds
    codes, keys = pd.factorize(observations["target_hash"].to_numpy())
    order = np.lexsort((seconds, codes))
    bounds = np.searchsorted(codes[order], np.arange(len(keys) + 1))
    limit = deadline.total_seconds() if deadline is not None else None
    statistics = []
    for k in range(len(keys)):
        rows = order[bounds[k] : bounds[k + 1]]
        low, high = _get_waits(seconds[rows], answers[rows], length)
        statistics.append(_get_wait_statistics(low, high, length, quantiles, limit))
    index = pd.Index(keys).get_indexer(result["target_hash"].to_numpy())
    empty = _get_wait_statistics(np.empty(0), np.empty(0), length, quantiles, limit)

    def column(name: str) -> np.ndarray:
        return np.array(
            [statistics[i][name] if i >= 0 else empty[name] for i in index], dtype=float
        )

    columns: dict[str, object] = {
        "responded": column("answered"),
        "response": pd.to_timedelta(column("response"), unit="s"),
        **{
            f"response p{round(100 * q)}": pd.to_timedelta(column(f"q{q}"), unit="s")
            for q in quantiles
        },
        "max_response": pd.to_timedelta(column("max"), unit="s"),
    }
    if deadline is not None:
        columns["within_deadline"] = column("within")
    return result.assign(**columns)
