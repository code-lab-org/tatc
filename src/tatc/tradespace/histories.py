"""
Many histories of missions in a replicate (simulated by
`Deployment.simulate_many`): their evaluation with a cache of the metrics
of their phases (see `PhaseCache`), their summaries as the records of
evaluations, with control variates, and the allocation of evaluations
between replicates and histories. The exact expectations of missions
without restoration are in the module `expectations`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .criteria import Requirement
from .mission import Capability, MissionHistory, summarize_phases
from .statistics import Probability, Quantile, StandardDeviation, Statistic


def evaluate_histories(
    histories: Sequence[MissionHistory],
    capability: Capability,
    *,
    cost_model: Any = None,
    growth: Any = 1.0,
    requirements: Sequence[Requirement] = (),
) -> pd.DataFrame:
    """
    Evaluates histories of a mission (such as those of a replicate, see
    `Deployment.simulate_many`): for each, the mean of each metric over the
    phases of the mission (`<metric> mission`, see `summarize_phases`), the
    time at which a phase first meets each requirement on a metric
    (`<metric> time`, infinite if none does), its measures (see
    `MissionHistory.metrics`), and its cost (see `CostEstimate.metrics`,
    with a cost model). The histories share the observations of the
    replicate, so that each set of operating satellites is evaluated once,
    by the capability of the satellites (see `Capability.evaluate`, with the
    cache of its metric).

    Args:
        histories (Sequence[MissionHistory]): The histories.
        capability (Capability): The capability of the satellites: their
            observations (of the replicate) and metric.
        cost_model (CostModel | None): A model of the cost of each history.
        growth (float | Mapping[str, float]): The growth of the cost.
        requirements (Sequence[Requirement]): Requirements on the metrics of
            the phases, if any.

    Returns:
        pandas.DataFrame: a row per history (with its likelihood ratio,
            `weight`, if the histories are of stressed failures)
    """
    rows = []
    weighted = any(history.weight != 1.0 for history in histories)
    for history in histories:
        phases = capability.evaluate(history.phases(capability.satellites))
        row = {
            **summarize_phases(phases, requirements=requirements),
            **history.metrics(),
        }
        if weighted:
            row["weight"] = history.weight
        if cost_model is not None:
            row.update(cost_model.estimate(history, growth).metrics())
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_histories(
    histories: pd.DataFrame,
    *,
    statistics: Sequence[Statistic] = (),
    controls: Mapping[str, float] | None = None,
) -> dict[str, float]:
    """
    Summarizes histories of a mission (see `evaluate_histories`) as the
    metrics of an evaluation: the mean of each column over the histories
    (by the same name), and `statistics` of columns over the histories, as
    those of `Criteria` are over the replicates, named after the histories:
    a `Quantile` (`<column> histories
    p<percent>`, such as `"accurate mission histories p10"` for the 0.1
    quantile, the risk of the mission), a `Probability` (`P(<column>
    <operator> <value>) histories`, such as `"P(accurate time <= 4)
    histories"`, the share of the histories that meet it), or a
    `StandardDeviation` (`<column> histories std`, as `allocate_histories`
    needs of a pilot study). The statistics over the histories of a
    replicate end in `histories`, apart from the statistics over the
    replicates (see `Criteria`), which are taken of these.

    Missing values (`NaN`, such as the `restoration years` of a history
    without restorations) are left out of the statistics of their column,
    but infinite ones are not: the mean of the time until a history meets a
    requirement (`<metric> time`, see `evaluate_histories`) is infinite if a
    history never meets it, and its quantiles, or the share of the
    histories that meet it by a time (such as `Probability("accurate time",
    "<=", 4)`), summarize it.

    With `controls`, the means are estimated with control variates: columns
    whose expected values are known (such as the `availability` of a
    mission, which simulations alone estimate precisely and cheaply), and
    which vary with the other columns from history to history. The mean of
    each other column is corrected by the deviation of the controls' means
    from their expected values, times the regression of the column on the
    controls over the histories, which removes the share of its variance
    that the controls explain (a share of `rho ** 2` for a correlation
    `rho`): fewer histories give the same precision. The other statistics
    are not corrected.

    Histories of stressed failures (see `Deployment.simulate_many`) are
    weighted by their likelihood ratios (`weight`), for every statistic, and
    the effective number of histories (`effective histories`, the square of
    the sum of the weights over the sum of their squares) measures how much
    a few weights dominate.

    The statistics are over the histories of a replicate, whose statistics
    over the replicates (such as their mean, see `summarize`) combine the
    risk of the mission with the uncertainty of the replicates' other
    parameters (such as their targets).

    Args:
        histories (pandas.DataFrame): A row per history.
        statistics (Sequence[Statistic]): Statistics of columns over the
            histories: `Quantile`, `Probability`, or `StandardDeviation`.
        controls (Mapping[str, float] | None): The expected values of
            columns that are control variates, if any.

    Returns:
        dict[str, float]: the summary

    Raises:
        TypeError: if a statistic is not one of these.
        ValueError: if a statistic is of a column that the histories do not
            have.
    """
    statistics = _history_statistics(statistics)
    numeric = histories.select_dtypes("number")
    named = [statistic.column for statistic in statistics]
    unknown = sorted(set(named) - set(numeric.columns))
    if unknown:
        raise ValueError(
            "statistics of columns that the histories do not have: "
            + f"{unknown} (they have {sorted(numeric.columns)})"
        )
    weights = np.ones(len(numeric))
    if "weight" in numeric.columns:
        weights = numeric.pop("weight").to_numpy(dtype=float)
    summary = {
        column: _weighted_mean(numeric[column].to_numpy(dtype=float), weights)
        for column in numeric.columns
    }
    if "weight" in histories.columns:
        summary["effective histories"] = float(weights.sum() ** 2 / np.sum(weights**2))
    if controls:
        summary.update(_controlled_means(numeric, controls, weights))
    for statistic in statistics:
        values = numeric[statistic.column].to_numpy(dtype=float)
        summary[_history_name(statistic)] = _weighted_statistic(
            statistic, values, weights
        )
    return summary


_HISTORY_STATISTICS = (Quantile, Probability, StandardDeviation)
"""The statistics of columns over histories (see `summarize_histories`)."""


def _history_statistics(statistics: Sequence[Statistic]) -> tuple[Statistic, ...]:
    """
    Gets statistics over histories (see `summarize_histories`), raising a
    `TypeError` if one is not a statistic of a column of a kind that
    summarizes histories.
    """
    if isinstance(statistics, Statistic):
        statistics = [statistics]
    statistics = tuple(statistics)
    for statistic in statistics:
        if type(statistic) not in _HISTORY_STATISTICS or not isinstance(
            statistic.column, str
        ):
            raise TypeError(
                "statistics over histories are a Quantile, Probability, or "
                + f"StandardDeviation of a column, not {statistic!r}"
            )
    return statistics


def _history_name(statistic: Statistic) -> str:
    """
    Gets the name of a statistic of a column over the histories of a
    mission (see `summarize_histories`): `<column> histories p<percent>` of a
    `Quantile`, `P(<column> <operator> <value>) histories` of a
    `Probability`, and `<column> histories std` of a `StandardDeviation`.
    """
    if isinstance(statistic, Quantile):
        return f"{statistic.column} histories p{100 * statistic.q:g}"
    if isinstance(statistic, Probability):
        return f"{statistic.name} histories"
    return f"{statistic.column} histories std"


def _weighted_statistic(
    statistic: Statistic, values: np.ndarray, weights: np.ndarray
) -> float:
    """
    Computes a statistic over histories (see `summarize_histories`) of their
    values, by their weights, leaving out missing values (`NaN`).
    """
    if isinstance(statistic, Quantile):
        return _weighted_quantile(values, weights, statistic.q)
    known = ~np.isnan(values)
    if isinstance(statistic, Probability):
        if not known.any() or weights[known].sum() <= 0:
            return math.nan
        requirement = Requirement(
            statistic.column, statistic.operator, statistic.threshold
        )
        met = requirement.met(values[known])
        return float(np.sum(weights[known] * met) / np.sum(weights[known]))
    mean = _weighted_mean(values, weights)
    variance = _weighted_mean((values - mean) ** 2, weights)
    count = known.sum()
    return float(np.sqrt(variance * count / (count - 1))) if count > 1 else math.nan


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    """
    The weighted mean of the known values (not `NaN`, but infinite ones),
    `NaN` without any.
    """
    known = ~np.isnan(values)
    if not known.any() or weights[known].sum() <= 0:
        return math.nan
    return float(np.sum(weights[known] * values[known]) / np.sum(weights[known]))


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """
    The quantile of the known values (not `NaN`, but infinite ones), by
    their weights: the least value whose cumulative weight reaches the share
    `q` of the total (as numpy's `inverted_cdf`), or numpy's linear quantile
    of finite values without weights.
    """
    known = ~np.isnan(values)
    values, weights = values[known], weights[known]
    if len(values) == 0:
        return math.nan
    if np.all(weights == weights[0]) and np.all(np.isfinite(values)):
        return float(np.quantile(values, q))
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order]) / weights.sum()
    position = min(int(np.searchsorted(cumulative, q - 1e-12)), len(values) - 1)
    return float(values[order][position])


def _controlled_means(
    numeric: pd.DataFrame, controls: Mapping[str, float], weights: np.ndarray
) -> dict:
    """
    Estimates the means of columns with control variates (see
    `summarize_histories`): each column's mean less the regression of the
    column on the controls applied to the deviations of the controls' means
    from their expected values, over the histories in which the column and
    the controls are known.
    """
    names = list(controls)
    missing = set(names) - set(numeric.columns)
    if missing:
        raise ValueError(f"no columns of controls {sorted(missing)}")
    expected = np.array([float(controls[name]) for name in names])
    means = {name: float(value) for name, value in zip(names, expected)}
    for column in numeric.columns:
        if column in controls:
            continue
        data = numeric[[column, *names]].to_numpy(dtype=float)
        known = np.all(np.isfinite(data), axis=1)
        data, w = data[known], weights[known]
        if len(data) <= len(names) + 1 or w.sum() <= 0:
            continue
        y, x = data[:, 0], data[:, 1:]
        y_mean, x_mean = np.average(y, weights=w), np.average(x, axis=0, weights=w)
        centered = x - x_mean
        if not np.any(centered):
            continue
        root = np.sqrt(w)[:, None]
        beta = np.linalg.lstsq(root * centered, root[:, 0] * (y - y_mean), rcond=None)[
            0
        ]
        means[column] = float(y_mean - (x_mean - expected) @ beta)
    return means


def allocate_histories(
    records: pd.DataFrame,
    metric: str,
    *,
    histories: int,
    replicate_seconds: float,
    history_seconds: float,
    standard_error: float | None = None,
    max_histories: int = 10_000,
) -> dict[str, float]:
    """
    Recommends how to divide the evaluations of a design between replicates
    (each a sample of its uncertain parameters other than failures, such as
    its targets, whose observations are expensive to compute) and histories
    of the failures of its mission in each replicate (cheap, with a
    `PhaseCache`), from a pilot study: the records of its replicates (with
    the statistic `StandardDeviation(metric)` of `Failures`, see
    `summarize_histories`), each the mean (`<metric>`) and standard
    deviation (`<metric> histories std`) of a metric over `histories`
    histories.

    The variance of a replicate's mean is `between + within / K` for `K`
    histories, with `within` the variance between the histories of a
    replicate and `between` the variance between replicates of their
    expected values, at a cost of `replicate_seconds + K * history_seconds`.
    The variance of the mean of the replicates per second of evaluation is
    least for `K = sqrt(replicate_seconds * within / (history_seconds *
    between))`: more histories where failures vary more than replicates, or
    where a history costs less than a replicate. With a `standard_error` of
    the mean, the replicates that achieve it at that number of histories.

    The variances are estimated from independent histories, so the pilot's
    must be (`sampling="random"`, see `Deployment.simulate_many`): the means
    of stratified histories vary less than those of independent ones, which
    would hide the variance between replicates (a warning is given if the
    means vary less than independent histories can). Studies with stratified
    histories then need fewer histories than recommended. With a standard
    error that the least replicates (two) achieve, the fewest histories that
    achieve it are recommended.

    Args:
        records (pandas.DataFrame): The records of the pilot replicates of a
            design.
        metric (str): The metric.
        histories (int): The histories of each pilot replicate.
        replicate_seconds (float): The cost of a replicate without its
            histories (such as computing the observations of its targets).
        history_seconds (float): The cost of a history.
        standard_error (float | None): A standard error of the mean of the
            metric over the replicates to achieve, if any.
        max_histories (int): The most histories to recommend.

    Returns:
        dict[str, float]: the recommended `histories` (and `replicates`, with
            a standard error), the standard deviations `between` and
            `within`, and the time of a replicate (`replicate seconds`) and of
            the evaluations (`seconds`, with a standard error)
    """
    if len(records) < 2:
        raise ValueError("a pilot needs at least two replicates")
    if histories < 2 or replicate_seconds <= 0 or history_seconds <= 0:
        raise ValueError("a pilot needs two or more histories and positive costs")
    means = records[metric].to_numpy(dtype=float)
    within = float(
        np.mean(records[f"{metric} histories std"].to_numpy(dtype=float) ** 2)
    )
    observed = float(np.var(means, ddof=1))
    if observed < 0.5 * within / histories:
        warnings.warn(
            "the replicates' means vary less than independent histories can: "
            + "pilots need independent histories (sampling='random')",
            stacklevel=2,
        )
    between = max(observed - within / histories, 0.0)
    if between == 0:
        count = max_histories if within > 0 else 1
    else:
        count = math.sqrt(replicate_seconds * within / (history_seconds * between))
    count = int(min(max(round(count), 1), max_histories))
    cost = replicate_seconds + count * history_seconds
    result = {
        "histories": float(count),
        "between": math.sqrt(between),
        "within": math.sqrt(within),
        "replicate seconds": cost,
    }
    if standard_error is not None:
        if standard_error <= 0:
            raise ValueError("standard_error must be positive")
        variance = between + within / count
        replicates = max(math.ceil(variance / standard_error**2), 2)
        if replicates == 2 and 2 * standard_error**2 > between:
            # the least replicates suffice: the fewest histories that do
            needed = math.ceil(within / (2 * standard_error**2 - between))
            count = int(min(max(needed, 1), count))
            cost = replicate_seconds + count * history_seconds
            result.update(histories=float(count), **{"replicate seconds": cost})
        result.update(replicates=float(replicates), seconds=replicates * cost)
    return result
