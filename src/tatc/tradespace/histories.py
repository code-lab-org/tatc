"""
Many histories of missions in a replicate (see `Deployment.simulate_many`):
their simulation (independent or stratified), evaluation with a cache of
the metrics of their phases (see `PhaseCache`), and summaries as the
records of evaluations, with control variates and the allocation of
evaluations between replicates and histories.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from ..schemas import Satellite
from .cache import get_satellite_key
from .criteria import _OPERATORS, Requirement
from .deployment import Deployment
from .failures import FailureModel, RestorationPolicy
from .mission import (
    MissionHistory,
    PhaseCache,
    _Simulation,
    evaluate_phases,
    summarize_phases,
)


def evaluate_missions(
    histories: Sequence[MissionHistory],
    satellites: Sequence[Satellite],
    observations: pd.DataFrame,
    metric: Callable[[pd.DataFrame], Mapping[str, float]] | PhaseCache,
    context: Any = None,
    cost_model: Any = None,
    growth: Any = 1.0,
) -> pd.DataFrame:
    """
    Evaluates histories of a mission (such as those of a replicate, see
    `Deployment.simulate_many`): for each, the mean of each metric over the
    phases of the mission (`<metric> mean`, see `summarize_phases`), its
    measures (see `MissionHistory.metrics`), and its cost (see
    `CostEstimate.metrics`, with a cost model). The histories share the
    observations of the replicate, so that a `PhaseCache` (as the metric)
    evaluates each set of operating satellites once.

    Args:
        histories (Sequence[MissionHistory]): The histories.
        satellites (Sequence[Satellite]): The satellites (see
            `MissionHistory.phases`).
        observations (pandas.DataFrame): The observations of every
            satellite, with their names (`satellite`).
        metric (Callable | PhaseCache): A function of the observations of a
            phase that gives its metrics, or a cache of one.
        context (Any): A hashable identifier of the observations, for a
            cache.
        cost_model (CostModel | None): A model of the cost of each history.
        growth (float | Mapping[str, float]): The growth of the cost.

    Returns:
        pandas.DataFrame: a row per history
    """
    rows = []
    keys = {s.name: get_satellite_key(s) for s in satellites}
    for history in histories:
        phases = history.phases(satellites)
        if isinstance(metric, PhaseCache):
            capability = metric.evaluate(phases, observations, keys, context)
        else:
            capability = evaluate_phases(phases, observations, metric)
        row = {**summarize_phases(capability), **history.metrics()}
        if cost_model is not None:
            row.update(cost_model.estimate(history, growth).metrics())
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_missions(
    missions: pd.DataFrame,
    quantiles: Mapping[str, float | Sequence[float]] | None = None,
    requirements: Sequence[Requirement] = (),
    controls: Mapping[str, float] | None = None,
    spread: bool = False,
) -> dict[str, float]:
    """
    Summarizes histories of a mission (see `evaluate_missions`) as the
    metrics of an evaluation: the mean of each column over the histories
    (by the same name), quantiles of columns over the histories (`<column>
    p<percent>`, such as `"accurate mean p10"` for the 0.1 quantile, the
    risk of the mission), and the share of the histories that meet each
    requirement on a column (`P(<column> <operator> <value>)`).

    With `controls`, the means are estimated with control variates: columns
    whose expected values are known (such as the `availability` of a
    mission, which simulations alone estimate precisely and cheaply, see
    `Deployment.expected_metrics`), and which vary with the other columns
    from history to history. The mean of each other column is corrected by
    the deviation of the controls' means from their expected values, times
    the regression of the column on the controls over the histories,
    which removes the share of its variance that the controls explain (a
    share of `rho ** 2` for a correlation `rho`): fewer histories give the
    same precision. The quantiles and probabilities are not corrected.

    The quantiles and probabilities are over the histories of a replicate,
    whose statistics over the replicates (such as their mean, see
    `aggregate`) combine the risk of the mission with the uncertainty of the
    replicates' other parameters (such as their targets).

    Args:
        missions (pandas.DataFrame): A row per history.
        quantiles (Mapping | None): The quantiles of columns, by column.
        requirements (Sequence[Requirement]): Requirements on columns.
        controls (Mapping[str, float] | None): The expected values of
            columns that are control variates, if any.
        spread (bool): `True` to add the standard deviation of each column
            over the histories (`<column> sd`), as for `allocate_histories`.

    Returns:
        dict[str, float]: the summary
    """
    numeric = missions.select_dtypes("number")
    summary = {column: float(numeric[column].mean()) for column in numeric.columns}
    if controls:
        summary.update(_controlled_means(numeric, controls))
    if spread:
        summary.update(
            {f"{column} sd": float(numeric[column].std()) for column in numeric.columns}
        )
    for column, values in (quantiles or {}).items():
        for q in [values] if isinstance(values, (int, float)) else values:
            summary[f"{column} p{100 * q:g}"] = float(numeric[column].quantile(q))
    for requirement in requirements:
        met = _OPERATORS[requirement.operator](
            missions[requirement.name].to_numpy(dtype=float), requirement.value
        )
        statement = f"{requirement.name} {requirement.operator} {requirement.value:g}"
        summary[f"P({statement})"] = float(np.mean(met))
    return summary


def _controlled_means(numeric: pd.DataFrame, controls: Mapping[str, float]) -> dict:
    """
    Estimates the means of columns with control variates (see
    `summarize_missions`): each column's mean less the regression of the
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
        data = data[np.all(np.isfinite(data), axis=1)]
        if len(data) <= len(names) + 1:
            continue
        y, x = data[:, 0], data[:, 1:]
        centered = x - x.mean(axis=0)
        if not np.any(centered):
            continue
        beta = np.linalg.lstsq(centered, y - y.mean(), rcond=None)[0]
        means[column] = float(y.mean() - (x.mean(axis=0) - expected) @ beta)
    return means


class _Sampler:
    """
    The uniform numbers of the failures of a history of a replicate,
    stratified over the histories of the replicate (see
    `Deployment.simulate_many`): a Latin hypercube in which, for each draw
    (each life and launch), history `k` of `K` takes a uniform number in its
    stratum of a permutation of the `K` strata of [0, 1). Each draw is keyed
    by its key, so that the same draw of the same history of a replicate
    takes the same number in every design (common random numbers).
    """

    def __init__(self, seed: int, index: int, histories: int):
        self.seed, self.index, self.histories = seed, index, histories

    def __call__(self, key: tuple, dimensions: int) -> np.ndarray:
        strata = np.random.default_rng([self.seed, 1, self.histories, *key])
        orders = [strata.permutation(self.histories) for _ in range(dimensions)]
        jitter = np.random.default_rng([self.seed, 2, self.index, *key])
        offsets = jitter.random(dimensions)
        return np.array(
            [
                (order[self.index] + offset) / self.histories
                for order, offset in zip(orders, offsets)
            ]
        )


SAMPLINGS = ("random", "stratified")
"""The samplings of the failures of histories (see `Deployment.simulate_many`)."""


def _simulate_many(
    deployment: Deployment,
    failures: FailureModel,
    policy: RestorationPolicy | None,
    seed: int,
    histories: int,
    sampling: str,
) -> list[MissionHistory]:
    """
    Simulates the histories of a replicate (see `Deployment.simulate_many`).
    """
    if sampling not in SAMPLINGS:
        raise ValueError(f"sampling must be one of {SAMPLINGS}, not {sampling!r}")
    if sampling == "random":
        return [
            deployment.simulate(failures, policy, (seed, k)) for k in range(histories)
        ]
    return [
        _Simulation(
            deployment,
            failures,
            policy or RestorationPolicy(),
            (seed, k),
            _Sampler(seed, k, histories),
        ).run()
        for k in range(histories)
    ]


def expected_metrics(
    deployment: Deployment,
    failures: FailureModel,
    policy: RestorationPolicy | None = None,
    histories: int = 1000,
    seed: int = 0,
    sampling: str = "stratified",
) -> dict[str, float]:
    """
    Estimates the expected measures of a mission (see
    `MissionHistory.metrics`, such as its `availability`) from many
    simulated histories, without evaluating their phases: precise and
    cheap, as the expected values of control variates (see
    `summarize_missions`).

    Args:
        deployment (Deployment): The deployment.
        failures (FailureModel): The failures.
        policy (RestorationPolicy | None): The policy for failed satellites.
        histories (int): The number of histories.
        seed (int): The seed of the histories.
        sampling (str): The sampling of their failures.

    Returns:
        dict[str, float]: the mean of each measure
    """
    rows = pd.DataFrame(
        [
            history.metrics()
            for history in _simulate_many(
                deployment, failures, policy, seed, histories, sampling
            )
        ]
    )
    return {column: float(rows[column].mean()) for column in rows.columns}


def allocate_histories(
    records: pd.DataFrame,
    metric: str,
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
    `PhaseCache`), from a pilot study: the records of its replicates (see
    `summarize_missions`, with `spread`), each the mean (`<metric>`) and
    standard deviation (`<metric> sd`) of a metric over `histories`
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
            `within`, and the cost of a replicate (`replicate cost`, seconds)
            and of the evaluations (`cost`, seconds, with a standard error)
    """
    if len(records) < 2:
        raise ValueError("a pilot needs at least two replicates")
    if histories < 2 or replicate_seconds <= 0 or history_seconds <= 0:
        raise ValueError("a pilot needs two or more histories and positive costs")
    means = records[metric].to_numpy(dtype=float)
    within = float(np.mean(records[f"{metric} sd"].to_numpy(dtype=float) ** 2))
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
        "replicate cost": cost,
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
            result.update(histories=float(count), **{"replicate cost": cost})
        result.update(replicates=float(replicates), cost=replicates * cost)
    return result
