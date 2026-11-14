"""
Many histories of missions in a replicate (see `Deployment.simulate_many`):
their simulation (independent or stratified), evaluation with a cache of
the metrics of their phases (see `PhaseCache`), and summaries as the
records of evaluations, with control variates and the allocation of
evaluations between replicates and histories; and the exact expectations
of missions without restoration (`exact_expectations`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import itertools
import math
import warnings
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from .criteria import Requirement
from .failures import FailureModel, RestorationPolicy
from .mission import (
    PHASE_COLUMNS,
    Capability,
    MissionHistory,
    Phase,
    _position_names,
    summarize_phases,
)
from .simulation import _Simulation

if TYPE_CHECKING:
    from .deployment import Deployment


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
    quantiles: Mapping[str, float | Sequence[float]] | None = None,
    requirements: Sequence[Requirement] = (),
    controls: Mapping[str, float] | None = None,
    spread: bool = False,
) -> dict[str, float]:
    """
    Summarizes histories of a mission (see `evaluate_histories`) as the
    metrics of an evaluation: the mean of each column over the histories
    (by the same name), quantiles of columns over the histories (`<column>
    histories p<percent>`, such as `"accurate mission histories p10"` for the
    0.1 quantile, the risk of the mission), and the share of the histories
    that meet each requirement on a column (`P(<column> <operator>
    <value>)`). The quantiles (and standard deviations) over the histories
    of a replicate are named apart from statistics over the replicates
    (such as a `Quantile`, see `Criteria`), which are taken of these.

    Missing values (`NaN`, such as the `restoration years` of a history
    without restorations) are left out of the statistics of their column,
    but infinite ones are not: the mean of the time until a history meets a
    requirement (`<metric> time`, see `evaluate_histories`) is infinite if a
    history never meets it, and its quantiles, or the share of the
    histories that meet it by a time (a requirement on it, such as
    `Requirement("accurate time", "<=", 4)`), summarize it.

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

    Histories of stressed failures (see `Deployment.simulate_many`) are
    weighted by their likelihood ratios (`weight`), for every statistic, and
    the effective number of histories (`effective histories`, the square of
    the sum of the weights over the sum of their squares) measures how much
    a few weights dominate.

    The quantiles and probabilities are over the histories of a replicate,
    whose statistics over the replicates (such as their mean, see
    `aggregate`) combine the risk of the mission with the uncertainty of the
    replicates' other parameters (such as their targets).

    Args:
        histories (pandas.DataFrame): A row per history.
        quantiles (Mapping | None): The quantiles of columns, by column.
        requirements (Sequence[Requirement]): Requirements on columns.
        controls (Mapping[str, float] | None): The expected values of
            columns that are control variates, if any.
        spread (bool): `True` to add the standard deviation of each column
            over the histories (`<column> histories std`), as for
            `allocate_histories`.

    Returns:
        dict[str, float]: the summary

    Raises:
        ValueError: if a quantile or requirement is of a column that the
            histories do not have.
    """
    numeric = histories.select_dtypes("number")
    named = [*(quantiles or {}), *(requirement.name for requirement in requirements)]
    unknown = sorted(set(named) - set(numeric.columns))
    if unknown:
        raise ValueError(
            "quantiles or requirements of columns that the histories do not have: "
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
    if spread:
        for column in numeric.columns:
            values = numeric[column].to_numpy(dtype=float)
            known = np.isfinite(values)
            mean = _weighted_mean(values, weights)
            variance = _weighted_mean((values - mean) ** 2, weights)
            count = known.sum()
            summary[f"{column} histories std"] = (
                float(np.sqrt(variance * count / (count - 1)))
                if count > 1
                else math.nan
            )
    for column, values in (quantiles or {}).items():
        for q in [values] if isinstance(values, (int, float)) else values:
            summary[f"{column} histories p{100 * q:g}"] = _weighted_quantile(
                numeric[column].to_numpy(dtype=float), weights, q
            )
    for requirement in requirements:
        met = requirement.met(histories[requirement.name].to_numpy(dtype=float))
        summary[f"P({requirement})"] = float(np.sum(weights * met) / np.sum(weights))
    return summary


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
    model: FailureModel,
    policy: RestorationPolicy | None,
    seed: int,
    histories: int,
    sampling: str,
    stress: float = 1.0,
    launch_stress: float | None = None,
) -> list[MissionHistory]:
    """
    Simulates the histories of a replicate (see `Deployment.simulate_many`).
    """
    if sampling not in SAMPLINGS:
        raise ValueError(f"sampling must be one of {SAMPLINGS}, not {sampling!r}")
    if stress <= 0 or (launch_stress is not None and launch_stress <= 0):
        raise ValueError("stress must be positive")
    launch_stress = stress if launch_stress is None else launch_stress
    return [
        _Simulation(
            deployment,
            model,
            policy or RestorationPolicy(),
            (seed, k),
            _Sampler(seed, k, histories) if sampling == "stratified" else None,
            stress,
            launch_stress,
        ).run()
        for k in range(histories)
    ]


def _expected_metrics(
    deployment: Deployment,
    model: FailureModel,
    policy: RestorationPolicy | None,
    seed: int,
    histories: int,
    sampling: str,
) -> dict[str, float]:
    """
    Estimates the expected measures of a mission from many simulated
    histories (see `Deployment.expected_metrics`).
    """
    rows = pd.DataFrame(
        [
            history.metrics()
            for history in _simulate_many(
                deployment, model, policy, seed, histories, sampling
            )
        ]
    )
    return {column: float(rows[column].mean()) for column in rows.columns}


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
    `PhaseCache`), from a pilot study: the records of its replicates (see
    `summarize_histories`, with `spread`), each the mean (`<metric>`) and
    standard deviation (`<metric> histories std`) of a metric over `histories`
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


def _launch_groups(
    deployment: Deployment, model: FailureModel, times: np.ndarray
) -> tuple[int, list[tuple]]:
    """
    Gets the launch groups of a deployment without restoration: the
    positions of each launch of the deployment (in order of the planes'
    positions), which its planned replacements carry again. For each group,
    at each time: the probability that each of its launches so far is the
    latest that succeeded (`weights`, launches by times), the survival of
    the satellites of each since its launch (`survival`, launches by times),
    the probability that none has succeeded (`none`), and whether the group
    has launched.
    """
    offsets = np.cumsum([0] + [plane.satellites for plane in deployment.planes])
    launches: dict[tuple[int, int], list] = {}
    for launch in sorted(deployment.launches, key=lambda l: l.time):
        launches.setdefault((launch.plane, launch.first), []).append(launch)
    p = model.launch
    groups = []
    for (plane, first), sequence in launches.items():
        built = deployment.planes[plane].satellites
        members = [
            int(offsets[plane]) + position
            for position in range(first, min(first + sequence[0].satellites, built))
        ]
        if not members:
            continue
        life = model.life_of(sequence[0].kind)
        flown = np.array([[t >= launch.time for t in times] for launch in sequence])
        # launches after each one, up to each time
        later = np.flip(np.cumsum(np.flip(flown, axis=0), axis=0), axis=0) - flown
        weights = np.where(flown, (1 - p) * p**later, 0.0)
        survival = np.array(
            [
                (
                    np.ones(len(times))
                    if life is None
                    else life.survival(np.maximum(times - launch.time, 0.0))
                )
                for launch in sequence
            ]
        )
        none = p ** flown.sum(axis=0)
        groups.append((members, weights, survival, none, flown.any(axis=0)))
    return int(offsets[-1]), groups


def _group_operating(group: tuple, count: int) -> np.ndarray:
    """
    Gets the probability at each time that exactly a given set of `count`
    satellites of a launch group operates (see `_set_weight`).
    """
    members, weights, survival, none, _ = group
    probability = np.sum(
        weights * survival**count * (1 - survival) ** (len(members) - count), axis=0
    )
    return probability + none if count == 0 else probability


def _set_weight(active: np.ndarray, groups: list) -> float:
    """
    Gets the share of the mission's time in which exactly the active
    positions operate, in expectation: at each time, the product over the
    launch groups of the probability that exactly the active satellites of
    the group operate, which (as they are of one type and launched together)
    depends only on their number `a` of the `n` of the group: the sum over
    its launches of the probability that the launch is the latest that
    succeeded times `S ** a * (1 - S) ** (n - a)`, for the survival `S` of
    its satellites, and the probability that none has succeeded if `a` is 0.
    """
    probability = 1.0
    for group in groups:
        probability = probability * _group_operating(group, int(active[group[0]].sum()))
    return float(np.mean(probability))


def _missing_levels(groups: list, count: int, tolerance: float) -> int:
    """
    Gets the most positions missing of the sets that exact expectations
    enumerate (see `exact_expectations`): the least number such that the
    sets of that many missing positions or fewer operate all but `tolerance`
    of the mission's time, in expectation. The distribution of the number
    missing at each time is the convolution over the launch groups of the
    distributions of theirs (each the probability of a set times the number
    of sets of its size), cheap to compute without enumerating the sets.
    """
    if not groups:
        return count
    distribution = np.zeros((count + 1, len(groups[0][3])))
    distribution[0] = 1.0
    for group in groups:
        size = len(group[0])
        combined = np.zeros_like(distribution)
        for gone in range(size + 1):
            share = math.comb(size, gone) * _group_operating(group, size - gone)
            combined[gone:] += share * distribution[: count + 1 - gone]
        distribution = combined
    covered = np.cumsum(distribution.mean(axis=1))
    reached = np.flatnonzero(covered >= 1 - tolerance - 1e-12)
    return int(reached[0]) if len(reached) > 0 else count


def _exact_sets(
    deployment: Deployment,
    model: FailureModel,
    tolerance: float = 1e-3,
    steps: int = 400,
) -> int:
    """
    Gets the number of sets of satellites that exact expectations of a
    deployment enumerate (see `exact_expectations`), without enumerating
    them: the sets of at most the number of positions missing that cover
    all but `tolerance` of the mission's time.
    """
    times = _times(deployment, steps)
    count, groups = _launch_groups(deployment, model, times)
    levels = _missing_levels(groups, count, tolerance)
    return sum(math.comb(count, missing) for missing in range(levels + 1))


def _times(deployment: Deployment, steps: int) -> np.ndarray:
    """The times at which exact expectations integrate probabilities."""
    begin, end = deployment.begin, deployment.end
    return begin + (np.arange(steps) + 0.5) * (end - begin) / steps


MAX_SETS = 200_000
"""The most sets of satellites that exact expectations enumerate by default
(see `exact_expectations`), as those of the failures of a study do (see
`Failures`)."""


def exact_expectations(
    deployment: Deployment,
    model: FailureModel,
    capability: Capability,
    *,
    tolerance: float = 1e-3,
    steps: int = 400,
    max_sets: int = MAX_SETS,
) -> dict[str, float]:
    """
    Computes the expected capability of a mission without restoration
    (failed satellites are not restored, but may be replaced as planned at
    the end of their design life) over its failures: the expectation of the
    mean of each metric over the phases of the mission (`<metric> mission`, as
    `summarize_phases` gives for one history) and of its availability
    (`availability`), exactly rather than from samples of histories (as
    `Deployment.expected_metrics` estimates the measures of missions).

    Without restoration, the satellites fail independently (but for the
    satellites of a launch, which a failed launch loses together), so that
    the probability that a set of satellites operates at a time is a product
    over the launches of their survival (see `Lifetime.survival`). The
    expected mean of a metric is the sum over the sets of satellites of the
    metric of each set (evaluated once, as with a `PhaseCache`) times the
    share of the mission's time in which exactly that set operates, in
    expectation (integrated over `steps` times). The sets are enumerated by
    the number of satellites missing until their shares cover all but
    `tolerance` of the mission's time (a number known before they are
    enumerated, from the distribution of the number missing), and evaluated
    in order of their shares until they do (`coverage`, the share covered,
    and `sets`, the sets evaluated); the expectation is over the sets
    evaluated.

    Without sampling error, an expectation estimates the mean of many
    histories with none of their variance (see `Deployment.simulate_many`),
    for missions of a few to a dozen or two satellites; for larger ones
    (whose sets number more than `max_sets`, which is known before any is
    enumerated), or with restoration (whose shared production and launches
    make failures dependent), simulate histories.

    Args:
        deployment (Deployment): The deployment.
        model (FailureModel): The failures.
        capability (Capability): The capability of the satellites (in the
            order of the planes, see `MissionHistory.phases`): their
            observations and metric (see `Capability.evaluate`).
        tolerance (float): The share of the mission's time that the sets
            evaluated may leave out.
        steps (int): The times at which the probabilities are integrated.
        max_sets (int): The most sets to enumerate.

    Returns:
        dict[str, float]: the expected metrics

    Raises:
        ValueError: if the deployment is staged (see `Stage`), whose launches
            depend on its rules, or has more sets to enumerate than
            `max_sets` (raised before any is enumerated).
    """
    # pylint: disable=too-many-locals
    if deployment.stages:
        raise ValueError(
            "exact expectations need a deployment that is not staged: "
            + "simulate histories of a staged deployment instead"
        )
    names = list(_position_names(deployment.planes, capability.satellites).values())
    begin, end = deployment.begin, deployment.end
    count, groups = _launch_groups(deployment, model, _times(deployment, steps))
    levels = _missing_levels(groups, count, tolerance)
    if sum(math.comb(count, missing) for missing in range(levels + 1)) > max_sets:
        raise ValueError("too many sets of satellites to enumerate: simulate histories")
    # availability: the expected years of operation of the positions, as a
    # share of those of the planned deployment
    operating = sum(
        len(members) * float(np.mean(np.sum(weights * survival, axis=0)))
        for members, weights, survival, _, _ in groups
    )
    planned = sum(
        len(members) * float(np.mean(launched)) for members, _, _, _, launched in groups
    )
    availability = operating / planned if planned > 0 else math.nan
    # enumerate the sets by the number of positions missing
    weights: list[tuple[float, tuple[int, ...]]] = []
    for missing in range(levels + 1):
        for gone in itertools.combinations(range(count), missing):
            active = np.ones(count, dtype=bool)
            active[list(gone)] = False
            weight = _set_weight(active, groups)
            if weight > 0:
                weights.append((weight, tuple(np.flatnonzero(active))))
    weights.sort(key=lambda item: -item[0])
    totals: dict[str, float] = {}
    evaluated, sets = 0.0, 0
    for weight, members in weights:
        phase = Phase(begin, end, tuple(names[k] for k in members))
        values = capability.evaluate([phase]).iloc[0]
        for name, value in values.items():
            if name in PHASE_COLUMNS:
                continue
            totals[f"{name} mission"] = totals.get(
                f"{name} mission", 0.0
            ) + weight * float(value)
        evaluated += weight
        sets += 1
        if evaluated >= 1 - tolerance:
            break
    result = {name: total / evaluated for name, total in totals.items()}
    result.update(availability=availability, coverage=evaluated, sets=float(sets))
    return result
