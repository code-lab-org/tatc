"""
Pareto efficiency of designs with several objectives: non-dominated
designs, Pareto ranks, and the hypervolume of a set of designs (and of the
designs of a search as it proceeds).

Designs are judged by criteria (see `Criteria`): objectives (columns of a
data frame, such as the statistics of a summary of a study, see `aggregate`)
to minimize or maximize, requirements that efficient designs meet, the
tolerances of the objectives, and the groups in which designs are judged
separately. A mapping of objectives to their directions (`"min"` or
`"max"`) is the criteria of these objectives alone. A design with a missing
objective (`NaN`, such as a statistic of a metric that none of its feasible
replicates gives) is not efficient, and dominates no other.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .criteria import (
    Criteria,
    Requirement,
    _as_criteria,
    _bootstrap_error,
    _check_summary,
    _combiners,
    _costs_of,
    _get_tolerances,
    _groups,
    _resamples,
    _resolve,
    _signs,
    _summary_groups,
    aggregate,
)
from .dominance import _dominated, _efficient


def _get_costs(frame: pd.DataFrame, objectives: Mapping[str, str]) -> np.ndarray:
    """
    Gets the objectives of a data frame as costs (to minimize): maximized
    objectives are negated (missing values are `NaN`).
    """
    missing = [column for column in objectives if column not in frame.columns]
    if missing:
        raise ValueError(f"no columns of the objectives {missing}")
    columns = []
    for column, direction in objectives.items():
        values = frame[column].to_numpy(dtype=float)
        columns.append(values if direction == "min" else -values)
    return np.column_stack(columns) if len(frame) > 0 else np.empty((0, len(columns)))


def _get_violation(
    frame: pd.DataFrame, requirements: Sequence[Requirement]
) -> np.ndarray:
    """
    Gets the total violation of requirements by the rows of a data frame,
    from the columns of their statistics (infinite where a statistic is
    missing).
    """
    violation = np.zeros(len(frame))
    for requirement in requirements:
        if requirement.name not in frame.columns:
            raise ValueError(
                f"no column of the statistic {requirement.name!r} of a requirement: "
                + "aggregate the results with the criteria (see aggregate)"
            )
        violation = violation + requirement.violation(
            frame[requirement.name].to_numpy(dtype=float)
        )
    return np.where(np.isnan(violation), np.inf, violation)


def pareto_front(
    summary: pd.DataFrame, criteria: Criteria | Mapping[str, str]
) -> pd.Series:
    """
    Finds the Pareto-efficient (non-dominated) rows of a summary of designs
    (a row per design, such as of `aggregate`, not the results of a study,
    which have a row per replicate): those that meet
    the requirements of the criteria (by the columns of their statistics)
    and that no other row that meets them dominates, being no worse in every
    objective and better in at least one, in each group of the criteria's
    `by`. Rows with the same objectives are all efficient if any is, and a
    row with a missing objective (`NaN`) is not. Summaries of criteria give
    whether each design is efficient already (`pareto`, see `aggregate`).

    With a tolerance of an objective (see `Criteria`), a row dominates
    another only if it would still dominate it if it were worse by the
    tolerance in that objective: a value in the objective's units, or a
    share of the other row's value (a percentage, such as `"5%"`), which
    suits objectives that range over orders of magnitude. Rows within the
    tolerances of the front (that no other row is clearly better than) are
    then efficient as well, which suits objectives estimated with an
    uncertainty (such as a few standard errors of their means, see
    `summarize`), or by a cheaper analysis (see `screen`).

    Args:
        summary (pandas.DataFrame): The rows (such as designs).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).

    Returns:
        pandas.Series: whether each row is efficient (boolean, with the
            summary's index)

    Raises:
        ValueError: if the summary has several rows of a design (in a group
            of the criteria's `by`), as the results of a study do.
    """
    criteria = _as_criteria(criteria)
    _check_summary(summary, criteria.by)
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerance)
    efficient = np.zeros(len(summary), dtype=bool)
    for _, positions in _groups(summary, criteria.by):
        group = summary.iloc[positions]
        meets = _get_violation(group, criteria.requirements) == 0
        costs = _get_costs(group, criteria.objectives)
        efficient[positions] = _efficient(costs, meets, tolerances)
    return pd.Series(efficient, index=summary.index, name="pareto")


def pareto_ranks(
    summary: pd.DataFrame, criteria: Criteria | Mapping[str, str]
) -> pd.Series:
    """
    Ranks the rows of a summary of designs (a row per design, see
    `pareto_front`) by non-dominated sorting, in each group
    of the criteria's `by`: rank 0 is the Pareto front (see `pareto_front`),
    rank 1 the front of the others, and so on, with the tolerances of the
    criteria. Rows that miss the requirements of the criteria are ranked
    after those that meet them, in order of their violation (the sum of the
    violations of the requirements, see `Requirement`), and rows with a
    missing objective (`NaN`) that meet them last.

    Args:
        summary (pandas.DataFrame): The rows (such as designs).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).

    Returns:
        pandas.Series: the rank of each row (with the summary's index)

    Raises:
        ValueError: if the summary has several rows of a design (in a group
            of the criteria's `by`), as the results of a study do.
    """
    criteria = _as_criteria(criteria)
    _check_summary(summary, criteria.by)
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerance)
    ranks = np.full(len(summary), -1, dtype=int)
    for _, positions in _groups(summary, criteria.by):
        group = summary.iloc[positions]
        violation = _get_violation(group, criteria.requirements)
        costs = _get_costs(group, criteria.objectives)
        known = ~np.isnan(costs).any(axis=1)
        meets = (violation == 0) & known
        remaining, costs = positions[meets], costs[meets]
        rank = 0
        while len(remaining) > 0:
            front = ~_dominated(costs, tolerances)
            ranks[remaining[front]] = rank
            remaining, costs = remaining[~front], costs[~front]
            rank += 1
        missing = violation != 0
        if missing.any():
            order = pd.Series(violation[missing]).rank(method="dense").to_numpy()
            ranks[positions[missing]] = rank + order.astype(int) - 1
            rank = int(ranks[positions[missing]].max()) + 1
        unknown = (violation == 0) & ~known
        ranks[positions[unknown]] = rank
    return pd.Series(ranks, index=summary.index, name="pareto_rank")


def _hypervolume(points: np.ndarray, reference: np.ndarray) -> float:
    """
    Computes the hypervolume (to minimize) of points that are all better than
    the reference in every dimension, by slicing the last dimension (exact,
    for few dimensions).
    """
    if len(points) == 0:
        return 0.0
    if points.shape[1] == 1:
        return float(reference[0] - points[:, 0].min())
    points = points[~_dominated(points)]
    order = np.argsort(points[:, -1], kind="stable")
    points = points[order]
    if points.shape[1] == 2:
        # sweep: by increasing second cost, the first cost decreases
        x = points[:, 0]
        upper = np.append(points[1:, 1], reference[1])
        return float(np.sum((reference[0] - x) * (upper - points[:, 1])))
    volume = 0.0
    for i in range(len(points)):
        upper = points[i + 1, -1] if i + 1 < len(points) else reference[-1]
        height = upper - points[i, -1]
        if height > 0:
            volume += height * _hypervolume(points[: i + 1, :-1], reference[:-1])
    return float(volume)


def hypervolume(
    summary: pd.DataFrame,
    criteria: Criteria | Mapping[str, str],
    reference: Mapping[str, float] | Sequence[float],
) -> float:
    """
    Computes the hypervolume of the rows of a summary of designs (a row per
    design, see `pareto_front`): the volume of the
    objective space that they dominate, up to a reference point that is
    worse than every row of interest in every objective. Rows that are not
    better than the reference in every objective, that miss the requirements
    of the criteria, or that have a missing objective add nothing (the
    tolerances of the criteria do not apply). A larger hypervolume is a better set of designs,
    which makes it a measure of the progress of a search, or of how close a
    search comes to the true Pareto front.

    The hypervolume is exact. Its computation slices the objectives in turn,
    which takes long with many objectives and many efficient rows.

    Args:
        summary (pandas.DataFrame): The rows (such as designs), of one group
            (criteria with `by` are not accepted).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).
        reference (Mapping[str, float] | Sequence[float]): The reference
            point: the value of each objective (by column, or in order).

    Returns:
        float: the hypervolume, in the product of the units of the objectives

    Raises:
        ValueError: if the summary has several rows of a design, as the
            results of a study do (see `hypervolume_progress` for those of a
            search).
    """
    criteria = _as_criteria(criteria)
    if criteria.by:
        raise ValueError(
            "a hypervolume is of the rows of one group: give criteria without `by`"
        )
    _check_summary(summary)
    objectives = criteria.objectives
    if isinstance(reference, Mapping):
        reference = [reference[column] for column in objectives]
    if len(reference) != len(objectives):
        raise ValueError("the reference point needs a value for each objective")
    bound = _signs(criteria) * np.asarray(reference, dtype=float)
    if len(summary) == 0:
        return 0.0
    summary = summary[_get_violation(summary, criteria.requirements) == 0]
    costs = _get_costs(summary, objectives)
    costs = costs[np.all(costs < bound, axis=1)]
    return _hypervolume(costs, bound)


def hypervolume_progress(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str],
    reference: Mapping[str, float] | Sequence[float],
    *,
    step: int = 10,
) -> pd.Series:
    """
    Gets the hypervolume of the designs evaluated by a search (such as an
    optimizer, see `optimize_optuna`) as it proceeds: of the first designs
    evaluated (in the order of their first records in the results of the
    search), judged by their statistics over their replicates (see
    `aggregate`), in steps of `step` designs (and of all of them).

    Args:
        results (pandas.DataFrame): The results of the search (see
            `Study.run`).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).
        reference (Mapping[str, float] | Sequence[float]): The reference
            point (see `hypervolume`).
        step (int): The number of designs of each step.

    Returns:
        pandas.Series: the hypervolume, by the number of designs evaluated
    """
    if step < 1:
        raise ValueError("step must be at least 1")
    criteria = _as_criteria(criteria)
    summary = aggregate(results, criteria).set_index("design")
    order = results["design"].drop_duplicates().tolist()
    counts = list(range(step, len(order) + 1, step))
    if order and (not counts or counts[-1] != len(order)):
        counts.append(len(order))
    return pd.Series(
        [
            hypervolume(summary[summary.index.isin(order[:n])], criteria, reference)
            for n in counts
        ],
        index=pd.Index(counts, name="designs"),
        name="hypervolume",
        dtype=float,
    )


def pareto_probability(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str],
    *,
    n_resamples: int = 500,
    seed: int = 0,
    design_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Estimates the probability that each design is Pareto-efficient by
    criteria (and meets their requirements, if any), given the uncertainty
    of its replicates, by a bootstrap: each resample draws replicates with
    replacement, the same replicates for every design (which keeps the
    correlation between designs of common random numbers, see
    `Design.realize`), computes the statistics of each design over the drawn
    replicates, and finds the front of the designs that meet the
    requirements (in each group of the criteria's `by`). The probability is
    the share of resamples in which a design is on that front.

    Designs that are efficient with a probability between 0 and 1 are those
    whose efficiency the replicates do not settle: more replicates (see
    `refine`) would settle it, unless designs differ by less than matters (as
    near saturation), which a tolerance of the objectives settles: designs
    within the tolerances of the front are efficient (see `pareto_front`).
    With few replicates, the bootstrap has few distinct resamples and
    understates the uncertainty. Designs may have different replicates (as
    after `refine`): replicates are drawn from all of them, and each design's
    statistics are those of its drawn replicates (or of all of its
    replicates, if none is drawn). With criteria `across` the scenarios of a
    study (see `Across`), a drawn replicate is drawn in every scenario, in
    which it has the same uncertain parameters.

    Args:
        results (pandas.DataFrame): The results of a study (see
            `Study.run`), with a record per replicate of each design (and
            group of the criteria's `by`). Infeasible results are ignored.
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (metrics whose means they are) and their directions
            (`"min"` or `"max"`).
        n_resamples (int): The number of resamples (at least 2).
        seed (int): The seed of the resamples.
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every
            replicate of a design (see `aggregate`).

    Returns:
        pandas.DataFrame: a row per design (and group): the summary of the
            criteria over every replicate (see `aggregate`: the values of the
            design, the statistics, the number of `replicates`, the
            `violation` of the requirements, whether the design meets them,
            `meets_requirements`, and whether it is efficient among those
            that do, `pareto`), the standard error of each objective
            (`<objective> se`, its standard deviation over the resamples),
            the probability that it is efficient (`pareto_probability`), and
            the probability that it meets the requirements
            (`requirements_probability`)
    """
    criteria = _as_criteria(criteria)
    resolved = _resolve(criteria)
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerance)
    combiners = None if criteria.across is None else _combiners(criteria)
    rng = np.random.default_rng(seed)
    frames = []
    for frame, table in _summary_groups(
        results,
        resolved,
        criteria.by,
        criteria.requirements,
        criteria.across,
        combiners,
        design_columns,
    ):
        values, violation = table.evaluate()
        counts = np.zeros(len(table.designs))
        met = np.zeros(len(table.designs))
        # the sums of the objectives (and of their squares) over the resamples
        sums = {name: np.zeros(len(table.designs)) for name in criteria.objectives}
        squares = {name: np.zeros(len(table.designs)) for name in criteria.objectives}
        for drawn in _resamples(len(table.replicates), n_resamples, rng):
            sample, sample_violation = table.evaluate(drawn)
            # designs without a drawn replicate keep the statistics of all of theirs
            for name in sample:
                missing = np.isnan(sample[name]) & ~np.isnan(values[name])
                sample[name][missing] = values[name][missing]
            missing = np.isnan(sample_violation)
            sample_violation[missing] = violation[missing]
            meets = sample_violation == 0
            counts += _efficient(_costs_of(sample, criteria), meets, tolerances)
            met += meets
            for name in criteria.objectives:
                sums[name] += sample[name]
                squares[name] += sample[name] ** 2
        for name, value in values.items():
            frame[name] = value
        for name in criteria.objectives:
            frame[f"{name} se"] = _bootstrap_error(
                sums[name], squares[name], n_resamples
            )
        meets = violation == 0
        frame["replicates"] = table.counts
        frame["violation"] = violation
        frame["meets_requirements"] = meets
        frame["pareto"] = _efficient(_costs_of(values, criteria), meets, tolerances)
        frame["pareto_probability"] = counts / n_resamples
        frame["requirements_probability"] = met / n_resamples
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)
