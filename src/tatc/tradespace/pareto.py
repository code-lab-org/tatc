"""
Pareto efficiency of designs with several objectives: non-dominated
designs, Pareto ranks, and the hypervolume of a set of designs (and of the
designs of a search as it proceeds).

Designs are judged by criteria (see `Criteria`): objectives (columns of a
data frame, such as the statistics of a summary of a study, see `aggregate`)
to minimize or maximize, requirements that efficient designs meet, the
tolerances of the objectives, and the groups in which designs are judged
separately. A mapping of objectives to their directions (`"min"` or
`"max"`) is the criteria of these objectives alone.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .criteria import (
    Criteria,
    Requirement,
    _across_groups,
    _as_criteria,
    _combiners,
    _get_tolerances,
    _groups,
    _resolve,
    _signs,
    _Statistics,
    aggregate,
)


def _get_costs(frame: pd.DataFrame, objectives: Mapping[str, str]) -> np.ndarray:
    """
    Gets the objectives of a data frame as costs (to minimize): maximized
    objectives are negated.
    """
    missing = [column for column in objectives if column not in frame.columns]
    if missing:
        raise ValueError(f"no columns of the objectives {missing}")
    columns = []
    for column, direction in objectives.items():
        values = frame[column].to_numpy(dtype=float)
        columns.append(values if direction == "min" else -values)
    costs = np.column_stack(columns) if len(frame) > 0 else np.empty((0, len(columns)))
    if np.isnan(costs).any():
        raise ValueError("objectives must not be NaN")
    return costs


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


_BLOCK_ELEMENTS = 2_000_000
"""The most elements of the temporary arrays of a block of comparisons."""

_CHUNK = 64
"""The points of each chunk of a sweep (see `_dominated_plain`)."""


def _dominated_by(
    points: np.ndarray, others: np.ndarray, tolerances: np.ndarray | None = None
) -> np.ndarray:
    """
    Finds the points dominated by any of others, even if the other is
    penalized by the tolerances (see `_dominated`), comparing them in blocks
    whose temporary arrays are limited to a few million elements.
    """
    dominated = np.zeros(len(points), dtype=bool)
    if len(points) == 0 or len(others) == 0:
        return dominated
    if tolerances is None:
        tolerances = np.zeros((2, points.shape[1]))
    absolute, relative = tolerances
    penalties = np.any(tolerances)
    block = max(1, _BLOCK_ELEMENTS // max(len(others) * points.shape[1], 1))
    for start in range(0, len(points), block):
        part = points[start : start + block]
        # [i, j]: other j dominates point i of the block
        penalized = others[None, :, :]
        if penalties:
            penalized = penalized + absolute + relative * np.abs(part[:, None, :])
        no_worse = np.all(penalized <= part[:, None, :], axis=2)
        better = np.any(penalized < part[:, None, :], axis=2)
        dominated[start : start + block] = np.any(no_worse & better, axis=1)
    return dominated


def _dominated_plain(costs: np.ndarray) -> np.ndarray:
    """
    Finds the points dominated by another (without tolerances) by a sweep in
    lexicographic order, in which a point can only be dominated by points
    before it, and is dominated by one if it is by any of the efficient
    points before it (dominance is transitive). The sweep compares each
    chunk of points with the efficient points before it, then those it does
    not dominate with each other.
    """
    n = len(costs)
    order = np.lexsort(costs.T[::-1])
    ordered = costs[order]
    dominated = np.zeros(n, dtype=bool)
    efficient = np.empty((0, costs.shape[1]))
    for start in range(0, n, _CHUNK):
        part = ordered[start : start + _CHUNK]
        beaten = _dominated_by(part, efficient)
        # a point that an earlier efficient point dominates dominates nothing
        # that it does not: compare only the others with each other
        survivors = part[~beaten]
        beaten[~beaten] = _dominated_by(survivors, survivors)
        dominated[order[start : start + _CHUNK]] = beaten
        efficient = np.vstack([efficient, part[~beaten]])
    return dominated


def _dominated(costs: np.ndarray, tolerances: np.ndarray | None = None) -> np.ndarray:
    """
    Finds the points (rows of costs) dominated by another, even if the other
    is penalized by the tolerances: no worse in every cost, and better in at
    least one, by at least its tolerance (absolute, plus a share of the
    magnitude of the cost of the point it would dominate, see
    `_get_tolerances`).

    A point that is dominated with tolerances is also dominated by a point
    that is efficient without them: the penalty depends only on the point
    dominated, so any point no worse than its dominator in every cost (as
    some efficient point is) dominates it too. Points are therefore compared
    only with the efficient points (see `_dominated_plain`), which are few.
    """
    plain = _dominated_plain(costs)
    if tolerances is None or not np.any(tolerances):
        return plain
    return _dominated_by(costs, costs[~plain], tolerances)


def pareto_front(
    frame: pd.DataFrame, criteria: Criteria | Mapping[str, str]
) -> pd.Series:
    """
    Finds the Pareto-efficient (non-dominated) rows of a data frame (such as
    the designs of a summary of a study, see `aggregate`): those that meet
    the requirements of the criteria (by the columns of their statistics)
    and that no other row that meets them dominates, being no worse in every
    objective and better in at least one, in each group of the criteria's
    `by`. Rows with the same objectives are all efficient if any is.

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
        frame (pandas.DataFrame): The rows (such as designs).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).

    Returns:
        pandas.Series: whether each row is efficient (boolean, with the
            frame's index)
    """
    criteria = _as_criteria(criteria)
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerance)
    efficient = np.zeros(len(frame), dtype=bool)
    for _, positions in _groups(frame, criteria.by):
        group = frame.iloc[positions]
        meets = _get_violation(group, criteria.requirements) == 0
        if meets.any():
            costs = _get_costs(group[meets], criteria.objectives)
            efficient[positions[meets]] = ~_dominated(costs, tolerances)
    return pd.Series(efficient, index=frame.index, name="pareto")


def pareto_ranks(
    frame: pd.DataFrame, criteria: Criteria | Mapping[str, str]
) -> pd.Series:
    """
    Ranks the rows of a data frame by non-dominated sorting, in each group
    of the criteria's `by`: rank 0 is the Pareto front (see `pareto_front`),
    rank 1 the front of the others, and so on, with the tolerances of the
    criteria. Rows that miss the requirements of the criteria are ranked
    after those that meet them, in order of their violation (the sum of the
    violations of the requirements, see `Requirement`).

    Args:
        frame (pandas.DataFrame): The rows (such as designs).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).

    Returns:
        pandas.Series: the rank of each row (with the frame's index)
    """
    criteria = _as_criteria(criteria)
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerance)
    ranks = np.full(len(frame), -1, dtype=int)
    for _, positions in _groups(frame, criteria.by):
        group = frame.iloc[positions]
        violation = _get_violation(group, criteria.requirements)
        meets = violation == 0
        costs = _get_costs(group[meets], criteria.objectives)
        remaining = positions[meets]
        rank = 0
        while len(remaining) > 0:
            front = ~_dominated(costs, tolerances)
            ranks[remaining[front]] = rank
            remaining, costs = remaining[~front], costs[~front]
            rank += 1
        if (~meets).any():
            order = pd.Series(violation[~meets]).rank(method="dense").to_numpy()
            ranks[positions[~meets]] = rank + order.astype(int) - 1
    return pd.Series(ranks, index=frame.index, name="pareto_rank")


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
    frame: pd.DataFrame,
    criteria: Criteria | Mapping[str, str],
    reference: Mapping[str, float] | Sequence[float],
) -> float:
    """
    Computes the hypervolume of the rows of a data frame: the volume of the
    objective space that they dominate, up to a reference point that is
    worse than every row of interest in every objective. Rows that are not
    better than the reference in every objective, or that miss the
    requirements of the criteria, add nothing (the tolerances of the
    criteria do not apply). A larger hypervolume is a better set of designs,
    which makes it a measure of the progress of a search, or of how close a
    search comes to the true Pareto front.

    The hypervolume is exact. Its computation slices the objectives in turn,
    which takes long with many objectives and many efficient rows.

    Args:
        frame (pandas.DataFrame): The rows (such as designs), of one group
            (criteria with `by` are not accepted).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).
        reference (Mapping[str, float] | Sequence[float]): The reference
            point: the value of each objective (by column, or in order).

    Returns:
        float: the hypervolume, in the product of the units of the objectives
    """
    criteria = _as_criteria(criteria)
    if criteria.by:
        raise ValueError(
            "a hypervolume is of the rows of one group: give criteria without `by`"
        )
    objectives = criteria.objectives
    if isinstance(reference, Mapping):
        reference = [reference[column] for column in objectives]
    if len(reference) != len(objectives):
        raise ValueError("the reference point needs a value for each objective")
    bound = _signs(criteria) * np.asarray(reference, dtype=float)
    if len(frame) == 0:
        return 0.0
    frame = frame[_get_violation(frame, criteria.requirements) == 0]
    costs = _get_costs(frame, objectives)
    costs = costs[np.all(costs < bound, axis=1)]
    return _hypervolume(costs, bound)


def hypervolume_progress(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str],
    reference: Mapping[str, float] | Sequence[float],
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
    n_resamples: int = 500,
    seed: int = 0,
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
        n_resamples (int): The number of resamples.
        seed (int): The seed of the resamples.

    Returns:
        pandas.DataFrame: a row per design (and group): the statistics of
            the criteria over every replicate, the standard error of each
            objective (`<objective> se`, its standard deviation over the
            resamples), the number of `replicates`,
            whether the design meets the requirements (`meets_requirements`)
            and is efficient among those that do (`pareto`) by these
            statistics, the probability that it is efficient
            (`pareto_probability`), and the probability that it meets the
            requirements (`requirements_probability`)
    """
    criteria = _as_criteria(criteria)
    resolved = _resolve(criteria)
    requirements = criteria.requirements
    signs = _signs(criteria)
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerance)
    feasible = results[results["feasible"]] if "feasible" in results else results
    combiners = None if criteria.across is None else _combiners(criteria)
    groups = _across_groups(results, criteria.across)
    rng = np.random.default_rng(seed)
    frames = []
    for key, positions in _groups(feasible, criteria.by):
        group = feasible.iloc[positions]
        table = _Statistics(
            group, resolved, requirements, criteria.across, combiners, groups
        )
        designs, replicates = table.designs, table.replicates
        values, violation = table.evaluate()

        def front(values, violation):
            costs = np.column_stack(
                [s * values[name] for s, name in zip(signs, criteria.objectives)]
            )
            ok = (violation == 0) & ~np.isnan(costs).any(axis=1)
            efficient = np.zeros(len(costs), dtype=bool)
            efficient[ok] = ~_dominated(costs[ok], tolerances)
            return efficient, ok

        efficient, meets = front(values, violation)
        counts = np.zeros(len(designs))
        met = np.zeros(len(designs))
        # the sums of the objectives (and of their squares) over the resamples
        sums = {name: np.zeros(len(designs)) for name in criteria.objectives}
        squares = {name: np.zeros(len(designs)) for name in criteria.objectives}
        for _ in range(n_resamples):
            drawn = rng.integers(0, len(replicates), len(replicates))
            sample, sample_violation = table.evaluate(drawn)
            # designs without a drawn replicate keep the statistics of all of theirs
            for name in sample:
                missing = np.isnan(sample[name]) & ~np.isnan(values[name])
                sample[name][missing] = values[name][missing]
            missing = np.isnan(sample_violation)
            sample_violation[missing] = violation[missing]
            sample_front, sample_meets = front(sample, sample_violation)
            counts += sample_front
            met += sample_meets
            for name in criteria.objectives:
                sums[name] += sample[name]
                squares[name] += sample[name] ** 2
        frame = pd.DataFrame({"design": designs})
        for name, value in zip(criteria.by, key):
            frame[name] = value
        for name, value in values.items():
            frame[name] = value
        for name in criteria.objectives:
            mean = sums[name] / n_resamples
            variance = np.maximum(squares[name] / n_resamples - mean**2, 0.0)
            frame[f"{name} se"] = np.sqrt(variance)
        frame["replicates"] = table.counts
        frame["meets_requirements"] = meets
        frame["pareto"] = efficient
        frame["pareto_probability"] = counts / n_resamples
        frame["requirements_probability"] = met / n_resamples
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)
