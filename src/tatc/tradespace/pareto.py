"""
Pareto efficiency of designs with several objectives: non-dominated
designs, Pareto ranks, and the hypervolume of a set of designs (and of the
designs of a search as it proceeds).

Designs are judged by criteria (see `Criteria`): objectives (columns of a
data frame, such as the statistics of a summary of a study, see `summarize`)
to minimize or maximize, requirements that efficient designs meet, the
tolerances of the objectives, and the groups in which designs are judged
separately. A mapping of objectives to their directions (`"min"` or
`"max"`) is the criteria of these objectives alone. A design with a missing
objective (`NaN`, such as a statistic of a metric that none of its feasible
replicates gives) is not efficient, and dominates no other. The probability
that a design is efficient, given the uncertainty of its replicates, is that
of `pareto_probability` (of the summaries of studies, see `summarize`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .criteria import (
    Criteria,
    Objectives,
    _as_criteria,
    _costs,
    _signs,
    _tolerances,
    _violation,
)
from .dominance import dominated, efficient
from .summary import _summarize
from .tables import check_summary, groups


def _get_costs(frame: pd.DataFrame, criteria: Criteria) -> np.ndarray:
    """
    Gets the objectives of criteria of a data frame (of their columns) as
    costs (to minimize, see `Criteria._costs`).
    """
    missing = [column for column in criteria.objectives if column not in frame]
    if missing:
        raise ValueError(f"no columns of the objectives {missing}")
    return _costs(criteria, frame)


def _get_violation(frame: pd.DataFrame, criteria: Criteria) -> np.ndarray:
    """
    Gets the total violation of the requirements of criteria by the rows of a
    data frame, from the columns of their statistics (infinite where a
    statistic is missing).
    """
    missing = [r.name for r in criteria.requirements if r.name not in frame.columns]
    if missing:
        raise ValueError(
            f"no columns of the statistics {missing} of requirements: summarize "
            + "the results with the criteria (see summarize)"
        )
    violation = _violation(criteria, frame, len(frame))
    return np.where(np.isnan(violation), np.inf, violation)


def pareto_front(summary: pd.DataFrame, criteria: Criteria | Objectives) -> pd.Series:
    """
    Finds the Pareto-efficient (non-dominated) rows of a summary of designs
    (a row per design, such as of `summarize`, not the results of a study,
    which have a row per replicate): those that meet
    the requirements of the criteria (by the columns of their statistics)
    and that no other row that meets them dominates, being no worse in every
    objective and better in at least one, in each group of the criteria's
    `by`. Rows with the same objectives are all efficient if any is, and a
    row with a missing objective (`NaN`) is not. Summaries of criteria give
    whether each design is efficient already (`pareto`, see `summarize`).

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
        criteria (Criteria | Objectives): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).

    Returns:
        pandas.Series: whether each row is efficient (boolean, with the
            summary's index)

    Raises:
        ValueError: if the summary has several rows of a design (in a group
            of the criteria's `by`), as the results of a study do.
    """
    criteria = _as_criteria(criteria)
    check_summary(summary, criteria.by)
    front = np.zeros(len(summary), dtype=bool)
    for _, positions in groups(summary, criteria.by):
        group = summary.iloc[positions]
        meets = _get_violation(group, criteria) == 0
        costs = _get_costs(group, criteria)
        front[positions] = efficient(costs, meets, _tolerances(criteria))
    return pd.Series(front, index=summary.index, name="pareto")


def pareto_ranks(summary: pd.DataFrame, criteria: Criteria | Objectives) -> pd.Series:
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
        criteria (Criteria | Objectives): The criteria, or the
            objectives (columns) and their directions (`"min"` or `"max"`).

    Returns:
        pandas.Series: the rank of each row (with the summary's index)

    Raises:
        ValueError: if the summary has several rows of a design (in a group
            of the criteria's `by`), as the results of a study do.
    """
    criteria = _as_criteria(criteria)
    check_summary(summary, criteria.by)
    tolerances = _tolerances(criteria)
    ranks = np.full(len(summary), -1, dtype=int)
    for _, positions in groups(summary, criteria.by):
        group = summary.iloc[positions]
        violation = _get_violation(group, criteria)
        costs = _get_costs(group, criteria)
        known = ~np.isnan(costs).any(axis=1)
        meets = (violation == 0) & known
        remaining, costs = positions[meets], costs[meets]
        rank = 0
        while len(remaining) > 0:
            front = ~dominated(costs, tolerances)
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
    points = points[~dominated(points)]
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
    criteria: Criteria | Objectives,
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
        criteria (Criteria | Objectives): The criteria, or the
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
    check_summary(summary)
    objectives = criteria.objectives
    if isinstance(reference, Mapping):
        reference = [reference[column] for column in objectives]
    if len(reference) != len(objectives):
        raise ValueError("the reference point needs a value for each objective")
    bound = _signs(criteria) * np.asarray(reference, dtype=float)
    if len(summary) == 0:
        return 0.0
    summary = summary[_get_violation(summary, criteria) == 0]
    costs = _get_costs(summary, criteria)
    costs = costs[np.all(costs < bound, axis=1)]
    return _hypervolume(costs, bound)


def hypervolume_progress(
    results: pd.DataFrame,
    criteria: Criteria | Objectives,
    reference: Mapping[str, float] | Sequence[float],
    *,
    step: int = 10,
) -> pd.Series:
    """
    Gets the hypervolume of the designs evaluated by a search (such as an
    optimizer, see `optimize_optuna`) as it proceeds: of the first designs
    evaluated (in the order of their first records in the results of the
    search), judged by their statistics over their replicates (see
    `summarize`), in steps of `step` designs (and of all of them).

    Args:
        results (pandas.DataFrame): The results of the search (see
            `Study.run`).
        criteria (Criteria | Objectives): The criteria, or the
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
    summary = _summarize(results, criteria, design_values=False)
    summary = summary.set_index("design")
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
