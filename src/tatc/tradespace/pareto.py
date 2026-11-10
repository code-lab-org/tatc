"""
Pareto efficiency of designs with several objectives: non-dominated
designs, Pareto ranks, and the hypervolume of a set of designs.

Objectives are given as a mapping of the columns of a data frame (such as a
summary of a study, see `summarize`) to their directions, `"min"` or
`"max"`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

Objectives = Mapping[str, str]
"""Objectives: columns and their directions (`"min"` or `"max"`)."""


def _get_costs(frame: pd.DataFrame, objectives: Objectives) -> np.ndarray:
    """
    Gets the objectives of a data frame as costs (to minimize): maximized
    objectives are negated.
    """
    if len(objectives) == 0:
        raise ValueError("at least one objective is required")
    columns = []
    for column, direction in objectives.items():
        if direction not in ("min", "max"):
            raise ValueError(
                f"direction of {column!r} must be 'min' or 'max', not {direction!r}"
            )
        values = frame[column].to_numpy(dtype=float)
        columns.append(values if direction == "min" else -values)
    costs = np.column_stack(columns)
    if np.isnan(costs).any():
        raise ValueError("objectives must not be NaN")
    return costs


def _get_tolerances(
    objectives: Objectives, tolerance: Mapping[str, float] | None
) -> np.ndarray:
    """
    Gets the tolerances of the objectives (zero if not given).
    """
    tolerance = tolerance or {}
    unknown = set(tolerance) - set(objectives)
    if unknown:
        raise ValueError(f"tolerances of columns that are not objectives: {unknown}")
    values = np.array([float(tolerance.get(column, 0.0)) for column in objectives])
    if np.any(values < 0):
        raise ValueError("tolerances must not be negative")
    return values


def _dominated(
    costs: np.ndarray, tolerances: np.ndarray | None = None, block: int = 512
) -> np.ndarray:
    """
    Finds the points (rows of costs) dominated by another, even if the other
    is penalized by the tolerances: no worse in every cost, and better in at
    least one, by at least its tolerance.
    """
    n = len(costs)
    tolerances = np.zeros(costs.shape[1]) if tolerances is None else tolerances
    dominated = np.zeros(n, dtype=bool)
    for start in range(0, n, block):
        part = costs[start : start + block]
        # [i, j]: point j dominates point i of the block
        penalized = costs[None, :, :] + tolerances
        no_worse = np.all(penalized <= part[:, None, :], axis=2)
        better = np.any(penalized < part[:, None, :], axis=2)
        dominated[start : start + block] = np.any(no_worse & better, axis=1)
    return dominated


def pareto_front(
    frame: pd.DataFrame,
    objectives: Objectives,
    tolerance: Mapping[str, float] | None = None,
) -> pd.Series:
    """
    Finds the Pareto-efficient (non-dominated) rows of a data frame: those
    that no other row dominates, being no worse in every objective and
    better in at least one. Rows with the same objectives are all efficient
    if any is.

    With a `tolerance` of an objective, a row dominates another only if it
    would still dominate it if it were worse by the tolerance in that
    objective. Rows within the tolerances of the front (that no other row is
    clearly better than) are then efficient as well, which suits objectives
    estimated with an uncertainty (such as a few standard errors of their
    means, see `summarize`), or by a cheaper analysis (see `screen`).

    Args:
        frame (pandas.DataFrame): The rows (such as designs).
        objectives (Mapping[str, str]): The objectives: columns and their
            directions (`"min"` or `"max"`).
        tolerance (Mapping[str, float] | None): The tolerances of objectives
            (by column, zero if not given).

    Returns:
        pandas.Series: whether each row is efficient (boolean, with the
            frame's index)
    """
    if len(frame) == 0:
        return pd.Series([], index=frame.index, dtype=bool)
    costs = _get_costs(frame, objectives)
    tolerances = _get_tolerances(objectives, tolerance)
    return pd.Series(~_dominated(costs, tolerances), index=frame.index, name="pareto")


def pareto_ranks(
    frame: pd.DataFrame,
    objectives: Objectives,
    tolerance: Mapping[str, float] | None = None,
) -> pd.Series:
    """
    Ranks the rows of a data frame by non-dominated sorting: rank 0 is the
    Pareto front, rank 1 the front of the others, and so on, optionally with
    tolerances of the objectives (see `pareto_front`).

    Args:
        frame (pandas.DataFrame): The rows (such as designs).
        objectives (Mapping[str, str]): The objectives: columns and their
            directions (`"min"` or `"max"`).
        tolerance (Mapping[str, float] | None): The tolerances of objectives
            (by column, zero if not given).

    Returns:
        pandas.Series: the rank of each row (with the frame's index)
    """
    ranks = np.full(len(frame), -1, dtype=int)
    if len(frame) > 0:
        costs = _get_costs(frame, objectives)
        tolerances = _get_tolerances(objectives, tolerance)
        remaining = np.arange(len(frame))
        rank = 0
        while len(remaining) > 0:
            front = ~_dominated(costs[remaining], tolerances)
            ranks[remaining[front]] = rank
            remaining = remaining[~front]
            rank += 1
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
    objectives: Objectives,
    reference: Mapping[str, float] | Sequence[float],
) -> float:
    """
    Computes the hypervolume of the rows of a data frame: the volume of the
    objective space that they dominate, up to a reference point that is
    worse than every row of interest in every objective. Rows that are not
    better than the reference in every objective add nothing. A larger
    hypervolume is a better set of designs, which makes it a measure of the
    progress of a search, or of how close a search comes to the true Pareto
    front.

    The hypervolume is exact. Its computation slices the objectives in turn,
    which takes long with many objectives and many efficient rows.

    Args:
        frame (pandas.DataFrame): The rows (such as designs).
        objectives (Mapping[str, str]): The objectives: columns and their
            directions (`"min"` or `"max"`).
        reference (Mapping[str, float] | Sequence[float]): The reference
            point: the value of each objective (by column, or in order).

    Returns:
        float: the hypervolume, in the product of the units of the objectives
    """
    if isinstance(reference, Mapping):
        reference = [reference[column] for column in objectives]
    if len(reference) != len(objectives):
        raise ValueError("the reference point needs a value for each objective")
    signs = np.array([1.0 if d == "min" else -1.0 for d in objectives.values()])
    bound = signs * np.asarray(reference, dtype=float)
    if len(frame) == 0:
        return 0.0
    costs = _get_costs(frame, objectives)
    costs = costs[np.all(costs < bound, axis=1)]
    return _hypervolume(costs, bound)


def pareto_probability(
    results: pd.DataFrame,
    objectives: Objectives,
    by: Sequence[str] = (),
    n_resamples: int = 500,
    seed: int = 0,
) -> pd.DataFrame:
    """
    Estimates the probability that each design is Pareto-efficient, given the
    uncertainty of the means of its objectives over replicates, by a
    bootstrap: each resample draws replicates with replacement, the same
    replicates for every design (which keeps the correlation between designs
    of common random numbers, see `Design.realize`), and finds the front of
    the means of the objectives over the drawn replicates. The probability
    is the share of resamples in which a design is efficient.

    Designs that are efficient with a probability between 0 and 1 are those
    whose efficiency the replicates do not settle: more replicates (or a
    tolerance, see `pareto_front`) would settle it. With few replicates,
    the bootstrap has few distinct resamples and understates the
    uncertainty.

    Args:
        results (pandas.DataFrame): The results of a study (see
            `Study.run`), with a record per replicate of each design (and
            group of `by`). Infeasible results are ignored.
        objectives (Mapping[str, str]): The objectives: columns and their
            directions (`"min"` or `"max"`).
        by (Sequence[str]): Columns that group records of an evaluation (such
            as a sensor), whose fronts are found separately.
        n_resamples (int): The number of resamples.
        seed (int): The seed of the resamples.

    Returns:
        pandas.DataFrame: a row per design (and group): the means of the
            objectives over every replicate, whether the design is efficient
            by these means (`pareto`), and the probability that it is
            (`pareto_probability`)
    """
    feasible = results[results["feasible"]] if "feasible" in results else results
    columns = list(objectives)
    signs = _get_costs(pd.DataFrame({column: [1.0] for column in columns}), objectives)[
        0
    ]
    rng = np.random.default_rng(seed)
    frames = []
    groups = feasible.groupby(list(by), sort=False) if by else [((), feasible)]
    for key, group in groups:
        # the costs of each objective, design (rows), and replicate (columns)
        tables = [
            group.pivot_table(index="design", columns="replicate", values=column)
            for column in columns
        ]
        designs, replicates = tables[0].index, tables[0].columns
        costs = np.stack(
            [
                sign * table.loc[designs, replicates].to_numpy(dtype=float)
                for sign, table in zip(signs, tables)
            ]
        )
        means = np.nanmean(costs, axis=2).T
        counts = np.zeros(len(designs))
        for _ in range(n_resamples):
            drawn = rng.integers(0, len(replicates), len(replicates))
            sample = np.nanmean(costs[:, :, drawn], axis=2).T
            # a design without a value in any drawn replicate is not ranked
            valid = ~np.isnan(sample).any(axis=1)
            front = np.zeros(len(designs), dtype=bool)
            front[valid] = ~_dominated(sample[valid])
            counts += front
        frame = pd.DataFrame({"design": designs})
        for name, value in zip(by, key if isinstance(key, tuple) else (key,)):
            frame[name] = value
        for i, column in enumerate(columns):
            frame[column] = signs[i] * means[:, i]
        frame["pareto"] = ~_dominated(means)
        frame["pareto_probability"] = counts / n_resamples
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)
