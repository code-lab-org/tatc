"""
Dominance between points of several costs (objectives to minimize, such as
the objectives of designs, see `pareto_front`), with tolerances within
which differences do not matter.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np

BLOCK_ELEMENTS = 2_000_000
"""The most elements of the temporary arrays of a block of comparisons."""

CHUNK = 64
"""The points of each chunk of a sweep (see `dominated_plain`)."""


def dominated_by(
    points: np.ndarray, others: np.ndarray, tolerances: np.ndarray | None = None
) -> np.ndarray:
    """
    Finds the points dominated by any of others, even if the other is
    penalized by the tolerances (see `dominated`), comparing them in blocks
    whose temporary arrays are limited to a few million elements.

    Args:
        points (numpy.ndarray): The points (a row each).
        others (numpy.ndarray): The others (a row each).
        tolerances (numpy.ndarray | None): The tolerances of the costs (see
            `dominated`), or `None` for none.

    Returns:
        numpy.ndarray: whether each point is dominated
    """
    beaten = np.zeros(len(points), dtype=bool)
    if len(points) == 0 or len(others) == 0:
        return beaten
    if tolerances is None:
        tolerances = np.zeros((2, points.shape[1]))
    absolute, relative = tolerances
    penalties = np.any(tolerances)
    block = max(1, BLOCK_ELEMENTS // max(len(others) * points.shape[1], 1))
    for start in range(0, len(points), block):
        part = points[start : start + block]
        # [i, j]: other j dominates point i of the block
        penalized = others[None, :, :]
        if penalties:
            penalized = penalized + absolute + relative * np.abs(part[:, None, :])
        no_worse = np.all(penalized <= part[:, None, :], axis=2)
        better = np.any(penalized < part[:, None, :], axis=2)
        beaten[start : start + block] = np.any(no_worse & better, axis=1)
    return beaten


def dominated_plain(costs: np.ndarray) -> np.ndarray:
    """
    Finds the points dominated by another (without tolerances) by a sweep in
    lexicographic order, in which a point can only be dominated by points
    before it, and is dominated by one if it is by any of the efficient
    points before it (dominance is transitive). The sweep compares each
    chunk of points with the efficient points before it, then those it does
    not dominate with each other.

    Args:
        costs (numpy.ndarray): The points (a row each).

    Returns:
        numpy.ndarray: whether each point is dominated
    """
    n = len(costs)
    order = np.lexsort(costs.T[::-1])
    ordered = costs[order]
    beaten_points = np.zeros(n, dtype=bool)
    front = np.empty((0, costs.shape[1]))
    for start in range(0, n, CHUNK):
        part = ordered[start : start + CHUNK]
        beaten = dominated_by(part, front)
        # a point that an earlier efficient point dominates dominates nothing
        # that it does not: compare only the others with each other
        survivors = part[~beaten]
        beaten[~beaten] = dominated_by(survivors, survivors)
        beaten_points[order[start : start + CHUNK]] = beaten
        front = np.vstack([front, part[~beaten]])
    return beaten_points


def dominated(costs: np.ndarray, tolerances: np.ndarray | None = None) -> np.ndarray:
    """
    Finds the points (rows of costs) dominated by another, even if the other
    is penalized by the tolerances: no worse in every cost, and better in at
    least one, by at least its tolerance (absolute, plus a share of the
    magnitude of the cost of the point it would dominate: the rows of
    `tolerances`, see `criteria._tolerances`).

    A point that is dominated with tolerances is also dominated by a point
    that is efficient without them: the penalty depends only on the point
    dominated, so any point no worse than its dominator in every cost (as
    some efficient point is) dominates it too. Points are therefore compared
    only with the efficient points (see `dominated_plain`), which are few.

    Args:
        costs (numpy.ndarray): The points (a row each).
        tolerances (numpy.ndarray | None): The absolute tolerances of the
            costs (first row) and their shares of the magnitudes of the costs
            (second row), or `None` for none.

    Returns:
        numpy.ndarray: whether each point is dominated
    """
    plain = dominated_plain(costs)
    if tolerances is None or not np.any(tolerances):
        return plain
    return dominated_by(costs, costs[~plain], tolerances)


def efficient(
    costs: np.ndarray, admissible: np.ndarray, tolerances: np.ndarray | None = None
) -> np.ndarray:
    """
    Finds the efficient points (rows of costs) among the admissible ones
    (such as the designs that meet requirements): those that no other
    admissible point dominates (see `dominated`). Points with a missing cost
    (`NaN`) are neither efficient nor dominate others.

    Args:
        costs (numpy.ndarray): The points (a row each).
        admissible (numpy.ndarray): Whether each point is admissible.
        tolerances (numpy.ndarray | None): The tolerances of the costs (see
            `dominated`), or `None` for none.

    Returns:
        numpy.ndarray: whether each point is efficient
    """
    known = np.asarray(admissible, dtype=bool)
    if len(costs) > 0:
        known = known & ~np.isnan(costs).any(axis=1)
    found = np.zeros(len(costs), dtype=bool)
    if known.any():
        found[known] = ~dominated(costs[known], tolerances)
    return found
