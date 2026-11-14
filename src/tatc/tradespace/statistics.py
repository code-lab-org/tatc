"""
Statistics of the metrics of designs over their replicates: such as a mean,
a standard deviation, a quantile, or the probability of meeting a threshold.

The outcomes of a design are uncertain if its replicates differ (see
`Uncertain`), and a statistic over the replicates states how a design is
judged (see `Criteria`): by its mean, by a pessimistic quantile, or by how
reliably it meets a threshold.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import operator as op
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

_OPERATORS = {">=": op.ge, ">": op.gt, "<=": op.le, "<": op.lt}


def _names(value: str | Sequence[str]) -> list[str]:
    """Gets names given as a name or a sequence of names."""
    return [value] if isinstance(value, str) else list(value)


def _is_derived(metric: Any) -> bool:
    """
    Whether a metric is computed from others in each replicate (such as a
    `Utility`): an object with the `columns` it needs and an `evaluate`
    method of their values.
    """
    return (
        not isinstance(metric, (str, Statistic))
        and hasattr(metric, "columns")
        and callable(getattr(metric, "evaluate", None))
    )


def _as_statistic(name: str, value: Any) -> Statistic:
    """
    Gets a statistic given by name: a statistic, or the mean of a metric
    given as a string or a metric computed from others (such as a
    `Utility`).
    """
    if isinstance(value, str) or _is_derived(value):
        return Mean(value)
    if not isinstance(value, Statistic):
        raise TypeError(f"statistic {name!r} is not a Statistic")
    return value


def _named_statistics(statistics: Any) -> dict[str, Statistic]:
    """
    Gets statistics by name: of a mapping of names to statistics (or metrics,
    whose means they are), or of a sequence of statistics (or metrics), named
    by their conventional names (see `Statistic.name`).
    """
    if statistics is None:
        return {}
    if isinstance(statistics, Mapping):
        return {name: _as_statistic(name, value) for name, value in statistics.items()}
    if isinstance(statistics, (str, Statistic)) or _is_derived(statistics):
        statistics = [statistics]
    named: dict[str, Statistic] = {}
    for item in statistics:
        statistic = _as_statistic(str(item), item)
        if statistic.name in named and named[statistic.name] != statistic:
            raise ValueError(f"statistics of the same name: {statistic.name!r}")
        named[statistic.name] = statistic
    return named


class Statistic:
    """
    Base class of statistics of a metric over the replicates of each design:
    a column of the results of a study, or a metric computed from others in
    each replicate (such as a `Utility`, whose mean is the expected utility).

    The column `feasible` of the results is the feasibility of each replicate
    (1 if it is feasible, 0 if not), whose statistics are over every
    replicate of a design, feasible or not (its mean is the share of its
    replicates that are feasible, see `aggregate`); the statistics of other
    metrics are over the feasible replicates.

    Args:
        column (str | Utility): The metric.
    """

    def __init__(self, column: Any):
        self.column = column

    @property
    def metric(self) -> str:
        """str: The name of the statistic's metric: its column, or `"utility"`
        of a metric computed from others."""
        return "utility" if _is_derived(self.column) else str(self.column)

    @property
    def name(self) -> str:
        """str: The conventional name of the statistic (the suffixes of which
        are those of `summarize` and `summarize_histories`): such as `cost`
        (a mean), `cost std`, `cost se`, `cost min`, `cost max`, `cost p70`
        (a quantile), or `P(cost <= 150)`, by which criteria name the
        statistics given without names."""
        return self.metric

    @property
    def columns(self) -> tuple[str, ...]:
        """tuple[str, ...]: The columns of the results that the statistic
        needs: its metric, or those its metric is computed from."""
        if _is_derived(self.column):
            return tuple(self.column.columns)
        return (self.column,)

    def values_of(self, matrices: Mapping[str, np.ndarray]) -> np.ndarray:
        """
        Gets the values of the statistic's metric in each replicate of each
        design, from those of the columns of the results.

        Args:
            matrices (Mapping[str, numpy.ndarray]): The values of each column
                (shape (designs, replicates)).

        Returns:
            numpy.ndarray: the values of the metric (shape (designs,
                replicates))
        """
        if _is_derived(self.column):
            values = {name: matrices[name] for name in self.columns}
            return np.asarray(self.column.evaluate(values), dtype=float)
        return matrices[self.column]

    def reduce(self, values: np.ndarray) -> np.ndarray:
        """
        Computes the statistic of each design (row) from its values in its
        replicates (columns), ignoring missing values (`NaN`).

        Args:
            values (numpy.ndarray): The values (shape (designs, replicates)).

        Returns:
            numpy.ndarray: the statistic of each design (`NaN` without values)
        """
        raise NotImplementedError

    def __eq__(self, other: object) -> bool:
        return type(self) is type(other) and vars(self) == vars(other)

    def __hash__(self) -> int:
        return hash((type(self), tuple(sorted(vars(self).items()))))

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.column!r})"


def _all_nan(values: np.ndarray) -> np.ndarray:
    return np.all(np.isnan(values), axis=1)


def _nan_reduce(function, values: np.ndarray) -> np.ndarray:
    """
    Applies a nan-aware reduction along replicates, giving `NaN` (without a
    warning) for designs without values.
    """
    result = np.full(values.shape[0], np.nan)
    valid = ~_all_nan(values)
    if np.any(valid):
        result[valid] = function(values[valid], axis=1)
    return result


class Mean(Statistic):
    """
    The mean of a metric over the replicates.
    """

    def reduce(self, values: np.ndarray) -> np.ndarray:
        return _nan_reduce(np.nanmean, values)


class StandardDeviation(Statistic):
    """
    The (sample) standard deviation of a metric over the replicates, a
    measure of how much a design's outcome depends on its uncertain
    parameters (`NaN` with fewer than two replicates).
    """

    @property
    def name(self) -> str:
        return f"{self.metric} std"

    def reduce(self, values: np.ndarray) -> np.ndarray:
        counts = np.sum(~np.isnan(values), axis=1)
        result = np.full(values.shape[0], np.nan)
        valid = counts > 1
        if np.any(valid):
            result[valid] = np.nanstd(values[valid], axis=1, ddof=1)
        return result


class StandardError(Statistic):
    """
    The standard error of the mean of a metric over the replicates: its
    standard deviation over the square root of the number of replicates, a
    measure of the uncertainty of the mean (`NaN` with fewer than two
    replicates).
    """

    @property
    def name(self) -> str:
        return f"{self.metric} se"

    def reduce(self, values: np.ndarray) -> np.ndarray:
        counts = np.sum(~np.isnan(values), axis=1)
        deviation = StandardDeviation(self.column).reduce(values)
        return deviation / np.sqrt(np.maximum(counts, 1))


class Minimum(Statistic):
    """
    The least value of a metric over the replicates (the worst case of a
    metric to maximize).
    """

    @property
    def name(self) -> str:
        return f"{self.metric} min"

    def reduce(self, values: np.ndarray) -> np.ndarray:
        return _nan_reduce(np.nanmin, values)


class Maximum(Statistic):
    """
    The greatest value of a metric over the replicates (the worst case of a
    metric to minimize).
    """

    @property
    def name(self) -> str:
        return f"{self.metric} max"

    def reduce(self, values: np.ndarray) -> np.ndarray:
        return _nan_reduce(np.nanmax, values)


class Quantile(Statistic):
    """
    A quantile of a metric over the replicates: for example, the 0.1
    quantile of a metric to maximize is a value that the design exceeds in
    90% of the replicates (by linear interpolation between them).

    Args:
        column (str | Utility): The metric.
        q (float): The quantile, between 0 and 1.
    """

    def __init__(self, column: str, q: float):
        super().__init__(column)
        if not 0 <= q <= 1:
            raise ValueError("q must be between 0 and 1")
        self.q = q

    def reduce(self, values: np.ndarray) -> np.ndarray:
        # numpy's nanquantile, vectorized: sort each row (missing values
        # last), and interpolate linearly between the values around the
        # quantile's position among each row's values
        ordered = np.sort(values, axis=1)
        counts = np.sum(~np.isnan(values), axis=1)
        position = self.q * np.maximum(counts - 1, 0)
        below = np.floor(position).astype(int)
        above = np.minimum(below + 1, np.maximum(counts - 1, 0))
        rows = np.arange(len(values))
        weight = position - below
        with np.errstate(invalid="ignore"):
            result = (1 - weight) * ordered[rows, below] + weight * ordered[rows, above]
        result[weight == 0] = ordered[rows, below][weight == 0]
        result[counts == 0] = np.nan
        return result

    @property
    def name(self) -> str:
        return f"{self.metric} p{100 * self.q:g}"

    def __repr__(self) -> str:
        return f"Quantile({self.column!r}, {self.q!r})"


class Probability(Statistic):
    """
    The probability that a metric meets a threshold: the share (0 to 1) of
    the replicates in which it does, such as the reliability of a design.

    Args:
        column (str | Utility): The metric.
        operator (str): The comparison: `">="`, `">"`, `"<="`, or `"<"`.
        threshold (float): The threshold.
    """

    def __init__(self, column: str, operator: str, threshold: float):
        super().__init__(column)
        if operator not in _OPERATORS:
            raise ValueError(f"operator must be one of {list(_OPERATORS)}")
        self.operator, self.threshold = operator, threshold

    def reduce(self, values: np.ndarray) -> np.ndarray:
        valid = ~np.isnan(values)
        with np.errstate(invalid="ignore"):
            met = _OPERATORS[self.operator](values, self.threshold) & valid
        counts = valid.sum(axis=1)
        result = np.full(values.shape[0], np.nan)
        result[counts > 0] = met.sum(axis=1)[counts > 0] / counts[counts > 0]
        return result

    @property
    def name(self) -> str:
        return f"P({self.metric} {self.operator} {self.threshold:g})"

    def __repr__(self) -> str:
        return f"Probability({self.column!r}, {self.operator!r}, {self.threshold!r})"
