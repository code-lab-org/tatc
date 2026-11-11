"""
Criteria of designs over their replicates: statistics of the metrics of
their evaluations (such as a mean, a quantile, or the probability of meeting
a threshold), and requirements on these statistics.

The outcomes of a design are uncertain if its replicates differ (see
`Uncertain`), and a statistic over the replicates states how a design is
judged: by its mean, by a pessimistic quantile, or by how reliably it meets
a threshold. A requirement on a statistic is an output constraint, which a
design meets or violates (by an amount) once it is evaluated, as opposed to
the constraints of a design space, which are known before.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import operator as op
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

_OPERATORS = {">=": op.ge, ">": op.gt, "<=": op.le, "<": op.lt}


class Statistic:
    """
    Base class of statistics of a metric (a column of the results of a study)
    over the replicates of each design.

    Args:
        column (str): The metric.
    """

    def __init__(self, column: str):
        self.column = column

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

    def reduce(self, values: np.ndarray) -> np.ndarray:
        counts = np.sum(~np.isnan(values), axis=1)
        result = np.full(values.shape[0], np.nan)
        valid = counts > 1
        if np.any(valid):
            result[valid] = np.nanstd(values[valid], axis=1, ddof=1)
        return result


class Minimum(Statistic):
    """
    The least value of a metric over the replicates (the worst case of a
    metric to maximize).
    """

    def reduce(self, values: np.ndarray) -> np.ndarray:
        return _nan_reduce(np.nanmin, values)


class Maximum(Statistic):
    """
    The greatest value of a metric over the replicates (the worst case of a
    metric to minimize).
    """

    def reduce(self, values: np.ndarray) -> np.ndarray:
        return _nan_reduce(np.nanmax, values)


class Quantile(Statistic):
    """
    A quantile of a metric over the replicates: for example, the 0.1
    quantile of a metric to maximize is a value that the design exceeds in
    90% of the replicates (by linear interpolation between them).

    Args:
        column (str): The metric.
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

    def __repr__(self) -> str:
        return f"Quantile({self.column!r}, {self.q!r})"


class Probability(Statistic):
    """
    The probability that a metric meets a threshold: the share (0 to 1) of
    the replicates in which it does, such as the reliability of a design.

    Args:
        column (str): The metric.
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

    def __repr__(self) -> str:
        return f"Probability({self.column!r}, {self.operator!r}, {self.threshold!r})"


@dataclass(frozen=True)
class Requirement:
    """
    A requirement on a statistic of a design (by its name in the statistics,
    see `aggregate`, or a metric, whose mean it then is): for example,
    `Requirement("sufficient p10", ">=", 80)`. A design that does not meet
    it violates it by the difference between the statistic and the value.

    Args:
        name (str): The statistic (or metric).
        operator (str): The comparison: `">="`, `">"`, `"<="`, or `"<"`.
        value (float): The required value.
    """

    name: str
    operator: str
    value: float

    def __post_init__(self):
        if self.operator not in _OPERATORS:
            raise ValueError(f"operator must be one of {list(_OPERATORS)}")

    def violation(self, statistic: np.ndarray) -> np.ndarray:
        """
        Computes the violations of the requirement: 0 where the statistic
        meets it, and otherwise how far it misses the value (`NaN` where the
        statistic is). A strict comparison that the statistic misses by
        equality is violated by a small amount.

        Args:
            statistic (numpy.ndarray): The statistic of each design.

        Returns:
            numpy.ndarray: the violations
        """
        statistic = np.asarray(statistic, dtype=float)
        if self.operator in (">=", ">"):
            gap = self.value - statistic
        else:
            gap = statistic - self.value
        with np.errstate(invalid="ignore"):
            met = _OPERATORS[self.operator](statistic, self.value)
        violation = np.where(met, 0.0, np.maximum(gap, 1e-9))
        return np.where(np.isnan(statistic), np.nan, violation)


def _get_statistics(
    names: Sequence[str],
    statistics: Mapping[str, Statistic | str] | None,
) -> dict[str, Statistic]:
    """
    Gets the statistics of names: those given (a metric given as a string is
    its mean), and the means of the metrics of names not given.
    """
    statistics = dict(statistics or {})
    resolved = {
        name: Mean(value) if isinstance(value, str) else value
        for name, value in statistics.items()
    }
    for name in names:
        if name not in resolved:
            resolved[name] = Mean(name)
    return resolved


def _get_matrices(
    group: pd.DataFrame, columns: Sequence[str]
) -> tuple[pd.Index, pd.Index, dict[str, np.ndarray]]:
    """
    Gets the values of metrics of each design (rows) in each replicate
    (columns) of the results of a group.
    """
    missing = [column for column in columns if column not in group.columns]
    if missing:
        raise ValueError(f"metrics that are not in the results: {missing}")
    designs = pd.Index(group["design"].drop_duplicates())
    replicates = pd.Index(sorted(group["replicate"].unique()))
    matrices = {}
    for column in columns:
        table = group.pivot_table(
            index="design", columns="replicate", values=column, dropna=False
        )
        matrices[column] = table.reindex(index=designs, columns=replicates).to_numpy(
            dtype=float
        )
    return designs, replicates, matrices


def _evaluate(
    statistics: Mapping[str, Statistic],
    requirements: Sequence[Requirement],
    matrices: Mapping[str, np.ndarray],
    columns: np.ndarray | None = None,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """
    Computes statistics (of the replicates `columns` of the matrices, or all)
    and the total violation of the requirements of each design.
    """
    values = {}
    for name, statistic in statistics.items():
        matrix = matrices[statistic.column]
        values[name] = statistic.reduce(
            matrix if columns is None else matrix[:, columns]
        )
    n = len(next(iter(matrices.values()))) if matrices else 0
    violation = np.zeros(n)
    for requirement in requirements:
        violation = violation + requirement.violation(values[requirement.name])
    return values, violation


def _check_names(
    statistics: Mapping[str, Statistic], requirements: Sequence[Requirement]
) -> None:
    for requirement in requirements:
        if requirement.name not in statistics:
            raise ValueError(f"no statistic named {requirement.name!r}")


def aggregate(
    results: pd.DataFrame,
    statistics: Mapping[str, Statistic | str] | Sequence[str],
    requirements: Sequence[Requirement] = (),
    by: Sequence[str] = (),
    design_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Aggregates the feasible results of a study (see `Study.run`) over the
    replicates of each design (and each group of the columns `by`, such as a
    sensor): the statistics, by name (a metric given as a string, or in a
    list of names, is its mean), the number of `replicates`, the `violation`
    of the requirements (their sum, 0 if the design meets them all), and
    whether the design meets them (`meets_requirements`), with the values of the
    design.

    Args:
        results (pandas.DataFrame): The results.
        statistics (Mapping[str, Statistic | str] | Sequence[str]): The
            statistics, by name, or the metrics whose means to compute.
        requirements (Sequence[Requirement]): Requirements on the statistics
            (or on the means of metrics, which are added to the statistics).
        by (Sequence[str]): Other columns that group records of an
            evaluation.
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every
            replicate of a design (except the metrics).

    Returns:
        pandas.DataFrame: a row per design (and group)
    """
    if not isinstance(statistics, Mapping):
        statistics = {name: Mean(name) for name in statistics}
    resolved = _get_statistics([r.name for r in requirements], statistics)
    _check_names(resolved, requirements)
    feasible = results[results["feasible"]] if "feasible" in results else results
    metrics = sorted({statistic.column for statistic in resolved.values()})
    if design_columns is None:
        excluded = {
            "design",
            "replicate",
            "feasible",
            "error",
            "seconds",
            *metrics,
            *by,
        }
        candidates = [c for c in feasible.columns if c not in excluded]
        constant = feasible.groupby("design")[candidates].nunique(dropna=False)
        design_columns = [c for c in candidates if (constant[c] <= 1).all()]
    design_columns = [c for c in design_columns if c not in resolved and c not in by]
    frames = []
    groups = feasible.groupby(list(by), sort=False) if by else [((), feasible)]
    for key, group in groups:
        designs, _, matrices = _get_matrices(group, metrics)
        values, violation = _evaluate(resolved, requirements, matrices)
        frame = pd.DataFrame({"design": designs})
        for name, value in zip(by, key if isinstance(key, tuple) else (key,)):
            frame[name] = value
        first = group.drop_duplicates("design").set_index("design")
        for column in design_columns:
            frame[column] = first.loc[designs, column].to_numpy()
        for name, value in values.items():
            frame[name] = value
        frame["replicates"] = (
            group.groupby("design")["replicate"].nunique().reindex(designs).to_numpy()
        )
        frame["violation"] = violation
        frame["meets_requirements"] = violation == 0
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=["design", *by, *resolved, "replicates"])
    return pd.concat(frames, ignore_index=True)
