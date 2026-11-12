"""
Criteria of designs over their replicates: statistics of the metrics of
their evaluations (such as a mean, a quantile, or the probability of meeting
a threshold), requirements on these statistics, and the criteria by which
designs are judged (objectives, requirements, and tolerances).

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
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from .study import _STATUS_COLUMNS

_OPERATORS = {">=": op.ge, ">": op.gt, "<=": op.le, "<": op.lt}

_DIRECTIONS = ("min", "max")


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


class Statistic:
    """
    Base class of statistics of a metric over the replicates of each design:
    a column of the results of a study, or a metric computed from others in
    each replicate (such as a `Utility`, whose mean is the expected utility).

    Args:
        column (str | Utility): The metric.
    """

    def __init__(self, column: Any):
        self.column = column

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

    def reduce(self, values: np.ndarray) -> np.ndarray:
        counts = np.sum(~np.isnan(values), axis=1)
        deviation = StandardDeviation(self.column).reduce(values)
        return deviation / np.sqrt(np.maximum(counts, 1))


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
    A requirement on a statistic of a design (by its name in the statistics
    of criteria, see `Criteria`, or a metric, whose mean it then is): for
    example, `Requirement("sufficient p10", ">=", 80)`.

    A design that does not meet it violates it by the difference between the
    statistic and the value, relative to a `scale`: by default the magnitude
    of the value (or 1 if the value is 0), so that the violations of
    requirements in different units are comparable, as they are added (see
    `aggregate`, `screen`, and the constraint-aware optimizers of
    `optimize_optuna` and `optimize_pymoo`).

    Args:
        name (str): The statistic (or metric).
        operator (str): The comparison: `">="`, `">"`, `"<="`, or `"<"`.
        value (float): The required value.
        scale (float | None): The scale of violations, in the units of the
            statistic, or `None` for the magnitude of the value (or 1 if it
            is 0).
    """

    name: str
    operator: str
    value: float
    scale: float | None = None

    def __post_init__(self):
        if self.operator not in _OPERATORS:
            raise ValueError(f"operator must be one of {list(_OPERATORS)}")
        if self.scale is not None and not self.scale > 0:
            raise ValueError("scale must be positive")

    def violation(self, statistic: np.ndarray) -> np.ndarray:
        """
        Computes the violations of the requirement: 0 where the statistic
        meets it, and otherwise how far it misses the value, relative to the
        scale (`NaN` where the statistic is). A strict comparison that the
        statistic misses by equality is violated by a small amount.

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
        scale = self.scale if self.scale is not None else (abs(self.value) or 1.0)
        with np.errstate(invalid="ignore"):
            met = _OPERATORS[self.operator](statistic, self.value)
        violation = np.where(met, 0.0, np.maximum(gap / scale, 1e-9))
        return np.where(np.isnan(statistic), np.nan, violation)


Tolerance = float | str
"""A tolerance of an objective: a value in its units, or a share of its value
(a string of a percentage, such as `"5%"`)."""


def _get_tolerances(
    objectives: Mapping[str, str], tolerance: Mapping[str, Tolerance]
) -> np.ndarray:
    """
    Gets the tolerances of the objectives (zero if not given): their absolute
    values (first row) and their shares of the values of the objectives
    (second row).
    """
    unknown = set(tolerance) - set(objectives)
    if unknown:
        raise ValueError(f"tolerances of names that are not objectives: {unknown}")
    values = np.zeros((2, len(objectives)))
    for i, name in enumerate(objectives):
        value = tolerance.get(name, 0.0)
        if isinstance(value, str):
            if not value.strip().endswith("%"):
                raise ValueError(
                    f"a tolerance string must be a percentage, not {value!r}"
                )
            values[1, i] = float(value.strip()[:-1]) / 100
        else:
            values[0, i] = float(value)
    if np.any(values < 0):
        raise ValueError("tolerances must not be negative")
    return values


@dataclass(frozen=True)
class Criteria:
    """
    The criteria by which designs are judged from the results of their
    replicates: objectives to minimize or maximize and requirements to
    meet, the statistics over the replicates that they name (see
    `Statistic`), the tolerances of the objectives within which differences
    do not matter, and the columns that group the records of an evaluation
    (such as a sensor), in each group of which designs are judged
    separately.

    Objectives and requirements name statistics: those of `statistics`, or
    otherwise the means of the metrics of their names. A design is
    efficient if it meets the requirements and no other design that meets
    them is as good in every objective and better in at least one, by more
    than the tolerances (see `pareto_front`).

    The functions that judge designs (`aggregate`, `pareto_front`,
    `pareto_ranks`, `hypervolume`, `pareto_probability`, `refine`,
    `confirm`, `screen`, and the searches of `optimize_optuna` and
    `optimize_pymoo`) take criteria, or a mapping of objectives to their
    directions, which are criteria of the means of metrics alone.

    Args:
        objectives (Mapping[str, str]): The objectives: statistics (or
            metrics) and their directions (`"min"` or `"max"`).
        statistics (Mapping[str, Statistic | str]): Statistics, by name, of
            objectives and requirements that are not means of the metrics of
            their names (a metric given as a string, or a `Utility`, is its
            mean), and others to compute with them (see `aggregate`).
        requirements (Sequence[Requirement]): Requirements on statistics.
        tolerance (Mapping[str, float | str]): The tolerances of objectives
            (by name, zero if not given): values in their units, or
            percentages of their values (such as `"5%"`).
        by (Sequence[str] | str): The columns that group the records of an
            evaluation (such as a sensor), in each group of which designs are
            judged separately.
    """

    objectives: Mapping[str, str]
    statistics: Mapping[str, Statistic | str] = field(default_factory=dict)
    requirements: Sequence[Requirement] = ()
    tolerance: Mapping[str, Tolerance] = field(default_factory=dict)
    by: Sequence[str] | str = ()

    def __post_init__(self):
        if not isinstance(self.objectives, Mapping):
            raise TypeError("objectives must be a mapping of names to directions")
        objectives = dict(self.objectives)
        if len(objectives) == 0:
            raise ValueError("at least one objective is required")
        for name, direction in objectives.items():
            if direction not in _DIRECTIONS:
                raise ValueError(
                    f"direction of {name!r} must be 'min' or 'max', not {direction!r}"
                )
        statistics = {
            name: Mean(value) if isinstance(value, str) or _is_derived(value) else value
            for name, value in dict(self.statistics).items()
        }
        for name, statistic in statistics.items():
            if not isinstance(statistic, Statistic):
                raise TypeError(f"statistic {name!r} is not a Statistic")
        requirements = tuple(self.requirements)
        for requirement in requirements:
            if not isinstance(requirement, Requirement):
                raise TypeError(f"{requirement!r} is not a Requirement")
        tolerance = dict(self.tolerance)
        _get_tolerances(objectives, tolerance)  # checks the tolerances
        by = (self.by,) if isinstance(self.by, str) else tuple(self.by)
        object.__setattr__(self, "objectives", objectives)
        object.__setattr__(self, "statistics", statistics)
        object.__setattr__(self, "requirements", requirements)
        object.__setattr__(self, "tolerance", tolerance)
        object.__setattr__(self, "by", by)

    @property
    def names(self) -> list[str]:
        """list[str]: The names of the statistics: of the objectives, of the
        requirements, and the others of `statistics`, in that order."""
        return list(
            dict.fromkeys(
                [
                    *self.objectives,
                    *(requirement.name for requirement in self.requirements),
                    *self.statistics,
                ]
            )
        )

    def statistic(self, name: str) -> Statistic:
        """
        Gets the statistic of a name: that of `statistics`, or the mean of
        the metric of the name.

        Args:
            name (str): The name.

        Returns:
            Statistic: the statistic
        """
        return self.statistics[name] if name in self.statistics else Mean(name)

    @property
    def metrics(self) -> list[str]:
        """list[str]: The metrics (columns of the results of a study) of
        the statistics."""
        return sorted(
            {column for name in self.names for column in self.statistic(name).columns}
        )


def _as_criteria(criteria: Criteria | Mapping[str, str]) -> Criteria:
    """
    Gets criteria, or the criteria of the means of objectives.
    """
    if isinstance(criteria, Criteria):
        return criteria
    if isinstance(criteria, Mapping):
        return Criteria(criteria)
    raise TypeError(
        "criteria must be Criteria or a mapping of objectives to their directions"
    )


def _resolve(criteria: Criteria) -> dict[str, Statistic]:
    """
    Gets the statistic of each name of criteria.
    """
    return {name: criteria.statistic(name) for name in criteria.names}


def _signs(criteria: Criteria) -> np.ndarray:
    """
    Gets the signs that make the objectives of criteria costs (to minimize).
    """
    return np.array(
        [
            1.0 if direction == "min" else -1.0
            for direction in criteria.objectives.values()
        ]
    )


def _groups(frame: pd.DataFrame, by: Sequence[str]) -> Iterator[tuple[tuple, Any]]:
    """
    Iterates over the groups of a data frame by columns (or the frame, if
    none): the values of the columns, and the positions of the group's rows.
    """
    if not by:
        yield (), np.arange(len(frame))
        return
    missing = [column for column in by if column not in frame.columns]
    if missing:
        raise ValueError(f"no columns {missing} to group by")
    for key, positions in frame.groupby(list(by), sort=False).indices.items():
        yield (key if isinstance(key, tuple) else (key,)), positions


def _get_matrices(
    group: pd.DataFrame, columns: Sequence[str]
) -> tuple[pd.Index, pd.Index, dict[str, np.ndarray]]:
    """
    Gets the values of metrics of each design (rows) in each replicate
    (columns) of the results of a group, which must have one record per
    evaluation (a design in a replicate).
    """
    missing = [column for column in columns if column not in group.columns]
    if missing:
        raise ValueError(f"metrics that are not in the results: {missing}")
    keys = ["design", "replicate"]
    if group.duplicated(keys).any():
        candidates = [
            column
            for column in group.columns
            if column not in (*_STATUS_COLUMNS, *columns)
        ]
        counts = group.groupby(keys)[candidates].nunique(dropna=False).max()
        varying = list(counts.index[counts > 1])
        raise ValueError(
            "the results have several records of an evaluation (a design in a "
            + "replicate): give the columns that separate them as `by` (of the "
            + "criteria)"
            + (f", such as {varying}" if varying else "")
        )
    designs = pd.Index(group["design"].drop_duplicates())
    replicates = pd.Index(sorted(group["replicate"].unique()))
    matrices = {}
    for column in columns:
        table = group.pivot(index="design", columns="replicate", values=column)
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
        matrix = statistic.values_of(matrices)
        values[name] = statistic.reduce(
            matrix if columns is None else matrix[:, columns]
        )
    n = len(next(iter(matrices.values()))) if matrices else 0
    violation = np.zeros(n)
    for requirement in requirements:
        violation = violation + requirement.violation(values[requirement.name])
    return values, violation


def _design_columns(results: pd.DataFrame, excluded: set[str]) -> list[str]:
    """
    Gets the columns of the values of designs: those (other than `excluded`)
    that are the same in every record of each design.
    """
    candidates = [column for column in results.columns if column not in excluded]
    if not candidates or len(results) == 0:
        return []
    constant = results.groupby("design")[candidates].nunique(dropna=False)
    return [column for column in candidates if (constant[column] <= 1).all()]


def aggregate(
    results: pd.DataFrame,
    statistics: Criteria | Mapping[str, Statistic | str] | Sequence[str],
    by: Sequence[str] | str = (),
    design_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Aggregates the feasible results of a study (see `Study.run`) over the
    replicates of each design (and each group of the columns `by`, such as a
    sensor): the statistics, by name (a metric given as a string, or in a
    list of names, is its mean), the number of `replicates`, and the values
    of the design.

    With criteria (see `Criteria`), the statistics that they name, in each
    group of their `by`, and the `violation` of their requirements (the sum
    of the violations, see `Requirement`, 0 if the design meets them all)
    and whether the design meets them (`meets_requirements`).

    The results must have one record per evaluation (a design in a
    replicate) in each group.

    Args:
        results (pandas.DataFrame): The results.
        statistics (Criteria | Mapping[str, Statistic | str] | Sequence[str]):
            Criteria, the statistics by name, or the metrics whose means to
            compute.
        by (Sequence[str] | str): Other columns that group records of an
            evaluation (those of the criteria, with criteria).
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every
            replicate of a design (except the metrics).

    Returns:
        pandas.DataFrame: a row per design (and group)
    """
    criteria = statistics if isinstance(statistics, Criteria) else None
    if criteria is not None:
        if by:
            raise ValueError("criteria group records by their own `by`")
        resolved, requirements, by = (
            _resolve(criteria),
            criteria.requirements,
            criteria.by,
        )
    else:
        if not isinstance(statistics, Mapping):
            statistics = {name: name for name in statistics}
        resolved = {
            name: Mean(value) if isinstance(value, str) or _is_derived(value) else value
            for name, value in statistics.items()
        }
        requirements = ()
        by = (by,) if isinstance(by, str) else tuple(by)
    feasible = results[results["feasible"]] if "feasible" in results else results
    metrics = sorted(
        {column for statistic in resolved.values() for column in statistic.columns}
    )
    if design_columns is None:
        design_columns = _design_columns(feasible, {*_STATUS_COLUMNS, *metrics, *by})
    design_columns = [c for c in design_columns if c not in resolved and c not in by]
    frames = []
    for key, positions in _groups(feasible, by):
        group = feasible.iloc[positions]
        designs, _, matrices = _get_matrices(group, metrics)
        values, violation = _evaluate(resolved, requirements, matrices)
        frame = pd.DataFrame({"design": designs})
        for name, value in zip(by, key):
            frame[name] = value
        first = group.drop_duplicates("design").set_index("design")
        for column in design_columns:
            frame[column] = first.loc[designs, column].to_numpy()
        for name, value in values.items():
            frame[name] = value
        frame["replicates"] = (
            group.groupby("design")["replicate"].nunique().reindex(designs).to_numpy()
        )
        if criteria is not None:
            frame["violation"] = violation
            frame["meets_requirements"] = violation == 0
        frames.append(frame)
    if not frames:
        judged = ["violation", "meets_requirements"] if criteria is not None else []
        return pd.DataFrame(columns=["design", *by, *resolved, "replicates", *judged])
    return pd.concat(frames, ignore_index=True)


def summarize(
    results: pd.DataFrame,
    metrics: Sequence[str],
    by: Sequence[str] | str = (),
    design_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Summarizes the metrics of the feasible results of a study over their
    replicates (see `aggregate`): for each design (and each group of the
    columns `by`, such as a sensor), the mean of each metric, its standard
    deviation (`<metric> std`) and standard error (`<metric> sem`) between
    replicates, and the number of replicates (`replicates`), with the values
    of the design.

    Args:
        results (pandas.DataFrame): The results (see `Study.run`).
        metrics (Sequence[str]): The columns of the metrics.
        by (Sequence[str] | str): Other columns that group records of an
            evaluation.
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every
            replicate of a design (except the metrics).

    Returns:
        pandas.DataFrame: the summary, a row per design (and group)
    """
    statistics: dict[str, Statistic] = {}
    for metric in metrics:
        statistics[metric] = Mean(metric)
        statistics[f"{metric} std"] = StandardDeviation(metric)
        statistics[f"{metric} sem"] = StandardError(metric)
    return aggregate(results, statistics, by=by, design_columns=design_columns)
