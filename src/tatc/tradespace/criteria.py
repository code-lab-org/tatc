"""
Criteria by which designs are judged from the statistics of the metrics of
their evaluations over their replicates (see `Statistic`): objectives to
minimize or maximize, requirements on statistics, the tolerances of the
objectives, and the rules that combine the statistics of a design across
groups of its results (such as scenarios).

A requirement on a statistic is an output constraint, which a design meets
or violates (by an amount) once it is evaluated, as opposed to the
constraints of a design space, which are known before.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from . import tables
from .statistics import _OPERATORS, Mean, Statistic, _named_statistics

_DIRECTIONS = ("min", "max")


@dataclass(frozen=True)
class Requirement:
    """
    A requirement on a statistic of a design (by its name in the statistics
    of criteria, see `Criteria`, or a metric, whose mean it then is): for
    example, `Requirement("sufficient p10", ">=", 80)`, or
    `Requirement("feasible", ">=", 1)`, which requires every replicate of a
    design to be feasible (`feasible` is the share of its replicates that
    are, see `aggregate`).

    A design that does not meet it violates it by the difference between the
    statistic and the value, relative to a `scale`: by default the magnitude
    of the value (or 1 if the value is 0), so that the violations of
    requirements in different units are comparable. Summaries add them (the
    `violation` of a design, see `aggregate`), by which the designs that miss
    requirements are ranked (see `pareto_ranks` and `screen`), and the
    constraint-aware optimizers of `optimize_optuna` and `optimize_pymoo`
    keep each as a constraint of its own.

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

    def met(self, statistic: np.ndarray) -> np.ndarray:
        """
        Gets whether values of the statistic meet the requirement (`False`
        where they are `NaN`).

        Args:
            statistic (numpy.ndarray): The values.

        Returns:
            numpy.ndarray: whether each meets it
        """
        values = np.asarray(statistic, dtype=float)
        with np.errstate(invalid="ignore"):
            return _OPERATORS[self.operator](values, self.value)

    def __str__(self) -> str:
        return f"{self.name} {self.operator} {self.value:g}"

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
        violation = np.where(self.met(statistic), 0.0, np.maximum(gap / scale, 1e-9))
        return np.where(np.isnan(statistic), np.nan, violation)


_RULES = ("worst", "mean")


@dataclass(frozen=True)
class Across:
    """
    A rule that combines the statistics of each design across the groups of
    a column of the results, such as the scenarios of a study (see `Study`),
    so that criteria (see `Criteria`) judge each design once, rather than in
    each group (as criteria `by` the column do): each statistic is computed
    over the replicates in each group, and then combined across the groups.

    * `"worst"`: a statistic is its worst value in any group: the largest of
      an objective to minimize (the smallest of one to maximize), and the
      worst of a statistic of requirements for them (the smallest for `">="`),
      so that a design meets a requirement only if it does in every group,
      as a design robust to the scenarios must.
    * `"mean"`: a statistic is its mean over the groups, weighted by their
      `weights` (equal if not given), as an expectation over scenarios of
      these probabilities.

    `rules` combine some statistics (by name) otherwise than `rule`: for
    example, objectives of the mean over scenarios, with a requirement in
    every scenario on a statistic of another name::

        Criteria(
            {"accurate": "max", "satellites": "min"},
            statistics={"accurate anywhere": "accurate"},
            requirements=[Requirement("accurate anywhere", ">=", 40)],
            across=Across("scenario", "mean", rules={"accurate anywhere": "worst"}),
        )

    The worst of a statistic needs its direction: that of the objective or
    the requirements that name it, which must agree; other statistics are
    combined by their means, unless a rule asks for their worst. The spreads
    that `summarize` adds are those of the combined statistics (see
    `summarize`). Requirements and tolerances apply to the combined
    statistics. A design is judged only if it has feasible results in every
    group of the results (in each group of the criteria's `by`). Statistics
    of the feasibility of replicates (such as `feasible`, the share of a
    design's replicates that are feasible, see `aggregate`) are not
    combined: they are over the evaluations of every group together.

    Combining the statistics of the groups (such as the worst of the means)
    differs from a statistic of the worst group of each replicate (the mean
    of the worst in each replicate), which, for a mean, is more pessimistic.

    Args:
        column (str): The column of the groups, such as `"scenario"`.
        rule (str): How statistics are combined: `"worst"` or `"mean"`.
        weights (Mapping[Any, float] | None): The weights of the mean of the
            groups (by value of the column), or `None` for equal weights:
            normalized over the groups of the results, each of which must
            have one.
        rules (Mapping[str, str]): The rules of some statistics (by name), in
            place of `rule`.
    """

    column: str = "scenario"
    rule: str = "worst"
    weights: Mapping[Any, float] | None = None
    rules: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self):
        if not isinstance(self.column, str) or not self.column:
            raise TypeError("the column must be a name (a non-empty string)")
        rules = dict(self.rules)
        for rule in [self.rule, *rules.values()]:
            if rule not in _RULES:
                raise ValueError(f"a rule must be one of {list(_RULES)}, not {rule!r}")
        weights = None if self.weights is None else dict(self.weights)
        if weights is not None and not all(
            np.isfinite(float(w)) and float(w) > 0 for w in weights.values()
        ):
            raise ValueError("the weights of groups must be positive")
        object.__setattr__(self, "rules", rules)
        object.__setattr__(self, "weights", weights)

    def weights_of(self, groups: Sequence[Any]) -> np.ndarray:
        """
        Gets the weights of the mean of groups, normalized.

        Args:
            groups (Sequence[Any]): The groups (values of the column).

        Returns:
            numpy.ndarray: the weight of each group, which sum to 1
        """
        if self.weights is None:
            return np.full(len(groups), 1.0 / max(len(groups), 1))
        missing = [group for group in groups if group not in self.weights]
        if missing:
            raise ValueError(f"groups of {self.column!r} without weights: {missing}")
        weights = np.array([float(self.weights[group]) for group in groups])
        return weights / weights.sum()


Tolerance = float | str
"""A tolerance of an objective: a value in its units, or a share of its value
(a string of a percentage, such as `"5%"`)."""


def _tolerance_array(
    objectives: Mapping[str, str], tolerances: Mapping[str, Tolerance]
) -> np.ndarray:
    """
    Gets the tolerances of the objectives (zero if not given): their absolute
    values (first row) and their shares of the values of the objectives
    (second row), checking them.
    """
    unknown = set(tolerances) - set(objectives)
    if unknown:
        raise ValueError(f"tolerances of names that are not objectives: {unknown}")
    values = np.zeros((2, len(objectives)))
    for i, name in enumerate(objectives):
        value = tolerances.get(name, 0.0)
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
    separately::

        criteria = Criteria(
            {"cost p70": "min", "coverage": "max"},
            statistics=[Quantile("cost", 0.7)],  # named cost p70
            requirements=[Requirement("cost p70", "<=", 150)],
            tolerances={"coverage": 1.0},
        )
        summary = aggregate(results, criteria)
        front = summary[summary.pareto]

    Objectives and requirements name statistics: those of `statistics`, or
    otherwise the means of the metrics of their names. A design is
    efficient if it meets the requirements and no other design that meets
    them is as good in every objective and better in at least one, by more
    than the tolerances (see `pareto_front`).

    Designs are judged by their feasible replicates. The share of a design's
    replicates that are feasible is the statistic `feasible` (see
    `aggregate`), which requirements may name, such as
    `Requirement("feasible", ">=", 1)` of designs feasible in every
    replicate. Designs without a feasible replicate are not judged.

    Designs evaluated in several groups of a column (such as the scenarios
    of a study, see `Study`) are judged in each group with criteria `by` the
    column, or once with criteria `across` it (see `Across`), which combine
    their statistics across the groups: their worst, or their mean.

    The functions that judge designs (`aggregate`, `summarize`,
    `pareto_front`, `pareto_ranks`, `hypervolume`, `hypervolume_progress`,
    `pareto_probability`, `robustness`, `refine`, `confirm`, `screen`, and the
    searches of `optimize_optuna` and `optimize_pymoo`) take criteria, or a
    mapping of objectives to their directions, which are criteria of the
    means of metrics alone.

    Args:
        objectives (Mapping[str, str]): The objectives: statistics (or
            metrics) and their directions (`"min"` or `"max"`).
        statistics (Mapping[str, Statistic | str] | Sequence): Statistics,
            by name, of objectives and requirements that are not means of
            the metrics of their names (a metric given as a string, or a
            `Utility`, is its mean), and others to compute with them (see
            `aggregate`), or a sequence of them, named by their conventional
            names (see `Statistic.name`, such as `cost p70` of
            `Quantile("cost", 0.7)`).
        requirements (Sequence[Requirement]): Requirements on statistics.
        tolerances (Mapping[str, float | str]): The tolerances of objectives
            (by name, zero if not given): values in their units, or
            percentages of their values (such as `"5%"`).
        by (Sequence[str] | str): The columns that group the records of an
            evaluation (such as a sensor), in each group of which designs are
            judged separately.
        across (Across | None): The rule that combines the statistics of
            each design across the groups of a column (such as the
            scenarios), if any.
    """

    objectives: Mapping[str, str]
    statistics: Mapping[str, Statistic | str] | Sequence[Statistic | str] = field(
        default_factory=dict
    )
    requirements: Sequence[Requirement] = ()
    tolerances: Mapping[str, Tolerance] = field(default_factory=dict)
    by: Sequence[str] | str = ()
    across: Across | None = None
    tolerance_array: np.ndarray = field(
        init=False, repr=False, compare=False, metadata={"fingerprint": False}
    )
    """numpy.ndarray: The tolerances of the objectives (zero if not given), in
    their order: their values in their units (first row), and their shares
    of the values of the objectives (second row)."""

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
        statistics = _named_statistics(self.statistics)
        requirements = tuple(self.requirements)
        for requirement in requirements:
            if not isinstance(requirement, Requirement):
                raise TypeError(f"{requirement!r} is not a Requirement")
        tolerances = dict(self.tolerances)
        tolerance_array = _tolerance_array(objectives, tolerances)
        by = tuple(tables.names(self.by))
        object.__setattr__(self, "objectives", objectives)
        object.__setattr__(self, "statistics", statistics)
        object.__setattr__(self, "requirements", requirements)
        object.__setattr__(self, "tolerances", tolerances)
        object.__setattr__(self, "tolerance_array", tolerance_array)
        object.__setattr__(self, "by", by)
        if self.across is not None:
            if not isinstance(self.across, Across):
                raise TypeError("across must be an Across (or None)")
            if self.across.column in by:
                raise ValueError(
                    f"criteria cannot both group and combine {self.across.column!r}"
                )
            unknown = sorted(set(self.across.rules) - set(self.names))
            if unknown:
                raise ValueError(f"rules of names that are not statistics: {unknown}")
            _ = self.combiners  # checks that each statistic has a worst, if needed

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

    @classmethod
    def of(cls, criteria: Criteria | Mapping[str, str]) -> Criteria:
        """
        Gets criteria given as criteria, or as a mapping of objectives to their
        directions (the criteria of the means of the metrics of their names).

        Args:
            criteria (Criteria | Mapping[str, str]): The criteria, or the
                objectives and their directions (`"min"` or `"max"`).

        Returns:
            Criteria: the criteria

        Raises:
            TypeError: if the criteria are neither.
        """
        if isinstance(criteria, Criteria):
            return criteria
        if isinstance(criteria, Mapping):
            return cls(criteria)
        raise TypeError(
            "criteria must be Criteria or a mapping of objectives to their directions"
        )

    @property
    def resolved(self) -> dict[str, Statistic]:
        """dict[str, Statistic]: The statistic of each name (see `names` and
        `statistic`)."""
        return {name: self.statistic(name) for name in self.names}

    @property
    def signs(self) -> np.ndarray:
        """numpy.ndarray: The sign of each objective that makes it a cost (to
        minimize): 1 for an objective to minimize, and -1 for one to
        maximize."""
        return np.array(
            [
                1.0 if direction == "min" else -1.0
                for direction in self.objectives.values()
            ]
        )

    def costs(self, values: Mapping[str, Any]) -> np.ndarray:
        """
        Gets the objectives as costs (to minimize, see `signs`), from the
        values of their statistics (by name, such as the columns of a
        summary): maximized objectives are negated.

        Args:
            values (Mapping[str, Any]): The values of the statistics of
                designs, by name.

        Returns:
            numpy.ndarray: the costs (a row per design, and a column per
                objective, `NaN` where a value is missing)
        """
        columns = [
            np.asarray(values[name], dtype=float) * sign
            for name, sign in zip(self.objectives, self.signs)
        ]
        count = len(columns[0]) if columns else 0
        return np.column_stack(columns) if count else np.empty((0, len(columns)))

    def violation(self, values: Mapping[str, Any], count: int) -> np.ndarray:
        """
        Gets the total violation of the requirements by designs (the sum of
        their violations, see `Requirement.violation`), from the values of
        their statistics (by name, such as the columns of a summary).

        Args:
            values (Mapping[str, Any]): The values of the statistics of
                designs, by name.
            count (int): The number of designs.

        Returns:
            numpy.ndarray: the violation of each design (0 if it meets every
                requirement, and `NaN` where a statistic of a requirement is)
        """
        total = np.zeros(count)
        for requirement in self.requirements:
            statistic = np.asarray(values[requirement.name], dtype=float)
            total = total + requirement.violation(statistic)
        return total

    @property
    def combiners(self) -> dict[str, str]:
        """dict[str, str]: How each statistic is combined across the groups of
        `across` (see `Across`), by name: by its largest (`"max"`) or smallest
        (`"min"`) value, its worst, or by its `"mean"` (empty without
        `across`)."""
        across = self.across
        if across is None:
            return {}
        worst: dict[str, set[str]] = {}
        for name, direction in self.objectives.items():
            worst.setdefault(name, set()).add("max" if direction == "min" else "min")
        for requirement in self.requirements:
            sense = "min" if requirement.operator in (">=", ">") else "max"
            worst.setdefault(requirement.name, set()).add(sense)
        combiners = {}
        for name in self.names:
            rule = across.rules.get(name, across.rule)
            if rule == "mean" or (name not in worst and name not in across.rules):
                # statistics that are neither objectives nor required (such as
                # others to compute with them) have no worst: their means
                combiners[name] = "mean"
            elif name not in worst:
                raise ValueError(
                    f"statistic {name!r} has no worst across {across.column!r}: it "
                    + "is neither an objective nor required (combine it by its "
                    + "mean, see Across)"
                )
            elif len(worst[name]) > 1:
                raise ValueError(
                    f"the objective and requirements of statistic {name!r} disagree "
                    + "on its worst: give it another name (in the statistics) for each"
                )
            else:
                combiners[name] = next(iter(worst[name]))
        return combiners
