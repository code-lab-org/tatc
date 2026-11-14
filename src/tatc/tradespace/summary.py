"""
Summaries of the results of trade studies: the statistics of each design over
its replicates (see `Statistic`), its judgment by criteria (see `Criteria`),
and the probability that it is efficient, given the uncertainty of its
replicates (by a bootstrap).

Designs are judged by their feasible replicates: the statistics of a metric
are over the replicates in which a design is feasible, and the share of its
replicates that are feasible is the statistic `feasible` (over every
replicate of the design), which requirements may name. Designs without a
feasible replicate are infeasible, and are left out of summaries.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import replace
from typing import Any

import numpy as np
import pandas as pd

from .criteria import (
    Across,
    Criteria,
    Requirement,
    _as_criteria,
    _combiners,
    _costs_of,
    _get_tolerances,
    _resolve,
    _violation_of,
)
from .dominance import _efficient
from .statistics import (
    Mean,
    StandardDeviation,
    StandardError,
    Statistic,
    _named_statistics,
)
from .study import _STATUS_COLUMNS

_FEASIBLE = "feasible"
"""The statistic (and column of summaries) of the share of a design's
replicates that are feasible, and the column of the results of a study of
whether each record is feasible."""

_JUDGMENT = ("violation", "meets_requirements", "pareto")
"""The columns of the judgment of designs by criteria in their summaries (see
`aggregate`)."""

_SEPARATE = (
    "give the columns that separate them as `by` (of the criteria), or combine "
    + "their groups `across` (see Across)"
)
"""How to separate several records of an evaluation for criteria."""

_SUMMARIZE = (
    "summarize the results of a study (a row per replicate) first, with a row "
    + "per design (see aggregate), or judge designs in groups (see Criteria's "
    + "`by`)"
)
"""How to judge the results of a study (a row per replicate) by criteria."""


def _resamples(count: int, n_resamples: int, seed: Any) -> np.ndarray:
    """
    Gets the replicates of each resample of a bootstrap: a row per resample of
    `count` replicates (their positions) drawn with replacement, from a seed
    (or a random number generator, which draws them).
    """
    if n_resamples < 2:
        raise ValueError("n_resamples must be at least 2")
    return np.random.default_rng(seed).integers(count, size=(n_resamples, count))


def _bootstrap_error(sums: np.ndarray, squares: np.ndarray, n: int) -> np.ndarray:
    """
    Gets the standard errors of statistics from the sums of their values (and
    of their squares) over `n` resamples of a bootstrap: their standard
    deviations over the resamples.
    """
    mean = sums / n
    return np.sqrt(np.maximum(squares - n * mean**2, 0.0) / max(n - 1, 1))


def _feasible_records(results: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Gets the feasible and the infeasible records of results (every record is
    feasible without a column `feasible`).
    """
    if _FEASIBLE not in results.columns:
        return results, results.iloc[:0]
    status = results[_FEASIBLE].astype(bool).to_numpy()
    return results[status], results[~status]


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


def _check_summary(
    summary: pd.DataFrame, by: Sequence[str] = (), hint: str = _SUMMARIZE
) -> None:
    """
    Checks that a data frame is a summary of designs (a row per design, in
    each group of the columns `by`), not the results of a study (a row per
    replicate), which the functions of summaries would judge as if each
    replicate were a design, raising a `ValueError` with a hint of what to do
    if not.
    """
    if "design" in summary.columns:
        keys = summary[["design", *(c for c in by if c in summary.columns)]]
    elif summary.index.name == "design" and not by:
        keys = summary.index.to_frame()
    else:
        return
    if keys.duplicated().any():
        raise ValueError(f"the data frame has several rows of a design: {hint}")


def _get_matrices(
    group: pd.DataFrame,
    columns: Sequence[str],
    designs: pd.Index | None = None,
    replicates: pd.Index | None = None,
    hint: str = _SEPARATE,
) -> tuple[pd.Index, pd.Index, dict[str, np.ndarray]]:
    """
    Gets the values of metrics of each design (rows) in each replicate
    (columns) of the results of a group, which must have one record per
    evaluation (a design in a replicate, or a `ValueError` with a hint of
    how to separate them): of the designs and replicates of the group, or of
    those given (`NaN` without a record).
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
            + f"replicate): {hint}"
            + (f", such as {varying}" if varying else "")
        )
    if designs is None:
        designs = pd.Index(group["design"].drop_duplicates())
    if replicates is None:
        replicates = pd.Index(sorted(group["replicate"].unique()))
    matrices = {}
    for column in columns:
        table = group.pivot(index="design", columns="replicate", values=column)
        matrices[column] = table.reindex(index=designs, columns=replicates).to_numpy(
            dtype=float
        )
    return designs, replicates, matrices


def _replicates_of(feasible: pd.DataFrame, infeasible: pd.DataFrame) -> pd.Index:
    """Gets the replicates of feasible and infeasible records, in order."""
    return pd.Index(sorted({*feasible["replicate"], *infeasible["replicate"]}))


def _feasibility(
    feasible: pd.DataFrame,
    infeasible: pd.DataFrame,
    designs: pd.Index,
    replicates: pd.Index,
) -> np.ndarray:
    """
    Gets the feasibility of each design (rows) in each replicate (columns): 1
    if it has feasible records, 0 if it has infeasible ones, and `NaN` if it
    has none.
    """
    matrix = np.full((len(designs), len(replicates)), np.nan)
    for records, flag in ((infeasible, 0.0), (feasible, 1.0)):
        rows = designs.get_indexer(records["design"])
        columns = replicates.get_indexer(records["replicate"])
        known = (rows >= 0) & (columns >= 0)
        matrix[rows[known], columns[known]] = flag
    return matrix


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
    return values, _violation_of(values, requirements, n)


def _across_groups(results: pd.DataFrame, across: Across | None) -> list | None:
    """
    Gets the groups of results (values of a column, in order) across which
    criteria combine statistics, if any.
    """
    if across is None:
        return None
    if across.column not in results.columns:
        raise ValueError(
            f"the results have no column {across.column!r} to combine across "
            + "(such as the scenarios of a study, see Study)"
        )
    return list(pd.unique(results[across.column]))


def _attributed(
    infeasible: pd.DataFrame, by: Sequence[str], key: Sequence[Any]
) -> pd.DataFrame:
    """
    Gets the infeasible records of a group of the columns `by` (of the values
    `key`): those of its values, or without a value (such as a sensor, which
    an infeasible evaluation does not give), which are of every group.
    """
    if not by or infeasible.empty:
        return infeasible
    mask = np.ones(len(infeasible), dtype=bool)
    for column, value in zip(by, key):
        values = infeasible[column]
        mask &= ((values == value) | values.isna()).to_numpy()
    return infeasible[mask]


class _Statistics:
    """
    The statistics of the designs of the feasible results of a group (of the
    criteria's `by`), over every replicate or some of them (see `evaluate`),
    from the values of their metrics in each replicate: in each group of a
    column, combined across the groups (see `Across`), if any. Designs
    without feasible results in every group are left out.

    The statistics of the feasibility of replicates (`feasible`, see
    `aggregate`) are over the replicates of the feasible and infeasible
    records of the designs, and, across groups, over those of every group
    together (they are not combined).

    Args:
        results (pandas.DataFrame): The feasible results of a group.
        statistics (Mapping[str, Statistic]): The statistics, by name.
        requirements (Sequence[Requirement]): The requirements.
        across (Across | None): The rule that combines the groups of a
            column, if any.
        combiners (Mapping[str, str] | None): How each statistic is combined
            (see `_combiners`), with a rule.
        groups (Sequence | None): The groups of the column (of every result),
            with a rule.
        infeasible (pandas.DataFrame | None): The infeasible results of the
            group, if any.
    """

    def __init__(
        self,
        results: pd.DataFrame,
        statistics: Mapping[str, Statistic],
        requirements: Sequence[Requirement] = (),
        across: Across | None = None,
        combiners: Mapping[str, str] | None = None,
        groups: Sequence | None = None,
        infeasible: pd.DataFrame | None = None,
    ):
        self.statistics, self.requirements = statistics, requirements
        self.combiners = combiners
        infeasible = results.iloc[:0] if infeasible is None else infeasible
        columns = {c for statistic in statistics.values() for c in statistic.columns}
        metrics = sorted(columns - {_FEASIBLE})
        # statistics of the feasibility alone are over every group together
        self.pooled = {
            name
            for name, statistic in statistics.items()
            if set(statistic.columns) == {_FEASIBLE}
        }
        if across is None:
            self.designs = pd.Index(results["design"].drop_duplicates())
            infeasible = infeasible[infeasible["design"].isin(self.designs)]
            self.replicates = _replicates_of(results, infeasible)
            matrices = _get_matrices(results, metrics, self.designs, self.replicates)[2]
            if _FEASIBLE in columns:
                matrices[_FEASIBLE] = _feasibility(
                    results, infeasible, self.designs, self.replicates
                )
            self.matrices = [matrices]
            self.counts = (
                results.groupby("design")["replicate"]
                .nunique()
                .reindex(self.designs)
                .to_numpy()
            )
            return
        combined = [name for name in statistics if name not in self.pooled]
        self.weights = (
            across.weights_of(groups)
            if any(combiners[name] == "mean" for name in combined)
            else None
        )
        parts = [
            results[(results[across.column] == group).to_numpy()] for group in groups
        ]
        complete = set.intersection(*(set(part["design"]) for part in parts))
        designs = results["design"].drop_duplicates()
        self.designs = pd.Index(designs[designs.isin(complete)])
        infeasible = infeasible[infeasible["design"].isin(self.designs)]
        self.replicates = _replicates_of(results, infeasible)
        self.matrices = []
        for group, part in zip(groups, parts):
            matrices = _get_matrices(part, metrics, self.designs, self.replicates)[2]
            if _FEASIBLE in columns:
                failed = infeasible[(infeasible[across.column] == group).to_numpy()]
                matrices[_FEASIBLE] = _feasibility(
                    part, failed, self.designs, self.replicates
                )
            self.matrices.append(matrices)
        self.counts = (
            results.groupby(["design", across.column])["replicate"]
            .nunique()
            .groupby(level="design")
            .min()
            .reindex(self.designs)
            .to_numpy()
        )

    def _pooled(self, statistic: Statistic, columns: np.ndarray | None) -> np.ndarray:
        """
        Computes a statistic over the replicates (`columns`, or all) of every
        group together.
        """
        matrix = np.hstack([statistic.values_of(m) for m in self.matrices])
        if columns is not None:
            n = len(self.replicates)
            matrix = matrix[
                :, np.concatenate([columns + k * n for k in range(len(self.matrices))])
            ]
        return statistic.reduce(matrix)

    def evaluate(
        self, columns: np.ndarray | None = None
    ) -> tuple[dict[str, np.ndarray], np.ndarray]:
        """
        Computes the statistics of each design (of the replicates `columns`
        of the matrices, or all) and the total violation of the requirements.
        """
        if self.combiners is None:
            return _evaluate(
                self.statistics, self.requirements, self.matrices[0], columns
            )
        combined = {
            name: statistic
            for name, statistic in self.statistics.items()
            if name not in self.pooled
        }
        parts = [_evaluate(combined, (), m, columns)[0] for m in self.matrices]
        values = {}
        for name, statistic in self.statistics.items():
            if name in self.pooled:
                values[name] = self._pooled(statistic, columns)
                continue
            stacked = np.vstack([part[name] for part in parts])
            how = self.combiners[name]
            if how == "max":
                values[name] = stacked.max(axis=0)
            elif how == "min":
                values[name] = stacked.min(axis=0)
            else:
                values[name] = self.weights @ stacked
        return values, _violation_of(values, self.requirements, len(self.designs))


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


def _value_columns(
    feasible: pd.DataFrame,
    statistics: Mapping[str, Statistic],
    by: Sequence[str],
    across: Across | None,
    design_columns: Sequence[str] | None,
) -> list[str]:
    """
    Gets the columns of the values of designs that a summary keeps: those
    given, or those that are the same in every record of each design (but
    the metrics of the statistics, the columns of groups, and status
    columns).
    """
    combined = [] if across is None else [across.column]
    if design_columns is None:
        metrics = {
            column for statistic in statistics.values() for column in statistic.columns
        }
        design_columns = _design_columns(
            feasible, {*_STATUS_COLUMNS, *metrics, *by, *combined}
        )
    return [
        c for c in design_columns if c not in statistics and c not in (*by, *combined)
    ]


def _with_values(
    frame: pd.DataFrame, group: pd.DataFrame, designs: Any, columns: Sequence[str]
) -> None:
    """
    Adds the values of designs (columns of the first record of each design of
    a group of results) to the rows of a summary of the designs.
    """
    first = group.drop_duplicates("design").set_index("design")
    for column in columns:
        frame[column] = first.loc[designs, column].to_numpy()


def _summary_groups(
    results: pd.DataFrame,
    statistics: Mapping[str, Statistic],
    by: Sequence[str],
    requirements: Sequence[Requirement] = (),
    across: Across | None = None,
    combiners: Mapping[str, str] | None = None,
    design_columns: Sequence[str] | None = None,
) -> Iterator[tuple[pd.DataFrame, _Statistics]]:
    """
    Iterates over the groups (of the columns `by`) of the feasible results of
    a study, for summaries of their designs: a data frame of the designs of
    each group (`design`, the values of `by`, and the values of the designs,
    see `aggregate`), and their statistics (see `_Statistics`). Results
    without a feasible record have no groups.
    """
    feasible, infeasible = _feasible_records(results)
    if feasible.empty:
        return
    groups = _across_groups(results, across)
    columns = _value_columns(feasible, statistics, by, across, design_columns)
    for key, positions in _groups(feasible, by):
        group = feasible.iloc[positions]
        failed = _attributed(infeasible, by, key)
        table = _Statistics(
            group, statistics, requirements, across, combiners, groups, failed
        )
        frame = pd.DataFrame({"design": table.designs})
        for name, value in zip(by, key):
            frame[name] = value
        _with_values(frame, group, table.designs, columns)
        yield frame, table


def _criteria_of(
    criteria: Criteria | Mapping[str, str] | Sequence[str] | str,
    statistics: Mapping[str, Statistic | str] | Sequence[str] | None,
    by: Sequence[str] | str,
) -> tuple[Criteria | None, dict[str, Statistic], tuple[str, ...]]:
    """
    Gets the criteria (if any) of the arguments of `aggregate`, the
    statistics to compute (by name), and the columns that group records.
    """
    by = (by,) if isinstance(by, str) else tuple(by)
    extra = _named_statistics(statistics)
    if isinstance(criteria, Mapping) and not set(criteria.values()) <= {"min", "max"}:
        raise ValueError(
            "a mapping of criteria gives the directions of objectives ('min' or "
            + "'max'): give statistics by name as `statistics`"
        )
    if isinstance(criteria, (Criteria, Mapping)):
        if by and isinstance(criteria, Criteria):
            raise ValueError(
                "criteria group records by their own `by`: give it to the criteria"
            )
        judged = _as_criteria(criteria)
        if by:
            judged = replace(judged, by=by)
        if extra:
            clashes = sorted(
                name
                for name in extra
                if name in judged.names and judged.statistic(name) != extra[name]
            )
            if clashes:
                raise ValueError(
                    f"statistics that the criteria name otherwise: {clashes}"
                )
            judged = replace(judged, statistics={**judged.statistics, **extra})
        return judged, _resolve(judged), judged.by
    metrics = [criteria] if isinstance(criteria, str) else list(criteria)
    resolved: dict[str, Statistic] = {name: Mean(name) for name in metrics}
    resolved.update(extra)
    if not resolved:
        raise ValueError(
            "give criteria, the metrics whose means to compute, or statistics"
        )
    return None, resolved, by


def _with_feasibility(statistics: Mapping[str, Statistic]) -> dict[str, Statistic]:
    """
    Gets statistics with the share of the replicates of each design that are
    feasible (`feasible`), unless a statistic has the name.
    """
    if _FEASIBLE in statistics:
        return dict(statistics)
    return {**statistics, _FEASIBLE: Mean(_FEASIBLE)}


def _judge(
    frame: pd.DataFrame,
    values: Mapping[str, np.ndarray],
    violation: np.ndarray,
    criteria: Criteria,
) -> None:
    """
    Adds the judgment of designs by criteria to the rows of a summary of
    them, from the values of their statistics and the violation of the
    criteria's requirements: the `violation`, whether each design meets the
    requirements (`meets_requirements`), and whether it is efficient among
    those that do (`pareto`).
    """
    meets = violation == 0
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerances)
    violations, requirements, efficient = _JUDGMENT
    frame[violations] = violation
    frame[requirements] = meets
    frame[efficient] = _efficient(_costs_of(values, criteria), meets, tolerances)


def _assemble(
    frame: pd.DataFrame,
    table: _Statistics,
    values: Mapping[str, np.ndarray],
    violation: np.ndarray,
    criteria: Criteria | None,
    extra: Mapping[str, np.ndarray] | None = None,
) -> pd.DataFrame:
    """
    Adds to the rows of the designs of a group (of `_summary_groups`) their
    statistics, `extra` columns (if any), the number of their feasible
    `replicates`, the share of their replicates that are feasible
    (`feasible`), and their judgment by criteria (if any, see `_judge`).
    """
    for name, value in values.items():
        if name != _FEASIBLE:
            frame[name] = value
    for name, value in (extra or {}).items():
        frame[name] = value
    frame["replicates"] = table.counts
    frame[_FEASIBLE] = values[_FEASIBLE]
    if criteria is not None:
        _judge(frame, values, violation, criteria)
    return frame


def _summaries(
    frames: Sequence[pd.DataFrame],
    by: Sequence[str],
    statistics: Mapping[str, Statistic],
    criteria: Criteria | None,
    extra: Sequence[str] = (),
    after: Sequence[str] = (),
) -> pd.DataFrame:
    """
    Gets a summary of the frames of its groups (see `_assemble`), or an empty
    one of its columns (without the values of designs) if there are none,
    such as of results without a feasible record: the statistics, `extra`
    columns, `replicates`, `feasible`, the judgment (with criteria), and the
    columns `after`.
    """
    if frames:
        return pd.concat(frames, ignore_index=True)
    names = [name for name in statistics if name != _FEASIBLE]
    dtypes = {
        "design": object,
        **{column: object for column in by},
        **{name: float for name in [*names, *extra]},
        "replicates": int,
        _FEASIBLE: float,
    }
    if criteria is not None:
        dtypes.update(zip(_JUDGMENT, (float, bool, bool)))
    dtypes.update({column: float for column in after})
    return pd.DataFrame(
        {column: pd.Series(dtype=dtype) for column, dtype in dtypes.items()}
    )


def aggregate(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str] | Sequence[str] | str = (),
    *,
    statistics: Mapping[str, Statistic | str] | Sequence | None = None,
    by: Sequence[str] | str = (),
    design_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Aggregates the results of a study (see `Study.run`, a row per replicate)
    into a summary of its designs (a row per design, and per group of the
    columns `by`, such as a sensor): statistics over the feasible replicates
    of each design, the number of these `replicates`, the share of its
    replicates that are feasible (`feasible`), and the values of the design
    (see `Design.flat`), from which `DesignSpace.designs_of` rebuilds the
    designs.

    The statistics are those of criteria (see `Criteria`, or a mapping of
    objectives to their directions, which are the means of the metrics of
    their names), or the means of metrics (a list of their names), and other
    `statistics`, by name (a statistic, or a metric given as a string, whose
    mean it is), or a list of statistics (or metrics), named by their
    conventional names (see `Statistic.name`)::

        aggregate(results, {"satellites": "min", "sufficient": "max"})
        aggregate(results, ["sufficient", "accurate"])
        aggregate(results, statistics=[Quantile("cost", 0.7)])  # cost p70
        aggregate(results, COST, statistics=["cost development"])

    With criteria, the summary also gives the `violation` of their
    requirements (the sum of the violations, see `Requirement`, 0 if the
    design meets them all), whether the design meets them
    (`meets_requirements`), and whether it is efficient among the designs
    that do (`pareto`, see `pareto_front`, in each group of the criteria's
    `by`; a design with a missing objective, `NaN`, is not). Criteria
    `across` a column combine the statistics of each design across its
    groups (see `Across`), of which `replicates` is the fewest of any group.

    Designs are judged by their feasible replicates. `feasible` is a
    statistic over every replicate of a design (the mean of the feasibility
    of each, 1 or 0, see `Statistic`), and across the groups of criteria
    `across` a column, over the replicates of every group together, which
    requirements may name: `Requirement("feasible", ">=", 1)` requires every
    replicate to be feasible. In groups of columns `by` that infeasible
    records lack (such as a sensor, which an infeasible evaluation does not
    give), the infeasible replicates of a design count in each of its groups.
    Designs without a feasible replicate (in a group) are infeasible, and are
    left out: a summary of results without a feasible record is empty.

    The results must have one record per evaluation (a design in a
    replicate) in each group (and in each group of the column of the
    criteria's `across`).

    Args:
        results (pandas.DataFrame): The results.
        criteria (Criteria | Mapping[str, str] | Sequence[str] | str): The
            criteria (or objectives and their directions), or the metrics
            whose means to compute.
        statistics (Mapping[str, Statistic | str] | Sequence | None): Other
            statistics, by name, or a sequence of statistics (or metrics,
            whose means they are), named by their conventional names.
        by (Sequence[str] | str): The columns that group records of an
            evaluation (those of the criteria, if they have any).
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every
            replicate of a design (except the metrics).

    Returns:
        pandas.DataFrame: a row per feasible design (and group)
    """
    judged, resolved, by = _criteria_of(criteria, statistics, by)
    statistics = _with_feasibility(resolved)
    across, combiners, requirements = None, None, ()
    if judged is not None:
        requirements, across = judged.requirements, judged.across
        combiners = None if across is None else _combiners(judged)
    frames = [
        _assemble(frame, table, *table.evaluate(), judged)
        for frame, table in _summary_groups(
            results, statistics, by, requirements, across, combiners, design_columns
        )
    ]
    return _summaries(frames, by, statistics, judged)


def summarize(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str] | Sequence[str] | str = (),
    *,
    statistics: Mapping[str, Statistic | str] | Sequence | None = None,
    by: Sequence[str] | str = (),
    design_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Summarizes the feasible results of a study over their replicates (see
    `aggregate`), with the spread of the means of metrics: for each design
    (and each group of the columns `by`, such as a sensor), the statistics
    of criteria (or the means of a list of metrics) and of other
    `statistics`, the standard deviation (`<name> std`) and standard error
    (`<name> se`) between replicates of each that is the mean of a metric
    (a column, or a metric computed from others, such as a `Utility`), the
    number of feasible replicates (`replicates`), the share of the replicates
    that are feasible (`feasible`), and the values of the design. Of criteria
    (or a mapping of objectives to their directions), the summary is theirs
    (see `aggregate`: the violation of their requirements, and whether each
    design is efficient, `pareto`)::

        summary = summarize(results, {"satellites": "min", "sufficient": "max"})
        front = summary[summary.pareto]
        summarize(results, ACCURACY, statistics=["sufficient"])  # and its mean

    Args:
        results (pandas.DataFrame): The results (see `Study.run`).
        criteria (Criteria | Mapping[str, str] | Sequence[str] | str): The
            criteria (or objectives and their directions), or the metrics
            whose means to summarize.
        statistics (Mapping[str, Statistic | str] | Sequence | None): Other
            statistics (see `aggregate`).
        by (Sequence[str] | str): Other columns that group records of an
            evaluation (those of criteria, if they have any).
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every
            replicate of a design (except the metrics).

    Returns:
        pandas.DataFrame: the summary, a row per feasible design (and group)
    """
    judged, resolved, by = _criteria_of(criteria, statistics, by)
    summarized: dict[str, Statistic] = {}
    spread: dict[str, Statistic] = {}
    for name, statistic in resolved.items():
        summarized[name] = statistic
        if type(statistic) is Mean and _FEASIBLE not in statistic.columns:
            spread[f"{name} std"] = StandardDeviation(statistic.column)
            spread[f"{name} se"] = StandardError(statistic.column)
            summarized.update(
                {
                    f"{name} std": spread[f"{name} std"],
                    f"{name} se": spread[f"{name} se"],
                }
            )
    if judged is None:
        # each mean of a metric followed by its spread
        return aggregate(
            results, statistics=summarized, by=by, design_columns=design_columns
        )
    return aggregate(results, judged, statistics=spread, design_columns=design_columns)


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

    Designs are judged by their feasible replicates, and the replicates in
    which a design is infeasible are drawn too: the statistics of a resample
    are those of its feasible replicates, and the share of its replicates
    that are feasible (`feasible`, see `aggregate`) is that of the drawn
    replicates.

    Args:
        results (pandas.DataFrame): The results of a study (see
            `Study.run`), with a record per replicate of each design (and
            group of the criteria's `by`).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (metrics whose means they are) and their directions
            (`"min"` or `"max"`).
        n_resamples (int): The number of resamples (at least 2).
        seed (int): The seed of the resamples.
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every
            replicate of a design (see `aggregate`).

    Returns:
        pandas.DataFrame: a row per feasible design (and group), empty
            without one: the summary of the criteria over every replicate
            (see `aggregate`: the values of the design, the statistics, the
            standard error of each objective, `<objective> se`, its standard
            deviation over the resamples, the number of feasible
            `replicates`, the share of the replicates that are `feasible`,
            the `violation` of the requirements, whether the design meets
            them, `meets_requirements`, and whether it is efficient among
            those that do, `pareto`), the probability that it is efficient
            (`pareto_probability`), and the probability that it meets the
            requirements (`requirements_probability`)
    """
    criteria = _as_criteria(criteria)
    statistics = _with_feasibility(_resolve(criteria))
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerances)
    combiners = None if criteria.across is None else _combiners(criteria)
    errors = [f"{name} se" for name in criteria.objectives]
    rng = np.random.default_rng(seed)
    frames = []
    for frame, table in _summary_groups(
        results,
        statistics,
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
        bootstrap = {
            error: _bootstrap_error(sums[name], squares[name], n_resamples)
            for name, error in zip(criteria.objectives, errors)
        }
        _assemble(frame, table, values, violation, criteria, bootstrap)
        frame["pareto_probability"] = counts / n_resamples
        frame["requirements_probability"] = met / n_resamples
        frames.append(frame)
    return _summaries(
        frames,
        criteria.by,
        statistics,
        criteria,
        extra=errors,
        after=["pareto_probability", "requirements_probability"],
    )
