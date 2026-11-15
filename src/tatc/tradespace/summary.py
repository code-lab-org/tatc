"""
Summaries of the results of trade studies: the statistics of each design over
its replicates (see `Statistic`) and their spreads, its judgment by criteria
(see `Criteria`), and the probability that it is efficient, given the
uncertainty of its replicates (by a bootstrap).

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
    _as_criteria,
    _combiners,
    _costs,
    _resolved,
    _statistic,
    _tolerances,
    _violation,
    _weights_of,
)
from .dominance import efficient
from .statistics import (
    Mean,
    Statistic,
    _named_statistics,
    _standard_deviation,
    _standard_error,
)
from .tables import (
    FEASIBLE,
    JUDGMENT,
    draw_resamples,
    feasible_records,
    groups,
    matrices,
    names,
    value_columns,
    with_values,
)

_NEGLIGIBLE = 1e-12
"""The share of the sum of the magnitudes of its terms below which a variance
combined from the covariances of groups is zero (rather than the error of
rounding of terms that cancel)."""


def _bootstrap_error(
    sums: np.ndarray, squares: np.ndarray, n: np.ndarray
) -> np.ndarray:
    """
    Gets the standard errors of statistics from the sums of their values (and
    of their squares) over the resamples of a bootstrap in which each is
    known (`n`, of each): their standard deviations over these resamples
    (`NaN` with fewer than two).
    """
    n = np.asarray(n, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sums / n
        variance = np.maximum(squares - n * mean**2, 0.0) / (n - 1)
    return np.where(n >= 2, np.sqrt(variance), np.nan)


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
    values: Mapping[str, np.ndarray],
    columns: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """
    Computes statistics of each design from the values of the columns of the
    results in its replicates (of the replicates `columns`, or all).
    """
    computed = {}
    for name, statistic in statistics.items():
        matrix = statistic.values_of(values)
        computed[name] = statistic.reduce(
            matrix if columns is None else matrix[:, columns]
        )
    return computed


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


def _covariances(
    a: np.ndarray, b: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Gets the sample covariance of the values of each design (rows) in two
    groups over the replicates (columns) in which both have values, and the
    numbers of these replicates and of the values of each group (the
    covariance is `NaN` with fewer than two shared replicates).
    """
    both = ~np.isnan(a) & ~np.isnan(b)
    shared = both.sum(axis=1)
    counts = (~np.isnan(a)).sum(axis=1), (~np.isnan(b)).sum(axis=1)
    x, y = np.where(both, a, 0.0), np.where(both, b, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mx = x.sum(axis=1) / shared
        my = y.sum(axis=1) / shared
        products = np.where(both, (x - mx[:, None]) * (y - my[:, None]), 0.0)
        covariance = products.sum(axis=1) / (shared - 1)
    covariance = np.where(shared > 1, covariance, np.nan)
    return covariance, shared, *counts


def _spread_of_mean(
    parts: Sequence[np.ndarray], weights: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the spread of the weighted mean of the means of the values of each
    design (rows) in groups (`parts`) over their replicates (columns), from
    the covariances of the groups over the replicates they share (none
    between groups without shared replicates): the standard deviation of the
    weighted sum of the values of a replicate, and the standard error of the
    weighted mean of the groups' means.
    """
    if len(parts) == 1:
        # the spread of the replicates of one group
        return _standard_deviation(parts[0]), _standard_error(parts[0])
    # the variance of the weighted sum, and the squared standard error, from
    # terms whose magnitudes tell rounding from a variance
    totals = np.zeros((2, len(parts[0])))
    magnitudes = np.zeros((2, len(parts[0])))
    for s, a in enumerate(parts):
        for t, b in enumerate(parts):
            covariance, shared, n_a, n_b = _covariances(a, b)
            if s != t:
                covariance = np.nan_to_num(covariance)
            weighted = weights[s] * weights[t] * covariance
            with np.errstate(invalid="ignore", divide="ignore"):
                terms = np.vstack([weighted, weighted * shared / (n_a * n_b)])
            totals += terms
            magnitudes += np.abs(terms)
    totals = np.where(totals <= _NEGLIGIBLE * magnitudes, 0.0, totals)
    return np.sqrt(totals[0]), np.sqrt(totals[1])


class _Statistics:
    """
    The statistics of the designs of the feasible results of a group (of the
    criteria's `by`), over every replicate or some of them (see `evaluate`),
    from the values of their metrics in each replicate: in each group of a
    column, combined across the groups (see `Across`), or, without criteria
    across a column, in one group of every result. Designs without feasible
    results in every group are left out.

    The statistics of the feasibility of replicates (`feasible`, see
    `summarize`) are over the replicates of the feasible and infeasible
    records of the designs, and, across groups, over those of every group
    together (they are not combined).

    Args:
        results (pandas.DataFrame): The feasible results of a group.
        statistics (Mapping[str, Statistic]): The statistics, by name.
        criteria (Criteria | None): The criteria (of the requirements, and
            of the rule that combines the groups of a column), if any.
        levels (Sequence | None): The groups of the column (of every result)
            that the criteria combine, if any.
        infeasible (pandas.DataFrame | None): The infeasible results of the
            group, if any.
    """

    def __init__(
        self,
        results: pd.DataFrame,
        statistics: Mapping[str, Statistic],
        criteria: Criteria | None = None,
        levels: Sequence | None = None,
        infeasible: pd.DataFrame | None = None,
    ):
        self.statistics, self.criteria = statistics, criteria
        across = None if criteria is None else criteria.across
        infeasible = results.iloc[:0] if infeasible is None else infeasible
        columns = {c for statistic in statistics.values() for c in statistic.columns}
        metrics = sorted(columns - {FEASIBLE})
        # statistics of the feasibility alone are over every group together
        self.pooled = {
            name
            for name, statistic in statistics.items()
            if set(statistic.columns) == {FEASIBLE}
        }
        if across is None:
            # one group of every result, whose statistics are its own
            column, self.combiners = None, {}
            parts, failures = [results], [infeasible]
            self.weights = np.ones(1)
        else:
            column, self.combiners = across.column, _combiners(criteria)
            parts = [results[(results[column] == group).to_numpy()] for group in levels]
            failures = [
                infeasible[(infeasible[column] == group).to_numpy()] for group in levels
            ]
            self.weights = _weights_of(across, levels)
        designs = results["design"].drop_duplicates()
        if len(parts) > 1:
            # the designs with feasible results in every group
            complete = set.intersection(*(set(part["design"]) for part in parts))
            designs = designs[designs.isin(complete)]
        self.designs = pd.Index(designs)
        failures = [failed[failed["design"].isin(self.designs)] for failed in failures]
        self.replicates = _replicates_of(results, pd.concat(failures))
        self.matrices = []
        for part, failed in zip(parts, failures):
            values = matrices(part, metrics, self.designs, self.replicates)[2]
            values[FEASIBLE] = _feasibility(part, failed, self.designs, self.replicates)
            self.matrices.append(values)
        # the fewest replicates of a design in any group
        counts = results.groupby(["design", *([] if column is None else [column])])[
            "replicate"
        ].nunique()
        if column is not None:
            counts = counts.groupby(level="design").min()
        self.counts = counts.reindex(self.designs).to_numpy()

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

    def violation(self, values: Mapping[str, np.ndarray]) -> np.ndarray:
        """Gets the total violation of the requirements of the criteria."""
        if self.criteria is None:
            return np.zeros(len(self.designs))
        return _violation(self.criteria, values, len(self.designs))

    def evaluate(
        self, columns: np.ndarray | None = None
    ) -> tuple[dict[str, np.ndarray], np.ndarray]:
        """
        Computes the statistics of each design (of the replicates `columns`
        of the matrices, or all) and the total violation of the requirements.
        """
        combined = {
            name: statistic
            for name, statistic in self.statistics.items()
            if name not in self.pooled
        }
        parts = [_evaluate(combined, m, columns) for m in self.matrices]
        values = {}
        for name, statistic in self.statistics.items():
            if name in self.pooled:
                values[name] = self._pooled(statistic, columns)
                continue
            if len(parts) == 1:
                values[name] = parts[0][name]
                continue
            stacked = np.vstack([part[name] for part in parts])
            how = self.combiners.get(name, "mean")
            if how == "max":
                values[name] = stacked.max(axis=0)
            elif how == "min":
                values[name] = stacked.min(axis=0)
            else:
                values[name] = self.weights @ stacked
        return values, self.violation(values)

    def judged(self, columns: np.ndarray) -> np.ndarray:
        """
        Gets whether each design has a feasible record among some replicates
        (`columns`) in every group, without which it is not judged by them.
        """
        found = np.ones(len(self.designs), dtype=bool)
        for values in self.matrices:
            found &= np.any(values[FEASIBLE][:, columns] == 1, axis=1)
        return found

    def recorded(self, columns: np.ndarray) -> np.ndarray:
        """
        Gets whether each design has a record (feasible or not) among some
        replicates (`columns`), in any group.
        """
        found = np.zeros(len(self.designs), dtype=bool)
        for values in self.matrices:
            found |= np.any(~np.isnan(values[FEASIBLE][:, columns]), axis=1)
        return found

    def spread(self, name: str, statistic: Statistic) -> tuple[np.ndarray, np.ndarray]:
        """
        Computes the spread of a mean of a metric (of a name) over the
        replicates of each design: its standard deviation and standard error
        (see `summarize`), and, across groups, those of the statistic of the
        groups that it combines into.
        """
        parts = [statistic.values_of(m) for m in self.matrices]
        how = self.combiners.get(name, "mean")
        if how in ("max", "min"):
            # the spread of the group of the worst value
            means = np.vstack([statistic.reduce(part) for part in parts])
            deviations = np.vstack([_standard_deviation(part) for part in parts])
            errors = np.vstack([_standard_error(part) for part in parts])
            filled = np.where(
                np.isnan(means), -np.inf if how == "max" else np.inf, means
            )
            chosen = filled.argmax(axis=0) if how == "max" else filled.argmin(axis=0)
            rows = np.arange(len(self.designs))
            missing = np.isnan(means).any(axis=0)
            return (
                np.where(missing, np.nan, deviations[chosen, rows]),
                np.where(missing, np.nan, errors[chosen, rows]),
            )
        return _spread_of_mean(parts, self.weights)


def _summary_groups(
    results: pd.DataFrame,
    statistics: Mapping[str, Statistic],
    by: Sequence[str],
    criteria: Criteria | None = None,
    design_columns: Sequence[str] | None = None,
) -> Iterator[tuple[pd.DataFrame, _Statistics]]:
    """
    Iterates over the groups (of the columns `by`) of the feasible results of
    a study, for summaries of their designs: a data frame of the designs of
    each group (`design`, the values of `by`, and the values of the designs,
    see `summarize`: those of `design_columns`, or, if `None`, those that are
    the same in every feasible record of a design), and their statistics
    (see `_Statistics`). Results without a feasible record have no groups.
    """
    feasible, infeasible = feasible_records(results)
    if feasible.empty:
        return
    across = None if criteria is None else criteria.across
    combined = _across_groups(results, across)
    grouped = [*by, *([] if across is None else [across.column])]
    columns = value_columns(results, statistics, grouped, design_columns)
    for key, positions in groups(feasible, by):
        group = feasible.iloc[positions]
        failed = _attributed(infeasible, by, key)
        table = _Statistics(group, statistics, criteria, combined, failed)
        frame = pd.DataFrame({"design": table.designs})
        for name, value in zip(by, key):
            frame[name] = value
        with_values(frame, group, table.designs, columns)
        yield frame, table


def _criteria_of(
    criteria: Criteria | Mapping[str, str] | None,
    statistics: Mapping[str, Statistic | str] | Sequence | None,
    by: Sequence[str] | str,
) -> tuple[Criteria | None, dict[str, Statistic], tuple[str, ...]]:
    """
    Gets the criteria (if any) of the arguments of `summarize`, the
    statistics to compute (by name), and the columns that group records.
    """
    by = tuple(names(by))
    extra = _named_statistics(statistics)
    if criteria is None:
        if not extra:
            raise ValueError(
                "give criteria, or the statistics (or metrics, whose means they "
                + "are) to summarize as `statistics`"
            )
        return None, extra, by
    if isinstance(criteria, Mapping) and not all(
        isinstance(value, str) for value in criteria.values()
    ):
        # statistics given as criteria (strings that are not directions, such
        # as "minimize", are misspelled directions, which Criteria rejects)
        raise ValueError(
            "a mapping of criteria gives the directions of objectives ('min' or "
            + "'max'): give statistics by name as `statistics`"
        )
    judged = _as_criteria(criteria)
    if by:
        if isinstance(criteria, Criteria):
            raise ValueError(
                "criteria group records by their own `by`: give it to the criteria"
            )
        judged = replace(judged, by=by)
    if extra:
        clashes = sorted(
            name
            for name in extra
            if name in judged.names and _statistic(judged, name) != extra[name]
        )
        if clashes:
            raise ValueError(f"statistics that the criteria name otherwise: {clashes}")
        judged = replace(judged, statistics={**judged.statistics, **extra})
    return judged, _resolved(judged), judged.by


def _with_feasibility(statistics: Mapping[str, Statistic]) -> dict[str, Statistic]:
    """
    Gets statistics with the share of the replicates of each design that are
    feasible (`feasible`), unless a statistic has the name.
    """
    if FEASIBLE in statistics:
        return dict(statistics)
    return {**statistics, FEASIBLE: Mean(FEASIBLE)}


def _spreads(statistics: Mapping[str, Statistic]) -> dict[str, list[str]]:
    """
    Gets the columns of the spreads of the means of metrics (a column, or a
    metric computed from others, such as a `Utility`) among statistics (by
    name): their standard deviations and standard errors (`<name> std` and
    `<name> se`), but those that statistics name.
    """
    return {
        name: [
            column
            for column in (f"{name} std", f"{name} se")
            if column not in statistics
        ]
        for name, statistic in statistics.items()
        if type(statistic) is Mean and FEASIBLE not in statistic.columns
    }


def _ordered(
    statistics: Mapping[str, Statistic], added: Mapping[str, Sequence[str]]
) -> list[str]:
    """Gets the columns of statistics, each followed by those `added` to it
    (such as its spreads)."""
    return [column for name in statistics for column in (name, *added.get(name, ()))]


def _with_spreads(
    values: dict[str, np.ndarray],
    table: _Statistics,
    statistics: Mapping[str, Statistic],
    spreads: Mapping[str, Sequence[str]],
) -> None:
    """Adds the spreads of means (see `_spreads`) to the values of a group."""
    for name, columns in spreads.items():
        deviation, error = table.spread(name, statistics[name])
        for column, value in ((f"{name} std", deviation), (f"{name} se", error)):
            if column in columns:
                values[column] = value


def _check_names(
    statistics: Mapping[str, Statistic],
    results: pd.DataFrame,
    added: Sequence[str] = (),
) -> None:
    """
    Checks that no statistic but the mean of a metric (of its name), and no
    column `added` to a summary (such as a standard error), has the name of
    a column of the results, which it would be taken for (such as a
    probability over the replicates for the share of the histories of a
    replicate that meet a requirement, see `summarize_histories`).
    """
    clashes = sorted(
        {
            *(
                name
                for name, statistic in statistics.items()
                if name in results.columns and statistic != Mean(name)
            ),
            *(name for name in added if name in results.columns),
        }
    )
    if clashes:
        raise ValueError(
            f"statistics with the names of columns of the results: {clashes} (a "
            + "statistic of a column's name is its mean): give them other names"
        )


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
    violations, requirements, front = JUDGMENT
    frame[violations] = violation
    frame[requirements] = meets
    frame[front] = efficient(_costs(criteria, values), meets, _tolerances(criteria))


def _assemble(
    frame: pd.DataFrame,
    table: _Statistics,
    values: Mapping[str, np.ndarray],
    violation: np.ndarray,
    criteria: Criteria | None,
) -> pd.DataFrame:
    """
    Adds to the rows of the designs of a group (of `_summary_groups`) their
    statistics (and the other columns of `values`, in order), the number of
    their feasible `replicates`, the share of their replicates that are
    feasible (`feasible`), and their judgment by criteria (if any, see
    `_judge`).
    """
    for name, value in values.items():
        if name != FEASIBLE:
            frame[name] = value
    frame["replicates"] = table.counts
    frame[FEASIBLE] = values[FEASIBLE]
    if criteria is not None:
        _judge(frame, values, violation, criteria)
    return frame


def _summaries(
    frames: Sequence[pd.DataFrame],
    by: Sequence[str],
    columns: Sequence[str],
    criteria: Criteria | None,
    after: Sequence[str] = (),
) -> pd.DataFrame:
    """
    Gets a summary of the frames of its groups (see `_assemble`), or an empty
    one of its columns (without the values of designs) if there are none,
    such as of results without a feasible record: the `columns` of
    statistics, `replicates`, `feasible`, the judgment (with criteria), and
    the columns `after`.
    """
    if frames:
        return pd.concat(frames, ignore_index=True)
    dtypes = {
        "design": object,
        **{column: object for column in by},
        **{name: float for name in columns if name != FEASIBLE},
        "replicates": int,
        FEASIBLE: float,
    }
    if criteria is not None:
        dtypes.update(zip(JUDGMENT, (float, bool, bool)))
    dtypes.update({column: float for column in after})
    return pd.DataFrame(
        {column: pd.Series(dtype=dtype) for column, dtype in dtypes.items()}
    )


def _without_constant(summary: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """
    Removes the columns of spreads (such as standard errors) of a summary
    that are zero (or `NaN`) for every design, such as of a metric that is
    the same in every replicate of each design (a summary without designs
    keeps them).
    """
    if summary.empty:
        return summary
    constant = [
        column
        for column in columns
        if column in summary.columns
        and (summary[column].fillna(0.0).to_numpy(dtype=float) == 0).all()
    ]
    return summary.drop(columns=constant)


def _summarize(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str] | None = None,
    *,
    statistics: Mapping[str, Statistic | str] | Sequence | None = None,
    by: Sequence[str] | str = (),
    spread: bool = False,
    design_values: bool = True,
) -> pd.DataFrame:
    """
    Summarizes results by the statistics of their designs (see `summarize`),
    as the steps of trade studies judge designs: without spreads by default
    (the summaries that steps give have them, as `summarize` does), and,
    without `design_values`, without the values of the designs (only their
    keys, `design`), which saves finding them.
    """
    judged, resolved, by = _criteria_of(criteria, statistics, by)
    statistics = _with_feasibility(resolved)
    spreads = _spreads(statistics) if spread else {}
    added = [column for columns in spreads.values() for column in columns]
    _check_names(statistics, results, added)
    order = _ordered(statistics, spreads)
    frames = []
    for frame, table in _summary_groups(
        results, statistics, by, judged, None if design_values else []
    ):
        values, violation = table.evaluate()
        _with_spreads(values, table, statistics, spreads)
        ordered = {name: values[name] for name in order}
        frames.append(_assemble(frame, table, ordered, violation, judged))
    summary = _summaries(frames, by, order, judged)
    return _without_constant(summary, added)


def summarize(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str] | None = None,
    *,
    statistics: Mapping[str, Statistic | str] | Sequence | None = None,
    by: Sequence[str] | str = (),
    spread: bool = True,
) -> pd.DataFrame:
    """
    Summarizes the results of a study (see `Study.run`, a row per replicate)
    with a row per design (and per group of the columns `by`, such as a
    sensor): statistics over the feasible replicates of each design, the
    spread of each that is the mean of a metric (with `spread`), the number
    of these `replicates`, the share of its replicates that are feasible
    (`feasible`), and the values of the design (see `Design.flat`), from
    which `DesignSpace.designs_of` rebuilds the designs.

    The statistics are those of criteria (see `Criteria`, or a mapping of
    objectives to their directions, which are the means of the metrics of
    their names), and other `statistics`: a sequence of statistics or
    metrics (a metric given as a string, or a `Utility`, is its mean),
    named by their conventional names (see `Statistic.name`), or a mapping
    of names to them::

        summarize(results, {"satellites": "min", "sufficient": "max"})
        summarize(results, statistics=["sufficient", "accurate"])
        summarize(results, statistics=[Quantile("cost", 0.7)])  # cost p70
        summarize(results, COST, statistics=["cost development"])

    A statistic of the name of a column of the results is the mean of the
    column: any other statistic of the name raises an error, as its values
    would be taken for those of the column.

    The spread of a mean of a metric (a column, or a metric computed from
    others, such as a `Utility`) follows it: the standard deviation
    (`<name> std`) and the standard error (`<name> se`) between the
    replicates. A spread that is zero for every design (of a metric that is
    the same in every replicate of each design, such as a number of
    satellites, or of designs with one replicate) is left out (but those of
    a summary without designs).

    With criteria, the summary also gives the `violation` of their
    requirements (the sum of the violations, see `Requirement`, 0 if the
    design meets them all), whether the design meets them
    (`meets_requirements`), and whether it is efficient among the designs
    that do (`pareto`, see `pareto_front`, in each group of the criteria's
    `by`; a design with a missing objective, `NaN`, is not)::

        summary = summarize(results, criteria)
        front = summary[summary.pareto]

    Criteria `across` the groups of a column (such as scenarios, see
    `Across`) combine the statistics of each design across its groups, of
    which `replicates` is the fewest of any group. The spread of a combined
    statistic is that of the worst group's value (of the worst), or, of a
    weighted mean, the standard deviation of the weighted sum of the values
    of a replicate in each group and the standard error of the weighted mean
    of the groups' means, from the covariances of the groups over the
    replicates they share (the same replicate in every group, in which it
    has the same uncertain parameters, see `Study`), so that what the groups
    share cancels. Groups without shared replicates are independent: the
    standard error is then the square root of the weighted sum of the
    squared standard errors of the groups (by the squared weights).

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
        criteria (Criteria | Mapping[str, str] | None): The criteria, or the
            objectives (metrics whose means they are) and their directions
            (`"min"` or `"max"`), if any.
        statistics (Mapping[str, Statistic | str] | Sequence | None): Other
            statistics: a sequence of statistics (or metrics, whose means
            they are), named by their conventional names, or a mapping of
            names to them.
        by (Sequence[str] | str): The columns that group records of an
            evaluation (those of the criteria, if they have any).
        spread (bool): `True` to give the spreads of means.

    Returns:
        pandas.DataFrame: a row per feasible design (and group)

    Raises:
        TypeError: if the criteria are neither criteria nor a mapping of
            directions (give metrics as `statistics`).
        ValueError: if a statistic other than the mean of a column has the
            column's name, or the results have several records of an
            evaluation (in a group).
    """
    return _summarize(results, criteria, statistics=statistics, by=by, spread=spread)


def pareto_probability(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str],
    *,
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
    after `refine`): replicates are drawn from all of them, each design's
    statistics are those of its drawn replicates, and a design none of whose
    replicates is drawn keeps the statistics of all of them. With criteria
    `across` the scenarios of a study (see `Across`), a drawn replicate is
    drawn in every scenario, in which it has the same uncertain parameters.

    Designs are judged by their feasible replicates, and the replicates in
    which a design is infeasible are drawn too: the statistics of a resample
    are those of the drawn replicates in which the design is feasible, and
    the share of its replicates that are feasible (`feasible`, see
    `summarize`) is that of the drawn replicates. A design whose drawn
    replicates are all infeasible (in a scenario, with criteria `across`
    scenarios) is not judged in the resample: it neither meets the
    requirements nor is efficient.

    Args:
        results (pandas.DataFrame): The results of a study (see
            `Study.run`), with a record per replicate of each design (and
            group of the criteria's `by`).
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (metrics whose means they are) and their directions
            (`"min"` or `"max"`).
        n_resamples (int): The number of resamples (at least 2).
        seed (int): The seed of the resamples.

    Returns:
        pandas.DataFrame: a row per feasible design (and group), empty
            without one: the summary of the criteria over every replicate
            (see `summarize`, with the spreads of means), in which each
            objective that is not a mean (such as a quantile) is followed by
            its standard error, its standard deviation over the resamples
            in which the design is judged (`<objective> se`, `NaN` with
            fewer than two, and left out if it is zero for every design, as
            `summarize` leaves out spreads), and the probability that the
            design is efficient (`pareto_probability`) and that it meets the
            requirements (`requirements_probability`)

    Raises:
        ValueError: if a statistic other than the mean of a column has the
            column's name, or the results have several records of an
            evaluation (in a group).
    """
    criteria = _as_criteria(criteria)
    statistics = _with_feasibility(_resolved(criteria))
    spreads = _spreads(statistics)
    # the standard errors of the objectives that are not means, by the bootstrap
    booted = [
        name
        for name in criteria.objectives
        if name not in spreads and f"{name} se" not in statistics
    ]
    added = {**spreads, **{name: [f"{name} se"] for name in booted}}
    columns = [column for spread in added.values() for column in spread]
    _check_names(statistics, results, columns)
    order = _ordered(statistics, added)
    tolerances = _tolerances(criteria)
    rng = np.random.default_rng(seed)
    frames = []
    for frame, table in _summary_groups(results, statistics, criteria.by, criteria):
        values, violation = table.evaluate()
        counts = np.zeros(len(table.designs))
        met = np.zeros(len(table.designs))
        # the sums of the objectives (and of their squares) over the resamples
        # in which they are known (in which each design is judged)
        sums = {name: np.zeros(len(table.designs)) for name in booted}
        squares = {name: np.zeros(len(table.designs)) for name in booted}
        known = {name: np.zeros(len(table.designs)) for name in booted}
        for drawn in draw_resamples(len(table.replicates), n_resamples, rng):
            sample, sample_violation = table.evaluate(drawn)
            # designs without a drawn replicate keep the statistics of all of
            # theirs, and those without a feasible one are not judged
            undrawn = ~table.recorded(drawn)
            for name in sample:
                sample[name][undrawn] = values[name][undrawn]
            sample_violation[undrawn] = violation[undrawn]
            judged = table.judged(drawn) | undrawn
            meets = (sample_violation == 0) & judged
            counts += efficient(_costs(criteria, sample), meets, tolerances)
            met += meets
            for name in booted:
                value = sample[name]
                found = ~np.isnan(value)
                sums[name] += np.where(found, value, 0.0)
                squares[name] += np.where(found, value, 0.0) ** 2
                known[name] += found
        computed = dict(values)
        _with_spreads(computed, table, statistics, spreads)
        for name in booted:
            computed[f"{name} se"] = _bootstrap_error(
                sums[name], squares[name], known[name]
            )
        ordered = {name: computed[name] for name in order}
        _assemble(frame, table, ordered, violation, criteria)
        frame["pareto_probability"] = counts / n_resamples
        frame["requirements_probability"] = met / n_resamples
        frames.append(frame)
    summary = _summaries(
        frames,
        criteria.by,
        order,
        criteria,
        after=["pareto_probability", "requirements_probability"],
    )
    return _without_constant(summary, columns)
