"""
Tables of trade studies: the results of a study (a row per evaluation, a
design in a replicate) and the summaries of its designs (a row per design):
their columns, their groups, the values of the metrics of designs in each
replicate, the values of designs, and checks that a table has the rows that
a function needs.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .space import STATUS_COLUMNS

FEASIBLE = "feasible"
"""The column of the results of a study of whether each record is feasible,
and the statistic (and column of summaries) of the share of a design's
replicates that are feasible (see `aggregate`)."""

JUDGMENT = ("violation", "meets_requirements", "pareto")
"""The columns of the judgment of designs by criteria in their summaries (see
`aggregate`)."""

SEPARATE = (
    "give the columns that separate them as `by` (of the criteria), or combine "
    + "their groups `across` (see Across)"
)
"""How to separate several records of an evaluation for criteria."""

SUMMARIZE = (
    "summarize the results of a study (a row per replicate) first, with a row "
    + "per design (see aggregate), or select the rows of one group (such as those "
    + "of a scenario)"
)
"""How to give each design one row, rather than the results of a study (a
row per replicate)."""


def names(value: str | Iterable[str]) -> list[str]:
    """
    Gets names given as a name or a sequence of names.

    Args:
        value (str | Iterable[str]): The name, or names.

    Returns:
        list[str]: the names
    """
    return [value] if isinstance(value, str) else list(value)


def groups(frame: pd.DataFrame, by: Sequence[str]) -> Iterator[tuple[tuple, Any]]:
    """
    Iterates over the groups of a data frame by columns (or the frame, if
    none): the values of the columns, and the positions of the group's rows.

    Args:
        frame (pandas.DataFrame): The data frame.
        by (Sequence[str]): The columns.

    Returns:
        Iterator[tuple[tuple, numpy.ndarray]]: the values and positions of
            each group

    Raises:
        ValueError: if the frame lacks a column.
    """
    if not by:
        yield (), np.arange(len(frame))
        return
    missing = [column for column in by if column not in frame.columns]
    if missing:
        raise ValueError(f"no columns {missing} to group by")
    for key, positions in frame.groupby(list(by), sort=False).indices.items():
        yield (key if isinstance(key, tuple) else (key,)), positions


def check_summary(
    summary: pd.DataFrame, by: Sequence[str] = (), hint: str = SUMMARIZE
) -> None:
    """
    Checks that a data frame is a summary of designs (a row per design, in
    each group of the columns `by`), not the results of a study (a row per
    replicate), which the functions of summaries would judge as if each
    replicate were a design.

    Args:
        summary (pandas.DataFrame): The data frame.
        by (Sequence[str]): The columns of groups, in each of which a design
            may have a row.
        hint (str): What to do if the data frame is not a summary.

    Raises:
        ValueError: if the data frame has several rows of a design (in a
            group), with the hint.
    """
    if "design" in summary.columns:
        keys = summary[["design", *(c for c in by if c in summary.columns)]]
    elif summary.index.name == "design" and not by:
        keys = summary.index.to_frame()
    else:
        return
    if keys.duplicated().any():
        raise ValueError(f"the data frame has several rows of a design: {hint}")


def check_records(
    results: pd.DataFrame,
    keys: Sequence[str] = ("design", "replicate"),
    hint: str = SEPARATE,
    metrics: Sequence[str] = (),
) -> None:
    """
    Checks that results have one record per evaluation (a design in a
    replicate, or a combination of the values of other `keys`).

    Args:
        results (pandas.DataFrame): The results.
        keys (Sequence[str]): The columns of an evaluation.
        hint (str): How to separate several records of an evaluation.
        metrics (Sequence[str]): The metrics, whose values may differ
            between the records of an evaluation without separating them.

    Raises:
        ValueError: if the results have several records of an evaluation,
            with the hint and the columns whose values differ between them
            (which may separate them).
    """
    keys = list(keys)
    if not results.duplicated(keys).any():
        return
    candidates = [
        column
        for column in results.columns
        if column not in (*STATUS_COLUMNS, *keys, *metrics) or column == "scenario"
    ]
    varying: list[str] = []
    if candidates:
        counts = results.groupby(keys)[candidates].nunique(dropna=False).max()
        varying = list(counts.index[counts > 1])
    raise ValueError(
        "the results have several records of an evaluation (a design in a "
        + f"replicate): {hint}"
        + (f", such as {varying}" if varying else "")
    )


def matrices(
    group: pd.DataFrame,
    columns: Sequence[str],
    designs: pd.Index | None = None,
    replicates: pd.Index | None = None,
    hint: str = SEPARATE,
) -> tuple[pd.Index, pd.Index, dict[str, np.ndarray]]:
    """
    Gets the values of metrics of each design (rows) in each replicate
    (columns) of the results of a group, which must have one record per
    evaluation (see `check_records`): of the designs and replicates of the
    group, or of those given (`NaN` without a record).

    Args:
        group (pandas.DataFrame): The results.
        columns (Sequence[str]): The metrics.
        designs (pandas.Index | None): The designs (keys), or `None` for
            those of the results, in order.
        replicates (pandas.Index | None): The replicates, or `None` for those
            of the results, sorted.
        hint (str): How to separate several records of an evaluation.

    Returns:
        tuple[pandas.Index, pandas.Index, dict[str, numpy.ndarray]]: the
            designs, the replicates, and the values of each metric (shape
            (designs, replicates))

    Raises:
        ValueError: if the results lack a metric, or have several records
            of an evaluation.
    """
    missing = [column for column in columns if column not in group.columns]
    if missing:
        raise ValueError(f"metrics that are not in the results: {missing}")
    check_records(group, hint=hint, metrics=columns)
    if designs is None:
        designs = pd.Index(group["design"].drop_duplicates())
    if replicates is None:
        replicates = pd.Index(sorted(group["replicate"].unique()))
    values = {}
    for column in columns:
        table = group.pivot(index="design", columns="replicate", values=column)
        values[column] = table.reindex(index=designs, columns=replicates).to_numpy(
            dtype=float
        )
    return designs, replicates, values


def feasible_records(results: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Gets the feasible and the infeasible records of results (every record is
    feasible without a column `feasible`).

    Args:
        results (pandas.DataFrame): The results.

    Returns:
        tuple[pandas.DataFrame, pandas.DataFrame]: the feasible records and
            the infeasible ones
    """
    if FEASIBLE not in results.columns:
        return results, results.iloc[:0]
    status = results[FEASIBLE].astype(bool).to_numpy()
    return results[status], results[~status]


def _constant_columns(results: pd.DataFrame, excluded: set[str]) -> list[str]:
    """
    Gets the columns (other than `excluded`) that are the same in every
    record of each design.
    """
    candidates = [column for column in results.columns if column not in excluded]
    if not candidates or len(results) == 0:
        return []
    constant = results.groupby("design")[candidates].nunique(dropna=False)
    return [column for column in candidates if (constant[column] <= 1).all()]


def value_columns(
    results: pd.DataFrame,
    statistics: Mapping[str, Any],
    by: Sequence[str],
    design_columns: Sequence[str] | None,
) -> list[str]:
    """
    Gets the columns of the values of designs that a summary of results
    keeps: those given, or those that are the same in every feasible record
    of each design (but the metrics of the statistics, the columns of groups,
    and status columns), without the names of statistics and the columns of
    groups.

    Args:
        results (pandas.DataFrame): The results.
        statistics (Mapping[str, Statistic]): The statistics of the summary,
            by name.
        by (Sequence[str]): The columns of groups (such as those of the
            criteria's `by`, and the column they combine `across`).
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every feasible
            record of a design.

    Returns:
        list[str]: the columns
    """
    if design_columns is None:
        metrics = {
            column for statistic in statistics.values() for column in statistic.columns
        }
        design_columns = _constant_columns(
            feasible_records(results)[0], {*STATUS_COLUMNS, *metrics, *by}
        )
    return [c for c in design_columns if c not in statistics and c not in by]


def with_values(
    frame: pd.DataFrame, results: pd.DataFrame, designs: Any, columns: Sequence[str]
) -> None:
    """
    Adds the values of designs (columns of the first record of each design of
    results) to the rows of the designs (in order) of a data frame.

    Args:
        frame (pandas.DataFrame): The data frame, a row per design.
        results (pandas.DataFrame): The results.
        designs (Any): The keys of the designs of the rows.
        columns (Sequence[str]): The columns of the values of designs.
    """
    first = results.drop_duplicates("design").set_index("design")
    for column in columns:
        frame[column] = first.loc[designs, column].to_numpy()


def design_values(
    results: pd.DataFrame, designs: Sequence[str], columns: Sequence[str]
) -> pd.DataFrame:
    """
    Gets the values of designs from their records in results (see
    `with_values`).

    Args:
        results (pandas.DataFrame): The results.
        designs (Sequence[str]): The keys of the designs.
        columns (Sequence[str]): The columns of the values of designs.

    Returns:
        pandas.DataFrame: a row per design: its key (`design`), and its values
    """
    frame = pd.DataFrame({"design": list(designs)})
    with_values(frame, results, list(designs), columns)
    return frame


def draw_resamples(count: int, n_resamples: int, seed: Any) -> np.ndarray:
    """
    Gets the replicates of each resample of a bootstrap: a row per resample of
    `count` replicates (their positions) drawn with replacement.

    Args:
        count (int): The number of replicates.
        n_resamples (int): The number of resamples (at least 2).
        seed (Any): The seed, or a random number generator, which draws them.

    Returns:
        numpy.ndarray: the positions of the replicates (shape (n_resamples,
            count))

    Raises:
        ValueError: if `n_resamples` is less than 2.
    """
    if n_resamples < 2:
        raise ValueError("n_resamples must be at least 2")
    return np.random.default_rng(seed).integers(count, size=(n_resamples, count))


def efficient_keys(summary: pd.DataFrame) -> list:
    """
    Gets the keys of the efficient designs of a summary of criteria (whose
    `pareto` is true, in any group), each once, in order.

    Args:
        summary (pandas.DataFrame): The summary (see `aggregate`).

    Returns:
        list: the keys
    """
    rows = summary["pareto"].to_numpy(dtype=bool)
    return list(pd.unique(summary.loc[rows, "design"]))
