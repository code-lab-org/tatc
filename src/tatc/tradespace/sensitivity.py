"""
Sensitivity of the metrics of designs to their design variables: which
decisions drive which outcomes, from the designs of a study (an enumeration,
or a sample such as a Latin hypercube).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def _levels(values: pd.Series, max_levels: int) -> pd.Series:
    """
    Gets the levels of a variable: its values, or (for a numeric variable
    with more than `max_levels` values) equal-count bins of them.
    """
    if pd.api.types.is_numeric_dtype(values) and values.nunique() > max_levels:
        return pd.qcut(values, max_levels, duplicates="drop").astype(str)
    return values


def main_effects(
    frame: pd.DataFrame, variable: str, metric: str, max_levels: int = 10
) -> pd.DataFrame:
    """
    Computes the main effects of a design variable on a metric: the mean (and
    standard deviation and count) of the metric over the designs at each
    level of the variable, among the designs in which the variable is active
    (not `NaN`, as a variable of a branch not chosen is).

    Args:
        frame (pandas.DataFrame): The designs, with their values (see
            `Design.flat`) and metrics (such as a summary of a study, see
            `summarize` or `aggregate`).
        variable (str): The column of the variable.
        metric (str): The column of the metric.
        max_levels (int): The most levels of a numeric variable, beyond which
            its values are grouped in bins of equal counts.

    Returns:
        pandas.DataFrame: a row per level: the `mean`, `std`, and `count` of
            the metric
    """
    active = frame.dropna(subset=[variable, metric])
    levels = _levels(active[variable], max_levels)
    return active.groupby(levels, sort=True)[metric].agg(["mean", "std", "count"])


def sensitivity(
    frame: pd.DataFrame,
    variables: Sequence[str],
    metrics: Sequence[str],
    max_levels: int = 10,
) -> pd.DataFrame:
    """
    Computes the sensitivity of metrics to design variables: for each
    variable and metric, the share of the variance of the metric explained
    by the levels of the variable (the correlation ratio, eta squared, of a
    one-way analysis of variance), among the designs in which the variable is
    active, and the share of the designs in which it is (`active`).

    A share near 1 means that the variable alone determines the metric; a
    share near 0, that the metric hardly depends on it. The shares of a full
    factorial enumeration are first-order effects (the shares of different
    variables do not add to 1 if variables interact); those of a sample
    estimate them. The shares of variables that are active in different
    designs (such as those of different branches of a choice) are of
    different designs, and compare how much each matters where it applies.

    Args:
        frame (pandas.DataFrame): The designs, with their values (see
            `Design.flat`) and metrics (such as a summary of a study, see
            `summarize` or `aggregate`).
        variables (Sequence[str]): The columns of the variables.
        metrics (Sequence[str]): The columns of the metrics.
        max_levels (int): The most levels of a numeric variable, beyond which
            its values are grouped in bins of equal counts.

    Returns:
        pandas.DataFrame: a row per variable: the share of each metric's
            variance it explains, and the share of designs in which it is
            `active`, sorted by the share of the first metric
    """
    rows = []
    for variable in variables:
        row: dict = {"variable": variable}
        active = frame[variable].notna()
        row["active"] = float(active.mean())
        for metric in metrics:
            subset = frame[active & frame[metric].notna()]
            values = subset[metric].astype(float)
            total = float(((values - values.mean()) ** 2).sum())
            if len(subset) < 2 or total == 0:
                row[metric] = np.nan
                continue
            levels = _levels(subset[variable], max_levels)
            means = values.groupby(levels).transform("mean")
            row[metric] = float(((means - values.mean()) ** 2).sum()) / total
        rows.append(row)
    result = pd.DataFrame(rows).set_index("variable")
    result = result[[*metrics, "active"]]
    return result.sort_values(metrics[0], ascending=False) if metrics else result
