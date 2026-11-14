"""
Scenarios: conditions under which designs are compared rather than averaged
over (see `Study`), such as the requirements of a mission or the field it
must sample, and the robustness of designs across them.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .criteria import Criteria, _as_criteria, _get_tolerances, _resolve, _signs
from .summary import _feasible_records, _value_columns, _with_values, aggregate


def _scales(criteria: Criteria, scale: Mapping[str, float] | None) -> np.ndarray:
    """
    Gets the scale of each objective of criteria: the given one, or its
    absolute tolerance (if any), or 1.
    """
    objectives = list(criteria.objectives)
    scale = dict(scale or {})
    unknown = sorted(set(scale) - set(objectives))
    if unknown:
        raise ValueError(f"scales of names that are not objectives: {unknown}")
    absolute = _get_tolerances(criteria.objectives, criteria.tolerances)[0]
    values = np.array(
        [
            float(scale[name]) if name in scale else (absolute[i] or 1.0)
            for i, name in enumerate(objectives)
        ]
    )
    if not np.all(values > 0):
        raise ValueError("the scales of the objectives must be positive")
    return values


def robustness(
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str],
    *,
    scale: Mapping[str, float] | None = None,
    design_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """
    Measures the robustness of designs across scenarios (see `Study`): how
    often each design is efficient, and how far it falls from the front of
    each scenario, from the results of a study in each scenario (judged by
    the criteria in each scenario, see `aggregate` and `pareto_front`). The
    scenarios are the groups of the column of the criteria's `by` (by
    default, `"scenario"`, as of a study with scenarios).

    * `pareto_trace`: the share of the scenarios in which the design is
      efficient (with the tolerances of the criteria, within which designs
      near the front are efficient, as in epoch-era analysis's fuzzy
      normalized Pareto trace). A design efficient in every scenario has a
      trace of 1.
    * `<scenario> regret`: in each scenario, the improvement in every
      objective (in units of its `scale`) that the design would need to match
      a design of the scenario's front (the additive epsilon indicator, as
      the fuzzy Pareto number of epoch-era analysis): 0 if it is efficient,
      and infinite if it misses the requirements of the criteria, is
      infeasible in every replicate, or has a missing objective (`NaN`, which
      is not efficient, see `pareto_front`). With one objective, it is the
      difference from the best design.
    * `max_regret`: the largest regret of the design over the scenarios, which
      the most robust designs minimize (minimax regret).

    With objectives of few levels (such as numbers of satellites), regret is
    at most a design's distance in them from the cheapest design of the
    front, whatever its other objectives, as dropping that much matches it:
    for the regret of the other objectives at a cost, require the cost (see
    `Requirement`) rather than minimizing it.

    Designs are compared with the others of the results, by the statistics of
    their replicates, whose uncertainty these measures do not show (see
    `pareto_probability`, with criteria `by` the scenario).

    Args:
        results (pandas.DataFrame): The results of a study with scenarios
            (see `Study.run`).
        criteria (Criteria | Mapping[str, str]): The criteria, by which
            designs are judged in each scenario (`by` none, for the column
            `"scenario"`, or the column of the scenarios), or the objectives
            (metrics whose means they are) and their directions (`"min"` or
            `"max"`).
        scale (Mapping[str, float] | None): The units of regret of
            objectives (by name): by default, the absolute tolerance of an
            objective (if it has one), or else 1 (its units).
        design_columns (Sequence[str] | None): The columns of the values of
            designs, or `None` for those that are the same in every feasible
            record of a design (see `aggregate`).

    Returns:
        pandas.DataFrame: a row per design (`design`): the values of the
            design (as a summary keeps them, see `aggregate`), each objective
            in each scenario (`<scenario> <objective>`), or once
            (`<objective>`) if it is a value of the designs (a column of
            `design_columns`, such as a number of satellites), its
            `pareto_trace`, `max_regret`, and regret in each scenario
            (`<scenario> regret`)

    Raises:
        ValueError: if the results have no scenarios, or the criteria group
            records by several columns or combine them across scenarios (see
            `Across`, by which searches find robust designs).
    """
    criteria = _as_criteria(criteria)
    if criteria.across is not None:
        raise ValueError(
            "robustness judges designs in each scenario: give criteria without "
            + "`across`"
        )
    if len(criteria.by) > 1:
        raise ValueError(
            "robustness judges designs in each scenario: give criteria `by` the "
            + f"column of the scenarios alone, not {list(criteria.by)}"
        )
    scenario = criteria.by[0] if criteria.by else "scenario"
    if scenario not in results.columns:
        raise ValueError(
            f"the results have no column {scenario!r} of scenarios (see Study, "
            + "or give the column as the criteria's `by`)"
        )
    judged = dataclasses.replace(criteria, by=(scenario,))
    objectives = list(criteria.objectives)
    # each objective in units of its scale, larger being better
    units = -_signs(criteria) / _scales(criteria, scale)
    summary = aggregate(results, judged, design_columns=[])
    efficient = summary["pareto"].to_numpy(dtype=bool)
    scenarios = list(pd.unique(results[scenario]))
    designs = list(pd.unique(results["design"]))
    regret = pd.DataFrame(np.inf, index=designs, columns=scenarios)
    trace = pd.Series(0.0, index=designs)
    for name in scenarios:
        rows = (summary[scenario] == name).to_numpy()
        group, front = summary[rows], efficient[rows]
        keys = group["design"].to_numpy()
        values = group[objectives].to_numpy(dtype=float) * units
        needed = np.full(len(group), np.inf)
        if front.any():
            # the improvement in every objective that each design needs to
            # match a design of the front: the least over the front of the
            # largest shortfall over the objectives
            shortfall = values[front][None, :, :] - values[:, None, :]
            needed = np.maximum(shortfall.max(axis=2).min(axis=1), 0.0)
        # a design with a missing objective is not efficient: infinite regret
        needed = np.where(np.isnan(values).any(axis=1), np.inf, needed)
        meets = group["meets_requirements"].to_numpy(dtype=bool)
        regret.loc[keys, name] = np.where(front, 0.0, np.where(meets, needed, np.inf))
        trace[keys[front]] += 1
    # the values of the designs (of their first records, as of a summary), and
    # the objectives that are values of the designs, which are the same in
    # every scenario
    statistics = _resolve(judged)
    feasible, _ = _feasible_records(results)
    every = _value_columns(feasible, {}, (scenario,), None, design_columns)
    metrics = {column for s in statistics.values() for column in s.columns}
    columns = [c for c in every if c not in metrics and c not in statistics]
    collapsed = {
        name
        for name in objectives
        if name in every and set(statistics[name].columns) <= set(every)
    }
    frame = pd.DataFrame({"design": designs})
    _with_values(frame, results, designs, columns)
    for name in objectives:
        table = summary.pivot(index="design", columns=scenario, values=name)
        table = table.reindex(index=designs, columns=scenarios)
        if name in collapsed:
            frame[name] = table.bfill(axis=1).iloc[:, 0].to_numpy()
        else:
            for group in scenarios:
                frame[f"{group} {name}"] = table[group].to_numpy()
    frame["pareto_trace"] = (trace / max(len(scenarios), 1)).to_numpy()
    frame["max_regret"] = regret.max(axis=1, skipna=False).to_numpy()
    for name in scenarios:
        frame[f"{name} regret"] = regret[name].to_numpy()
    return frame
