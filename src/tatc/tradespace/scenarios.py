"""
Scenarios: conditions under which designs are compared rather than averaged
over (see `Study`), such as the requirements of a mission or the field it
must sample, and the robustness of designs across them.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping

import numpy as np
import pandas as pd

from .criteria import Criteria, _as_criteria, _get_tolerances, aggregate
from .pareto import pareto_front


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
    absolute = _get_tolerances(criteria.objectives, criteria.tolerance)[0]
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
    scenario: str = "scenario",
    scale: Mapping[str, float] | None = None,
) -> pd.DataFrame:
    """
    Measures the robustness of designs across scenarios (see `Study`): how
    often each design is efficient, and how far it falls from the front of
    each scenario, from the results of a study in each scenario (judged by
    the criteria in each scenario, see `aggregate` and `pareto_front`).

    * `pareto_trace`: the share of the scenarios in which the design is
      efficient (with the tolerances of the criteria, within which designs
      near the front are efficient, as in epoch-era analysis's fuzzy
      normalized Pareto trace). A design efficient in every scenario has a
      trace of 1.
    * `<scenario> regret`: in each scenario, the improvement in every
      objective (in units of its `scale`) that the design would need to match
      a design of the scenario's front (the additive epsilon indicator, as
      the fuzzy Pareto number of epoch-era analysis): 0 if it is efficient,
      and infinite if it misses the requirements of the criteria (or is
      infeasible in every replicate). With one objective, it is the
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
            designs are judged in each scenario (`by` none, or the scenario),
            or the objectives (metrics whose means they are) and their
            directions (`"min"` or `"max"`).
        scenario (str): The column of the scenarios.
        scale (Mapping[str, float] | None): The units of regret of
            objectives (by name): by default, the absolute tolerance of an
            objective (if it has one), or else 1 (its units).

    Returns:
        pandas.DataFrame: a row per design (`design`): its `pareto_trace`,
            `max_regret`, and regret in each scenario (`<scenario> regret`)

    Raises:
        ValueError: if the results have no scenarios, or the criteria group
            records by other columns or combine them across scenarios (see
            `Across`, by which searches find robust designs).
    """
    criteria = _as_criteria(criteria)
    if criteria.across is not None:
        raise ValueError(
            "robustness judges designs in each scenario: give criteria without "
            + "`across`"
        )
    if criteria.by not in ((), (scenario,)):
        raise ValueError(
            f"criteria judge designs by the scenario ({scenario!r}) alone: "
            + f"give criteria without other columns `by`, not {list(criteria.by)}"
        )
    if scenario not in results.columns:
        raise ValueError(
            f"the results have no column {scenario!r} of scenarios (see Study)"
        )
    judged = dataclasses.replace(criteria, by=(scenario,))
    objectives = list(criteria.objectives)
    # each objective in units of its scale, larger being better
    signs = np.array(
        [1.0 if d == "max" else -1.0 for d in criteria.objectives.values()]
    )
    units = signs / _scales(criteria, scale)
    summary = aggregate(results, judged)
    efficient = pareto_front(summary, judged).to_numpy()
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
        meets = group["meets_requirements"].to_numpy(dtype=bool)
        regret.loc[keys, name] = np.where(front, 0.0, np.where(meets, needed, np.inf))
        trace[keys[front]] += 1
    frame = pd.DataFrame(
        {
            "design": designs,
            "pareto_trace": (trace / max(len(scenarios), 1)).to_numpy(),
            "max_regret": regret.max(axis=1).to_numpy(),
        }
    )
    for name in scenarios:
        frame[f"{name} regret"] = regret[name].to_numpy()
    return frame
