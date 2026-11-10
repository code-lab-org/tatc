"""
Screening of designs with analyses of increasing fidelity: cheap analyses
(such as shorter periods, fewer targets, or fewer replicates) rank many
designs, and only the best are analyzed with more expensive ones.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .criteria import Requirement, Statistic, _get_statistics, aggregate
from .pareto import pareto_ranks
from .space import Design
from .study import Study


def _rank(
    summary: pd.DataFrame,
    objectives: Mapping[str, str],
    tolerance: Mapping[str, float] | None,
) -> pd.Series:
    """
    Ranks designs: those that meet the requirements by non-dominated sorting
    of their objectives, then the others by their violations.
    """
    ranks = pd.Series(-1, index=summary.index, name="pareto_rank")
    meets = summary["meets_requirements"]
    if meets.any():
        ranks[meets] = pareto_ranks(summary[meets], objectives, tolerance)
    if (~meets).any():
        start = int(ranks[meets].max()) + 1 if meets.any() else 0
        violation = summary.loc[~meets, "violation"].fillna(np.inf)
        ranks[~meets] = start + violation.rank(method="dense").astype(int) - 1
    return ranks


@dataclass
class Screening:
    """
    The outcome of a screening (see `screen`).
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of every stage (see `Study.run`), with
    the `stage` (0 for the first)."""
    summaries: list[pd.DataFrame] = field(default_factory=list)
    """list[pandas.DataFrame]: For each stage, the summary of its designs
    (see `summarize`), with their Pareto rank (`pareto_rank`, -1 if
    infeasible) and whether they were `kept` for the next stage."""
    designs: list[Design] = field(default_factory=list)
    """list[Design]: The designs analyzed in the last stage."""

    @property
    def final(self) -> pd.DataFrame:
        """pandas.DataFrame: The summary of the last stage."""
        return self.summaries[-1]


def _per_stage(value, count: int, name: str) -> list:
    """
    Gets a value for each of `count` transitions between stages.
    """
    if isinstance(value, Sequence) and not isinstance(value, str):
        if len(value) != count:
            raise ValueError(f"{name} needs a value for each of {count} transitions")
        return list(value)
    return [value] * count


def screen(
    designs: Sequence[Design],
    studies: Sequence[Study],
    objectives: Mapping[str, str],
    keep: float | int | Sequence[float | int] = 0.25,
    tolerance: Mapping[str, float] | Sequence[Mapping[str, float] | None] | None = None,
    statistics: Mapping[str, Statistic | str] | None = None,
    requirements: Sequence[Requirement] = (),
) -> Screening:
    """
    Screens designs with studies of increasing fidelity (and cost): each
    study (stage) evaluates the designs kept by the one before, ranks them by
    non-dominated sorting of their objectives (see `pareto_ranks`): the means
    of metrics, or other statistics over their replicates (see `aggregate`),
    and keeps whole fronts, best first, until it keeps at least `keep`
    designs (a number, or a share of the feasible designs).

    Designs that miss a requirement are ranked after those that meet them
    all, in order of their violation (see `aggregate`), rather than dropped:
    a cheap analysis may miss a requirement that the full analysis meets, so
    designs that nearly meet the requirements are kept if `keep` reaches
    them.

    A cheap analysis ranks designs less exactly than an expensive one, so a
    design it ranks just behind another may be better: keep enough designs,
    or give a `tolerance` of each objective (in the units of the stage's
    analysis), within which differences do not separate fronts (see
    `pareto_front`), so that nearly efficient designs are kept. Compare the
    final front with one found without screening on a sample of designs to
    check that screening does not lose efficient designs.

    Args:
        designs (Sequence[Design]): The designs to screen.
        studies (Sequence[Study]): The studies of each stage, in order of
            increasing fidelity.
        objectives (Mapping[str, str]): The objectives: statistics (or
            metrics) and their directions (`"min"` or `"max"`).
        keep (float | int | Sequence): The designs kept by each stage but the
            last: a number (int), or a share of its feasible designs (float),
            or a value for each stage.
        tolerance (Mapping | Sequence | None): The tolerances of objectives
            in each stage but the last (see `pareto_front`), or the same for
            each stage.
        statistics (Mapping[str, Statistic | str] | None): Statistics, by
            name, of objectives and requirements that are not means.
        requirements (Sequence[Requirement]): Requirements on statistics.

    Returns:
        Screening: the results and summaries of each stage
    """
    if len(studies) == 0:
        raise ValueError("at least one study is required")
    transitions = len(studies) - 1
    keeps = _per_stage(keep, transitions, "keep")
    if tolerance is None or isinstance(tolerance, Mapping):
        tolerances = [tolerance] * transitions
    else:
        tolerances = _per_stage(tolerance, transitions, "tolerance")
    names = list(objectives) + [requirement.name for requirement in requirements]
    resolved = _get_statistics(names, statistics)
    current = list({design.key: design for design in designs}.values())
    screening = Screening(results=pd.DataFrame())
    frames = []
    for stage, study in enumerate(studies):
        results = study.run(current).assign(stage=stage)
        frames.append(results)
        summary = aggregate(results, resolved, requirements, design_columns=[])
        summary = summary.dropna(subset=list(objectives)).reset_index(drop=True)
        summary["pareto_rank"] = _rank(
            summary, objectives, tolerances[stage] if stage < transitions else None
        )
        if stage < transitions:
            amount = keeps[stage]
            count = (
                math.ceil(amount * len(summary))
                if isinstance(amount, float)
                else min(int(amount), len(summary))
            )
            ranks = summary["pareto_rank"].sort_values(kind="stable")
            last = ranks.iloc[count - 1] if count > 0 else -1
            summary["kept"] = summary["pareto_rank"] <= last
            kept = set(summary.loc[summary["kept"], "design"])
            current = [design for design in current if design.key in kept]
        else:
            summary["kept"] = True
        screening.summaries.append(summary.assign(stage=stage))
    screening.results = pd.concat(frames, ignore_index=True)
    screening.designs = current
    return screening
