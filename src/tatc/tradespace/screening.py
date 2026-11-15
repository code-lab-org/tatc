"""
Screening of designs with analyses of increasing fidelity: cheap analyses
(such as shorter periods, fewer targets, or fewer replicates) rank many
designs, and only the best are analyzed with more expensive ones.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
import numbers
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

from .criteria import Criteria
from .pareto import pareto_ranks
from .space import Design, _designs
from .study import Study
from .summary import aggregate


@dataclass
class Screening:
    """
    The outcome of a screening (see `screen`): as the outcomes of other steps
    of a trade study (see `Refinement`), its `results` and the `summary` of
    its designs (those of its last stage, the analysis of the highest
    fidelity), the `designs` it selects for the next step (those that
    survive every stage), and its `criteria` (those of the last stage),
    which refine and confirm its designs as for any study::

        screening = screen([quick, full], criteria, designs, keep=0.1)
        refinement = refine(full, criteria, screening.results, screening.designs)

    The results and summaries of every stage are kept as well (`history` and
    `summaries`).
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of the last stage (see `Study.run`)."""
    criteria: Criteria
    """Criteria: The criteria of the last stage."""
    history: pd.DataFrame = field(default_factory=pd.DataFrame)
    """pandas.DataFrame: The results of every stage, with the `stage` (0 for
    the first), such as for the cost of the screening."""
    summaries: list[pd.DataFrame] = field(default_factory=list)
    """list[pandas.DataFrame]: For each stage, the summary of its designs
    (see `aggregate`, with whether each is efficient by the stage's
    analysis, `pareto`), with their Pareto rank (`pareto_rank`, see
    `pareto_ranks`), whether they were `kept` for the next stage, and the
    `stage`."""
    designs: list[Design] = field(default_factory=list)
    """list[Design]: The designs that survive every stage (those the last
    stage analyzes), for the next step."""

    @property
    def summary(self) -> pd.DataFrame:
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


def _check_keep(value) -> None:
    """
    Checks a number of designs to keep (an integer of at least 1) or a share
    of them (a float greater than 0 and at most 1).
    """
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise TypeError(f"keep must be a number or a share of designs, not {value!r}")
    if isinstance(value, numbers.Integral):
        if value < 1:
            raise ValueError(f"a number of designs to keep must be at least 1: {value}")
    elif not 0 < value <= 1:
        raise ValueError(
            f"a share of designs to keep must be greater than 0 and at most 1: {value}"
        )


def screen(
    studies: Sequence[Study],
    criteria: Criteria | Mapping[str, str] | Sequence[Criteria | Mapping[str, str]],
    designs: Sequence[Design] | pd.DataFrame,
    *,
    keep: float | int | Sequence[float | int] = 0.25,
) -> Screening:
    """
    Screens designs with studies of increasing fidelity (and cost): each
    study (stage) evaluates the designs kept by the one before, ranks them by
    non-dominated sorting of the objectives of the criteria (see
    `pareto_ranks`): the means of metrics, or other statistics over their
    replicates (see `Criteria`), and keeps whole fronts, best first, until
    it keeps at least `keep` designs (a number, or a share of the feasible
    designs).

    Designs that miss a requirement are ranked after those that meet them
    all, in order of their violation, and designs with a missing objective
    last (see `pareto_ranks`), rather than dropped:
    a cheap analysis may miss a requirement that the full analysis meets, so
    designs that nearly meet the requirements are kept if `keep` reaches
    them.

    A cheap analysis ranks designs less exactly than an expensive one, so a
    design it ranks just behind another may be better: keep enough designs,
    or give the stages of cheap analyses criteria of their own, with a
    tolerance of each objective (in the units of the stage's analysis)
    within which differences do not separate fronts (see `pareto_front`),
    so that nearly efficient designs are kept::

        cheap = Criteria(objectives, tolerances={"sufficient": 5.0})
        screening = screen([quick, full], [cheap, objectives], designs, keep=0.1)

    Compare the final front with one found without screening on a sample of
    designs to check that screening does not lose efficient designs.

    The screening's `results` are those of its last stage, which the next
    steps continue (such as `refine` and `confirm`), and its `designs` are
    those that survive every stage (that the last stage analyzes); its
    `history` gives the results of every stage (with the `stage`), such as
    for the cost of the screening.

    Args:
        studies (Sequence[Study]): The studies of each stage, in order of
            increasing fidelity.
        criteria (Criteria | Mapping[str, str] | Sequence): The criteria,
            without `by` (see `Criteria`), or the objectives (metrics whose
            means they are) and their directions (`"min"` or `"max"`): the
            same for every stage, or those of each stage.
        designs (Sequence[Design] | pandas.DataFrame): The designs to
            screen, or a data frame of their values (see
            `DesignSpace.designs_of`).
        keep (float | int | Sequence): The designs kept by each stage but the
            last: a number (an int of at least 1), or a share of its feasible
            designs (a float greater than 0 and at most 1), or a value for
            each stage.

    Returns:
        Screening: the results and summary of the last stage, the designs
            that survive every stage, and the results and summaries of each
            stage
    """
    if len(studies) == 0:
        raise ValueError("at least one study is required")
    transitions = len(studies) - 1
    keeps = _per_stage(keep, transitions, "keep")
    for amount in keeps:
        _check_keep(amount)
    if isinstance(criteria, (Criteria, Mapping)):
        stages = [Criteria.of(criteria)] * len(studies)
    else:
        stages = [
            Criteria.of(c) for c in _per_stage(criteria, len(studies), "criteria")
        ]
    for criteria in stages:
        if criteria.by:
            raise ValueError(
                "screening ranks designs by one record per evaluation: "
                + "give criteria without `by`"
            )
    current = _designs(studies[0].space, designs)
    frames, summaries = [], []
    for stage, study in enumerate(studies):
        criteria = stages[stage]
        results = study.run(current)
        frames.append(results.assign(stage=stage))
        summary = aggregate(results, criteria)
        summary["pareto_rank"] = pareto_ranks(summary, criteria)
        if stage < transitions:
            amount = keeps[stage]
            count = (
                min(int(amount), len(summary))
                if isinstance(amount, numbers.Integral)
                else math.ceil(amount * len(summary))
            )
            ranks = summary["pareto_rank"].sort_values(kind="stable")
            last = ranks.iloc[count - 1] if count > 0 else -1
            summary["kept"] = summary["pareto_rank"] <= last
            kept = set(summary.loc[summary["kept"], "design"])
            current = [design for design in current if design.key in kept]
        else:
            summary["kept"] = True
        summaries.append(summary.assign(stage=stage))
    return Screening(
        results=results,
        criteria=criteria,
        history=pd.concat(frames, ignore_index=True),
        summaries=summaries,
        designs=current,
    )
