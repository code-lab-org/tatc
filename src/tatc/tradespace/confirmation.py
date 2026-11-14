"""
Confirmation of selected designs with replicates that did not select them,
which removes the optimism of selection (the winner's curse).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .criteria import Criteria, _as_criteria
from .space import Design, _designs
from .study import Study
from .summary import _FEASIBLE, _JUDGMENT, aggregate


@dataclass
class Confirmation:
    """
    The outcome of a confirmation (see `confirm`): as the outcomes of other
    steps of a trade study (see `Refinement`), its `results`, the `summary`
    of its designs, the `designs` it selects for a next step (those it
    confirmed), and its `criteria`. Its results are only those of the new
    replicates, which judge the designs without the optimism of their
    selection: the designs still efficient by them are those of
    `summary.pareto`.
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of the new replicates (see `Study.run`),
    which do not include those that selected the designs."""
    summary: pd.DataFrame
    """pandas.DataFrame: A row per design (and group): the values of the
    design, the statistics of the criteria by the new replicates, those of
    the replicates that selected the design (`<statistic> selected`), the
    optimism of the selection in each objective (`<objective> optimism`: by
    how much the selecting replicates made the design look better than the
    new replicates do, in its units, positive if they did), the number of
    new feasible `replicates`, the share of the new replicates that are
    `feasible`, the `violation` of the requirements, and whether the design
    meets them (`meets_requirements`) and is efficient among the confirmed
    designs that do (`pareto`), by the new replicates."""
    criteria: Criteria
    """Criteria: The criteria of the confirmation."""
    designs: list[Design] = field(default_factory=list)
    """list[Design]: The confirmed designs (every design evaluated with new
    replicates)."""

    @property
    def optimism(self) -> pd.Series:
        """pandas.Series: The optimism of selection: the mean over the designs
        of their optimism in each objective of the criteria (see
        `summary`), in its units, positive if the selecting replicates made
        them look better."""
        return pd.Series(
            {
                name: float(np.nanmean(self.summary[f"{name} optimism"]))
                for name in self.criteria.objectives
            },
            dtype=float,
        )


def confirm(
    study: Study,
    criteria: Criteria | Mapping[str, str],
    results: pd.DataFrame,
    designs: Sequence[Design | str] | pd.DataFrame | None = None,
    *,
    replicates: int | None = None,
) -> Confirmation:
    """
    Confirms selected designs (by default, the Pareto front of the results,
    or the designs that another step selects, such as the efficient designs
    of a refinement or a search, see `Refinement.designs`) with new
    replicates.

    Designs selected for their outcomes in some replicates are selected in
    part for the luck of those replicates: among many designs, the best
    looking are likely to have drawn favorable uncertain parameters, so
    their statistics in the replicates that selected them are optimistic
    (the winner's curse), even if each design's statistics are unbiased
    before selection. Replicates that played no part in the selection give
    unbiased statistics of the selected designs. They are the replicates
    after every replicate of the `results` (of any design), the same for
    every design, so that common random numbers still pair the designs (see
    `Design.realize`).

    Args:
        study (Study): The study that evaluates designs.
        criteria (Criteria | Mapping[str, str]): The criteria (see
            `Criteria`), or the objectives (metrics whose means they are)
            and their directions (`"min"` or `"max"`).
        results (pandas.DataFrame): Every result that informed the selection
            (see `Study.run`), such as those of an enumeration, refinement,
            or search.
        designs (Sequence[Design | str] | pandas.DataFrame | None): The
            selected designs (designs, their keys, or a data frame of them,
            such as rows of a summary, or the `designs` of another step),
            or `None` for those efficient by
            the criteria over every replicate of the results (see
            `aggregate` and `pareto_front`).
        replicates (int | None): The number of new replicates, or `None` for
            the study's number.

    Returns:
        Confirmation: the new results, the summary of each design, and the
            designs
    """
    criteria = _as_criteria(criteria)
    count = study.replicates if replicates is None else replicates
    if count < 1:
        raise ValueError("replicates must be at least 1")
    if designs is None:
        summary = aggregate(results, criteria, design_columns=[])
        designs = summary[summary["pareto"]]
    designs = _designs(study.space, designs, results)
    keys = {design.key for design in designs}
    start = int(results["replicate"].max()) + 1 if len(results) > 0 else 0
    fresh = study.run(designs, replicate_indices=range(start, start + count))
    names = [name for name in criteria.names if name != _FEASIBLE]
    columns = ["design", *criteria.by]
    confirmed = aggregate(fresh, criteria)
    selected = aggregate(
        results[results["design"].isin(keys)], criteria, design_columns=[]
    )[[*columns, *names]]
    summary = confirmed.merge(
        selected.rename(columns={name: f"{name} selected" for name in names}),
        on=columns,
        how="left",
    )
    for name, direction in criteria.objectives.items():
        sign = 1.0 if direction == "max" else -1.0
        summary[f"{name} optimism"] = sign * (
            summary[f"{name} selected"] - summary[name]
        )
    status = ["replicates", _FEASIBLE, *_JUDGMENT]
    values = [c for c in confirmed.columns if c not in {*columns, *names, *status}]
    summary = summary[
        [
            *columns,
            *values,
            *(c for name in names for c in (name, f"{name} selected")),
            *(f"{name} optimism" for name in criteria.objectives),
            *status,
        ]
    ]
    return Confirmation(
        results=fresh, summary=summary, criteria=criteria, designs=designs
    )
