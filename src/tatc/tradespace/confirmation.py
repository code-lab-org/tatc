"""
Confirmation of selected designs with replicates that did not select them,
which removes the optimism of selection (the winner's curse).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .criteria import Criteria, _as_criteria, aggregate
from .pareto import pareto_front
from .space import Design
from .study import Study


@dataclass
class Confirmation:
    """
    The outcome of a confirmation (see `confirm`).
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of the new replicates (see `Study.run`)."""
    summary: pd.DataFrame
    """pandas.DataFrame: A row per design (and group): the statistics of
    the criteria by the new replicates, those of the replicates that
    selected the design (`<statistic> selected`), the number of new
    `replicates`, and whether the design meets the requirements
    (`meets_requirements`) and is efficient among the confirmed designs that
    do (`pareto`), by the new replicates."""
    criteria: Criteria
    """Criteria: The criteria of the confirmation."""

    def optimism(self) -> pd.Series:
        """
        Gets the optimism of selection: the mean amount by which the selecting
        replicates made the designs look better than the new replicates do,
        in each objective of the criteria (positive if they looked better).

        Returns:
            pandas.Series: the optimism of each objective, in its units
        """
        return pd.Series(
            {
                name: float(
                    np.nanmean(
                        (self.summary[f"{name} selected"] - self.summary[name])
                        * (1.0 if direction == "max" else -1.0)
                    )
                )
                for name, direction in self.criteria.objectives.items()
            }
        )


def confirm(
    study: Study,
    designs: Sequence[Design],
    results: pd.DataFrame,
    criteria: Criteria | Mapping[str, str],
    replicates: int | None = None,
) -> Confirmation:
    """
    Confirms selected designs (such as a Pareto front, or the best designs
    of an optimizer) with new replicates.

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
        designs (Sequence[Design]): The selected designs.
        results (pandas.DataFrame): Every result that informed the selection
            (see `Study.run`), such as those of an enumeration, refinement,
            or search.
        criteria (Criteria | Mapping[str, str]): The criteria (see
            `Criteria`), or the objectives (metrics whose means they are)
            and their directions (`"min"` or `"max"`).
        replicates (int | None): The number of new replicates, or `None` for
            the study's number.

    Returns:
        Confirmation: the new results and the summary of each design
    """
    criteria = _as_criteria(criteria)
    count = study.replicates if replicates is None else replicates
    if count < 1:
        raise ValueError("replicates must be at least 1")
    designs = list({design.key: design for design in designs}.values())
    keys = {design.key for design in designs}
    start = int(results["replicate"].max()) + 1 if len(results) > 0 else 0
    fresh = study.run(designs, replicates=range(start, start + count))
    names = criteria.names
    columns = ["design", *criteria.by]
    selected = aggregate(
        results[results["design"].isin(keys)], criteria, design_columns=[]
    )
    confirmed = aggregate(fresh, criteria, design_columns=[])
    summary = confirmed[[*columns, *names, "replicates", "meets_requirements"]]
    summary = summary.merge(
        selected[[*columns, *names]].rename(
            columns={name: f"{name} selected" for name in names}
        ),
        on=columns,
        how="left",
    )
    summary = summary[
        [
            *columns,
            *(c for name in names for c in (name, f"{name} selected")),
            "replicates",
            "meets_requirements",
        ]
    ]
    # the front of the designs whose objectives are known, by the new replicates
    known = summary.dropna(subset=list(criteria.objectives))
    summary["pareto"] = False
    summary.loc[known.index, "pareto"] = pareto_front(known, criteria).to_numpy()
    return Confirmation(results=fresh, summary=summary, criteria=criteria)
