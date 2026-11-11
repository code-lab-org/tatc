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

from .criteria import Requirement, Statistic, _get_statistics, aggregate
from .pareto import Tolerance, pareto_front
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
    the new replicates, those of the replicates that selected the design
    (`<statistic> selected`), the number of new `replicates`, and whether
    the design meets the requirements (`meets_requirements`) and is
    efficient among the confirmed designs that do (`pareto`), by the new
    replicates."""

    def optimism(self, objectives: Mapping[str, str]) -> pd.Series:
        """
        Gets the optimism of selection: the mean amount by which the selecting
        replicates made the designs look better than the new replicates do,
        in each objective (positive if they looked better).

        Args:
            objectives (Mapping[str, str]): The objectives and their
                directions (`"min"` or `"max"`).

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
                for name, direction in objectives.items()
            }
        )


def confirm(
    study: Study,
    designs: Sequence[Design],
    results: pd.DataFrame,
    objectives: Mapping[str, str],
    replicates: int | None = None,
    statistics: Mapping[str, Statistic | str] | None = None,
    requirements: Sequence[Requirement] = (),
    by: Sequence[str] = (),
    tolerance: Mapping[str, Tolerance] | None = None,
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
        objectives (Mapping[str, str]): The objectives: statistics (or
            metrics) and their directions (`"min"` or `"max"`).
        replicates (int | None): The number of new replicates, or `None` for
            the study's number.
        statistics (Mapping[str, Statistic | str] | None): Statistics, by
            name, of objectives and requirements that are not means.
        requirements (Sequence[Requirement]): Requirements on statistics.
        by (Sequence[str]): Columns that group records of an evaluation (such
            as a sensor), whose fronts are found separately.
        tolerance (Mapping[str, float | str] | None): The tolerances of
            objectives (see `pareto_front`).

    Returns:
        Confirmation: the new results and the summary of each design
    """
    count = study.replicates if replicates is None else replicates
    if count < 1:
        raise ValueError("replicates must be at least 1")
    designs = list({design.key: design for design in designs}.values())
    keys = {design.key for design in designs}
    start = int(results["replicate"].max()) + 1 if len(results) > 0 else 0
    fresh = study.run(designs, replicates=range(start, start + count))
    names = list(dict.fromkeys([*objectives, *(r.name for r in requirements)]))
    resolved = _get_statistics(names, statistics)
    resolved = {name: resolved[name] for name in names}
    columns = ["design", *by]
    selected = aggregate(
        results[results["design"].isin(keys)],
        resolved,
        requirements,
        by=by,
        design_columns=[],
    )
    confirmed = aggregate(fresh, resolved, requirements, by=by, design_columns=[])
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
    summary["pareto"] = False
    groups = summary.groupby(list(by), sort=False) if by else [((), summary)]
    for _, group in groups:
        candidates = group[group["meets_requirements"]].dropna(subset=list(objectives))
        if len(candidates) > 0:
            efficient = pareto_front(candidates, objectives, tolerance)
            summary.loc[efficient.index, "pareto"] = efficient
    return Confirmation(results=fresh, summary=summary)
