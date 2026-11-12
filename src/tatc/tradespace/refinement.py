"""
Adaptive replication: replicates added only to the designs whose Pareto
efficiency their replicates do not settle, within a budget of evaluations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .criteria import Criteria, _as_criteria, _get_tolerances
from .pareto import pareto_probability
from .space import Design, _designs
from .study import Study


@dataclass
class Refinement:
    """
    The outcome of an adaptive replication (see `refine`): as the outcomes of
    other steps of a trade study (see `Confirmation`, `Screening`, and
    `Search`), its `results`, the `summary` of its designs, and the `designs`
    themselves, for the next step::

        refinement = refine(study, criteria, results)
        confirmation = confirm(study, criteria, refinement.results)
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of every evaluation, those given and
    those added (see `Study.run`)."""
    summary: pd.DataFrame
    """pandas.DataFrame: A row per design (and group of the criteria's `by`):
    the probability that it is efficient after the last round (see
    `pareto_probability`, with the values of the design), and whether its
    efficiency is `uncertain` (by the bounds of `refine`, or with fewer than
    its least replicates) and, if so, whether it is `indifferent` (its
    objectives are known to within the precision of their tolerances)."""
    history: pd.DataFrame = field(default_factory=pd.DataFrame)
    """pandas.DataFrame: For each round, the number of designs whose
    efficiency is uncertain (`uncertain`), of these, those that are
    indifferent (`indifferent`), the number to which replicates were added
    (`refined`), and the evaluations added (`evaluations`)."""
    designs: list[Design] = field(default_factory=list)
    """list[Design]: The designs that could be refined."""

    @property
    def uncertain(self) -> pd.DataFrame:
        """pandas.DataFrame: The rows of `summary` whose efficiency is still
        uncertain, and that more replicates would settle (those that are not
        indifferent)."""
        summary = self.summary
        return summary[summary["uncertain"] & ~summary["indifferent"]]

    @property
    def indifferent(self) -> pd.DataFrame:
        """pandas.DataFrame: The rows of `summary` whose efficiency is
        uncertain, but whose objectives are known to within the precision of
        their tolerances, so that whether they are within the tolerances of
        the front does not matter."""
        return self.summary[self.summary["indifferent"]]


def _precise(
    probability: pd.DataFrame, criteria: Criteria, precision: float | None
) -> np.ndarray:
    """
    Gets whether each design's objectives are known to within the precision
    of their tolerances: the standard error of each objective that varies
    between replicates is less than `precision` times its tolerance (at the
    design's value). Designs whose varying objectives lack a tolerance are
    not precise.
    """
    precise = np.zeros(len(probability), dtype=bool)
    if precision is None or not criteria.tolerance:
        return precise
    tolerances = _get_tolerances(criteria.objectives, criteria.tolerance)
    precise[:] = True
    for i, name in enumerate(criteria.objectives):
        error = probability[f"{name} se"].to_numpy(dtype=float)
        value = np.abs(probability[name].to_numpy(dtype=float))
        tolerance = tolerances[0, i] + tolerances[1, i] * value
        varies = np.nan_to_num(error) > 0
        precise &= ~varies | (error < precision * tolerance)
    return precise


def refine(
    study: Study,
    criteria: Criteria | Mapping[str, str],
    results: pd.DataFrame,
    designs: Sequence[Design | str] | pd.DataFrame | None = None,
    *,
    bounds: tuple[float, float] = (0.05, 0.95),
    step: int | None = None,
    max_replicates: int | None = None,
    min_replicates: int = 5,
    max_rounds: int = 10,
    precision: float | None = 0.5,
    budget: int | None = None,
    n_resamples: int = 500,
    seed: int = 0,
) -> Refinement:
    """
    Adds replicates to the designs whose Pareto efficiency (by criteria, see
    `Criteria`) is uncertain, until it is settled: in each round, the
    probability that each design is efficient (and meets the requirements)
    is estimated by a bootstrap (see
    `pareto_probability`), and each design whose probability is strictly
    between the `bounds` gets `step` more replicates, up to
    `max_replicates`. Rounds end when no design is uncertain (or can be
    refined), when the `budget` of evaluations is spent, or after
    `max_rounds`. If a round would spend more than is left of the budget, the
    designs with the fewest replicates (below `min_replicates`) are refined
    first, and then those whose probability is nearest 0.5.

    A bootstrap of few replicates has few distinct resamples, and can find a
    design surely inefficient (or efficient) by the luck of its replicates:
    with two replicates, designs as good as the front can be efficient in
    none of the resamples. Designs with fewer than `min_replicates` are
    therefore uncertain whatever their probability, and are refined first.

    Replicates are added in order (the next indices after a design's
    replicates), so that designs refined together share their uncertain
    parameters (with common random numbers, see `Design.realize`). Designs
    that are surely efficient or surely not keep their replicates, which
    saves the evaluations that uniform replication would spend on them.

    Designs that differ by less than any number of replicates could
    separate (such as designs near saturation) stay uncertain: give the
    criteria a tolerance of the objectives within which differences do not
    matter, so that such designs are all efficient (see
    `pareto_probability`). A design whose distance from the front is about
    its tolerance stays uncertain too, and more replicates settle which side
    of the tolerance it falls on only slowly, though it does not matter: a
    design is `indifferent`, and not refined, once the standard error of each
    objective that varies between its replicates (see `pareto_probability`)
    is less than `precision` times the objective's tolerance (at the
    design's value), which saves the evaluations that would otherwise go to
    such designs until `max_replicates`. A design whose varying objectives
    lack a tolerance is never indifferent. A warning is given if designs
    remain uncertain (but not indifferent) without a tolerance.

    Args:
        study (Study): The study that evaluates designs.
        criteria (Criteria | Mapping[str, str]): The criteria, or the
            objectives (metrics whose means they are) and their directions
            (`"min"` or `"max"`). A design is uncertain if it is in any
            group of the criteria's `by`.
        results (pandas.DataFrame): The results of the designs (see
            `Study.run`), such as those of an enumeration or a search.
        designs (Sequence[Design | str] | pandas.DataFrame | None): The
            designs that may be refined (designs, their keys, or a data
            frame of them, such as rows of a summary), or `None` for every
            design of the results (rebuilt from them, see
            `DesignSpace.designs_of`).
        bounds (tuple[float, float]): The probabilities between which a
            design is uncertain.
        step (int | None): The replicates added to an uncertain design in a
            round, or `None` for the study's number of replicates.
        max_replicates (int | None): The most replicates of a design, if
            limited.
        min_replicates (int): The fewest replicates with which a design's
            efficiency can be settled.
        max_rounds (int): The most rounds.
        precision (float | None): The standard error of an objective, as a
            share of its tolerance, below which an uncertain design is
            indifferent, or `None` to refine uncertain designs until
            `max_replicates`.
        budget (int | None): The most evaluations to add (a design in a
            replicate, which a study with scenarios evaluates in each), if
            limited.
        n_resamples (int): The number of resamples of each bootstrap.
        seed (int): The seed of the bootstraps.

    Returns:
        Refinement: the results, the final probabilities, the history, and
            the designs
    """
    criteria = _as_criteria(criteria)
    low, high = bounds
    if not 0 <= low < high <= 1:
        raise ValueError("bounds must satisfy 0 <= low < high <= 1")
    step = study.replicates if step is None else step
    if step < 1:
        raise ValueError("step must be at least 1")
    if max_replicates is not None and max_replicates < min_replicates:
        raise ValueError("max_replicates must be at least min_replicates")
    if precision is not None and precision <= 0:
        raise ValueError("precision must be positive (or None)")
    if budget is not None and budget < 0:
        raise ValueError("budget must not be negative")
    refined = _designs(study.space, results if designs is None else designs, results)
    by_key = {design.key: design for design in refined}
    frames = [results]
    history = []
    spent = 0

    def estimate(round_: int) -> pd.DataFrame:
        probability = pareto_probability(
            pd.concat(frames, ignore_index=True),
            criteria,
            n_resamples=n_resamples,
            seed=seed + round_,
        )
        few = probability["replicates"] < min_replicates
        probability["uncertain"] = (
            (probability["pareto_probability"] > low)
            & (probability["pareto_probability"] < high)
        ) | few
        probability["indifferent"] = (
            probability["uncertain"] & ~few & _precise(probability, criteria, precision)
        )
        return probability

    for round_ in range(max_rounds):
        probability = estimate(round_)
        grouped = probability.groupby("design")
        uncertain = grouped["uncertain"].any()
        # a design is indifferent if it is in every group in which it is uncertain
        settling = (
            (probability["uncertain"] & ~probability["indifferent"])
            .groupby(probability["design"])
            .any()
        )
        indifferent = uncertain & ~settling
        combined = pd.concat(frames, ignore_index=True)
        following = combined.groupby("design")["replicate"].max() + 1
        # the designs to refine, the fewest replicates and the least settled first
        distance = (probability["pareto_probability"] - 0.5).abs()
        priority = pd.DataFrame(
            {
                "few": grouped["replicates"].min() < min_replicates,
                "distance": distance.groupby(probability["design"]).min(),
                "replicates": grouped["replicates"].min(),
            }
        )
        candidates = priority[settling.reindex(priority.index, fill_value=False)]
        candidates = candidates.sort_values(
            ["few", "distance", "replicates"], ascending=[False, True, True]
        )
        left = None if budget is None else budget - spent
        added_to = {}
        for key in candidates.index:
            if key not in by_key:
                continue
            count = (
                step
                if max_replicates is None
                else min(step, max_replicates - int(following[key]))
            )
            if left is not None:
                count = min(count, left)
            if count <= 0:
                continue
            added_to[key] = count
            if left is not None:
                left -= count
                if left == 0:
                    break
        record = {
            "round": round_,
            "uncertain": int(uncertain.sum()),
            "indifferent": int(indifferent.sum()),
            "refined": len(added_to),
            "evaluations": 0,
        }
        if not added_to:
            history.append(record)
            break
        # designs with the same next replicates are evaluated together
        plan = pd.DataFrame(
            {
                "start": following[list(added_to)].astype(int),
                "count": list(added_to.values()),
            },
            index=list(added_to),
        )
        for (start, count), keys in plan.groupby(["start", "count"]):
            new = study.run(
                [by_key[key] for key in keys.index],
                replicates=range(int(start), int(start) + int(count)),
            )
            frames.append(new)
            record["evaluations"] += new.groupby(["design", "replicate"]).ngroups
        spent += record["evaluations"]
        history.append(record)
    else:
        probability = estimate(max_rounds)
    remaining = probability.loc[
        probability["uncertain"] & ~probability["indifferent"], "design"
    ].nunique()
    if remaining > 0 and not criteria.tolerance:
        warnings.warn(
            f"{remaining} designs remain uncertain: designs that differ by less "
            + "than replicates separate stay uncertain without a tolerance of "
            + "the objectives (see Criteria)",
            stacklevel=2,
        )
    return Refinement(
        results=pd.concat(frames, ignore_index=True),
        summary=probability,
        history=pd.DataFrame(history),
        designs=refined,
    )
