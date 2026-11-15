"""
Adaptive replication: replicates added only to the designs whose Pareto
efficiency their replicates do not settle, within a budget of evaluations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .alerts import _warn
from .criteria import Criteria, _as_criteria, _tolerances
from .space import Design, _designs
from .study import Study
from .summary import pareto_probability
from .tables import efficient_keys


@dataclass
class Refinement:
    """
    The outcome of an adaptive replication (see `refine`): as the outcomes of
    other steps of a trade study (see `Confirmation`, `Screening`, and
    `Search`), its `results`, the `summary` of its designs, the `designs` it
    selects for the next step, and its `criteria`::

        refinement = refine(study, criteria, results)
        confirmation = confirm(study, criteria, refinement.results, refinement.designs)

    The designs a refinement selects are those efficient by the criteria over
    every replicate after its last round (`pareto` of its `summary`), which a
    confirmation checks with new replicates (see `confirm`, which selects
    them by default): designs whose efficiency is still uncertain may be
    efficient by more replicates too, those that more replicates would
    settle::

        summary = refinement.summary
        unsettled = summary[summary.uncertain & ~summary.indifferent]
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of every evaluation, those given and
    those added (see `Study.run`)."""
    summary: pd.DataFrame
    """pandas.DataFrame: A row per design (and group of the criteria's `by`):
    the probability that it is efficient after the last round (see
    `pareto_probability`, with the values of the design), whether its
    efficiency is `uncertain` (its probability is strictly between the
    bounds of `refine`, or it has fewer feasible replicates than the least),
    and whether an uncertain design is `indifferent` (its objectives are
    known to within the precision of their tolerances, so that whether it is
    within the tolerances of the front does not matter, and it is not
    refined). Designs that more replicates would settle are those uncertain
    but not indifferent."""
    designs: list[Design]
    """list[Design]: The designs efficient after the last round (`pareto` of
    `summary`, in any group of the criteria's `by`), for the next step."""
    criteria: Criteria
    """Criteria: The criteria of the refinement."""
    rounds: pd.DataFrame = field(default_factory=pd.DataFrame)
    """pandas.DataFrame: A row per round (`round`, from 0): the number of
    designs whose efficiency is uncertain (`uncertain`, in any group), of
    these, those that are indifferent wherever they are uncertain
    (`indifferent`), the number to which replicates were added (`refined`),
    and the evaluations added (`evaluations`, a design in a replicate in
    each scenario of a study with scenarios)."""
    candidates: list[Design] = field(default_factory=list)
    """list[Design]: The designs that could be refined (those given to
    `refine`)."""

    @property
    def evaluations(self) -> int:
        """int: The number of evaluations added over every round (a design in
        a replicate in a scenario: a study with scenarios evaluates each
        replicate in each scenario)."""
        if "evaluations" not in self.rounds:
            return 0
        return int(self.rounds["evaluations"].sum())


def _precise(
    probability: pd.DataFrame, criteria: Criteria, precision: float | None
) -> np.ndarray:
    """
    Gets whether each design's objectives are known to within the precision
    of their tolerances: the standard error of each objective that varies
    between replicates is less than `precision` times its tolerance (at the
    design's value). Designs whose varying objectives lack a tolerance are
    not precise. An objective without a standard error (which
    `pareto_probability` leaves out where it is zero for every design) does
    not vary.
    """
    precise = np.zeros(len(probability), dtype=bool)
    if precision is None or not criteria.tolerances:
        return precise
    tolerances = _tolerances(criteria)
    precise[:] = True
    for i, name in enumerate(criteria.objectives):
        if f"{name} se" not in probability:
            continue
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
    first, and then those whose probability is nearest 0.5. A study with
    scenarios evaluates each replicate in each scenario, which spends an
    evaluation of the budget in each.

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

    Designs are judged by their feasible replicates (see `summarize`): a
    design has fewer than `min_replicates` if fewer of its replicates are
    feasible, and `max_replicates` limits all of its replicates. Designs
    without a feasible replicate are not refined (and a refinement of results
    without a feasible record adds none).

    The refinement selects the designs efficient after its last round for
    the next step (`Refinement.designs`, such as of `confirm`), and keeps
    those that could be refined (`Refinement.candidates`), its rounds
    (`Refinement.rounds`), and the number of evaluations it added
    (`Refinement.evaluations`).

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
            frame of them, such as rows of a summary, or the `designs` of
            another step), or `None` for every design of the results
            (rebuilt from them, see `DesignSpace.designs_of`).
        bounds (tuple[float, float]): The probabilities between which a
            design is uncertain.
        step (int | None): The replicates added to an uncertain design in a
            round, or `None` for the study's number of replicates.
        max_replicates (int | None): The most replicates of a design, if
            limited.
        min_replicates (int): The fewest feasible replicates with which a
            design's efficiency can be settled.
        max_rounds (int): The most rounds.
        precision (float | None): The standard error of an objective, as a
            share of its tolerance, below which an uncertain design is
            indifferent, or `None` to refine uncertain designs until
            `max_replicates`.
        budget (int | None): The most evaluations to add (a design in a
            replicate in a scenario: a replicate of a study with scenarios
            takes an evaluation in each), if limited.
        n_resamples (int): The number of resamples of each bootstrap.
        seed (int): The seed of the bootstraps.

    Returns:
        Refinement: the results, the final probabilities, the rounds, the
            efficient designs, and the candidates
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
    # the evaluations of a replicate of a design: one in each scenario
    cost = len(study.scenarios or [None])
    combined, rounds, spent = results, [], 0

    def estimate(round_: int) -> pd.DataFrame:
        probability = pareto_probability(
            combined, criteria, n_resamples=n_resamples, seed=seed + round_
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
        uncertain, indifferent, settling = _settling(probability)
        following = combined.groupby("design")["replicate"].max() + 1
        added_to = _allocate(
            probability,
            settling,
            following,
            list(by_key),
            step=step,
            min_replicates=min_replicates,
            max_replicates=max_replicates,
            budget=None if budget is None else budget - spent,
            cost=cost,
        )
        record = {
            "round": round_,
            "uncertain": int(uncertain.sum()),
            "indifferent": int(indifferent.sum()),
            "refined": len(added_to),
            "evaluations": 0,
        }
        if not added_to:
            rounds.append(record)
            break
        added = _evaluate_added(study, by_key, following, added_to)
        combined = pd.concat([combined, *added], ignore_index=True)
        record["evaluations"] = sum(_count_evaluations(new) for new in added)
        spent += record["evaluations"]
        rounds.append(record)
    else:
        probability = estimate(max_rounds)
    remaining = probability.loc[
        probability["uncertain"] & ~probability["indifferent"], "design"
    ].nunique()
    if remaining > 0 and not criteria.tolerances:
        _warn(
            f"{remaining} designs remain uncertain: designs that differ by less "
            + "than replicates separate stay uncertain without a tolerance of "
            + "the objectives (see Criteria)"
        )
    selected = [by_key.get(key, key) for key in efficient_keys(probability)]
    return Refinement(
        results=combined,
        summary=probability,
        designs=_designs(study.space, selected, combined),
        criteria=criteria,
        rounds=pd.DataFrame(rounds),
        candidates=refined,
    )


def _settling(probability: pd.DataFrame) -> tuple[pd.Series, pd.Series, pd.Series]:
    """
    Gets, by design, whether its efficiency is uncertain (in any group of the
    criteria's `by`), indifferent (uncertain only where its objectives are
    known to within the precision of their tolerances), and settling (more
    replicates would settle it), for a round of `refine`.
    """
    uncertain = probability.groupby("design")["uncertain"].any()
    settling = (
        (probability["uncertain"] & ~probability["indifferent"])
        .groupby(probability["design"])
        .any()
    )
    return uncertain, uncertain & ~settling, settling


def _allocate(
    probability: pd.DataFrame,
    settling: pd.Series,
    following: pd.Series,
    refinable: list[str],
    *,
    step: int,
    min_replicates: int,
    max_replicates: int | None,
    budget: int | None,
    cost: int = 1,
) -> dict[str, int]:
    """
    Gets the replicates to add to designs whose efficiency more replicates
    would settle (by key), for a round of `refine`: `step` each (up to
    `max_replicates`), the fewest replicates and the least settled first
    (nearest a probability of 0.5), within a budget of evaluations (if any),
    of which each replicate spends `cost` (one in each scenario).
    """
    grouped = probability.groupby("design")
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
    left: int | None = budget
    added: dict[str, int] = {}
    allowed = set(refinable)
    for key in candidates.index:
        if key not in allowed:
            continue
        count = step
        if max_replicates is not None:
            count = min(count, max_replicates - int(following[key]))
        if left is not None:
            count = min(count, left // cost)
        if count <= 0:
            continue
        added[key] = count
        if left is not None:
            left -= count * cost
            if left < cost:
                break
    return added


def _count_evaluations(results: pd.DataFrame) -> int:
    """
    Counts the evaluations of results: a design in a replicate (in a
    scenario, with scenarios), whatever the number of records of each.
    """
    keys = ["design", "replicate", *(["scenario"] if "scenario" in results else [])]
    return results.groupby(keys).ngroups


def _evaluate_added(
    study: Study,
    designs: dict[str, Design],
    following: pd.Series,
    added: dict[str, int],
) -> list[pd.DataFrame]:
    """
    Evaluates the replicates added to designs (by key, after the replicates
    each has): designs with the same new replicates together.
    """
    plan = pd.DataFrame(
        {"start": following[list(added)].astype(int), "count": list(added.values())},
        index=list(added),
    )
    return [
        study.run(
            [designs[key] for key in keys.index],
            replicate_indices=range(int(start), int(start) + int(count)),
        )
        for (start, count), keys in plan.groupby(["start", "count"])
    ]
