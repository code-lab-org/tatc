"""
Adaptive replication: replicates added only to the designs whose Pareto
efficiency their replicates do not settle.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import pandas as pd

from .criteria import Requirement, Statistic
from .pareto import pareto_probability
from .space import Design
from .study import Study


@dataclass
class Refinement:
    """
    The outcome of an adaptive replication (see `refine`).
    """

    results: pd.DataFrame
    """pandas.DataFrame: The results of every evaluation, those given and
    those added (see `Study.run`)."""
    probability: pd.DataFrame
    """pandas.DataFrame: The probability that each design is efficient after
    the last round (see `pareto_probability`)."""
    history: pd.DataFrame = field(default_factory=pd.DataFrame)
    """pandas.DataFrame: For each round, the number of designs whose
    efficiency is uncertain (`uncertain`), the number to which replicates
    were added (`refined`), and the evaluations added (`evaluations`)."""

    @property
    def uncertain(self) -> pd.DataFrame:
        """pandas.DataFrame: The rows of `probability` whose efficiency is
        still uncertain (as of the bounds of `refine`)."""
        return self.probability[self.probability["uncertain"]]


def refine(
    study: Study,
    designs: Sequence[Design],
    results: pd.DataFrame,
    objectives: Mapping[str, str],
    statistics: Mapping[str, Statistic | str] | None = None,
    requirements: Sequence[Requirement] = (),
    by: Sequence[str] = (),
    bounds: tuple[float, float] = (0.05, 0.95),
    step: int | None = None,
    max_replicates: int | None = None,
    max_rounds: int = 10,
    n_resamples: int = 500,
    seed: int = 0,
    tolerance: Mapping[str, float] | None = None,
) -> Refinement:
    """
    Adds replicates to the designs whose Pareto efficiency is uncertain,
    until it is settled: in each round, the probability that each design is
    efficient (and meets the requirements) is estimated by a bootstrap (see
    `pareto_probability`), and each design whose probability is strictly
    between the `bounds` gets `step` more replicates, up to
    `max_replicates`. Rounds end when no design is uncertain (or can be
    refined), or after `max_rounds`.

    Replicates are added in order (the next indices after a design's
    replicates), so that designs refined together share their uncertain
    parameters (with common random numbers, see `Design.realize`). Designs
    that are surely efficient or surely not keep their replicates, which
    saves the evaluations that uniform replication would spend on them.

    Designs that differ by less than any number of replicates could
    separate (such as designs near saturation) stay uncertain: give a
    `tolerance` of the objectives within which differences do not matter,
    so that such designs are all efficient (see `pareto_probability`).

    Args:
        study (Study): The study that evaluates designs.
        designs (Sequence[Design]): The designs (those of the results).
        results (pandas.DataFrame): The results of the designs (see
            `Study.run`).
        objectives (Mapping[str, str]): The objectives: statistics (or
            metrics) and their directions (`"min"` or `"max"`).
        statistics (Mapping[str, Statistic | str] | None): Statistics, by
            name, of objectives and requirements that are not means.
        requirements (Sequence[Requirement]): Requirements on statistics.
        by (Sequence[str]): Columns that group records of an evaluation (such
            as a sensor); a design is uncertain if it is in any group.
        bounds (tuple[float, float]): The probabilities between which a
            design is uncertain.
        step (int | None): The replicates added to an uncertain design in a
            round, or `None` for the study's number of replicates.
        max_replicates (int | None): The most replicates of a design, if
            limited.
        max_rounds (int): The most rounds.
        n_resamples (int): The number of resamples of each bootstrap.
        seed (int): The seed of the bootstraps.
        tolerance (Mapping[str, float] | None): The tolerances of objectives
            (by name, zero if not given, see `pareto_front`).

    Returns:
        Refinement: the results, the final probabilities, and the history
    """
    low, high = bounds
    if not 0 <= low < high <= 1:
        raise ValueError("bounds must satisfy 0 <= low < high <= 1")
    step = study.replicates if step is None else step
    if step < 1:
        raise ValueError("step must be at least 1")
    by_key = {design.key: design for design in designs}
    frames = [results]
    history = []

    def estimate(round_: int) -> pd.DataFrame:
        probability = pareto_probability(
            pd.concat(frames, ignore_index=True),
            objectives,
            by=by,
            n_resamples=n_resamples,
            seed=seed + round_,
            statistics=statistics,
            requirements=requirements,
            tolerance=tolerance,
        )
        probability["uncertain"] = (probability["pareto_probability"] > low) & (
            probability["pareto_probability"] < high
        )
        return probability

    for round_ in range(max_rounds):
        probability = estimate(round_)
        uncertain = probability.groupby("design")["uncertain"].any()
        uncertain = list(uncertain[uncertain].index)
        combined = pd.concat(frames, ignore_index=True)
        # the next replicate of each uncertain design
        following = combined.groupby("design")["replicate"].max() + 1
        refined = [
            key
            for key in uncertain
            if key in by_key
            and (max_replicates is None or following[key] + step <= max_replicates)
        ]
        if not refined:
            history.append(
                {
                    "round": round_,
                    "uncertain": len(uncertain),
                    "refined": 0,
                    "evaluations": 0,
                }
            )
            break
        added = 0
        # designs with the same next replicate are evaluated together
        for start, keys in pd.Series(following[refined]).groupby(following[refined]):
            new = study.run(
                [by_key[key] for key in keys.index],
                replicates=range(int(start), int(start) + step),
            )
            frames.append(new)
            added += new.groupby(["design", "replicate"]).ngroups
        history.append(
            {
                "round": round_,
                "uncertain": len(uncertain),
                "refined": len(refined),
                "evaluations": added,
            }
        )
    else:
        probability = estimate(max_rounds)
    return Refinement(
        results=pd.concat(frames, ignore_index=True),
        probability=probability,
        history=pd.DataFrame(history),
    )
