"""
Exact expectations of the capability of missions without restoration over
their failures (`exact_expectations`), rather than from samples of their
histories (see `Deployment.simulate_many`), as the failures of a study
take them where they apply (see `Failures`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import itertools
import math
from typing import TYPE_CHECKING

import numpy as np

from .failures import FailureModel
from .mission import _PHASE_COLUMNS, Capability, _position_names

if TYPE_CHECKING:
    from .deployment import Deployment

MAX_SETS = 4_096
"""The most sets of satellites that exact expectations enumerate by default
(see `exact_expectations`), as those of the failures of a study do (see
`Failures`): those of about 12 satellites, whose metrics they evaluate
(nearly every set) in about 35 times as many evaluations as 20 histories."""


class ExactExpectationsError(ValueError):
    """
    Raised if exact expectations do not apply to a mission (see
    `exact_expectations`): of a staged deployment, of a deployment whose
    policy restores failures, or of more sets of satellites than they
    enumerate. Its histories are simulated instead.
    """


def _launch_groups(
    deployment: Deployment, model: FailureModel, times: np.ndarray
) -> tuple[int, list[tuple]]:
    """
    Gets the launch groups of a deployment without restoration: the
    positions of each launch of the deployment (in order of the planes'
    positions), which its planned replacements carry again. For each group,
    at each time: the probability that each of its launches so far is the
    latest that succeeded (`weights`, launches by times), the survival of
    the satellites of each since its launch (`survival`, launches by times),
    the probability that none has succeeded (`none`), and whether the group
    has launched.
    """
    offsets = np.cumsum([0] + [plane.satellites for plane in deployment.planes])
    launches: dict[tuple[int, int], list] = {}
    for launch in sorted(deployment.launches, key=lambda l: l.time):
        launches.setdefault((launch.plane, launch.first), []).append(launch)
    p = model.launch
    groups = []
    for (plane, first), sequence in launches.items():
        built = deployment.planes[plane].satellites
        members = [
            int(offsets[plane]) + position
            for position in range(first, min(first + sequence[0].satellites, built))
        ]
        if not members:
            continue
        life = model.life_of(sequence[0].kind)
        flown = np.array([[t >= launch.time for t in times] for launch in sequence])
        # launches after each one, up to each time
        later = np.flip(np.cumsum(np.flip(flown, axis=0), axis=0), axis=0) - flown
        weights = np.where(flown, (1 - p) * p**later, 0.0)
        survival = np.array(
            [
                (
                    np.ones(len(times))
                    if life is None
                    else life.survival(np.maximum(times - launch.time, 0.0))
                )
                for launch in sequence
            ]
        )
        none = p ** flown.sum(axis=0)
        groups.append((members, weights, survival, none, flown.any(axis=0)))
    return int(offsets[-1]), groups


def _group_operating(group: tuple, count: int) -> np.ndarray:
    """
    Gets the probability at each time that exactly a given set of `count`
    satellites of a launch group operates (see `_set_weight`).
    """
    members, weights, survival, none, _ = group
    probability = np.sum(
        weights * survival**count * (1 - survival) ** (len(members) - count), axis=0
    )
    return probability + none if count == 0 else probability


def _set_weight(active: np.ndarray, groups: list) -> float:
    """
    Gets the share of the mission's time in which exactly the active
    positions operate, in expectation: at each time, the product over the
    launch groups of the probability that exactly the active satellites of
    the group operate, which (as they are of one type and launched together)
    depends only on their number `a` of the `n` of the group: the sum over
    its launches of the probability that the launch is the latest that
    succeeded times `S ** a * (1 - S) ** (n - a)`, for the survival `S` of
    its satellites, and the probability that none has succeeded if `a` is 0.
    """
    probability = 1.0
    for group in groups:
        probability = probability * _group_operating(group, int(active[group[0]].sum()))
    return float(np.mean(probability))


def _missing_levels(groups: list, count: int, tolerance: float) -> int:
    """
    Gets the most positions missing of the sets that exact expectations
    enumerate (see `exact_expectations`): the least number such that the
    sets of that many missing positions or fewer operate all but `tolerance`
    of the mission's time, in expectation. The distribution of the number
    missing at each time is the convolution over the launch groups of the
    distributions of theirs (each the probability of a set times the number
    of sets of its size), cheap to compute without enumerating the sets.
    """
    if not groups:
        return count
    distribution = np.zeros((count + 1, len(groups[0][3])))
    distribution[0] = 1.0
    for group in groups:
        size = len(group[0])
        combined = np.zeros_like(distribution)
        for gone in range(size + 1):
            share = math.comb(size, gone) * _group_operating(group, size - gone)
            combined[gone:] += share * distribution[: count + 1 - gone]
        distribution = combined
    covered = np.cumsum(distribution.mean(axis=1))
    reached = np.flatnonzero(covered >= 1 - tolerance - 1e-12)
    return int(reached[0]) if len(reached) > 0 else count


def _sets(count: int, levels: int) -> int:
    """The sets of `count` positions of which at most `levels` are missing."""
    return sum(math.comb(count, missing) for missing in range(levels + 1))


def _exact_sets(
    deployment: Deployment,
    model: FailureModel,
    tolerance: float = 1e-3,
    steps: int = 400,
) -> int:
    """
    Gets the number of sets of satellites that exact expectations of a
    deployment enumerate (see `exact_expectations`), without enumerating
    them: the sets of at most the number of positions missing that cover
    all but `tolerance` of the mission's time.
    """
    count, groups = _launch_groups(deployment, model, _times(deployment, steps))
    return _sets(count, _missing_levels(groups, count, tolerance))


def _times(deployment: Deployment, steps: int) -> np.ndarray:
    """The times at which exact expectations integrate probabilities."""
    begin, end = deployment.begin, deployment.end
    return begin + (np.arange(steps) + 0.5) * (end - begin) / steps


def exact_expectations(
    deployment: Deployment,
    model: FailureModel,
    capability: Capability,
    *,
    tolerance: float = 1e-3,
    steps: int = 400,
    max_sets: int = MAX_SETS,
) -> dict[str, float]:
    """
    Computes the expected capability of a mission without restoration (of a
    deployment whose policy does not restore failed satellites, see
    `RestorationPolicy`, though they may be replaced as planned at the end
    of their design life) over its failures: the expectation of the mean of
    each metric over the phases of the mission (`<metric> mission`, as
    `summarize_phases` gives for one history) and of its availability
    (`availability`), exactly rather than from samples of histories (see
    `Deployment.simulate_many`).

    Without restoration, the satellites fail independently (but for the
    satellites of a launch, which a failed launch loses together), so that
    the probability that a set of satellites operates at a time is a product
    over the launches of their survival (see `Lifetime.survival`). The
    expected mean of a metric is the sum over the sets of satellites of the
    metric of each set (evaluated once, as with a `PhaseCache`) times the
    share of the mission's time in which exactly that set operates, in
    expectation (integrated over `steps` times). The sets are enumerated by
    the number of satellites missing until their shares cover all but
    `tolerance` of the mission's time (a number known before they are
    enumerated, from the distribution of the number missing), and evaluated
    in order of their shares until they do (`coverage`, the share covered,
    and `sets`, the sets evaluated); the expectation is over the sets
    evaluated.

    Without sampling error, an expectation estimates the mean of many
    histories with none of their variance (see `Deployment.simulate_many`),
    for missions of a few to about a dozen satellites, but evaluates the
    metrics of nearly every set (about `2 ** n` of `n` satellites, against
    the few sets that histories evaluate); for larger ones (whose sets
    number more than `max_sets`, which is known before any is enumerated),
    or with restoration (whose shared production and launches make failures
    dependent), simulate histories.

    Args:
        deployment (Deployment): The deployment, whose policy does not
            restore failures.
        model (FailureModel): The failures.
        capability (Capability): The capability of the satellites (in the
            order of the planes, see `MissionHistory.phases`): their
            observations and metric (see `Capability.evaluate`).
        tolerance (float): The share of the mission's time that the sets
            evaluated may leave out.
        steps (int): The times at which the probabilities are integrated.
        max_sets (int): The most sets to enumerate.

    Returns:
        dict[str, float]: the expected metrics

    Raises:
        ExactExpectationsError: if the deployment is staged (see `Stage`),
            whose launches depend on its rules, its policy restores failures,
            or it has more sets to enumerate than `max_sets` (raised before
            any is enumerated).
    """
    # pylint: disable=too-many-locals
    if deployment.stages:
        raise ExactExpectationsError(
            "exact expectations need a deployment that is not staged: "
            + "simulate histories of a staged deployment instead"
        )
    if deployment.policy.restore:
        raise ExactExpectationsError(
            "exact expectations need a deployment whose policy does not restore "
            + "failures (RestorationPolicy(restore=False)): simulate histories instead"
        )
    names = list(_position_names(deployment.planes, capability.names).values())
    count, groups = _launch_groups(deployment, model, _times(deployment, steps))
    levels = _missing_levels(groups, count, tolerance)
    if _sets(count, levels) > max_sets:
        raise ExactExpectationsError(
            "too many sets of satellites to enumerate: simulate histories"
        )
    # availability: the expected years of operation of the positions, as a
    # share of those of the planned deployment
    operating = sum(
        len(members) * float(np.mean(np.sum(weights * survival, axis=0)))
        for members, weights, survival, _, _ in groups
    )
    planned = sum(
        len(members) * float(np.mean(launched)) for members, _, _, _, launched in groups
    )
    availability = operating / planned if planned > 0 else math.nan
    # enumerate the sets by the number of positions missing
    weights: list[tuple[float, tuple[int, ...]]] = []
    for missing in range(levels + 1):
        for gone in itertools.combinations(range(count), missing):
            active = np.ones(count, dtype=bool)
            active[list(gone)] = False
            weight = _set_weight(active, groups)
            if weight > 0:
                weights.append((weight, tuple(np.flatnonzero(active))))
    weights.sort(key=lambda item: -item[0])
    totals: dict[str, float] = {}
    evaluated, sets = 0.0, 0
    for weight, members in weights:
        # the metrics of the set, as of a phase in which it operates
        for name, value in capability(names[k] for k in members).items():
            if name in _PHASE_COLUMNS:
                continue
            totals[f"{name} mission"] = totals.get(
                f"{name} mission", 0.0
            ) + weight * float(value)
        evaluated += weight
        sets += 1
        if evaluated >= 1 - tolerance:
            break
    result = {name: total / evaluated for name, total in totals.items()}
    result.update(availability=availability, coverage=evaluated, sets=float(sets))
    return result
