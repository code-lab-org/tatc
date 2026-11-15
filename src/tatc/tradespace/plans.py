"""
The missions of the designs of an observation study (see `ObservationStudy`):
the plan of each design's mission (`MissionPlan`: its deployment, the policy
that restores the failures of its satellites, and its costs, from the
values of the design), and the failures of its satellites (`Failures`),
from which an evaluation plans the mission of a replicate of a design (see
`Mission`) and summarizes it into its record.

The parts refer to the values of designs as `tatc.tradespace.references`
describes.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
import numbers
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .cost import CostModel
from .criteria import Requirement
from .deployment import DeploymentSchedule, Stage
from .expectations import ExactExpectationsError, exact_expectations
from .failures import FailureModel, RestorationPolicy
from .histories import _history_statistics, evaluate_histories, summarize_histories
from .mission import Capability, Mission
from .references import Reference, configure, value_of
from .simulation import _check_sampling
from .space import Design
from .statistics import Statistic
from .values import ChoiceValues, Values

_NUMBERS = (numbers.Real, Mapping)
"""The values of fields of numbers: a number, or a mapping (such as of
each type of satellite to its number)."""


def _resolve(design: Design, value: Any, kinds: type | tuple[type, ...] = ()) -> Any:
    """
    Gets the value of a field of a plan for a design: the value given, if it
    is `None` or of the kinds of the field's values, or else the value of the
    design that it refers to (see `value_of`).
    """
    if value is None or isinstance(value, kinds):
        return value
    return value_of(design, value)


def _policy(value: Any) -> RestorationPolicy | None:
    """
    Gets the policy that a value gives: a policy, or the default policy
    configured from the values of a group (a mapping, or a choice of
    groups, see `configure`), or else `None`.
    """
    if isinstance(value, RestorationPolicy):
        return value
    if isinstance(value, ChoiceValues):
        if value.value is not None and not isinstance(value.value, Values):
            return None  # a choice of other values
    elif not isinstance(value, Mapping):
        return None
    return configure(RestorationPolicy(), value)


@dataclasses.dataclass(frozen=True)
class MissionPlan:
    """
    The deployment of the satellites of each design (see
    `DeploymentSchedule.plan`), the policy that restores the failures of its
    satellites (see `RestorationPolicy`), and its costs (see `CostModel`).

    Each part is given as a value, or as a reference to a value of the
    design (see `Reference`: its name or path, or a function of the realized
    design that computes it), which the plan reads for each design: a part
    that is not a value of its kind (a `DeploymentSchedule`, a `CostModel`, a
    `RestorationPolicy`, a sequence of `Stage`, or numbers) is a reference.
    The type of the satellites (`kind`) and the scores of cost risk
    (`scores`) are always references. A function sizes a part by the values
    of the design with relationships of the study's own (TAT-C includes
    none, see the module `tatc.tradespace.cost`): for example, the costs of
    a type of satellite from the mass of its payload, or the capacity of
    launches to the design's altitude::

        def schedule(d):
            return dataclasses.replace(SCHEDULE, launch_mass_capacity=capacity(d.altitude))

        PLAN = MissionPlan(schedule, COST_MODEL, kind="payload", masses=lambda d: mass(d.payload))

    The deployment of each design is planned with its policy (see
    `DeploymentSchedule.plan`), which gives its spares: in orbit in each
    plane (`RestorationPolicy.orbit_spares`) and on the ground, which are
    produced (and cost) whether or not the study has failures (see
    `Failures`). By default, the policy of a design is its value `policy`,
    if it has one that configures a policy (such as a group of
    `policy_space`, whose values configure the default policy, see
    `configure`), and otherwise the default policy, which restores failed
    satellites by producing new ones. A policy may also be given, or a
    reference to another value of the design (a policy, or a group that
    configures one), or a function that gives each design's policy (such as
    one that configures a base policy)::

        space = space.extend(policy=policy_space(orbit_spares=Integer(0, 1), procurement_years=0.5))
        PLAN = MissionPlan(SCHEDULE, COST_MODEL)  # the policy of each design: its values `policy`
        PLAN = MissionPlan(SCHEDULE, COST_MODEL, policy=lambda d: configure(BASE_POLICY, d.policy))

    Launches limited in mass, or priced by mass (see `CostModel.estimate`),
    take the masses of the satellites (`masses`).

    The growth of the costs of the components of the cost model (see
    `Component`) is drawn from scores that are values of the design (such as
    a group of `risk_space`), shared by designs that share components::

        space = space.extend(cost_risk=risk_space(COST_MODEL))
        PLAN = MissionPlan(SCHEDULE, COST_MODEL, kind="payload", scores="cost_risk")

    Args:
        schedule (DeploymentSchedule | Reference): The schedule of the
            deployment.
        cost_model (CostModel | Reference | None): The costs, or `None`.
        kind (Reference | None): The type of the satellites (for schedules,
            failures, and costs of several types).
        policy (RestorationPolicy | Reference | None): The policy for the
            failed satellites of each design, and its spares (or the values
            of a group that configure the default policy), or `None` for the
            design's value `policy`, if it has one that configures a policy,
            and otherwise the default policy.
        growth (float | Reference | Mapping): The growth of the costs (see
            `CostModel.estimate`): a factor, or a mapping of each element
            (such as `"development"`) to a factor, or to a mapping of each
            type of satellite to one (each a number or a reference).
        stages (Sequence[Stage] | Reference): The stages of the deployment
            (see `Stage`), such as a function of the design that gives them
            by a decision whether to stage it and the values its rules read.
            Their rules may also read the capability of the satellites
            operating at a decision, by the study's metrics (see
            `Decision.capability`).
        masses (float | Mapping | Reference | None): The mass (kg) of the
            satellites (or of each type, a mapping), for launches limited (or
            priced) in mass, or `None`.
        scores (Reference | None): The scores from which the growth of the
            components of the cost model is drawn (see `CostModel.draw`),
            such as a group of `risk_space`, or `None` for the expected
            costs of the components.
    """

    schedule: DeploymentSchedule | Reference
    cost_model: CostModel | Reference | None = None
    _: dataclasses.KW_ONLY
    kind: Reference | None = None
    policy: RestorationPolicy | Reference | None = None
    growth: float | Reference | Mapping[str, Any] = 1.0
    stages: Sequence[Stage] | Reference = ()
    masses: float | Mapping[Any, float] | Reference | None = None
    scores: Reference | None = None

    def _policy_of(self, design: Design) -> RestorationPolicy:
        """
        Gets the policy for the failed satellites of a design, and its
        spares, raising a `TypeError` if the policy it refers to is not a
        policy or a group that configures one.
        """
        if self.policy is None:
            named = _policy(design["policy"]) if "policy" in design else None
            return RestorationPolicy() if named is None else named
        policy = _policy(_resolve(design, self.policy, RestorationPolicy))
        if policy is None:
            raise TypeError(
                "the policy of a design is a RestorationPolicy, or values that "
                + f"configure one, not {value_of(design, self.policy)!r}"
            )
        return policy

    def _cost_model_of(self, design: Design) -> CostModel | None:
        """
        Gets the cost model of a design, with the growth of its components
        drawn from the design's scores (if any, see `CostModel.draw`).
        """
        model = _resolve(design, self.cost_model, CostModel)
        if model is None or self.scores is None:
            return model
        return model.draw(value_of(design, self.scores))

    def _growth_of(self, design: Design) -> Any:
        """
        Gets the growth of the costs of a design: a factor, or a mapping of
        each element to its factor (or to a mapping of each type of
        satellite to its factor).
        """
        growth = _resolve(design, self.growth, _NUMBERS)
        if not isinstance(growth, Mapping):
            return growth
        resolved: dict[str, Any] = {}
        for element, value in growth.items():
            value = _resolve(design, value, _NUMBERS)
            if isinstance(value, Mapping):
                value = {k: _resolve(design, g, numbers.Real) for k, g in value.items()}
            resolved[element] = value
        return resolved

    def plan_mission(
        self,
        design: Design,
        layout: Any,
        capability: Capability,
        *,
        failures: Failures | None = None,
        requirements: Sequence[Requirement] | None = None,
    ) -> Mission:
        """
        Plans the mission of a realized design: its deployment (with its
        policy and spares, and whose stages read the capability of the
        satellites operating at their decisions) and its costs.

        Args:
            design (Design): The (realized) design.
            layout (Any): The layout of its satellites (see `layout_of`).
            capability (Capability): The capability of its satellites.
            failures (Failures | None): The failures of its satellites, if
                any.
            requirements (Sequence[Requirement] | None): The requirements on
                the metrics of its phases (see `Mission.requirements`), or
                `None` not to summarize the phases of a mission without
                failures.

        Returns:
            Mission: the mission
        """
        stages = _resolve(design, self.stages, (list, tuple))
        deployment = _resolve(design, self.schedule, DeploymentSchedule).plan(
            layout,
            kind=value_of(design, self.kind),
            policy=self._policy_of(design),
            masses=_resolve(design, self.masses, _NUMBERS),
            stages=tuple(stages or ()),
            capability=capability,
        )
        return Mission(
            design,
            deployment,
            self._cost_model_of(design),
            self._growth_of(design),
            failures,
            None if requirements is None else tuple(requirements),
        )


@dataclasses.dataclass(frozen=True)
class Failures:
    """
    The failures of the satellites of each design: histories of its mission
    (see `Deployment.simulate_many`), drawn from a `seed` of each replicate
    and restored by the policy of its deployment (see `MissionPlan.policy`:
    by default by producing new satellites; with
    `RestorationPolicy(restore=False)`, not restored), summarized over the
    histories (see `summarize_histories`): the mean of each metric over the
    phases of the mission and the histories (and the time until a phase
    first meets each requirement of the study on its phases, see
    `ObservationStudy`), `statistics` of them over the histories (such as a
    `Quantile` of the cost), and their costs::

        statistics = [Quantile("accurate mission", 0.1), Quantile("cost", 0.7)]
        FAILURES = Failures(FAILURE_MODEL, seed="failure_seed", histories=20, statistics=statistics)

    The time until a history meets a requirement is infinite if it never
    does, so that its mean is infinite if any history never meets it, and
    its quantiles, or the share of the histories that meet it by a time (a
    `Probability`), summarize it.

    Without restoration, the means of the metrics and the availability are
    the exact expectations over the failures (see `exact_expectations`) if
    they apply: if the deployment is not staged (see `Stage`, whose rules
    make its launches depend on the history) and has at most `MAX_SETS`
    (200,000) sets of satellites to enumerate, which is known before they
    are enumerated; otherwise (such as for constellations of about 20
    satellites or more) they are the means of the histories.

    With a `stress`, the failures of satellites and launches are drawn from
    stressed distributions for importance sampling (see
    `Deployment.simulate_many`), each history weighted by its likelihood
    ratio. With `controls`, the means over the histories are estimated with
    control variates (see `summarize_histories`): measures of the missions
    (see `MissionHistory.metrics`, such as `availability`) whose expected
    values are estimated from `control_histories` other histories of the
    replicate, simulated without evaluating their phases.

    Args:
        model (FailureModel): The failures.
        seed (Reference): The value of the design that seeds the failures of
            its replicate (see `value_of`), such as of `Uncertain.seed`.
        histories (int): The histories of each replicate.
        sampling (str): The sampling of the histories (`"stratified"` or
            `"random"`, see `Deployment.simulate_many`).
        statistics (Sequence[Statistic]): Statistics of the columns of the
            histories over the histories: `Quantile`, `Probability`, or
            `StandardDeviation` (see `summarize_histories`).
        stress (float): The factor of the hazards of stressed failures, and
            of the probability that a launch fails (1 for none).
        controls (Sequence[str]): The measures of the missions that are
            control variates, if any.
        control_histories (int): The histories that estimate the expected
            values of the controls.

    Raises:
        TypeError: if the statistics are not of these kinds.
        ValueError: if the numbers of histories are less than 1, the
            sampling is unknown, or the stress is not positive.
    """

    model: FailureModel
    _: dataclasses.KW_ONLY
    seed: Reference
    histories: int = 1
    sampling: str = "stratified"
    statistics: Sequence[Statistic] = ()
    stress: float = 1.0
    controls: Sequence[str] = ()
    control_histories: int = 1000

    def __post_init__(self):
        if isinstance(self.controls, str):
            raise TypeError("controls are a sequence of names of measures")
        object.__setattr__(self, "controls", tuple(self.controls))
        object.__setattr__(self, "statistics", _history_statistics(self.statistics))
        if self.seed is None:
            raise ValueError("failures need a seed: a value of the design")
        if self.control_histories < 1:
            raise ValueError("control_histories must be at least 1")
        _check_sampling(self.histories, self.sampling, self.stress)

    def _seed_of(self, design: Design) -> Any:
        """
        Gets the seed of the failures of a (realized) design's replicate,
        raising a `ValueError` if it is `None` (of a design that is not
        realized).
        """
        seed = value_of(design, self.seed)
        if seed is None:
            raise ValueError(
                "the seed of the failures of a design is None: failures need "
                + "a realized design (with a seed, such as of Uncertain.seed)"
            )
        return seed

    def _evaluate(self, mission: Mission) -> pd.DataFrame:
        """
        Evaluates the histories of a mission (see `evaluate_histories`): for
        each, the mean of each metric over its phases, the times until they
        meet the requirements of the mission, the measures of its mission,
        and its costs.
        """
        histories = mission.deployment.simulate_many(
            self.model,
            self._seed_of(mission.design),
            histories=self.histories,
            sampling=self.sampling,
            stress=self.stress,
        )
        return evaluate_histories(
            histories,
            mission.capability,
            cost_model=mission.cost_model,
            growth=mission.growth,
            requirements=mission.requirements or (),
        )

    def _summarize(self, mission: Mission) -> dict[str, float]:
        """
        Summarizes the histories of a mission as the metrics of its
        evaluation (see `summarize_histories`): with the controls (if any),
        and with exact expectations of its means if they apply.
        """
        seed = self._seed_of(mission.design)
        controls = None
        if self.controls:
            # other histories of the replicate (by a seed of its seed)
            other = int(np.random.default_rng([int(seed), 7]).integers(2**31))
            # pylint: disable-next=protected-access
            expected = mission.deployment._expected_metrics(
                self.model,
                other,
                histories=self.control_histories,
                sampling=self.sampling,
            )
            unknown = sorted(set(self.controls) - set(expected))
            if unknown:
                raise ValueError(
                    f"controls that are not measures of missions: {unknown} "
                    + f"(they are {sorted(expected)})"
                )
            controls = {name: expected[name] for name in self.controls}
        record = summarize_histories(
            self._evaluate(mission), statistics=self.statistics, controls=controls
        )
        try:
            exact = exact_expectations(
                mission.deployment, self.model, mission.capability
            )
        except ExactExpectationsError:
            exact = {}  # the means of the histories instead
        record.update(
            {
                key: value
                for key, value in exact.items()
                if key.endswith(" mission") or key == "availability"
            }
        )
        return record
