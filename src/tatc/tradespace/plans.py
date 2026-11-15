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
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from .cost import CostModel
from .criteria import Requirement
from .deployment import DeploymentSchedule, Stage
from .expectations import ExactExpectationsError, exact_expectations
from .failures import FailureModel, RestorationPolicy
from .histories import evaluate_histories, summarize_histories
from .mission import Capability, Mission, MissionHistory
from .references import Reference, _number, configure, value_of
from .simulation import SAMPLINGS
from .space import Design
from .values import ChoiceValues


@dataclasses.dataclass(frozen=True)
class MissionPlan:
    """
    The deployment of the satellites of each design (see
    `DeploymentSchedule.plan`), the policy that restores the failures of its
    satellites (see `RestorationPolicy`), and its costs (see `CostModel`),
    from the values of the design (see `value_of`).

    The schedule, the policy, and the cost model may be functions of the
    (realized) design, which size them by its values with relationships of
    the study's own (TAT-C includes none, see the module
    `tatc.tradespace.cost`): for example, the costs of a type of satellite
    from the mass of its payload, or the capacity of launches to the
    design's altitude::

        def schedule(d):
            return dataclasses.replace(SCHEDULE, launch_mass_capacity=capacity(d.altitude))

        PLAN = MissionPlan(schedule, COST_MODEL, kind="payload")

    The deployment of each design is planned with its policy (see
    `DeploymentSchedule.plan`), which gives its spares: in orbit in each
    plane (`RestorationPolicy.orbit_spares`) and on the ground, which are
    produced (and cost) whether or not the study has failures (see
    `Failures`). A policy may be a decision of each design,
    such as a group of `policy_space` whose values configure the default
    policy (see `configure`), or a function of the design that gives its
    policy (such as one that configures a base policy)::

        space = space.extend(policy=policy_space(orbit_spares=Integer(0, 1), procurement_years=0.5))
        PLAN = MissionPlan(SCHEDULE, COST_MODEL, policy="policy")
        PLAN = MissionPlan(SCHEDULE, COST_MODEL, policy=lambda d: configure(BASE_POLICY, d.policy))

    Launches limited in mass take the masses of the satellites from the cost
    model (`SatelliteCost.mass`, of each type), unless `masses` are given,
    and launches priced by mass take the same masses (see
    `CostModel.estimate`).

    The growth of the costs of the components of the cost model (see
    `Component`) is drawn from scores that are values of the design (such as
    a group of `risk_space`), shared by designs that share components::

        space = space.extend(cost_risk=risk_space(COST_MODEL))
        PLAN = MissionPlan(SCHEDULE, COST_MODEL, kind="payload", scores="cost_risk")

    Args:
        schedule (DeploymentSchedule | Callable): The schedule of the
            deployment, or a function of the design that gives it.
        cost_model (CostModel | Callable | None): The costs, a function of the
            design that gives them, or `None`.
        kind (Reference | None): The type of the satellites (for schedules,
            failures, and costs of several types), as a value of the design.
        policy (RestorationPolicy | Reference): The policy for the failed
            satellites of each design, and its spares (by default, failed
            satellites are restored by producing new ones), or a value of
            the design (see `value_of`) that is its policy, or the values of
            a group (such as of `policy_space`) that configure the default
            policy (see `configure`).
        growth (float | Reference | Mapping): The growth of the costs (see
            `CostModel.estimate`): a factor, a value of the design, or a
            mapping of each element (such as `"development"`) to a factor or
            a value of the design, or to a mapping of each type of satellite
            to one.
        stages (Sequence[Stage] | Callable): The stages of the deployment
            (see `Stage`), or a function of the (realized) design that gives
            them, such as of a decision whether to stage it and of the values
            its rules read. Their rules may also read the capability of the
            satellites operating at a decision, by the study's metrics (see
            `Decision.capability`).
        masses (float | Mapping | Reference | None): The mass (kg) of the
            satellites (or of each type, a mapping), or a value of the
            design, for launches limited (or priced) in mass, or `None` for
            the masses of the cost model.
        scores (Reference | None): The scores from which the growth of the
            components of the cost model is drawn (see `CostModel.draw`), as
            a value of the design (such as a group of `risk_space`), or
            `None` for the expected costs of the components.
    """

    schedule: DeploymentSchedule | Callable[[Design], DeploymentSchedule]
    cost_model: CostModel | Callable[[Design], CostModel | None] | None = None
    _: dataclasses.KW_ONLY
    kind: Reference | None = None
    policy: RestorationPolicy | Reference = RestorationPolicy()
    growth: float | Reference | Mapping[str, float | Reference | Mapping[Any, Any]] = (
        1.0
    )
    stages: Sequence[Stage] | Callable[[Design], Sequence[Stage]] = ()
    masses: float | Mapping[Any, float] | Reference | None = None
    scores: Reference | None = None

    def schedule_of(self, design: Design) -> DeploymentSchedule:
        """
        Gets the schedule of the deployment of a design.

        Args:
            design (Design): The (realized) design.

        Returns:
            DeploymentSchedule: the schedule
        """
        if isinstance(self.schedule, DeploymentSchedule):
            return self.schedule
        return self.schedule(design)

    def policy_of(self, design: Design) -> RestorationPolicy:
        """
        Gets the policy for the failed satellites of a design, and its spares.

        Args:
            design (Design): The (realized) design.

        Returns:
            RestorationPolicy: the policy
        """
        if isinstance(self.policy, RestorationPolicy):
            return self.policy
        policy = value_of(design, self.policy)
        if isinstance(policy, RestorationPolicy):
            return policy
        if isinstance(policy, (Mapping, ChoiceValues)):
            return configure(RestorationPolicy(), policy)
        raise TypeError(
            "the policy of a design is a RestorationPolicy, or values that "
            + f"configure one, not {policy!r}"
        )

    def cost_model_of(self, design: Design) -> CostModel | None:
        """
        Gets the cost model of a design, with the growth of its components
        drawn from the design's scores (if any, see `CostModel.draw`).

        Args:
            design (Design): The (realized) design.

        Returns:
            CostModel | None: the cost model, if any
        """
        model = self.cost_model
        if model is not None and not isinstance(model, CostModel):
            model = model(design)
        if model is None or self.scores is None:
            return model
        return model.draw(value_of(design, self.scores))

    def masses_of(self, design: Design, cost_model: CostModel | None = None) -> Any:
        """
        Gets the masses of the satellites of a design, for launches limited in
        mass: the given masses (or the value of the design they name), or
        those of the types of satellites of its cost model, if any.

        Args:
            design (Design): The (realized) design.
            cost_model (CostModel | None): The design's cost model.

        Returns:
            Any: a mass, a mapping of types to masses, or `None`
        """
        if self.masses is not None:
            if isinstance(self.masses, (int, float, Mapping)):
                return self.masses
            return value_of(design, self.masses)
        if cost_model is None:
            return None
        satellites = cost_model.satellites
        if isinstance(satellites, Mapping):
            masses = {kind: cost.mass for kind, cost in satellites.items()}
            return None if all(m is None for m in masses.values()) else masses
        return satellites.mass

    def stages_of(self, design: Design) -> tuple:
        """
        Gets the stages of the deployment of a design.

        Args:
            design (Design): The (realized) design.

        Returns:
            tuple[Stage, ...]: the stages
        """
        stages = self.stages(design) if callable(self.stages) else self.stages
        return tuple(stages or ())

    def growth_of(self, design: Design) -> Any:
        """
        Gets the growth of the costs of a design.

        Args:
            design (Design): The (realized) design.

        Returns:
            Any: a factor, or a mapping of each element to its factor (or to
                a mapping of each type of satellite to its factor)
        """
        if not isinstance(self.growth, Mapping):
            return _number(design, self.growth)
        growth: dict[str, Any] = {}
        for element, value in self.growth.items():
            if isinstance(value, Mapping):
                growth[element] = {k: _number(design, g) for k, g in value.items()}
            else:
                growth[element] = _number(design, value)
        return growth

    def plan_mission(
        self,
        design: Design,
        layout: Any,
        capability: Capability,
        failures: Failures | None = None,
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

        Returns:
            Mission: the mission
        """
        policy = self.policy_of(design)
        cost_model = self.cost_model_of(design)
        deployment = self.schedule_of(design).plan(
            layout,
            kind=value_of(design, self.kind),
            policy=policy,
            masses=self.masses_of(design, cost_model),
            stages=self.stages_of(design),
        )
        deployment.capability = capability
        return Mission(
            design, deployment, cost_model, self.growth_of(design), capability, failures
        )


@dataclasses.dataclass(frozen=True)
class Failures:
    """
    The failures of the satellites of each design: histories of its mission
    (see `Deployment.simulate_many`), restored by the policy of its
    deployment (see `MissionPlan.policy`: by default by producing new
    satellites; with `RestorationPolicy(restore=False)`, not restored),
    summarized over the
    histories (see `summarize_histories`): the mean of each metric over the
    phases of the mission and the histories, `quantiles` and `requirements`
    of them, and their costs. With requirements on the metrics of the
    `phases`, each history also gives the time until a phase first meets
    each (`<metric> time`, see `summarize_phases`), summarized over the
    histories as the other columns: its quantiles, or the share of the
    histories that meet it by a time (a requirement on it), say more than
    its mean, which is infinite if a history never meets it.

    Without restoration, the means of the metrics and the availability are
    the exact expectations over the failures (see `exact_expectations`) if
    `exact`, if the histories have failures (with a `seed`), and if
    exact expectations apply: if the deployment is not staged (see `Stage`,
    whose rules make its launches depend on the history) and has at most
    `MAX_SETS` (200,000) sets of satellites to enumerate, which is known
    before they are enumerated; otherwise (such as for constellations of
    about 20 satellites or more) they are the means of the histories.

    With a `stress`, the failures are drawn from stressed distributions for
    importance sampling (see `Deployment.simulate_many`), each history
    weighted by its likelihood ratio. With `controls`, the means over the
    histories are estimated with control variates (see
    `summarize_histories`): measures of the missions (see
    `MissionHistory.metrics`, such as `availability`) whose expected values
    are estimated from `control_histories` other histories of the replicate,
    simulated without evaluating their phases (see
    `Deployment.expected_metrics`). With `spread`, the summaries add the
    standard deviation of each column over the histories (`<column>
    histories std`), as `allocate_histories` needs of a pilot study.

    Args:
        model (FailureModel): The failures.
        histories (int): The histories of each replicate.
        seed (Reference | None): The value of the design that seeds the
            failures (see `value_of`, such as of `Uncertain.seed`), or
            `None` for the one history without failures.
        sampling (str): The sampling of the histories (`"stratified"` or
            `"random"`, see `Deployment.simulate_many`).
        phases (Sequence[Requirement]): Requirements on the metrics of the
            phases of each history, whose times it gives (`<metric> time`).
        quantiles (Mapping | None): Quantiles of the columns of the
            histories (see `summarize_histories`).
        requirements (Sequence[Requirement]): Requirements on the columns of
            the histories (see `summarize_histories`), such as on a time.
        exact (bool): `True` for exact expectations without restoration.
        stress (float): The factor of the hazards of stressed failures (1
            for none).
        launch_stress (float | None): The factor of the probability that a
            launch fails, or `None` for the stress.
        controls (Sequence[str]): The measures of the missions that are
            control variates, if any.
        control_histories (int): The histories that estimate the expected
            values of the controls.
        spread (bool): `True` to add the standard deviations of the columns
            over the histories.

    Raises:
        ValueError: if the numbers of histories are less than 1, the
            sampling is unknown, or a stress is not positive.
    """

    model: FailureModel
    _: dataclasses.KW_ONLY
    histories: int = 1
    seed: Reference | None = None
    sampling: str = "stratified"
    phases: Sequence[Requirement] = ()
    quantiles: Mapping[str, float | Sequence[float]] | None = None
    requirements: Sequence[Requirement] = ()
    exact: bool = True
    stress: float = 1.0
    launch_stress: float | None = None
    controls: Sequence[str] = ()
    control_histories: int = 1000
    spread: bool = False

    def __post_init__(self):
        if isinstance(self.controls, str):
            raise TypeError("controls are a sequence of names of measures")
        object.__setattr__(self, "controls", tuple(self.controls))
        object.__setattr__(self, "phases", tuple(self.phases))
        object.__setattr__(self, "requirements", tuple(self.requirements))
        if not all(
            isinstance(requirement, Requirement)
            for requirement in (*self.phases, *self.requirements)
        ):
            raise TypeError("phases and requirements are sequences of Requirements")
        if self.histories < 1 or self.control_histories < 1:
            raise ValueError("histories and control_histories must be at least 1")
        if self.sampling not in SAMPLINGS:
            raise ValueError(
                f"sampling must be one of {SAMPLINGS}, not {self.sampling!r}"
            )
        if self.stress <= 0 or (
            self.launch_stress is not None and self.launch_stress <= 0
        ):
            raise ValueError("stress and launch_stress must be positive")

    def seed_of(self, design: Design) -> Any:
        """
        Gets the seed of the failures of a (realized) design's replicate.

        Args:
            design (Design): The (realized) design.

        Returns:
            Any: the seed, or `None` for the one history without failures
        """
        return value_of(design, self.seed)

    def simulate(self, mission: Mission) -> list[MissionHistory]:
        """
        Simulates the histories of a mission (see `Deployment.simulate_many`),
        restored by the policy of its deployment.

        Args:
            mission (Mission): The mission.

        Returns:
            list[MissionHistory]: the histories
        """
        return mission.deployment.simulate_many(
            self.model,
            self.seed_of(mission.design),
            histories=self.histories,
            sampling=self.sampling,
            stress=self.stress,
            launch_stress=self.launch_stress,
        )

    def evaluate(self, mission: Mission) -> pd.DataFrame:
        """
        Evaluates the histories of a mission (see `simulate` and
        `evaluate_histories`): for each, the mean of each metric over its
        phases, the times until they meet the requirements of the `phases`,
        the measures of its mission, and its costs.

        Args:
            mission (Mission): The mission.

        Returns:
            pandas.DataFrame: a row per history
        """
        return evaluate_histories(
            self.simulate(mission),
            mission.capability,
            cost_model=mission.cost_model,
            growth=mission.growth,
            requirements=self.phases,
        )

    def summarize(self, mission: Mission) -> dict[str, float]:
        """
        Summarizes the histories of a mission as the metrics of its
        evaluation (see `evaluate` and `summarize_histories`): with the
        controls (if any), and with exact expectations of its means if they
        apply (see `Failures`).

        Args:
            mission (Mission): The mission.

        Returns:
            dict[str, float]: the summary
        """
        seed = self.seed_of(mission.design)
        controls = None
        if self.controls and seed is not None:
            # other histories of the replicate (by a seed of its seed)
            other = int(np.random.default_rng([int(seed), 7]).integers(2**31))
            expected = mission.deployment.expected_metrics(
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
            self.evaluate(mission),
            quantiles=self.quantiles,
            requirements=self.requirements,
            controls=controls,
            spread=self.spread,
        )
        if self.exact and seed is not None:
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
