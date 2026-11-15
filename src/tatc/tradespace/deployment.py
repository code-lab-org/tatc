"""
Deployment of constellations: when the satellites of a design are developed,
produced, launched, and replaced, which dates their costs (see `CostModel`)
and sets the satellites that observe in each phase of the mission, for
measures of capability over the mission rather than of the full
constellation alone.

A `DeploymentSchedule` holds the rates and durations of a study (the
mission, the design life of satellites, the development and production of
each type of satellite, and the launch cadence and capacity), from which
each design's `Deployment` follows: a design of more satellites and planes
takes longer to produce and launch. A deployment is planned with a policy
that restores the failures of its satellites or not (see
`RestorationPolicy`), and its spares, and simulates a mission
(`Deployment.simulate`) with the failures of its satellites and launches
(see `FailureModel`), drawn for each replicate: the `MissionHistory` of
what is produced and launched, which a `CostModel` prices, and of the
satellites that operate. Their changes divide the mission into phases
(`MissionHistory.phases`) of constant sets of operating satellites, each
evaluated over the study's representative period (`evaluate_phases`, which
states what these quasi-static phases assume) and combined over the mission
(`summarize_phases`), such as the time-weighted
mean of a metric and the time until it first meets a requirement.

A deployment may be staged (see `Stage`): planes that are produced and
launched only if a rule decides to, at times during the mission, from what
it knows then (see `Decision`, which includes the capability of the
satellites in orbit), such as a second plane that is launched only if the
mission turns out to need it. Whether and how to stage a deployment
is a decision of a design, whose flexibility (the cost saved where a stage
is not needed, against the capability lost by deferring it where it is) the
analyses of a study compare with deployments that are not staged.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
import numbers
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from ..schemas import Satellite
from .constellations import layout_of
from .events import (
    Launch,
    PerType,
    PlannedPlane,
    Production,
    _Manifest,
    per_type,
    per_type_values,
)
from .failures import FailureModel, RestorationPolicy
from .mission import Capability, Decision, MissionHistory, Phase
from .simulation import _Simulation, simulate_histories
from .values import Values


@dataclass(frozen=True)
class Stage:
    """
    A stage of a deployment (see `DeploymentSchedule.plan`): planes that are
    deployed only if (and when) a rule decides to. At each of its `times`
    (years after the first launch) until it is deployed, the rule is given
    what is known then (see `Decision`) and returns whether to deploy the
    stage. A deployed stage's satellites are produced on their type's line
    `lead_years` after the decision (after the work before them), and
    launched at the next launch slots (in as many launches as they need, see
    `DeploymentSchedule`); they are then replaced at the end of their design
    life and restored if they fail, as the satellites of other planes. A
    stage that is not deployed is neither produced nor launched, and its
    satellites never operate; but a type of satellite whose planes are all
    staged is still developed, and costs its development (see
    `CostModel.estimate`), whether or not its stages are deployed: the cost
    of the option to deploy them.

    A rule that depends on the values of a design (such as an uncertain
    parameter or a value of a scenario, known by the time of the decision)
    reads them from the design, as stages made for each design (see
    `MissionPlan`) can: for example, `Stage([1], 1.0, lambda decision:
    design.decorrelation_time < timedelta(hours=6))`. A rule may instead
    read what the satellites in orbit achieve (see `Decision.capability`),
    as a mission that does not know such values measures them: for example,
    `Stage([1], 1.0, lambda decision: decision.capability["accurate"] <
    50)`, which deploys the plane if the satellites operating then sample
    fewer than half the points accurately.

    Args:
        planes (Sequence[int] | int): The planes of the stage (their indices,
            in the order of the layouts' planes, such as integers of numpy),
            produced in this order (on their types' lines) and launched in
            the order they are ready.
        times (float | Sequence[float]): The times of the decisions, in years
            after the first launch (a number, or a sequence or array of
            them).
        rule (Callable | None): A function of a `Decision` that returns
            whether to deploy the stage, or `None` to deploy it at the first
            time (a planned delay).
        lead_years (float): The time from a decision to the start of the
            production of the stage's satellites (such as to procure
            long-lead items).
    """

    planes: Sequence[int] | int
    times: float | Sequence[float] = 0.0
    rule: Callable[[Decision], bool] | None = None
    lead_years: float = 0.0

    def __post_init__(self):
        planes = self.planes
        if isinstance(planes, numbers.Integral) or np.ndim(planes) == 0:
            planes = (planes,)
        planes = tuple(int(plane) for plane in planes)
        if not planes:
            raise ValueError("a stage needs planes")
        if len(set(planes)) != len(planes) or min(planes) < 0:
            raise ValueError("the planes of a stage must be distinct indices")
        times = self.times
        if isinstance(times, numbers.Real) or np.ndim(times) == 0:
            times = (times,)
        times = tuple(sorted(float(time) for time in times))
        if not times or times[0] < 0:
            raise ValueError("a stage needs times of decisions, not negative")
        if self.lead_years < 0:
            raise ValueError("lead_years must not be negative")
        object.__setattr__(self, "planes", planes)
        object.__setattr__(self, "times", times)


@dataclass
class Deployment:
    """
    The planned deployment of a constellation: its development, production,
    launches, and planned replacements (see `DeploymentSchedule.plan`),
    without failures, and the policy with which it restores the failures of
    its satellites (see `simulate`).
    """

    schedule: DeploymentSchedule
    """DeploymentSchedule: The schedule of the deployment."""
    planes: list[PlannedPlane] = field(default_factory=list)
    """list[PlannedPlane]: The planes, in the order of the layouts' planes:
    the type of their satellites, the number of satellites of the layout,
    and of spares in orbit."""
    development: dict[Any, float] = field(default_factory=dict)
    """dict: The end of the development of each type of satellite (years)."""
    production: list[Production] = field(default_factory=list)
    """list[Production]: The production of the satellites, in order (of a
    sustained constellation, also of replacements launched after the end of
    the mission, see `DeploymentSchedule`)."""
    launches: list[Launch] = field(default_factory=list)
    """list[Launch]: The launches (and planned replacements) of the mission,
    in order."""
    masses: dict[Any, float | None] = field(default_factory=dict)
    """dict: The mass (kg) of each type of satellite, if given."""
    plane_values: list[Any] = field(default_factory=list)
    """list: The values of each plane in the design (see `Plane`), in the
    order of `planes`, as given to functions of the schedule of a plane."""
    stages: tuple[Stage, ...] = ()
    """tuple[Stage, ...]: The stages deployed only if their rules decide to
    (whose planes the production and launches above do not include)."""
    policy: RestorationPolicy = field(default_factory=RestorationPolicy)
    """RestorationPolicy: The policy for failed satellites, whose spares in
    orbit the planes above include (and whose spares on the ground its
    simulations produce, see `simulate`)."""
    capability: Capability | None = field(default=None, repr=False, compare=False)
    """Capability | None: The capability of the satellites, which the rules
    of stages may read at their decisions (see `Decision.capability`) in
    every simulation of the deployment, as the evaluations of an
    `ObservationStudy` with a `MissionPlan` give it."""
    _planned: tuple | None = field(default=None, repr=False, compare=False)

    @property
    def begin(self) -> float:
        """float: The first launch, which begins the mission (years)."""
        return min(launch.time for launch in self.launches)

    @property
    def end(self) -> float:
        """float: The end of the mission (years)."""
        return self.begin + self.schedule.years

    @property
    def deployed(self) -> float:
        """float: The time by which every satellite of the first generation
        is launched (years): of the planes that are not staged, and of the
        stages that their rules deploy, in the history of the mission without
        failures (see `planned` and `MissionHistory.deployed`, which is
        infinite if the deployment does not finish within the mission)."""
        return self.planned().deployed

    @property
    def units(self) -> dict[Any, int]:
        """dict: The satellites of each type that the plan produces: of the
        planes that are not staged (with their spares in orbit) and their
        planned replacements (and, of a sustained constellation, those in
        production at the end of the mission), but not the spares on the
        ground of the policy, the planes of stages, or restorations, which
        its simulations produce (see `planned`, whose production a
        `CostModel` prices)."""
        units: dict[Any, int] = {}
        for batch in self.production:
            if not self.schedule.sustained or batch.begin < self.end:
                units[batch.kind] = units.get(batch.kind, 0) + batch.satellites
        return units

    def simulate(
        self,
        model: FailureModel | None = None,
        seed: int | Sequence[int] | None = None,
    ) -> MissionHistory:
        """
        Simulates the mission: the deployment, with the failures of its
        satellites and launches (drawn from a seed of a replicate, such as an
        uncertain parameter, see `Uncertain.seed`) and the restoration of
        failed satellites by its policy (see `policy`, whose spares in orbit
        the plan produces and launches with their planes, and whose spares
        on the ground the simulation produces after the planes of the plan).
        Without a seed, nothing fails (but the spares are still produced).
        To simulate another policy, plan the deployment with it (see
        `DeploymentSchedule.plan`).

        A satellite's life is drawn when it is launched (or takes the place
        of a failed one), from uniform numbers of a random number generator
        seeded by the seed and its plane, position, and generation; the
        success of a planned launch (of the deployment, or a replacement)
        from one seeded by its plane, the position of its first satellite,
        and its generation; and that of a launch of restorations or spares
        from one seeded by its plane and its order among them: with common
        random numbers, the satellite at the same position of the same plane
        fails at the same age in the same replicate of every design. Spares
        in orbit do not age or fail until they take a failed satellite's
        place, which they then hold (and which planned replacements replace).

        A planned replacement (at the end of a design life) whose launch
        fails is not flown again: the satellite it would replace operates
        until it fails (and is then restored, by a policy that restores
        failures).

        Args:
            model (FailureModel | None): The failures, or `None` for none.
            seed (int | Sequence[int] | None): The seed of the replicate (or
                a sequence of integers, such as the seed and the index of a
                history, see `simulate_many`), or `None` for no failures.

        Returns:
            MissionHistory: the history of the mission
        """
        return _Simulation(self, model or FailureModel(), seed).run()

    def planned(self) -> MissionHistory:
        """
        Gets the history of the mission without failures (see `simulate`),
        simulated once (for its policy and capability): its deployment, with
        the stages that their rules deploy and the spares of its policy,
        whose costs a `CostModel` estimates (see `CostModel.estimate`).

        Returns:
            MissionHistory: the history
        """
        key = (self.policy, self.capability)
        if self._planned is None or any(
            a is not b for a, b in zip(self._planned[0], key)
        ):
            self._planned = (key, self.simulate())
        return self._planned[1]

    def simulate_many(
        self,
        model: FailureModel | None = None,
        seed: int | None = None,
        *,
        histories: int = 1,
        sampling: str = "stratified",
        stress: float = 1.0,
        launch_stress: float | None = None,
    ) -> list[MissionHistory]:
        """
        Simulates several histories of the mission in a replicate (see
        `simulate`), history `k` drawn from the seed and `k`, so that, with
        common random numbers, the same history of the same replicate draws
        the same failures in every design. Without a seed, the one history
        without failures.

        The failures of the histories are stratified by default
        (`"stratified"`: a Latin hypercube over the histories, in which each
        life and launch of each history draws from its own stratum of its
        distribution, see `Lifetime.from_uniform`), which spreads the
        histories over the distributions of the failures and estimates the
        means of the measures of the mission more precisely than independent
        histories (`"random"`), several times so for missions of a few
        satellites or many, at no cost.

        With a `stress`, the failures are drawn from stressed distributions
        for importance sampling, which samples rare outcomes (such as the
        loss of several satellites of a reliable constellation) more often:
        the hazards of the lives of the satellites times the stress, and the
        probability that a launch fails times the `launch_stress` (by
        default the stress, at most 0.5). Each history is weighted by its
        likelihood ratio (`MissionHistory.weight`, see `Lifetime.stressed`),
        which `summarize_histories` applies to every statistic. A stress
        that samples the outcomes of interest about as often as the others
        estimates their probabilities (and the quantiles of their tails)
        with far fewer histories; too large a stress makes a few weights
        dominate (see the effective number of histories).

        Many histories of the failures of a replicate, whose sample of
        targets (the expensive part of an evaluation) they share, measure the
        risk of a mission at little cost (see `evaluate_histories` and
        `PhaseCache`).

        Args:
            model (FailureModel | None): The failures, or `None` for none.
            seed (int | None): The seed of the replicate, or `None` for the
                one history without failures.
            histories (int): The number of histories.
            sampling (str): `"stratified"` or `"random"`.
            stress (float): The factor of the hazards of stressed failures.
            launch_stress (float | None): The factor of the probability that
                a launch fails, or `None` for the stress.

        Returns:
            list[MissionHistory]: the histories
        """
        if histories < 1:
            raise ValueError("histories must be at least 1")
        if seed is None:
            return [self.simulate(model, None)]
        return simulate_histories(
            self,
            model or FailureModel(),
            seed,
            histories=histories,
            sampling=sampling,
            stress=stress,
            launch_stress=launch_stress,
        )

    def expected_metrics(
        self,
        model: FailureModel | None = None,
        seed: int = 0,
        *,
        histories: int = 1000,
        sampling: str = "stratified",
    ) -> dict[str, float]:
        """
        Estimates the expected measures of the mission (see
        `MissionHistory.metrics`, such as its `availability`) from many
        simulated histories (see `simulate_many`), without evaluating their
        phases: precise and cheap, as the expected values of control
        variates (see `summarize_histories`, and the `controls` of
        `Failures`).

        Args:
            model (FailureModel | None): The failures, or `None` for none.
            seed (int): The seed of the histories.
            histories (int): The number of histories.
            sampling (str): The sampling of their failures.

        Returns:
            dict[str, float]: the mean of each measure
        """
        simulated = self.simulate_many(
            model, seed, histories=histories, sampling=sampling
        )
        rows = pd.DataFrame([history.metrics() for history in simulated])
        return {column: float(rows[column].mean()) for column in rows.columns}

    def phases(self, satellites: Sequence[Satellite | str]) -> list[Phase]:
        """
        Gets the phases of the mission without failures (see
        `MissionHistory.phases`), each evaluated over the same window (see
        `evaluate_phases`).

        Args:
            satellites (Sequence[Satellite | str]): The satellites (or their
                names), in the order of the planes.

        Returns:
            list[Phase]: the phases, in order
        """
        return self.planned().phases(satellites)


@dataclass(frozen=True)
class DeploymentSchedule:
    """
    A schedule of the deployment of constellations, from rates and
    durations, so that each design's deployment follows from its satellites
    and planes (see `plan`). Times are in years from the start of
    development.

    * Each type of satellite is developed over `development_years` (one for
      every type, or one for each).
    * Its satellites are then produced on one line at `production_rate`
      satellites a year (or at once, if `None`), plane by plane, in the
      order of the layouts' planes (or of `order`), with the spares in orbit
      of the deployment's policy. The line also produces, after the work
      before them, the spares on the ground of the policy, the planes of
      stages that are deployed, and the satellites that restore failed ones
      (see `RestorationPolicy` and `Stage`).
    * Each plane is launched when its satellites are ready, in as many
      launches as its satellites need (given a `launch_capacity` in number
      or mass, or one per plane), at most one every `launch_interval` years,
      in the order the planes are ready.
    * The mission lasts `years` from the first launch. Satellites are
      replaced at the end of their `design_life` (one for every type, or one
      for each; if before the end of the mission): those at the positions of
      the layout that each launch carries, not its spares in orbit, which do
      not age (see `Deployment.simulate`). Each replacement is produced at
      the production rate just before it is needed (planned apart from the
      line, so that replacements do not wait for each other) and launched
      at the next launch slot then (keeping the launch interval), and is
      replaced in turn a design life after its launch. A design life does
      not retire a satellite: one that is not replaced (as near the end of
      the mission, or if its replacement's launch fails) operates until it
      fails (by the wear-out of its life, see `FailureModel`).
    * A `sustained` constellation (such as of an operational program, whose
      mission is the period of a study) is replaced after the end of the
      mission as before it: the replacements needed after the end are
      produced (their production begun before the end, or within the
      development years of their type after it, for the development of a
      later block of production, see `SatelliteCost.block_size`), though
      not launched in the mission, and its costs are those incurred by the
      end (see `CostModel.estimate`), which count the units in production
      at the end by the share produced.

    A schedule that is itself a decision (such as the launch cadence) is a
    design variable, with the schedule made for each design (such as with
    `dataclasses.replace`).

    Args:
        years (float): The duration of the mission, from the first launch.
        design_life (float | Mapping | None): The design life of the
            satellites (of each type), or `None` if they last the mission (or
            fail, see `Deployment.simulate`).
        development_years (float | Mapping): The duration of development of
            each type of satellite, before production.
        production_rate (float | Mapping | None): The satellites of each type
            produced a year, or `None` for production that takes no time.
        launch_interval (float): The least time between launches.
        launch_capacity (int | None): The most satellites of a launch, if
            limited.
        launch_mass_capacity (float | Callable | None): The most mass (kg) of
            a launch, if limited (which requires the masses of the
            satellites): a value, or a function of the index and values of a
            plane (see `Plane`, as for `order`), such as for planes at
            different inclinations, to which a launcher lifts different
            masses.
        order (Callable | None): A function of the index and values of a
            plane (`None` for the planes of a Walker layout) that gives the
            key by which its type's planes are produced and launched, if
            not in the order of the layouts.
        sustained (bool): Whether the constellation is replaced after the
            end of the mission, whose costs are then those incurred by the
            end.
    """

    years: float = 1.0
    design_life: PerType | None = None
    development_years: PerType = 0.0
    production_rate: PerType | None = None
    launch_interval: float = 0.0
    launch_capacity: int | None = None
    launch_mass_capacity: float | Callable[[int, Any], float] | None = None
    order: Callable[[int, Any], Any] | None = None
    sustained: bool = False

    def __post_init__(self):
        if self.years <= 0:
            raise ValueError("years must be positive")
        if any(
            life is not None and life <= 0 for life in per_type_values(self.design_life)
        ):
            raise ValueError("design_life must be positive")
        if self.launch_interval < 0:
            raise ValueError("launch_interval must not be negative")
        if self.launch_capacity is not None and self.launch_capacity < 1:
            raise ValueError("launch_capacity must be at least 1")
        if any(
            rate is not None and rate <= 0
            for rate in per_type_values(self.production_rate)
        ):
            raise ValueError("production_rate must be positive")
        if any(years < 0 for years in per_type_values(self.development_years)):
            raise ValueError("development_years must not be negative")

    def design_life_of(self, kind: Any = None) -> float | None:
        """
        Gets the design life of a type of satellite.

        Args:
            kind (Any): The type of the satellites (see `plan`).

        Returns:
            float | None: the design life (years), or `None` if they last
                the mission
        """
        return per_type(self.design_life, kind, "design life")

    def production_years(self, kind: Any, satellites: int) -> float:
        """
        Gets the time to produce satellites of a type on its line.

        Args:
            kind (Any): The type of the satellites.
            satellites (int): The number of satellites.

        Returns:
            float: the years (0 for production that takes no time)
        """
        rate = per_type(self.production_rate, kind, "production rate")
        return 0.0 if rate is None else satellites / rate

    def mass_capacity_of(self, plane: tuple[int, Any] | None = None) -> float | None:
        """
        Gets the most mass (kg) of a launch to a plane, if limited.

        Args:
            plane (tuple[int, Any] | None): The index and values of the plane,
                for a capacity that is a function of the plane.

        Returns:
            float | None: the capacity, or `None` if not limited
        """
        capacity = self.launch_mass_capacity
        if capacity is None or not callable(capacity):
            return capacity
        if plane is None:
            raise ValueError("a launch capacity of each plane needs the plane")
        return float(capacity(*plane))

    def next_slot(self, booked: Sequence[float], ready: float) -> float:
        """
        Gets the earliest time of a launch at or after a time that keeps the
        launch interval from launches booked before (as plans and
        simulations book every launch: of deployments, replacements,
        restorations, and spares).

        Args:
            booked (Sequence[float]): The times of the launches booked, in
                order.
            ready (float): The time when the launch is ready.

        Returns:
            float: the time of the launch
        """
        interval = self.launch_interval
        time = ready
        if interval > 0:
            for other in booked:
                if other - interval < time < other + interval:
                    time = other + interval
        return time

    def per_launch(
        self, mass: float | None = None, plane: tuple[int, Any] | None = None
    ) -> float:
        """
        Gets the most satellites of a launch (to a plane), by the capacity of
        launches in number and in mass.

        Args:
            mass (float | None): The mass (kg) of each satellite.
            plane (tuple[int, Any] | None): The index and values of the plane,
                for a capacity in mass that is a function of the plane.

        Returns:
            float: the number of satellites (infinite if not limited)

        Raises:
            ValueError: if a capacity in mass is not given the mass, or a
                satellite exceeds it.
        """
        capacity = math.inf if self.launch_capacity is None else self.launch_capacity
        mass_capacity = self.mass_capacity_of(plane)
        if mass_capacity is not None:
            if mass is None:
                raise ValueError("a launch capacity in mass requires masses")
            fits = math.floor(mass_capacity / mass)
            if fits < 1:
                raise ValueError("a satellite exceeds the mass capacity of a launch")
            capacity = min(capacity, fits)
        return capacity

    def split(
        self,
        satellites: int,
        mass: float | None = None,
        plane: tuple[int, Any] | None = None,
    ) -> list[int]:
        """
        Gets the satellites of each launch of a plane: as few launches as
        carry them (see `per_launch`, at least one), sharing them evenly.

        Args:
            satellites (int): The number of satellites.
            mass (float | None): The mass (kg) of each satellite.
            plane (tuple[int, Any] | None): The index and values of the plane,
                for a capacity in mass that is a function of the plane.

        Returns:
            list[int]: the satellites of each launch
        """
        capacity = self.per_launch(mass, plane)
        count = max(1, math.ceil(satellites / capacity)) if capacity < math.inf else 1
        return [
            satellites // count + (1 if j < satellites % count else 0)
            for j in range(count)
        ]

    def launches_of(
        self,
        satellites: int,
        mass: float | None = None,
        plane: tuple[int, Any] | None = None,
    ) -> int:
        """
        Gets the launches of a plane of satellites (see `split`).

        Args:
            satellites (int): The number of satellites.
            mass (float | None): The mass (kg) of each satellite.
            plane (tuple[int, Any] | None): The index and values of the plane,
                for a capacity in mass that is a function of the plane.

        Returns:
            int: the number of launches
        """
        return len(self.split(satellites, mass, plane))

    def plan(
        self,
        layout: Any,
        *,
        kind: Any = None,
        policy: RestorationPolicy | None = None,
        masses: PerType | None = None,
        stages: Sequence[Stage] = (),
    ) -> Deployment:
        """
        Plans the deployment of a constellation: the layout (or layouts) of
        its satellites, of one type (`kind`), or a mapping of types to
        layouts (such as the shells of a mixed constellation), with the
        policy that restores the failures of its satellites (see
        `Deployment.simulate`): its spares in orbit in each plane are
        produced and launched with the plane, and wait in orbit to take the
        place of failed satellites, and its spares on the ground are
        produced after the planes (see `RestorationPolicy`). The planes of
        `stages` are deployed only if (and when) their rules decide to, in a
        simulation of the mission (see `Stage` and `Deployment.simulate`), and
        the others as planned.

        Args:
            layout (Any): The layout (or layouts) of the satellites, or a
                mapping of types of satellites to layouts.
            kind (Any): The type of the satellites of a layout.
            policy (RestorationPolicy | None): The policy for failed
                satellites, and its spares, or `None` to restore them by
                producing new ones (see `RestorationPolicy`, and `MissionPlan`,
                which plans the deployment of each design with its policy).
            masses (float | Mapping | None): The mass (kg) of the satellites
                (of each type), for a launch capacity in mass (and the price
                of launches by mass, see `CostModel.estimate`).
            stages (Sequence[Stage]): The stages deployed later, if their
                rules decide to.

        Returns:
            Deployment: the deployment

        Raises:
            ValueError: if a plane is in several stages or is not a plane of
                the layouts, if no satellites are planned, or if the spares of
                the policy lack a type of the layouts.
        """
        policy = RestorationPolicy() if policy is None else policy
        stages = tuple(stages)
        staged = [plane for stage in stages for plane in stage.planes]
        if len(set(staged)) != len(staged):
            raise ValueError("a plane is in several stages")
        # a mapping of types to layouts (the values of a layout are a mapping too)
        if isinstance(layout, Mapping) and not isinstance(layout, Values):
            fleets = dict(layout)
        else:
            fleets = {kind: layout}
        deployment = Deployment(schedule=self, stages=stages, policy=policy)
        produced = []  # the planes that are not staged, in their order
        for name, item in fleets.items():
            developed = float(per_type(self.development_years, name, "development"))
            spares = policy.orbit_spares_of(name)
            policy.ground_spares_of(name)  # (raises if it lacks the type)
            deployment.development[name] = developed
            deployment.masses[name] = (
                None if masses is None else per_type(masses, name, "mass")
            )
            listed = []
            for plane in layout_of(item).planes():
                index = len(deployment.planes)
                deployment.planes.append(
                    PlannedPlane(
                        name, plane.satellites, spares if plane.satellites > 0 else 0
                    )
                )
                deployment.plane_values.append(plane.values)
                if plane.satellites > 0 and index not in staged:
                    listed.append((index, plane))
            if self.order is not None:
                listed.sort(key=lambda item: self.order(item[0], item[1].values))
            produced.extend(index for index, _ in listed)
        unknown = sorted(set(staged) - set(range(len(deployment.planes))))
        if unknown:
            raise ValueError(
                f"stages of planes that the layouts do not have: {unknown}"
            )
        if not produced:
            raise ValueError("a deployment requires satellites (that are not staged)")
        # each type's planes produced on its line from the end of its development
        manifest = _Manifest(deployment, deployment.production)
        deployment.launches.extend(manifest.deploy(produced, 0.0))
        deployment.production.sort(key=lambda batch: batch.begin)
        deployment.launches.sort(key=lambda launch: (launch.time, launch.plane))
        return deployment
