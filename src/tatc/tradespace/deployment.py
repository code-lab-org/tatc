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
takes longer to produce and launch. A deployment simulates a mission
(`Deployment.simulate`) with the failures of its satellites and launches
(see `FailureModel`) and a policy that restores them or not (see
`RestorationPolicy`), drawn for each replicate: the `MissionHistory` of
what is produced and launched, which a `CostModel` prices, and of the
satellites that operate. Their changes divide the mission into phases
(`MissionHistory.phases`) of constant sets of operating satellites, each
evaluated over the study's representative period (`evaluate_phases`, which
states what these quasi-static phases assume) and combined over the mission
(`summarize_phases`), such as the time-weighted
mean of a metric and the time until it first meets a requirement.

A deployment may be staged (see `Stage`): planes that are produced and
launched only if a rule decides to, at times during the mission, from what
it knows then (see `Decision`), such as a second plane that is launched only
if the mission turns out to need it. Whether and how to stage a deployment
is a decision of a design, whose flexibility (the cost saved where a stage
is not needed, against the capability lost by deferring it where it is) the
analyses of a study compare with deployments that are not staged.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from ..schemas import Satellite
from .constellations import layout_of
from .events import Launch, PerType, Production, _per_type
from .failures import FailureModel, RestorationPolicy
from .histories import _simulate_many
from .histories import expected_metrics as _expected_metrics
from .mission import Decision, MissionHistory, Phase, _Simulation
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
    satellites never operate.

    A rule that depends on the values of a design (such as an uncertain
    parameter or a value of a scenario, known by the time of the decision)
    reads them from the design, as stages made for each design (see
    `MissionPlan`) can: for example, `Stage([1], 1.0, lambda decision:
    design.decorrelation_time < timedelta(hours=6))`.

    Args:
        planes (Sequence[int] | int): The planes of the stage (their indices,
            in the order of the layouts' planes), produced and launched in
            this order.
        times (float | Sequence[float]): The times of the decisions, in years
            after the first launch.
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
        planes = (self.planes,) if isinstance(self.planes, int) else self.planes
        planes = tuple(int(plane) for plane in planes)
        if not planes:
            raise ValueError("a stage needs planes")
        if len(set(planes)) != len(planes) or min(planes) < 0:
            raise ValueError("the planes of a stage must be distinct indices")
        times = self.times if isinstance(self.times, Sequence) else (self.times,)
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
    without failures.
    """

    schedule: DeploymentSchedule
    """DeploymentSchedule: The schedule of the deployment."""
    planes: list[tuple[Any, int, int]] = field(default_factory=list)
    """list: The planes, in the order of the layouts' planes: the type of
    their satellites, the number of satellites of the layout, and of
    spares."""
    development: dict[Any, float] = field(default_factory=dict)
    """dict: The end of the development of each type of satellite (years)."""
    production: list[Production] = field(default_factory=list)
    """list[Production]: The production of the satellites, in order."""
    launches: list[Launch] = field(default_factory=list)
    """list[Launch]: The launches (and planned replacements), in order."""
    masses: dict[Any, float | None] = field(default_factory=dict)
    """dict: The mass (kg) of each type of satellite, if given."""
    stages: tuple[Stage, ...] = ()
    """tuple[Stage, ...]: The stages deployed only if their rules decide to
    (whose planes the production and launches above do not include)."""
    _planned: MissionHistory | None = field(default=None, repr=False, compare=False)

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
        is launched (years)."""
        return max(
            launch.time for launch in self.launches if launch.purpose == "deployment"
        )

    @property
    def units(self) -> dict[Any, int]:
        """dict: The satellites produced of each type (including replacements
        and spares)."""
        units: dict[Any, int] = {}
        for batch in self.production:
            units[batch.kind] = units.get(batch.kind, 0) + batch.satellites
        return units

    def simulate(
        self,
        failures: FailureModel | None = None,
        policy: RestorationPolicy | None = None,
        seed: int | Sequence[int] | None = None,
    ) -> MissionHistory:
        """
        Simulates the mission: the deployment, with the failures of its
        satellites and launches (drawn from a seed of a replicate, such as an
        uncertain parameter, see `Uncertain.seed`) and the restoration of
        failed satellites by a policy. Without a seed, nothing fails (but a
        policy's spares on the ground are still produced).

        A satellite's life is drawn when it is launched (or takes the place
        of a failed one), from a random number generator seeded by the seed
        and its plane, position, and generation, and the success of a launch
        from one seeded by its plane and order: with common random numbers,
        the satellite at the same position of the same plane fails at the
        same age in the same replicate of every design. Spares in orbit do
        not age or fail until they take a failed satellite's place.

        A planned replacement (at the end of a design life) whose launch
        fails is not flown again: the satellite it would replace operates
        until it fails (and is then restored, by a policy that restores
        failures).

        Args:
            failures (FailureModel | None): The failures, or `None` for none.
            policy (RestorationPolicy | None): The policy for failed
                satellites, or `None` to restore them by producing new ones.
            seed (int | Sequence[int] | None): The seed of the replicate (or
                a sequence of integers, such as the seed and the index of a
                history, see `simulate_many`), or `None` for no failures.

        Returns:
            MissionHistory: the history of the mission
        """
        return _Simulation(
            self, failures or FailureModel(), policy or RestorationPolicy(), seed
        ).run()

    def planned(self) -> MissionHistory:
        """
        Gets the history of the mission without failures (see `simulate`),
        simulated once.

        Returns:
            MissionHistory: the history
        """
        if self._planned is None:
            self._planned = self.simulate()
        return self._planned

    def simulate_many(
        self,
        failures: FailureModel,
        policy: RestorationPolicy | None,
        seed: int | None,
        histories: int,
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
        which `summarize_missions` applies to every statistic. A stress
        that samples the outcomes of interest about as often as the others
        estimates their probabilities (and the quantiles of their tails)
        with far fewer histories; too large a stress makes a few weights
        dominate (see the effective number of histories).

        Many histories of the failures of a replicate, whose sample of
        targets (the expensive part of an evaluation) they share, measure the
        risk of a mission at little cost (see `evaluate_missions` and
        `PhaseCache`).

        Args:
            failures (FailureModel): The failures.
            policy (RestorationPolicy | None): The policy for failed
                satellites.
            seed (int | None): The seed of the replicate.
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
            return [self.simulate(failures, policy, None)]
        return _simulate_many(
            self, failures, policy, seed, histories, sampling, stress, launch_stress
        )

    def expected_metrics(
        self,
        failures: FailureModel,
        policy: RestorationPolicy | None = None,
        histories: int = 1000,
        seed: int = 0,
        sampling: str = "stratified",
    ) -> dict[str, float]:
        """
        Estimates the expected measures of the mission (see
        `MissionHistory.metrics`, such as its `availability`) from many
        simulated histories (see `expected_metrics`).

        Args:
            failures (FailureModel): The failures.
            policy (RestorationPolicy | None): The policy for failed
                satellites.
            histories (int): The number of histories.
            seed (int): The seed of the histories.
            sampling (str): The sampling of their failures.

        Returns:
            dict[str, float]: the mean of each measure
        """
        return _expected_metrics(self, failures, policy, histories, seed, sampling)

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
      order of the layouts' planes (or of `order`).
    * Each plane is launched when its satellites are ready, in as many
      launches as its satellites need (given a `launch_capacity` in number
      or mass, or one per plane), at most one every `launch_interval` years,
      in the order the planes are ready.
    * The mission lasts `years` from the first launch. Satellites are
      replaced at the end of their `design_life` (if before the end of the
      mission), each replacement produced at the production rate just
      before it is needed and launched then. A design life does not retire
      a satellite: one that is not replaced (as near the end of the mission,
      or if its replacement's launch fails) operates until it fails (by the
      wear-out of its life, see `FailureModel`).

    A schedule that is itself a decision (such as the launch cadence) is a
    design variable, with the schedule made for each design (such as with
    `dataclasses.replace`).

    Args:
        years (float): The duration of the mission, from the first launch.
        design_life (float | None): The design life of the satellites, or
            `None` if they last the mission (or fail, see
            `Deployment.simulate`).
        development_years (float | Mapping): The duration of development of
            each type of satellite, before production.
        production_rate (float | Mapping | None): The satellites of each type
            produced a year, or `None` for production that takes no time.
        launch_interval (float): The least time between launches.
        launch_capacity (int | None): The most satellites of a launch, if
            limited.
        launch_mass_capacity (float | None): The most mass (kg) of a launch,
            if limited (which requires the masses of the satellites).
        order (Callable | None): A function of the index and values of a
            plane (`None` for the planes of a Walker layout) that gives the
            key by which its type's planes are produced and launched, if
            not in the order of the layouts.
    """

    years: float = 1.0
    design_life: float | None = None
    development_years: PerType = 0.0
    production_rate: PerType | None = None
    launch_interval: float = 0.0
    launch_capacity: int | None = None
    launch_mass_capacity: float | None = None
    order: Callable[[int, Any], Any] | None = None

    def __post_init__(self):
        if self.years <= 0:
            raise ValueError("years must be positive")
        if self.design_life is not None and self.design_life <= 0:
            raise ValueError("design_life must be positive")
        if self.launch_interval < 0:
            raise ValueError("launch_interval must not be negative")
        if self.launch_capacity is not None and self.launch_capacity < 1:
            raise ValueError("launch_capacity must be at least 1")
        rates = self.production_rate
        values = rates.values() if isinstance(rates, Mapping) else [rates]
        if any(rate is not None and rate <= 0 for rate in values):
            raise ValueError("production_rate must be positive")
        years = self.development_years
        values = years.values() if isinstance(years, Mapping) else [years]
        if any(value < 0 for value in values):
            raise ValueError("development_years must not be negative")

    def launches_of(self, satellites: int, mass: float | None = None) -> int:
        """
        Gets the launches of a plane of satellites.

        Args:
            satellites (int): The number of satellites.
            mass (float | None): The mass (kg) of each satellite.

        Returns:
            int: the number of launches
        """
        count = 1
        if self.launch_capacity is not None:
            count = max(count, math.ceil(satellites / self.launch_capacity))
        if self.launch_mass_capacity is not None:
            if mass is None:
                raise ValueError("a launch capacity in mass requires masses")
            per_launch = math.floor(self.launch_mass_capacity / mass)
            if per_launch < 1:
                raise ValueError("a satellite exceeds the mass capacity of a launch")
            count = max(count, math.ceil(satellites / per_launch))
        return count

    def plan(
        self,
        layout: Any,
        kind: Any = None,
        spares: int = 0,
        masses: PerType | None = None,
        stages: Sequence[Stage] = (),
    ) -> Deployment:
        """
        Plans the deployment of a constellation: the layout (or layouts) of
        its satellites, of one type (`kind`), or a mapping of types to
        layouts (such as the shells of a mixed constellation). The planes of
        `stages` are deployed only if (and when) their rules decide to, in a
        simulation of the mission (see `Stage` and `Deployment.simulate`), and
        the others as planned.

        Args:
            layout (Any): The layout (or layouts) of the satellites, or a
                mapping of types of satellites to layouts.
            kind (Any): The type of the satellites of a layout.
            spares (int): The spares in orbit in each plane, produced and
                launched with it, which wait in orbit to take the place of
                failed satellites (see `RestorationPolicy`).
            masses (float | Mapping | None): The mass (kg) of the satellites
                (of each type), for a launch capacity in mass.
            stages (Sequence[Stage]): The stages deployed later, if their
                rules decide to.

        Returns:
            Deployment: the deployment

        Raises:
            ValueError: if a plane is in several stages or is not a plane of
                the layouts, or if no satellites are planned.
        """
        # pylint: disable=too-many-locals
        if spares < 0:
            raise ValueError("spares must not be negative")
        stages = tuple(stages)
        staged = [plane for stage in stages for plane in stage.planes]
        if len(set(staged)) != len(staged):
            raise ValueError("a plane is in several stages")
        # a mapping of types to layouts (the values of a layout are a mapping too)
        if isinstance(layout, Mapping) and not isinstance(layout, Values):
            fleets = dict(layout)
        else:
            fleets = {kind: layout}
        deployment = Deployment(schedule=self, stages=stages)
        planes = []  # (ready, index, kind, satellites, mass)
        for name, item in fleets.items():
            developed = float(_per_type(self.development_years, name, "development"))
            rate = _per_type(self.production_rate, name, "production rate")
            mass = None if masses is None else _per_type(masses, name, "mass")
            deployment.development[name] = developed
            deployment.masses[name] = mass
            listed = []
            for plane in layout_of(item).planes():
                size, values = plane.satellites, plane.values
                index = len(deployment.planes)
                deployment.planes.append((name, size, spares if size > 0 else 0))
                if size > 0 and index not in staged:
                    listed.append((index, size + spares, values))
            if self.order is not None:
                listed.sort(key=lambda plane: self.order(plane[0], plane[2]))
            line = developed
            for index, size, _ in listed:
                ready = line + (0.0 if rate is None else size / rate)
                deployment.production.append(Production(name, line, ready, size))
                planes.append((ready, index, name, size, mass))
                line = ready
        unknown = sorted(set(staged) - set(range(len(deployment.planes))))
        if unknown:
            raise ValueError(
                f"stages of planes that the layouts do not have: {unknown}"
            )
        if not planes:
            raise ValueError("a deployment requires satellites (that are not staged)")
        cadence = -math.inf
        for ready, index, name, size, mass in sorted(planes, key=lambda p: p[0]):
            count = self.launches_of(size, mass)
            first = max(ready, cadence)
            position = 0
            for j in range(count):
                carried = size // count + (1 if j < size % count else 0)
                deployment.launches.append(
                    Launch(
                        first + j * self.launch_interval, index, name, position, carried
                    )
                )
                position += carried
            cadence = first + count * self.launch_interval
        end = deployment.end
        if self.design_life is not None:
            for launch in list(deployment.launches):
                rate = _per_type(self.production_rate, launch.kind, "production rate")
                need, generation = launch.time + self.design_life, 1
                while need < end - 1e-9:
                    produced = need - (
                        0.0 if rate is None else launch.satellites / rate
                    )
                    deployment.production.append(
                        Production(
                            launch.kind,
                            produced,
                            need,
                            launch.satellites,
                            "replacement",
                        )
                    )
                    deployment.launches.append(
                        Launch(
                            need,
                            launch.plane,
                            launch.kind,
                            launch.first,
                            launch.satellites,
                            generation,
                            "replacement",
                        )
                    )
                    need, generation = need + self.design_life, generation + 1
        deployment.production.sort(key=lambda batch: batch.begin)
        deployment.launches.sort(key=lambda launch: (launch.time, launch.plane))
        return deployment
