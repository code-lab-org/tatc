"""
Studies of the observations of constellations: the targets of each replicate,
the observations of each satellite (cached), and the metrics of a design's
observations, its deployment and costs, the histories of its failures, and
the phases of its mission.

The parts refer to the values of designs by name (such as `"points_seed"`
for `d.points_seed`, or a path of `Design.flat`, such as
`"layout.walker.planes"`), or by a function of the (realized) design.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
import os
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timedelta
from functools import partial
from typing import Any

import pandas as pd

from ..analysis import collect_observations
from ..schemas import Satellite
from .cache import SatelliteCache
from .constellations import count_planes
from .cost import CostModel
from .criteria import Requirement
from .deployment import DeploymentSchedule
from .failures import FailureModel, RestorationPolicy
from .histories import evaluate_missions, expected_missions, summarize_missions
from .mission import Capability, PhaseCache, _given_design, summarize_phases
from .space import Design, DesignSpace, _scenario_of
from .study import Study, _takes_design
from .values import configure

_MAX_TARGETS = 32
"""The most sets of targets (of seeds) kept in memory."""


def value_of(design: Design, reference: Any) -> Any:
    """
    Gets a value of a design that a reference names: a string names a value
    (`d.<name>`, or a path of `Design.flat`, such as `"layout.walker.planes"`),
    a function computes it from the design, and `None` is `None`.

    Args:
        design (Design): The (realized) design.
        reference (Any): The reference: a string, a function, or `None`.

    Returns:
        Any: the value
    """
    if reference is None:
        return None
    if callable(reference):
        return reference(design)
    if not isinstance(reference, str):
        raise TypeError(f"a reference to a value of a design, not {reference!r}")
    if "." in reference or "[" in reference:
        values = design.flat()
        if reference not in values:
            raise KeyError(f"no value {reference!r} of the design")
        return values[reference]
    return getattr(design, reference)


@dataclasses.dataclass(frozen=True, eq=False, init=False)
class Targets:
    """
    The targets of observations (such as points): fixed (a data frame), or
    generated in each replicate by a function (such as
    `tatc.generation.generate_points_rotated_lattice`) given the `arguments`
    and `keywords`, and a seed that is a value of the design (as its keyword
    `seed`, such as of `Uncertain.seed`), so that each replicate samples
    other targets. Generated targets are kept in memory by seed.

    Args:
        source (pandas.DataFrame | Callable): The targets, or a function that
            generates them.
        *arguments (Any): The arguments of the function.
        seed (Any): The value of the design that seeds the function (see
            `value_of`), or `None` (to generate the same targets).
        **keywords (Any): The keyword arguments of the function.
    """

    source: Any
    arguments: tuple
    keywords: dict
    seed: Any

    def __init__(self, source: Any, *arguments: Any, seed: Any = None, **keywords: Any):
        if isinstance(source, pd.DataFrame):
            if arguments or keywords or seed is not None:
                raise ValueError("fixed targets take no arguments or seed")
        elif not callable(source):
            raise TypeError(
                "targets are a data frame or a function that generates them"
            )
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "arguments", tuple(arguments))
        object.__setattr__(self, "keywords", dict(keywords))
        object.__setattr__(self, "seed", seed)
        object.__setattr__(self, "_memory", {})

    def __getstate__(self) -> dict:
        # copies sent to parallel workers start with empty memories
        state = self.__dict__.copy()
        state["_memory"] = {}
        return state

    def seed_of(self, design: Design | None) -> Any:
        """
        Gets the seed of the targets of a design (`None` if they are fixed).

        Args:
            design (Design | None): The (realized) design.

        Returns:
            Any: the seed
        """
        if self.seed is None:
            return None
        if design is None:
            raise ValueError(
                "targets generated from a seed of the design need a design"
            )
        return value_of(design, self.seed)

    def of(self, design: Design | None = None) -> pd.DataFrame:
        """
        Gets the targets of a design (of its replicate).

        Args:
            design (Design | None): The (realized) design, or `None` for fixed
                targets.

        Returns:
            pandas.DataFrame: the targets
        """
        if isinstance(self.source, pd.DataFrame):
            return self.source
        seed = self.seed_of(design)
        memory = self._memory
        if seed not in memory:
            if len(memory) >= _MAX_TARGETS:
                memory.pop(next(iter(memory)))
            keywords = (
                self.keywords if self.seed is None else {**self.keywords, "seed": seed}
            )
            memory[seed] = self.source(*self.arguments, **keywords)
        return memory[seed]


def _observe(
    analysis: Callable,
    targets: Any,
    start: datetime,
    end: datetime,
    satellite: Satellite,
) -> pd.DataFrame:
    """The observations of the targets by a satellite (see `Observations`)."""
    return analysis(targets, satellite, start, end)


@dataclasses.dataclass(frozen=True, eq=False)
class Observations:
    """
    The observations of targets by satellites over a period, by an analysis
    of one satellite (such as `tatc.analysis.collect_observations`, called
    as `analysis(targets, satellite, start, end)`), cached for each satellite
    (see `SatelliteCache`) and each set of targets (of a seed, see
    `Targets`), in memory, or on disk in a `directory` (shared by the
    parallel workers of a study, and between sessions). The results on disk
    are keyed by the fingerprint of the analysis, which describes the
    targets by their content: the `namespace` need only identify data that
    the analysis reads (such as from files) and the code of other modules it
    calls (other than TAT-C's).

    The period is one window that represents the mission (see
    `evaluate_phases`): its start may be a value of the design (see
    `value_of`), such as an uncertain epoch drawn for each replicate over the
    mission (`Uncertain`), so that the replicates average over the geometry
    of the mission (such as the drift of the relative nodes of shells of
    different inclinations or altitudes) and its seasons. With common random
    numbers, the designs of a replicate share its epoch. The observations of
    each epoch (and set of targets) are cached apart.

    Args:
        targets (Targets | pandas.DataFrame): The targets.
        start (datetime | Any): The start of the period, or a value of the
            design that gives it (see `value_of`).
        duration (timedelta): The duration of the period.
        analysis (Callable): The analysis of one satellite.
        directory (str | os.PathLike | None): A directory of the cache on
            disk, or `None` to cache in memory only.
        namespace (str): An identifier of what the fingerprint of the
            analysis cannot see.
    """

    targets: Targets | pd.DataFrame
    start: datetime | Any
    duration: timedelta
    analysis: Callable = collect_observations
    directory: str | os.PathLike | None = dataclasses.field(
        default=None, metadata={"fingerprint": False}
    )
    namespace: str = ""

    def __post_init__(self):
        if isinstance(self.targets, pd.DataFrame):
            object.__setattr__(self, "targets", Targets(self.targets))
        object.__setattr__(self, "_caches", {})

    def __getstate__(self) -> dict:
        # copies sent to parallel workers start without caches
        state = self.__dict__.copy()
        state["_caches"] = {}
        return state

    @property
    def end(self) -> datetime:
        """datetime: The end of the period (of a fixed start, see `end_of`)."""
        if not isinstance(self.start, datetime):
            raise ValueError("the start is a value of the design: see end_of")
        return self.start + self.duration

    def start_of(self, design: Design | None = None) -> datetime:
        """
        Gets the start of the period of a design (of its replicate).

        Args:
            design (Design | None): The (realized) design, or `None` for a
                fixed start.

        Returns:
            datetime: the start
        """
        if isinstance(self.start, datetime):
            return self.start
        if design is None:
            raise ValueError("a start that is a value of the design needs a design")
        start = value_of(design, self.start)
        if not isinstance(start, datetime):
            raise TypeError(f"the start of the period is not a datetime: {start!r}")
        return start

    def end_of(self, design: Design | None = None) -> datetime:
        """
        Gets the end of the period of a design (of its replicate).

        Args:
            design (Design | None): The (realized) design, or `None` for a
                fixed start.

        Returns:
            datetime: the end
        """
        return self.start_of(design) + self.duration

    def context_of(self, design: Design | None = None) -> tuple:
        """
        Gets a key of the observations of a design: the seed of its targets
        and the period (such as the context of a `PhaseCache`).

        Args:
            design (Design | None): The (realized) design.

        Returns:
            tuple: the key
        """
        return (self.targets.seed_of(design), self.start_of(design), self.duration)

    def cache_of(self, design: Design | None = None) -> SatelliteCache:
        """
        Gets the cache of the observations of the targets of a design (by the
        seed of its targets and its period), in a namespace of the targets'
        seed and the period (and the `namespace`),
        so that a cache warns of results of other definitions (see
        `SatelliteCache`) of the same targets and period.

        Args:
            design (Design | None): The (realized) design.

        Returns:
            SatelliteCache: the cache
        """
        seed, start = self.targets.seed_of(design), self.start_of(design)
        if (seed, start) not in self._caches:
            if len(self._caches) >= _MAX_TARGETS:
                self._caches.pop(next(iter(self._caches)))
            function = partial(
                _observe,
                self.analysis,
                self.targets.of(design),
                start,
                start + self.duration,
            )
            context = (
                f"targets of seed {seed}, from {start.isoformat()} for {self.duration}"
            )
            namespace = f"{self.namespace}: {context}" if self.namespace else context
            self._caches[seed, start] = SatelliteCache(
                function, self.directory, namespace
            )
        return self._caches[seed, start]

    def __call__(
        self, satellites: Sequence[Satellite], design: Design | None = None
    ) -> pd.DataFrame:
        """
        Gets the observations of the targets of a design by satellites.

        Args:
            satellites (Sequence[Satellite]): The satellites.
            design (Design | None): The (realized) design.

        Returns:
            pandas.DataFrame: the observations (empty without satellites)
        """
        return self.cache_of(design)(satellites)


@dataclasses.dataclass(frozen=True)
class MissionPlan:
    """
    The deployment of the satellites of each design (see
    `DeploymentSchedule.plan`) and its costs (see `CostModel`), from the
    values of the design (see `value_of`).

    The schedule and the cost model may be functions of the (realized)
    design, which size them by its values with relationships of the study's
    own (TAT-C includes none, see the module `tatc.tradespace.cost`): for
    example, the costs of a type of satellite from the mass of its payload,
    or the capacity of launches to the design's altitude::

        def schedule(d):
            return dataclasses.replace(SCHEDULE, launch_mass_capacity=capacity(d.altitude))

        PLAN = MissionPlan(schedule, COST_MODEL, kind="payload")

    Launches limited in mass take the masses of the satellites from the cost
    model (`SatelliteCost.mass`, of each type), unless `masses` are given.

    Args:
        schedule (DeploymentSchedule | Callable): The schedule of the
            deployment, or a function of the design that gives it.
        cost_model (CostModel | Callable | None): The costs, a function of the
            design that gives them, or `None`.
        kind (Any): The kind of the satellites (for schedules and costs of
            several kinds), as a value of the design.
        spares (Any): The spares in orbit in each plane, or a value of the
            design (replaced by the `orbit_spares` of a design's policy, see
            `Failures`).
        growth (Any): The growth of the costs (see `CostModel.estimate`): a
            factor, a value of the design, or a mapping of each element (such
            as `"development"`) to a factor or a value of the design.
        stages (Sequence[Stage] | Callable): The stages of the deployment
            (see `Stage`), or a function of the (realized) design that gives
            them, such as of a decision whether to stage it and of the values
            its rules read. Their rules may also read the capability of the
            satellites operating at a decision, by the study's metrics (see
            `Decision.capability`).
        masses (Any): The mass (kg) of the satellites (or of each type, a
            mapping), or a value of the design, for launches limited in mass,
            or `None` for the masses of the cost model.
    """

    schedule: DeploymentSchedule | Callable[[Design], DeploymentSchedule]
    cost_model: CostModel | Callable[[Design], CostModel] | None = None
    kind: Any = None
    spares: Any = 0
    growth: Any = 1.0
    stages: Any = ()
    masses: Any = None

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

    def cost_model_of(self, design: Design) -> CostModel | None:
        """
        Gets the cost model of a design.

        Args:
            design (Design): The (realized) design.

        Returns:
            CostModel | None: the cost model, if any
        """
        if self.cost_model is None or isinstance(self.cost_model, CostModel):
            return self.cost_model
        return self.cost_model(design)

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
            Any: a factor, or a mapping of each element to its factor
        """
        if isinstance(self.growth, Mapping):
            return {element: _number(design, g) for element, g in self.growth.items()}
        return _number(design, self.growth)


def _number(design: Design, value: Any) -> Any:
    """A number, or a value of a design that a reference names."""
    return value if isinstance(value, (int, float)) else value_of(design, value)


@dataclasses.dataclass(frozen=True)
class Failures:
    """
    The failures of the satellites of each design: histories of its mission
    (see `Deployment.simulate_many`), restored by a policy (or by the policy
    of each design of a design space with a value `policy`, see
    `policy_space`, which configures the policy, see `configure`, and the
    spares in orbit), summarized over the histories (see
    `summarize_missions`): the mean of each metric over the phases of the
    mission and the histories, `quantiles` and `requirements` of them, and
    their costs. Without restoration, the means of the metrics and the
    availability are the exact expectations over the failures (see
    `expected_missions`) if `exact`.

    Args:
        model (FailureModel): The failures.
        policy (RestorationPolicy | None): The policy (`None` not to restore).
        histories (int): The histories of each replicate.
        seed (Any): The value of the design that seeds the failures (see
            `value_of`, such as of `Uncertain.seed`), or `None` for the one
            history without failures.
        sampling (str): The sampling of the histories (see
            `Deployment.simulate_many`).
        quantiles (Mapping | None): Quantiles of the metrics of the
            histories (see `summarize_missions`).
        requirements (Sequence[Requirement]): Requirements on the metrics of
            the histories (see `summarize_missions`).
        exact (bool): `True` for exact expectations without restoration.
    """

    model: FailureModel
    policy: RestorationPolicy | None = None
    histories: int = 1
    seed: Any = None
    sampling: str = "stratified"
    quantiles: Mapping[str, float | Sequence[float]] | None = None
    requirements: Sequence[Requirement] = ()
    exact: bool = True


@dataclasses.dataclass(frozen=True, eq=False)
class _Evaluation:
    """
    The evaluation of the satellites of a design (see `ObservationStudy`).
    """

    observations: Observations
    metrics: Callable[[pd.DataFrame], Mapping[str, Any]]
    layout: Any = None
    plan: MissionPlan | None = None
    failures: Failures | None = None
    phases: tuple[Requirement, ...] | None = None

    def __post_init__(self):
        object.__setattr__(self, "_phase_cache", PhaseCache(self.metrics))

    def _policy(self, design: Design) -> tuple[RestorationPolicy | None, Any]:
        """The policy of a design, and its spares in orbit in each plane."""
        base = None if self.failures is None else self.failures.policy
        spares = _number(design, self.plan.spares)
        decided = design.values.get("policy")
        if decided is None:
            return base, spares
        return configure(base or RestorationPolicy(), decided), decided.orbit_spares

    def __call__(self, satellites: Sequence[Satellite], design: Design) -> dict:
        observations = self.observations(satellites, design)
        record: dict[str, Any] = {"satellites": len(satellites)}
        layout = value_of(design, self.layout)
        if layout is not None:
            record["planes"] = count_planes(layout)
        # the metrics of phases are cached by the observations and, if they
        # are given the design, by its scenario
        context = self.observations.context_of(design)
        if _takes_design(self.metrics):
            context = (context, _scenario_of(design))
        if self.plan is not None:
            policy, spares = self._policy(design)
            kind = value_of(design, self.plan.kind)
            cost_model = self.plan.cost_model_of(design)
            deployment = self.plan.schedule_of(design).plan(
                layout,
                kind=kind,
                spares=spares,
                masses=self.plan.masses_of(design, cost_model),
                stages=self.plan.stages_of(design),
            )
            # the capability that the rules of stages read at their decisions
            deployment.capability = Capability(
                satellites, observations, self._phase_cache, context, design
            )
            growth = self.plan.growth_of(design)
            if cost_model is not None:
                record.update(cost_model.estimate(deployment, growth).metrics())
        metrics = _given_design(self.metrics, design)
        record.update(metrics(observations))
        if self.failures is not None:
            histories = deployment.simulate_many(
                self.failures.model,
                policy,
                value_of(design, self.failures.seed),
                self.failures.histories,
                sampling=self.failures.sampling,
            )
            missions = evaluate_missions(
                histories,
                satellites,
                observations,
                self._phase_cache,
                context=context,
                cost_model=cost_model,
                growth=growth,
                design=design,
            )
            record.update(
                summarize_missions(
                    missions,
                    quantiles=self.failures.quantiles,
                    requirements=self.failures.requirements,
                )
            )
            # exact expectations are of deployments that are not staged
            if (
                self.failures.exact
                and policy is not None
                and not policy.restore
                and not deployment.stages
            ):
                exact = expected_missions(
                    deployment,
                    self.failures.model,
                    satellites,
                    observations,
                    self._phase_cache,
                    context=context,
                    design=design,
                )
                record.update(
                    {
                        key: value
                        for key, value in exact.items()
                        if key.endswith(" mean") or key == "availability"
                    }
                )
        elif self.phases is not None:
            capability = self._phase_cache.evaluate(
                deployment.phases(satellites), observations, satellites, context, design
            )
            record.update(summarize_phases(capability, requirements=self.phases))
        return record


class ObservationStudy(Study):
    """
    A trade study of the observations of constellations: each replicate of a
    design builds its satellites, observes the targets of the replicate (see
    `Observations`), and records the number of `satellites` (and of `planes`,
    with a `layout`), the `metrics` of the observations, and, with a mission
    `plan`, the costs of its deployment, and either the summaries of the
    histories of its `failures` (see `Failures`) or of the `phases` of its
    deployment (see `evaluate_phases` and `summarize_phases`, with
    requirements, if any).

    The metrics are a function of the observations of any set of satellites
    (which may be empty), so that failures and phases evaluate them for the
    satellites that operate in each phase, over the one window of the
    observations (see `evaluate_phases` for what this assumes, and
    `Observations` for an epoch drawn for each replicate). A metrics function
    with a parameter named `design` is also given the realized design, such
    as for the window of an epoch drawn for each replicate or for values that
    scenarios set (see `Study`): the metrics of phases are cached by the
    satellites that operate, the targets and window of the observations, and
    the scenario, so such a function must depend on the design only through
    these. As a `Study`, it runs as any other (see `Study.run`), such as in
    parallel, with a store of its records, and in searches (such as
    `refine`, `confirm`, `screen`, and the optimizers). Its fingerprint (see
    `Study.fingerprint`) describes its parts, including targets by their
    content (or the recipe of generated targets), so that its store needs no
    namespace but for what the fingerprint cannot see.

    Args:
        space (DesignSpace): The design space.
        satellites (Callable): A function of a realized design that builds
            its satellites (such as with `build_constellation`).
        observations (Observations): The observations of the targets.
        metrics (Callable): A function of observations (and of the realized
            design, if it has a parameter named `design`) that gives their
            metrics (a mapping).
        layout (Any): The value of the design that is its layout (see
            `value_of` and `layout_of`), if any, which a `plan` requires.
        plan (MissionPlan | None): The deployment and costs, if any.
        failures (Failures | None): The failures, if any, which require a
            plan.
        phases (Sequence[Requirement] | bool | None): `True` (or requirements
            on the metrics) to summarize the phases of the deployment, which
            require a plan.
        **options (Any): The settings of the study (see `Study`), such as
            `replicates`, `seed`, `n_jobs`, and `store`.
    """

    def __init__(
        self,
        space: DesignSpace,
        satellites: Callable[[Design], Sequence[Satellite]],
        observations: Observations,
        metrics: Callable[[pd.DataFrame], Mapping[str, Any]],
        layout: Any = None,
        plan: MissionPlan | None = None,
        failures: Failures | None = None,
        phases: Sequence[Requirement] | bool | None = None,
        **options: Any,
    ):
        # the requirements on the phases, if they are summarized
        phases = None if phases is None or phases is False else phases
        phases = () if phases is True else None if phases is None else tuple(phases)
        if plan is None and (failures is not None or phases is not None):
            raise ValueError("failures and phases require a mission plan")
        if plan is not None and layout is None:
            raise ValueError("a mission plan requires a layout")
        if failures is not None and phases is not None:
            raise ValueError(
                "failures summarize the phases of their histories: give one or the other"
            )
        self.observations, self.metrics, self.layout = observations, metrics, layout
        self.plan, self.failures, self.phases = plan, failures, phases
        evaluation = _Evaluation(observations, metrics, layout, plan, failures, phases)
        super().__init__(space, satellites, evaluation, **options)

    def satellites_of(self, design: Design, replicate: int = 0) -> list[Satellite]:
        """
        Gets the satellites of a replicate of a design.

        Args:
            design (Design): The design.
            replicate (int): The replicate.

        Returns:
            list[Satellite]: the satellites
        """
        return list(self.build(self._realize(design, replicate)))

    def targets_of(self, design: Design, replicate: int = 0) -> pd.DataFrame:
        """
        Gets the targets of a replicate of a design.

        Args:
            design (Design): The design.
            replicate (int): The replicate.

        Returns:
            pandas.DataFrame: the targets
        """
        return self.observations.targets.of(self._realize(design, replicate))

    def observations_of(self, design: Design, replicate: int = 0) -> pd.DataFrame:
        """
        Gets the observations of a replicate of a design (such as for analyses
        other than its metrics).

        Args:
            design (Design): The design.
            replicate (int): The replicate.

        Returns:
            pandas.DataFrame: the observations
        """
        realized = self._realize(design, replicate)
        return self.observations(self.build(realized), realized)

    def _realize(self, design: Design, replicate: int) -> Design:
        return design.realize(replicate, self.seed, self.common_random_numbers)
