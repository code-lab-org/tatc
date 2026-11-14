"""
Studies of the observations of constellations: the targets of each replicate,
the observations of each satellite (cached), and the metrics of a design's
observations and, with a mission plan (see the module `plans`), its
deployment and costs, the histories of its failures, and the phases of its
mission.

The parts refer to the values of designs by a `Reference` (see `value_of`):
a name (such as `"points_seed"` for `d.points_seed`, or a path of
`Design.flat`, such as `"layout.walker.planes"`), or a function of the
(realized) design (such as `lambda d: d.layout`). A string given where a
reference is expected is always a reference, never a value itself: give a
function that returns a string (or a number, where numbers are references)
to give a value::

    MissionPlan(SCHEDULE, COST_MODEL, kind="field_of_regard")  # d.field_of_regard
    MissionPlan(SCHEDULE, COST_MODEL, kind=lambda d: "imager")  # always "imager"

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
import json
import os
import warnings
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timedelta
from functools import partial
from pathlib import Path
from typing import Any

import pandas as pd

from ..analysis import collect_observations
from ..schemas import Satellite
from .cache import SatelliteCache, _write
from .constellations import count_planes
from .criteria import Requirement
from .mission import Capability, PhaseCache, _given_design
from .plans import Failures, Mission, MissionPlan
from .references import Reference, value_of
from .space import Design, DesignSpace
from .store import _decode, _encode
from .study import Study, _Case, _takes_design
from .walks import _generic

_MAX_TARGETS = 32
"""The most sets of targets (of seeds) kept in memory."""


@dataclasses.dataclass(frozen=True, eq=False, init=False)
class Targets:
    """
    The targets of observations (such as points): fixed (a data frame), or
    generated in each replicate by a function (such as
    `tatc.generation.generate_points_fibonacci_lattice`) given the `arguments`
    and `keywords`, and a seed that is a value of the design (as its keyword
    `seed`, such as of `Uncertain.seed`), so that each replicate samples
    other targets. Generated targets are kept in memory by seed.

    Args:
        source (pandas.DataFrame | Callable): The targets, or a function that
            generates them.
        *arguments (Any): The arguments of the function.
        seed (Reference | None): The value of the design that seeds the
            function (see `value_of`), or `None` (to generate the same
            targets).
        **keywords (Any): The keyword arguments of the function.
    """

    source: Any
    arguments: tuple
    keywords: dict
    seed: Reference | None

    def __init__(
        self,
        source: Any,
        *arguments: Any,
        seed: Reference | None = None,
        **keywords: Any,
    ):
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
        return self._of_seed(self.seed_of(design))

    def _of_seed(self, seed: Any) -> pd.DataFrame:
        """Gets the targets of a seed (generated, or kept in memory)."""
        if isinstance(self.source, pd.DataFrame):
            return self.source
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
        start (datetime | Reference): The start of the period, or a value of
            the design that gives it (see `value_of`).
        duration (timedelta): The duration of the period.
        analysis (Callable): The analysis of one satellite.
        directory (str | os.PathLike | None): A directory of the cache on
            disk, or `None` to cache in memory only.
        namespace (str): An identifier of what the fingerprint of the
            analysis cannot see.

    Results of the analysis with other code (such as of an earlier version of
    TAT-C) are not reused, and a study warns of them (once, see
    `SatelliteCache`); `prune` removes them.
    """

    targets: Targets | pd.DataFrame
    start: datetime | Reference
    duration: timedelta
    analysis: Callable = collect_observations
    _: dataclasses.KW_ONLY
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
        seed and the period (and the `namespace`), so that the cache warns of
        results of the analysis with other code (see `SatelliteCache`) of the
        same targets and period.

        Args:
            design (Design | None): The (realized) design.

        Returns:
            SatelliteCache: the cache
        """
        return self._cache(self.targets.seed_of(design), self.start_of(design))

    def _namespace(self, seed: Any, start: datetime) -> str:
        """Gets the namespace of the cache of the targets of a seed and a start."""
        context = (
            f"targets of seed {seed}, from {start.isoformat()} for {self.duration}"
        )
        return f"{self.namespace}: {context}" if self.namespace else context

    def _cache(self, seed: Any, start: datetime) -> SatelliteCache:
        """
        Gets the cache of the observations of the targets of a seed from a
        start (see `cache_of`), with the seed and start in its namespace's
        folder (`context.json`), from which `prune` finds it.
        """
        if (seed, start) not in self._caches:
            if len(self._caches) >= _MAX_TARGETS:
                self._caches.pop(next(iter(self._caches)))
            function = partial(
                _observe,
                self.analysis,
                self.targets._of_seed(seed),  # pylint: disable=protected-access
                start,
                start + self.duration,
            )
            cache = SatelliteCache(
                function, self.directory, self._namespace(seed, start)
            )
            if self.directory is not None:
                context = json.dumps(_encode([seed, start])).encode()
                _write(cache.folder.parent / "context.json", context)
            self._caches[seed, start] = cache
        return self._caches[seed, start]

    def prune(self) -> int:
        """
        Removes the results on disk of the analysis with other code (see
        `SatelliteCache.prune`) from the caches of every set of targets (of
        each seed and start) of these observations in the directory, such as
        after a study warns of them. The results of other observations that
        share the directory (such as of other targets) are kept.

        Returns:
            int: the number of results removed (0 without a directory)
        """
        if self.directory is None or not Path(self.directory).is_dir():
            return 0
        removed = 0
        for folder in Path(self.directory).iterdir():
            try:
                seed, start = _decode(json.loads((folder / "context.json").read_text()))
                namespace = (folder / "namespace.txt").read_text(encoding="utf-8")
            except (OSError, ValueError, TypeError):
                continue
            if namespace != self._namespace(seed, start):
                continue  # of other observations
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                removed += self._cache(seed, start).prune()
        return removed

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


def _hashable(value: Any) -> Any:
    """A value, or its representation if it is not hashable."""
    try:
        hash(value)
    except TypeError:
        return repr(value)
    return value


def _uncertain_values(design: Design) -> tuple:
    """
    Gets the values of the uncertain parameters of a realized design, by
    path (see `DesignSpace.uncertain_paths`), as a key.
    """
    paths = set(design.space.uncertain_paths())
    return tuple(
        (path, _hashable(value))
        for path, value in design.flat().items()
        if _generic(path) in paths
    )


@dataclasses.dataclass(frozen=True, eq=False)
class _Evaluation:
    """
    The evaluation of the satellites of a design (see `ObservationStudy`):
    the metrics of the observations of every satellite and, with a plan, of
    the mission (see `Mission`): its deployment and costs, and either the
    phases of its deployment or the histories of its failures.
    """

    observations: Observations
    metrics: Callable[[pd.DataFrame], Mapping[str, Any]]
    layout: Reference | None = None
    plan: MissionPlan | None = None
    failures: Failures | None = None
    phases: tuple[Requirement, ...] | None = None
    _phase_cache: PhaseCache = dataclasses.field(
        init=False, repr=False, metadata={"fingerprint": False}
    )

    def __post_init__(self):
        object.__setattr__(self, "_phase_cache", PhaseCache(self.metrics))

    def __call__(self, satellites: Sequence[Satellite], design: Design) -> dict:
        observations = self.observations(satellites, design)
        record: dict[str, Any] = {"satellites": len(satellites)}
        layout = value_of(design, self.layout)
        if layout is not None:
            record["planes"] = count_planes(layout)
        mission = None
        if self.plan is not None:
            mission = self.mission_of(satellites, observations, design)
            record.update(mission.metrics())
        record.update(_given_design(self.metrics, design)(observations))
        if mission is not None:
            record.update(mission.summarize(self.phases))
        return record

    def context_of(self, design: Design) -> tuple:
        """
        Gets the context of the metrics of the phases of a design in the
        cache (see `PhaseCache`): its observations (see
        `Observations.context_of`) and, if the metrics are given the design,
        its scenario and the values of its uncertain parameters, which such
        metrics may read (designs of a replicate share them, with common
        random numbers).
        """
        context = self.observations.context_of(design)
        if _takes_design(self.metrics):
            context = (context, design.scenario, _uncertain_values(design))
        return context

    def mission_of(
        self,
        satellites: Sequence[Satellite],
        observations: pd.DataFrame,
        design: Design,
    ) -> Mission:
        """
        Plans the mission of a realized design (see `MissionPlan.mission_of`),
        with the capability of its satellites (and the cache of the metrics
        of their phases).
        """
        capability = Capability(
            satellites,
            observations,
            self._phase_cache,
            context=self.context_of(design),
            design=design,
        )
        return self.plan.mission_of(
            design, value_of(design, self.layout), capability, self.failures
        )


class ObservationStudy(Study):
    """
    A trade study of the observations of constellations: each replicate of a
    design builds its satellites, observes the targets of the replicate (see
    `Observations`), and records the metrics of the observations and, with a
    mission `plan`, of the mission (see `Mission`): its deployment and costs,
    and either the `phases` of its deployment (see `evaluate_phases` and
    `summarize_phases`, with requirements, if any) or the histories of its
    `failures` (see `Failures`)::

        study = ObservationStudy(
            space,
            build,  # the satellites of a design
            Observations(Targets(generate_points_fibonacci_lattice, count=100, rotate=True, seed="seed"), START, DAY),
            sample,  # the metrics of observations
            plan=MissionPlan(SCHEDULE, COST_MODEL),
            phases=True,
        )

    Its results have the columns of every study (see `Study.run`) and:

    * `satellites`, and `planes` (with a `layout`, by default the value
      `layout` of a design space that has one): their numbers.
    * the metrics of the observations of every satellite (by their names).
    * with a plan: `deployed`, the years from the start of development to the
      last launch of the first generation of the satellites in the history
      of the mission without failures (see `Mission.planned`), including the
      stages that their rules deploy then (see `Stage`); with stages, `stages
      deployed`, the number of stages that their rules deploy then (with
      failures, the mean over the histories); and, with a cost model, the
      costs of the deployment (`cost`, `cost development`, and so on, see
      `CostEstimate.metrics`; with failures, the means of the costs of the
      histories).
    * with phases: `<metric> mission`, the mean of each metric over the
      phases of the mission as its satellites are deployed (from the first
      launch), and `<metric> time`, the time at which a phase first meets
      each requirement on a metric (see `summarize_phases`), in years from
      the start of development, as `deployed` (infinite if none does).
    * with failures: the means over the histories of the failures of each
      replicate of `<metric> mission` (the mean over the phases of each
      history), of the measures of the missions (`availability`, `failures`,
      `restorations`, and so on, see `MissionHistory.metrics`), of their
      costs, and of the times of each history that meets the requirements of
      the `phases` of the failures (`<metric> time`); and the `quantiles` of
      columns over the histories (`<column> histories p<percent>`) and the
      shares of the histories that meet `requirements` (`P(<requirement>)`),
      see `summarize_histories`.

    Statistics over the replicates of these columns (see `Criteria`) judge
    designs: such as the mean over the replicates of `accurate mission`, or a
    quantile over them of `cost histories p70`.

    The metrics are a function of the observations of any set of satellites
    (which may be empty), so that failures and phases evaluate them for the
    satellites that operate in each phase, over the one window of the
    observations (see `evaluate_phases` for what this assumes, and
    `Observations` for an epoch drawn for each replicate). A metrics function
    with a parameter named `design` is also given the realized design, such
    as for the window of an epoch drawn for each replicate or for values that
    scenarios set (see `Study`): the metrics of phases are cached by the
    satellites that operate, the targets and window of the observations, the
    scenario, and the values of the uncertain parameters of the replicate
    (which the designs of a replicate share, with common random numbers), so
    such a function must depend on the design only through these. As a
    `Study`, it runs as any other (see `Study.run`), such as in parallel,
    with a store of its records, and in searches (such as `refine`,
    `confirm`, `screen`, and the optimizers). Its fingerprint (see
    `Study.fingerprint`) describes its parts, including targets by their
    content (or the recipe of generated targets), so that its store needs no
    namespace but for what the fingerprint cannot see.

    The parts of the evaluation of a design are also given apart, for
    analyses other than its metrics: its satellites, targets, and
    observations (`satellites_of`, `targets_of`, and `observations_of`), and
    its mission (`mission_of`): its planned deployment, the metrics of its
    phases, its cost estimate, and the histories of its failures.

    Args:
        space (DesignSpace): The design space.
        build (Callable): A function of a realized design that builds its
            satellites (such as with `build_constellation`).
        observations (Observations): The observations of the targets.
        metrics (Callable): A function of observations (and of the realized
            design, if it has a parameter named `design`) that gives their
            metrics (a mapping).
        layout (Reference | None): The value of the design that is its
            layout (see `value_of` and `layout_of`), which a `plan`
            requires, or `None` for the value `layout` of a design space
            that has one (and none otherwise).
        plan (MissionPlan | None): The deployment, policy, and costs, if any.
        failures (Failures | None): The failures, if any, which require a
            plan.
        phases (Sequence[Requirement] | bool | None): `True` (or requirements
            on the metrics) to summarize the phases of the deployment, which
            require a plan (with failures, the requirements on the phases of
            the histories are those of `Failures`).
        **options (Any): The settings of the study (see `Study`), such as
            `replicates`, `seed`, `n_jobs`, and `store`.

    Raises:
        ValueError: if failures or phases are given without a plan, a plan
            without a layout, failures with phases, or the quantiles or
            requirements of failures are of costs without a cost model.
    """

    def __init__(
        self,
        space: DesignSpace,
        build: Callable[[Design], Sequence[Satellite]],
        observations: Observations,
        metrics: Callable[[pd.DataFrame], Mapping[str, Any]],
        *,
        layout: Reference | None = None,
        plan: MissionPlan | None = None,
        failures: Failures | None = None,
        phases: Sequence[Requirement] | bool | None = None,
        **options: Any,
    ):
        # the requirements on the phases, if they are summarized
        phases = None if phases is None or phases is False else phases
        phases = () if phases is True else None if phases is None else tuple(phases)
        if layout is None and "layout" in space.children:
            layout = "layout"
        if plan is None and (failures is not None or phases is not None):
            raise ValueError("failures and phases require a mission plan")
        if plan is not None and layout is None:
            raise ValueError("a mission plan requires a layout")
        if failures is not None and phases is not None:
            raise ValueError(
                "failures summarize the phases of their histories: give the "
                + "requirements on them to Failures (phases), not to the study"
            )
        if failures is not None and plan.cost_model is None:
            named = [
                *(failures.quantiles or {}),
                *(r.name for r in failures.requirements),
            ]
            costs = [
                name for name in named if name == "cost" or name.startswith("cost ")
            ]
            if costs:
                raise ValueError(
                    f"quantiles or requirements of costs without a cost model: {costs}"
                )
        self.observations, self.metrics, self.layout = observations, metrics, layout
        self.plan, self.failures, self.phases = plan, failures, phases
        evaluation = _Evaluation(observations, metrics, layout, plan, failures, phases)
        super().__init__(space, build, evaluation, **options)

    def _arguments(self) -> dict[str, Any]:
        return {
            "space": self.space,
            "build": self.build,
            "observations": self.observations,
            "metrics": self.metrics,
            "layout": self.layout,
            "plan": self.plan,
            "failures": self.failures,
            "phases": self.phases,
            **self._settings(),
        }

    def _prepare(
        self,
        design: Design,
        replicate: int,
        shift: Mapping[str, int] | None,
        case: _Case,
    ) -> None:
        # the cache of the observations of the first evaluation, in this
        # process, so that it warns of results of other code here (once, see
        # `SatelliteCache`), rather than in a worker process
        realized = self._realized(design, replicate, shift, case)
        self.observations.cache_of(realized)

    def satellites_of(
        self, design: Design, replicate: int = 0, scenario: str | None = None
    ) -> list[Satellite]:
        """
        Gets the satellites of a replicate of a design (see `Study.realize`).

        Args:
            design (Design): The design.
            replicate (int): The replicate.
            scenario (str | None): A scenario of the study (by name), if any.

        Returns:
            list[Satellite]: the satellites
        """
        return list(self.build(self.realize(design, replicate, scenario)))

    def targets_of(
        self, design: Design, replicate: int = 0, scenario: str | None = None
    ) -> pd.DataFrame:
        """
        Gets the targets of a replicate of a design (see `Study.realize`).

        Args:
            design (Design): The design.
            replicate (int): The replicate.
            scenario (str | None): A scenario of the study (by name), if any.

        Returns:
            pandas.DataFrame: the targets
        """
        return self.observations.targets.of(self.realize(design, replicate, scenario))

    def observations_of(
        self, design: Design, replicate: int = 0, scenario: str | None = None
    ) -> pd.DataFrame:
        """
        Gets the observations of a replicate of a design (see `Study.realize`),
        such as for analyses other than its metrics.

        Args:
            design (Design): The design.
            replicate (int): The replicate.
            scenario (str | None): A scenario of the study (by name), if any.

        Returns:
            pandas.DataFrame: the observations
        """
        realized = self.realize(design, replicate, scenario)
        return self.observations(self.build(realized), realized)

    def mission_of(
        self, design: Design, replicate: int = 0, scenario: str | None = None
    ) -> Mission:
        """
        Gets the mission of a replicate of a design, as its evaluation plans
        it (see `MissionPlan.mission_of`), for analyses other than its
        metrics: its planned deployment (`Mission.deployment`, with the
        plan's schedule, type, policy and spares, masses, and stages, and the
        capability of its satellites, which the rules of stages read), the
        metrics of its phases (`Mission.phases`, which a study with `phases`
        summarizes), its cost estimate (`Mission.estimate`, the costs of a
        study without failures), and the histories of its failures
        (`Mission.histories`, a row per history, which a study with failures
        summarizes). The design is realized, built, and observed once for
        all of them.

        Args:
            design (Design): The design.
            replicate (int): The replicate.
            scenario (str | None): A scenario of the study (by name), if any.

        Returns:
            Mission: the mission

        Raises:
            ValueError: if the study has no mission plan.
        """
        if self.plan is None:
            raise ValueError("the study has no mission plan")
        realized = self.realize(design, replicate, scenario)
        satellites = list(self.build(realized))
        observations = self.observations(satellites, realized)
        return self.evaluate.mission_of(satellites, observations, realized)

    def prune(self) -> int:
        """
        Removes the evaluations of the study's store by this study with other
        code (see `Study.prune`), and the results of its observations' caches
        with other code (see `Observations.prune`, which returns their number
        when called alone).

        Returns:
            int: the number of evaluations removed (as `Study.prune`, not
                counting the results of the observations' caches)
        """
        evaluations = super().prune()
        self.observations.prune()
        return evaluations
