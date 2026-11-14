"""
Simulated missions: the history of a deployment with failures and
restorations (see `Deployment.simulate`, whose simulation is in the module
`simulation`), the decisions of the rules of its stages, the phases of a
mission in which constant sets of satellites operate, and their evaluation
over the study's representative period and combination over the mission.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import cached_property, partial
from typing import TYPE_CHECKING, Any, NamedTuple

import numpy as np
import pandas as pd

from ..schemas import Satellite
from .cache import get_satellite_key
from .criteria import Requirement
from .events import Launch, Production
from .study import _takes_design

if TYPE_CHECKING:
    from .deployment import Deployment, DeploymentSchedule


@dataclass(frozen=True)
class Phase:
    """
    A phase of a mission, in which a set of satellites operates (see
    `MissionHistory.phases`).
    """

    start: float
    """float: The start (years from the start of development)."""
    end: float
    """float: The end (years from the start of development)."""
    satellites: tuple[str, ...]
    """tuple[str, ...]: The names of the satellites that operate."""

    @property
    def duration(self) -> float:
        """float: The duration (years)."""
        return self.end - self.start


PHASE_COLUMNS = ("start", "end", "duration", "satellites")
"""The columns of the evaluations of phases that are not metrics (see
`evaluate_phases`)."""


def _names(satellites: Sequence[Satellite | str]) -> list[str]:
    return [s.name if isinstance(s, Satellite) else str(s) for s in satellites]


def _position_names(
    planes: Sequence[tuple[Any, int, int]], satellites: Sequence[Satellite | str]
) -> dict[tuple[int, int], str]:
    """
    Gets the name of the satellite at each position of each plane (by plane
    and position) of satellites in the order of the planes, without spares
    (see `MissionHistory.phases`), raising a `ValueError` if their numbers
    differ.
    """
    names = _names(satellites)
    positions = [
        (plane, position)
        for plane, (_, built, _) in enumerate(planes)
        for position in range(built)
    ]
    if len(names) != len(positions):
        raise ValueError(
            f"expected {len(positions)} satellites for the planes, not {len(names)}"
        )
    return dict(zip(positions, names))


def _satellite_keys(
    satellites: Sequence[Satellite | str] | Mapping[str, str],
) -> dict[str, str]:
    """
    Gets the keys of satellites in a cache of their metrics (see
    `PhaseCache`), by name: of their content (see `get_satellite_key`), or
    their names, if given as names (or the keys given, by name).
    """
    if isinstance(satellites, Mapping):
        return dict(satellites)
    return {
        name: get_satellite_key(s) if isinstance(s, Satellite) else name
        for name, s in zip(_names(satellites), satellites)
    }


class Failure(NamedTuple):
    """A failure of a satellite (see `MissionHistory.failures`)."""

    time: float
    """float: The time of the failure (years from the start of development)."""
    plane: int
    """int: The plane (its index)."""
    position: int
    """int: The position in the plane."""
    cause: str
    """str: `"life"`, or `"launch"` (lost in a launch)."""


class Restoration(NamedTuple):
    """A restoration of a failed satellite (see `MissionHistory.restorations`)."""

    plane: int
    """int: The plane (its index)."""
    position: int
    """int: The position in the plane."""
    failed: float
    """float: The time of the failure (years from the start of development)."""
    restored: float
    """float: The time of the restoration."""
    source: str
    """str: `"orbit"` (a spare in orbit), `"ground"` (a spare on the ground),
    or `"production"` (a satellite produced when needed)."""


@dataclass
class MissionHistory:
    """
    The history of a simulated mission (see `Deployment.simulate`): what is
    produced and launched (as a deployment, which a `CostModel` prices), when
    satellites fail and are restored, and when each operates.
    """

    deployment: Deployment
    """Deployment: The planned deployment."""
    production: list[Production] = field(default_factory=list)
    """list[Production]: The production of satellites, in order."""
    launches: list[Launch] = field(default_factory=list)
    """list[Launch]: The launches, in order."""
    operating: dict[tuple[int, int], list[tuple[float, float]]] = field(
        default_factory=dict
    )
    """dict: The intervals in which the satellite at each position of each
    plane (by plane and position) operates."""
    parked: list[tuple[float, float, int]] = field(default_factory=list)
    """list: The intervals in which spares wait in orbit, and their number."""
    failures: list[Failure] = field(default_factory=list)
    """list[Failure]: The failures: their time, plane, position, and cause
    (`"life"` or `"launch"`)."""
    restorations: list[Restoration] = field(default_factory=list)
    """list[Restoration]: The restorations: their plane, position, the time
    of the failure and of the restoration, and the source (`"orbit"`,
    `"ground"`, or `"production"`)."""
    weight: float = 1.0
    """float: The likelihood ratio of the history, drawn from stressed
    failures (see `Deployment.simulate_many`), or 1."""
    stages: list[float | None] = field(default_factory=list)
    """list: When each stage of the deployment was decided (years from the
    start of development, see `Stage`), or `None` if it was not."""

    @property
    def schedule(self) -> DeploymentSchedule:
        """DeploymentSchedule: The schedule of the deployment."""
        return self.deployment.schedule

    @property
    def planes(self) -> list[tuple[Any, int, int]]:
        """list: The planes of the deployment."""
        return self.deployment.planes

    @property
    def development(self) -> dict[Any, float]:
        """dict: The end of the development of each type of satellite."""
        return self.deployment.development

    @property
    def masses(self) -> dict[Any, float | None]:
        """dict: The mass of each type of satellite, if given."""
        return self.deployment.masses

    @property
    def begin(self) -> float:
        """float: The beginning of the mission (its planned first launch)."""
        return self.deployment.begin

    @property
    def end(self) -> float:
        """float: The end of the mission."""
        return self.deployment.end

    @property
    def deployed(self) -> float:
        """float: The time by which the satellites of the first generation
        are launched (years), whether their launches succeed or not: of the
        planes that are not staged, and of the stages that this history
        deployed (see `stages`)."""
        return max(
            launch.time for launch in self.launches if launch.purpose == "deployment"
        )

    def in_orbit(self) -> list[tuple[float, float, int]]:
        """
        Gets the intervals in which satellites are in orbit (operating, or
        spares waiting), for the cost of their operations.

        Returns:
            list: the intervals, and the number of satellites in each
        """
        spans = [
            (start, stop, 1)
            for intervals in self.operating.values()
            for start, stop in intervals
        ]
        return spans + list(self.parked)

    def phases(self, satellites: Sequence[Satellite | str]) -> list[Phase]:
        """
        Gets the phases of the mission, from its beginning to its end, in
        which constant sets of satellites operate.

        The satellites are those built for the layouts of the deployment
        (see `build_constellation`), in the order of their planes (and of
        several layouts, in the order given to `DeploymentSchedule.plan`),
        without spares: the satellite at each position of each plane is the
        one (and its replacements) that operates there.

        Args:
            satellites (Sequence[Satellite | str]): The satellites (or their
                names), in the order of the planes.

        Returns:
            list[Phase]: the phases, in order
        """
        named = _position_names(self.planes, satellites)
        names = list(named.values())
        intervals = {
            named[position]: spans for position, spans in self.operating.items()
        }
        begin, end = self.begin, self.end
        events = {begin, end}
        for spans in intervals.values():
            for start, stop in spans:
                events.update((start, stop))
        times = sorted(time for time in events if begin <= time <= end)
        phases: list[Phase] = []
        for start, stop in zip(times[:-1], times[1:]):
            if stop - start <= 1e-12:
                continue
            middle = (start + stop) / 2
            active = tuple(
                name
                for name in names
                if any(a <= middle < b for a, b in intervals.get(name, []))
            )
            if phases and phases[-1].satellites == active:
                phases[-1] = Phase(phases[-1].start, stop, active)
            else:
                phases.append(Phase(start, stop, active))
        return phases

    def metrics(self) -> dict[str, float]:
        """
        Gets measures of the mission as metrics of an evaluation: the
        failures of satellites (`failures`, including those lost in
        launches), the launches that failed (`launch failures`), the
        restorations (`restorations`), the spares used (`spares used`), the
        mean time from failure to restoration (`restoration years`, `NaN`
        without restorations), and the availability of the satellites
        (`availability`: their years of operation, as a share of those
        possible from the launches that deployed them to the end of the
        mission, whether these succeeded or not, and of the stages that this
        history deployed).

        Returns:
            dict[str, float]: the metrics
        """
        operated = sum(b - a for spans in self.operating.values() for a, b in spans)
        built = [b for _, b, _ in self.planes]
        possible = sum(
            max(self.end - launch.time, 0.0)
            * max(
                min(launch.first + launch.satellites, built[launch.plane])
                - launch.first,
                0,
            )
            for launch in self.launches
            if launch.purpose == "deployment"
        )
        delays = [r.restored - r.failed for r in self.restorations]
        metrics = {
            "failures": float(len(self.failures)),
            "launch failures": float(sum(launch.failed for launch in self.launches)),
            "restorations": float(len(self.restorations)),
            "spares used": float(
                sum(r.source != "production" for r in self.restorations)
            ),
            "restoration years": float(np.mean(delays)) if delays else math.nan,
            "availability": operated / possible if possible > 0 else math.nan,
        }
        if self.deployment.stages:
            metrics["stages deployed"] = float(
                sum(time is not None for time in self.stages)
            )
        return metrics


@dataclass(frozen=True)
class Decision:
    """
    What the rule of a stage of a deployment knows when it decides whether to
    deploy the stage (see `Stage`): the satellites operating and the history
    so far, and the capability of the satellites operating (see
    `capability`), such as for a stage deployed only if the satellites in
    orbit fall short of a requirement. Its time is from the start of
    development, as every time of a mission (see `since_launch` for the
    time after the first launch, as of the times of the stage).
    """

    time: float
    """float: The time of the decision (years from the start of
    development)."""
    stage: int
    """int: The stage (its index)."""
    satellites_operating: int
    """int: The number of satellites operating."""
    history: MissionHistory
    """MissionHistory: The history of the mission so far: its production,
    launches, failures, restorations, and the periods of operation that
    have ended."""
    stages_deployed: tuple[int, ...]
    """tuple[int, ...]: The stages deployed so far (their indices)."""
    _measure: Callable[[], dict[str, float]] | None = field(
        default=None, repr=False, compare=False
    )

    @property
    def since_launch(self) -> float:
        """float: The time of the decision after the first launch (years), as
        the times of stages are given (see `Stage`)."""
        return self.time - self.history.begin

    @cached_property
    def capability(self) -> dict[str, float]:
        """dict[str, float]: The metrics of the satellites operating at the
        decision, as of a phase in which they operate (see `Capability`):
        the capability of the satellites in orbit then, evaluated only if a
        rule reads it. It is the capability that the evaluation computes,
        as if the mission measured it exactly and at once: a rule that
        decides from an estimate adds its error (such as from an uncertain
        parameter of the design), and one that waits for a measurement
        decides later.

        Raises:
            ValueError: if the deployment has no capability (see
                `Deployment.capability`, which the evaluations of an
                `ObservationStudy` with a `MissionPlan` give)."""
        if self._measure is None:
            raise ValueError(
                "the capability of the satellites at a decision needs a "
                + "deployment with a capability (see Deployment.capability), "
                + "as the evaluations of an ObservationStudy with a plan give"
            )
        return dict(self._measure())


def evaluate_phases(
    phases: Sequence[Phase],
    observations: pd.DataFrame,
    metric: Callable[..., Mapping[str, float]],
    *,
    design: Any = None,
) -> pd.DataFrame:
    """
    Evaluates the phases of a mission: a metric of the observations of the
    satellites that operate in each phase, over the study's representative
    period (the observations of the full constellation, by satellite, such
    as from a `SatelliteCache`, filtered to those of each phase). Phases of
    the same satellites are evaluated once.

    The phases are quasi-static: each is evaluated over the same window, as
    if the geometry of the constellation (and the season) were that of the
    window throughout the mission. This misses the drift of the relative
    nodes of shells of different inclinations or altitudes (which precess at
    different rates), the drift of the local times of planes that are not
    sun-synchronous, and seasons (for metrics that depend on the sun, or on
    fields that vary with the seasons). Replicates average over these if the
    epoch of the window is an uncertain parameter drawn over the mission (see
    `Observations`), as are the relative nodes of shells (such as the nodes
    of `planes_space`). A constellation whose planes share an altitude and
    inclination keeps its geometry, so that only seasons remain.

    Args:
        phases (Sequence[Phase]): The phases (see `MissionHistory.phases`).
        observations (pandas.DataFrame): The observations of every
            satellite, with their names (`satellite`).
        metric (Callable): A function of the observations of a phase (which
            may be empty) that gives its metrics (and of the design, if it
            has a parameter named `design`, as the metrics of a study).
        design (Design | None): The (realized) design of the mission, for a
            metric with a parameter named `design`.

    Returns:
        pandas.DataFrame: a row per phase: its `start`, `end`, `duration`,
            number of `satellites`, and metrics
    """
    given = _given_design(metric, design)
    return _evaluate_phases(phases, observations, lambda subset, _: given(subset()))


def summarize_phases(
    phases: pd.DataFrame,
    *,
    metrics: Sequence[str] | None = None,
    requirements: Sequence[Requirement] = (),
) -> dict[str, float]:
    """
    Summarizes the capability of a mission over its phases (see
    `evaluate_phases`): the mean of each metric weighted by the durations of
    the phases (`<metric> mission`, over the mission, from its first launch),
    and, for each requirement on a metric, the time (years from the start of
    development) at which a phase first meets it (`<metric> time`, infinite
    if none does).

    Args:
        phases (pandas.DataFrame): The metrics of each phase (see
            `evaluate_phases`).
        metrics (Sequence[str] | None): The metrics, or `None` for every
            column other than those of the phases (see `PHASE_COLUMNS`).
        requirements (Sequence[Requirement]): Requirements on metrics.

    Returns:
        dict[str, float]: the summary, as metrics of an evaluation
    """
    if metrics is None:
        metrics = [column for column in phases.columns if column not in PHASE_COLUMNS]
    weights = phases["duration"].to_numpy(dtype=float)
    summary = {}
    for name in metrics:
        values = phases[name].to_numpy(dtype=float)
        summary[f"{name} mission"] = float(np.sum(values * weights) / np.sum(weights))
    for requirement in requirements:
        met = requirement.met(phases[requirement.name].to_numpy(dtype=float))
        starts = phases["start"].to_numpy(dtype=float)[met]
        summary[f"{requirement.name} time"] = (
            float(starts[0]) if len(starts) > 0 else math.inf
        )
    return summary


def _evaluate_phases(
    phases: Sequence[Phase],
    observations: pd.DataFrame,
    metric: Callable[[Callable[[], pd.DataFrame], frozenset], Mapping[str, float]],
) -> pd.DataFrame:
    """
    Evaluates phases with a metric of a function that gets their
    observations (only if needed) and the names of their satellites (see
    `evaluate_phases`), each set of satellites once.
    """
    rows, evaluated = [], {}
    for phase in phases:
        names = frozenset(phase.satellites)
        if names not in evaluated:
            subset = partial(_subset, observations, names)
            evaluated[names] = dict(metric(subset, names))
        rows.append(
            {
                "start": phase.start,
                "end": phase.end,
                "duration": phase.duration,
                "satellites": len(phase.satellites),
                **evaluated[names],
            }
        )
    return pd.DataFrame(rows)


def _subset(observations: pd.DataFrame, names: frozenset) -> pd.DataFrame:
    """The observations of the satellites of names (none if not by satellite)."""
    if "satellite" in observations.columns:
        return observations[observations["satellite"].isin(names)]
    return observations.iloc[0:0]


def _given_design(metric: Callable, design: Any) -> Callable:
    """
    Gets a metric of observations alone: the metric, given the design if it
    has a parameter named `design` (see `Study`).
    """
    return partial(metric, design=design) if _takes_design(metric) else metric


def _as_cache(metric: Callable[..., Mapping[str, float]] | PhaseCache) -> PhaseCache:
    """Gets a cache of a metric (see `PhaseCache`), or the cache given."""
    return metric if isinstance(metric, PhaseCache) else PhaseCache(metric)


class PhaseCache:
    """
    A cache of a metric of the phases of missions (see `evaluate_phases`),
    keyed by a context that identifies the observations (such as the seed of
    a replicate's sample of targets and the analysis period) and by the
    content of the satellites that operate (see `get_satellite_key`), so that
    the phases of missions with the same satellites (of other histories,
    policies, or designs, in the same replicate) are evaluated once.

    The metric must depend only on the observations (of the one window of
    the observations, see `evaluate_phases`), not on the names of the
    satellites. A metric with a parameter named `design` is given the design
    of the mission (see `evaluate`), on which it must depend only through
    what the context identifies (such as a scenario). Copies sent to the
    parallel workers of a study start with an empty memory, which the
    evaluations of a chunk share (consecutive evaluations of the same
    replicate, see `Study`).

    Args:
        metric (Callable): A function of the observations of a phase (which
            may be empty) that gives its metrics.
    """

    def __init__(self, metric: Callable[..., Mapping[str, float]]):
        self.metric = metric
        self.hits = 0
        """int: The phases taken from the cache."""
        self.misses = 0
        """int: The phases evaluated."""
        self._memory: dict = {}

    def __getstate__(self) -> dict:
        state = self.__dict__.copy()
        state["_memory"] = {}
        return state

    def evaluate(
        self,
        phases: Sequence[Phase],
        satellites: Sequence[Satellite | str] | Mapping[str, str],
        observations: pd.DataFrame,
        *,
        context: Any = None,
        design: Any = None,
    ) -> pd.DataFrame:
        """
        Evaluates the phases of a mission (see `evaluate_phases`), taking
        the metrics of the satellites of a phase from the cache if they are
        in it.

        Args:
            phases (Sequence[Phase]): The phases.
            satellites (Sequence[Satellite | str] | Mapping[str, str]): The
                satellites, which key the cache by their content (or by
                their names, if given as names, or by the keys given of their
                names, see `get_satellite_key`).
            observations (pandas.DataFrame): The observations of every
                satellite, with their names (`satellite`).
            context (Any): A hashable identifier of the observations (such
                as the seed of their targets and their period).
            design (Design | None): The (realized) design of the mission,
                for a metric with a parameter named `design`.

        Returns:
            pandas.DataFrame: a row per phase: its `start`, `end`,
                `duration`, number of `satellites`, and metrics
        """
        keys = _satellite_keys(satellites)

        def cached(subset: Callable[[], pd.DataFrame], names: frozenset) -> dict:
            del subset  # the cache filters the observations, if needed
            return self.metrics_of(
                {name: keys[name] for name in names},
                observations,
                context=context,
                design=design,
            )

        return _evaluate_phases(phases, observations, cached)

    def metrics_of(
        self,
        keys: Mapping[str, str],
        observations: pd.DataFrame,
        *,
        context: Any = None,
        design: Any = None,
    ) -> dict[str, float]:
        """
        Gets the metrics of satellites that operate together, from the cache
        if they are in it (keyed by the context and the keys of the
        satellites), or else of their observations (and kept).

        Args:
            keys (Mapping[str, str]): The keys of the satellites in the cache
                (such as of their content, see `get_satellite_key`), by
                their names in the observations.
            observations (pandas.DataFrame): The observations of every
                satellite, with their names (`satellite`), filtered to those
                of the satellites only if they are not in the cache.
            context (Any): A hashable identifier of the observations (such
                as the seed of their targets and their period).
            design (Design | None): The (realized) design of the mission,
                for a metric with a parameter named `design`.

        Returns:
            dict[str, float]: the metrics
        """
        key = (context, frozenset(keys.values()))
        if key in self._memory:
            self.hits += 1
        else:
            self.misses += 1
            metric = _given_design(self.metric, design)
            self._memory[key] = dict(metric(_subset(observations, frozenset(keys))))
        return dict(self._memory[key])

    def clear(self) -> None:
        """
        Clears the cache and its statistics.
        """
        self._memory.clear()
        self.hits = self.misses = 0

    def __repr__(self) -> str:
        return (
            f"PhaseCache({len(self._memory)} results, "
            + f"{self.hits} hits, {self.misses} misses)"
        )


class Capability:
    """
    The capability of the satellites that operate together at a time of a
    mission: the metrics of their observations, as of a phase in which they
    operate (over the study's representative period, see
    `evaluate_phases`), which the rules of stages read at their decisions
    (see `Decision.capability` and `Deployment.capability`). The
    evaluations of an `ObservationStudy` with a `MissionPlan` give their
    deployments the capability of their satellites, with the cache of their
    phases: the satellites operating at a decision are those of a phase of
    the mission, which the evaluation evaluates anyway.

    A capability also evaluates the phases of the mission's histories (see
    `evaluate`, `evaluate_histories`, and `exact_expectations`), with its
    cache.

    Args:
        satellites (Sequence[Satellite | str]): The satellites (or their
            names), in the order of the planes (see `MissionHistory.phases`).
        observations (pandas.DataFrame): The observations of every
            satellite, with their names (`satellite`).
        metric (Callable | PhaseCache): A function of the observations of
            satellites (which may be empty) that gives their metrics, or a
            cache of one (which keys satellites by their content, or by the
            names given for them).
        context (Any): A hashable identifier of the observations, for a
            cache (see `PhaseCache`).
        design (Design | None): The (realized) design of the mission, for a
            metric with a parameter named `design`.
    """

    def __init__(
        self,
        satellites: Sequence[Satellite | str],
        observations: pd.DataFrame,
        metric: Callable[..., Mapping[str, float]] | PhaseCache,
        *,
        context: Any = None,
        design: Any = None,
    ):
        self.satellites = list(satellites)
        """list: The satellites (or their names), in the order of the planes."""
        self.observations = observations
        """pandas.DataFrame: The observations of every satellite."""
        self.metric = metric
        """Callable | PhaseCache: The metric, or a cache of one."""
        self.context = context
        """Any: The identifier of the observations, for a cache."""
        self.design = design
        """Design | None: The (realized) design of the mission."""
        self._names = _names(self.satellites)

    @cached_property
    def _keys(self) -> dict[str, str]:
        """The keys of the satellites in a cache, by name."""
        return _satellite_keys(self.satellites)

    @cached_property
    def _cache(self) -> PhaseCache:
        """The cache of the metric (that given, or one of its own)."""
        return _as_cache(self.metric)

    def evaluate(self, phases: Sequence[Phase]) -> pd.DataFrame:
        """
        Evaluates the phases of a mission of the satellites (see
        `evaluate_phases`): the metrics of the satellites that operate in
        each, each set evaluated once (by the cache, with the context and the
        design, see `PhaseCache.evaluate`).

        Args:
            phases (Sequence[Phase]): The phases.

        Returns:
            pandas.DataFrame: a row per phase: its `start`, `end`,
                `duration`, number of `satellites`, and metrics
        """
        return self._cache.evaluate(
            phases,
            self._keys,
            self.observations,
            context=self.context,
            design=self.design,
        )

    def __call__(self, names: Iterable[str]) -> dict[str, float]:
        """
        Gets the metrics of satellites that operate together.

        Args:
            names (Iterable[str]): The names of the satellites.

        Returns:
            dict[str, float]: the metrics
        """
        names = frozenset(names)
        unknown = names - set(self._names)
        if unknown:
            raise ValueError(f"satellites that are not given: {sorted(unknown)}")
        return self._cache.metrics_of(
            {name: self._keys[name] for name in names},
            self.observations,
            context=self.context,
            design=self.design,
        )

    def at(
        self,
        planes: Sequence[tuple[Any, int, int]],
        positions: Iterable[tuple[int, int]],
    ) -> dict[str, float]:
        """
        Gets the metrics of the satellites at positions of the planes of a
        deployment (see `Deployment.planes`).

        Args:
            planes (Sequence[tuple[Any, int, int]]): The planes: the type of
                their satellites, and the number of satellites and of spares.
            positions (Iterable[tuple[int, int]]): The plane and position of
                each satellite.

        Returns:
            dict[str, float]: the metrics
        """
        named = _position_names(planes, self._names)
        return self(named[plane, position] for plane, position in positions)
