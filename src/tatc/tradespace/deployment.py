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
takes longer to produce and launch. Its launches divide the mission into
phases (`Deployment.phases`) of constant sets of operating satellites, each
evaluated over the study's representative period (`evaluate_phases`) and
combined over the mission (`summarize_phases`), such as the time-weighted
mean of a metric and the time until it first meets a requirement.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from ..schemas import Satellite
from .constellations import PLANES, SSO_PLANES, WALKER, _kind, _layout, _layouts
from .criteria import _OPERATORS, Requirement
from .values import Values

PerType = float | Mapping[Any, float]
"""A value for every type of satellite, or for each (by its name)."""


def _per_type(value: Any, kind: Any, name: str) -> Any:
    """Gets the value of a type of satellite."""
    if isinstance(value, Mapping):
        if kind not in value:
            raise ValueError(f"no {name} of satellites of kind {kind!r}")
        return value[kind]
    return value


def _planes(layout: Any) -> list[tuple[int, Any]]:
    """
    Gets the planes of a layout (or several): the number of satellites of
    each, and its values (or `None` for the planes of a Walker layout).
    """
    several = _layouts(layout)
    if several is not None:
        return [plane for item in several for plane in _planes(item)]
    layout = _layout(layout)
    kind = _kind(layout)
    if kind == WALKER:
        return [(int(layout.per_plane), None)] * int(layout.planes)
    if kind in (SSO_PLANES, PLANES):
        return [(int(plane.satellites), plane) for plane in layout.planes]
    raise ValueError(f"unknown layout {kind!r}")


@dataclass(frozen=True)
class Launch:
    """
    A launch of satellites of a plane (see `Deployment`).
    """

    time: float
    """float: The time (years from the start of development)."""
    plane: int
    """int: The plane (its index, in the order of the layouts' planes)."""
    kind: Any
    """Any: The type of its satellites."""
    first: int
    """int: The position in its plane of its first satellite."""
    satellites: int
    """int: The number of satellites (including spares)."""
    generation: int = 0
    """int: The generation (0 for the first, others replacements)."""


@dataclass(frozen=True)
class Production:
    """
    The production of satellites of a type (see `Deployment`).
    """

    kind: Any
    """Any: The type of the satellites."""
    begin: float
    """float: The start of their production (years)."""
    end: float
    """float: The end of their production (years)."""
    satellites: int
    """int: The number of satellites."""


@dataclass(frozen=True)
class Phase:
    """
    A phase of a mission, in which a set of satellites operates (see
    `Deployment.phases`).
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


@dataclass
class Deployment:
    """
    The deployment of a constellation: its development, production,
    launches, and replacements (see `DeploymentSchedule.plan`).
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
    """list[Launch]: The launches (and replacements), in order."""

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
        return max(launch.time for launch in self.launches if launch.generation == 0)

    @property
    def units(self) -> dict[Any, int]:
        """dict: The satellites produced of each type (including replacements
        and spares)."""
        units: dict[Any, int] = {}
        for batch in self.production:
            units[batch.kind] = units.get(batch.kind, 0) + batch.satellites
        return units

    def phases(
        self,
        satellites: Sequence[Satellite | str],
        failures: Mapping[str, float] | None = None,
    ) -> list[Phase]:
        """
        Gets the phases of the mission, from the first launch to its end, in
        which constant sets of satellites operate: each satellite from its
        launch, except, if it fails (at a time given by `failures`, see
        `draw_failures`), from its failure to the launch of its replacement
        (or the end of the mission).

        The satellites are those built for the layouts of the deployment
        (see `build_constellation`), in the order of their planes (and of
        several layouts, in the order given to `DeploymentSchedule.plan`),
        without spares. The satellites of a plane of several launches
        operate as each is launched, in order.

        Args:
            satellites (Sequence[Satellite | str]): The satellites (or their
                names), in the order of the planes.
            failures (Mapping[str, float] | None): The time (years from the
                start of development) at which satellites fail, by name.

        Returns:
            list[Phase]: the phases, in order
        """
        names = [s.name if isinstance(s, Satellite) else str(s) for s in satellites]
        expected = sum(built for _, built, _ in self.planes)
        if len(names) != expected:
            raise ValueError(
                f"expected {expected} satellites for the planes, not {len(names)}"
            )
        offsets = np.cumsum([0] + [built for _, built, _ in self.planes])
        begin, end = self.begin, self.end
        failures = failures or {}
        intervals: dict[str, list[tuple[float, float]]] = {}
        for launch in self.launches:
            if launch.generation != 0:
                continue
            built = self.planes[launch.plane][1]
            for position in range(
                launch.first, min(launch.first + launch.satellites, built)
            ):
                name = names[offsets[launch.plane] + position]
                start = launch.time
                failure = failures.get(name)
                if failure is None or failure >= end or failure < start:
                    intervals[name] = [(start, end)]
                    continue
                replaced = [
                    time for time in self._replacements(launch) if failure < time < end
                ]
                after = [(replaced[0], end)] if replaced else []
                intervals[name] = [(start, failure)] + after
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

    def _replacements(self, launch: Launch) -> list[float]:
        """Gets the times of the replacements of the satellites of a launch."""
        return sorted(
            other.time
            for other in self.launches
            if other.generation > 0
            and other.plane == launch.plane
            and other.first == launch.first
        )


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
      before it is needed and launched then.

    A schedule that is itself a decision (such as the launch cadence) is a
    design variable, with the schedule made for each design (such as with
    `dataclasses.replace`).

    Args:
        years (float): The duration of the mission, from the first launch.
        design_life (float | None): The design life of the satellites, or
            `None` if they last the mission.
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
    ) -> Deployment:
        """
        Plans the deployment of a constellation: the layout (or layouts) of
        its satellites, of one type (`kind`), or a mapping of types to
        layouts (such as the shells of a mixed constellation).

        Args:
            layout (Any): The layout (or layouts) of the satellites, or a
                mapping of types of satellites to layouts.
            kind (Any): The type of the satellites of a layout.
            spares (int): The spares in orbit in each plane, produced and
                launched with it.
            masses (float | Mapping | None): The mass (kg) of the satellites
                (of each type), for a launch capacity in mass.

        Returns:
            Deployment: the deployment
        """
        # pylint: disable=too-many-locals
        if spares < 0:
            raise ValueError("spares must not be negative")
        # a mapping of types to layouts (the values of a layout are a mapping too)
        if isinstance(layout, Mapping) and not isinstance(layout, Values):
            fleets = dict(layout)
        else:
            fleets = {kind: layout}
        deployment = Deployment(schedule=self)
        planes = []  # (ready, index, kind, satellites)
        for name, item in fleets.items():
            developed = float(_per_type(self.development_years, name, "development"))
            rate = _per_type(self.production_rate, name, "production rate")
            mass = None if masses is None else _per_type(masses, name, "mass")
            deployment.development[name] = developed
            listed = []
            for size, values in _planes(item):
                index = len(deployment.planes)
                deployment.planes.append((name, size, spares if size > 0 else 0))
                if size > 0:
                    listed.append((index, size + spares, values))
            if self.order is not None:
                listed.sort(key=lambda plane: self.order(plane[0], plane[2]))
            line = developed
            for index, size, _ in listed:
                ready = line + (0.0 if rate is None else size / rate)
                deployment.production.append(Production(name, line, ready, size))
                planes.append((ready, index, name, size, mass))
                line = ready
        if not planes:
            raise ValueError("a deployment requires satellites")
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
                        Production(launch.kind, produced, need, launch.satellites)
                    )
                    deployment.launches.append(
                        Launch(
                            need,
                            launch.plane,
                            launch.kind,
                            launch.first,
                            launch.satellites,
                            generation,
                        )
                    )
                    need, generation = need + self.design_life, generation + 1
        deployment.production.sort(key=lambda batch: batch.begin)
        deployment.launches.sort(key=lambda launch: (launch.time, launch.plane))
        return deployment


def draw_failures(
    deployment: Deployment,
    satellites: Sequence[Satellite | str],
    annual_rate: float,
    seed: int | None,
) -> dict[str, float]:
    """
    Draws the times at which satellites fail (see `Deployment.phases`), each
    after an exponentially distributed life from its launch (at a constant
    `annual_rate` of failure), from an uncertain seed of a replicate (see
    `Uncertain.seed`). The life of a satellite depends only on the seed and
    its position, so that with common random numbers, the satellite at the
    same position has the same life in the same replicate of every design.

    Args:
        deployment (Deployment): The deployment.
        satellites (Sequence[Satellite | str]): The satellites (or their
            names), in the order of the planes (see `Deployment.phases`).
        annual_rate (float): The rate of failure (per year).
        seed (int | None): The seed of the replicate, or `None` for no
            failures.

    Returns:
        dict[str, float]: the time of failure of each satellite that fails
            before the end of the mission
    """
    if annual_rate < 0:
        raise ValueError("annual_rate must not be negative")
    if seed is None or annual_rate == 0:
        return {}
    launched = {}
    for phase in deployment.phases(satellites):
        for name in phase.satellites:
            launched.setdefault(name, phase.start)
    names = [s.name if isinstance(s, Satellite) else str(s) for s in satellites]
    failures = {}
    for k, name in enumerate(names):
        life = np.random.default_rng([seed, k]).exponential(1 / annual_rate)
        time = launched[name] + life
        if time < deployment.end:
            failures[name] = float(time)
    return failures


def evaluate_phases(
    phases: Sequence[Phase],
    observations: pd.DataFrame,
    metric: Callable[[pd.DataFrame], Mapping[str, float]],
) -> pd.DataFrame:
    """
    Evaluates the phases of a mission: a metric of the observations of the
    satellites that operate in each phase, over the study's representative
    period (the observations of the full constellation, by satellite, such
    as from a `SatelliteCache`, filtered to those of each phase). Phases of
    the same satellites are evaluated once.

    Args:
        phases (Sequence[Phase]): The phases (see `Deployment.phases`).
        observations (pandas.DataFrame): The observations of every
            satellite, with their names (`satellite`).
        metric (Callable): A function of the observations of a phase (which
            may be empty) that gives its metrics.

    Returns:
        pandas.DataFrame: a row per phase: its `start`, `end`, `duration`,
            number of `satellites`, and metrics
    """
    rows, evaluated = [], {}
    for phase in phases:
        key = frozenset(phase.satellites)
        if key not in evaluated:
            if "satellite" in observations.columns:
                subset = observations[observations["satellite"].isin(key)]
            else:
                subset = observations.iloc[0:0]
            evaluated[key] = dict(metric(subset))
        rows.append(
            {
                "start": phase.start,
                "end": phase.end,
                "duration": phase.duration,
                "satellites": len(phase.satellites),
                **evaluated[key],
            }
        )
    return pd.DataFrame(rows)


def summarize_phases(
    capability: pd.DataFrame,
    metrics: Sequence[str] | None = None,
    requirements: Sequence[Requirement] = (),
) -> dict[str, float]:
    """
    Summarizes the capability of a mission over its phases (see
    `evaluate_phases`): the mean of each metric weighted by the durations of
    the phases (`<metric> mean`, over the mission, from its first launch),
    and, for each requirement on a metric, the time (years from the start of
    development) at which a phase first meets it (`<metric> time`, infinite
    if none does).

    Args:
        capability (pandas.DataFrame): The metrics of each phase.
        metrics (Sequence[str] | None): The metrics, or `None` for every
            column other than those of the phases.
        requirements (Sequence[Requirement]): Requirements on metrics.

    Returns:
        dict[str, float]: the summary, as metrics of an evaluation
    """
    if metrics is None:
        metrics = [
            column
            for column in capability.columns
            if column not in ("start", "end", "duration", "satellites")
        ]
    weights = capability["duration"].to_numpy(dtype=float)
    summary = {}
    for name in metrics:
        values = capability[name].to_numpy(dtype=float)
        summary[f"{name} mean"] = float(np.sum(values * weights) / np.sum(weights))
    for requirement in requirements:
        met = _OPERATORS[requirement.operator](
            capability[requirement.name].to_numpy(dtype=float), requirement.value
        )
        starts = capability["start"].to_numpy(dtype=float)[met]
        summary[f"{requirement.name} time"] = (
            float(starts[0]) if len(starts) > 0 else math.inf
        )
    return summary
