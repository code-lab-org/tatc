"""
Events of the deployments of constellations: the production and launches of
satellites (see `Deployment` and `MissionHistory`), their booking on the
production lines and at the launch slots of a schedule, which plans and
simulations share, and values given for each type of satellite.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import bisect
import heapq
import itertools
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, NamedTuple

if TYPE_CHECKING:
    from .deployment import Deployment, DeploymentSchedule

_PerType = float | Mapping[Any, float]
"""A value for every type of satellite, or for each (a mapping by type)."""


def _per_type(value: Any, kind: Any, name: str) -> Any:
    """
    Gets the value of a type of satellite of a value for every type or a
    mapping of each (see `_PerType`), raising a `ValueError` if the mapping
    lacks the type (`name` names the value in the error).
    """
    if isinstance(value, Mapping):
        if kind not in value:
            raise ValueError(f"no {name} of satellites of type {kind!r}")
        return value[kind]
    return value


def _per_type_values(value: Any) -> list:
    """
    Gets the values of every type of satellite of a value for every type or a
    mapping of each (see `_PerType`), such as to validate them.
    """
    return list(value.values()) if isinstance(value, Mapping) else [value]


def _production_years(schedule: DeploymentSchedule, kind: Any, count: int) -> float:
    """
    Gets the time to produce satellites of a type on its line (0 for
    production that takes no time).
    """
    rate = _per_type(schedule.production_rate, kind, "production rate")
    return 0.0 if rate is None else count / rate


def _next_slot(
    schedule: DeploymentSchedule, booked: Sequence[float], ready: float
) -> float:
    """
    Gets the earliest time of a launch at or after a time that keeps the
    launch interval from the launches booked before (in order), as plans and
    simulations book every launch: of deployments, replacements,
    restorations, and spares.
    """
    interval = schedule.launch_interval
    time = ready
    if interval > 0:
        for other in booked:
            if other - interval < time < other + interval:
                time = other + interval
    return time


def _per_launch(
    schedule: DeploymentSchedule,
    mass: float | None = None,
    plane: tuple[int, Any] | None = None,
) -> float:
    """
    Gets the most satellites of a launch (to a plane: its index and values,
    for a capacity in mass that is a function of the plane), by the capacity
    of launches in number and in mass (of satellites of a mass, kg), infinite
    if not limited, raising a `ValueError` if a capacity in mass lacks the
    mass (or the plane) or a satellite exceeds it.
    """
    capacity = (
        math.inf if schedule.launch_capacity is None else schedule.launch_capacity
    )
    mass_capacity = schedule.launch_mass_capacity
    if callable(mass_capacity):
        if plane is None:
            raise ValueError("a launch capacity of each plane needs the plane")
        mass_capacity = float(mass_capacity(*plane))
    if mass_capacity is not None:
        if mass is None:
            raise ValueError("a launch capacity in mass requires masses")
        fits = math.floor(mass_capacity / mass)
        if fits < 1:
            raise ValueError("a satellite exceeds the mass capacity of a launch")
        capacity = min(capacity, fits)
    return capacity


def _split(
    schedule: DeploymentSchedule,
    satellites: int,
    mass: float | None = None,
    plane: tuple[int, Any] | None = None,
) -> list[int]:
    """
    Gets the satellites of each launch of a plane: as few launches as carry
    them (see `_per_launch`, at least one), sharing them evenly.
    """
    capacity = _per_launch(schedule, mass, plane)
    count = max(1, math.ceil(satellites / capacity)) if capacity < math.inf else 1
    return [
        satellites // count + (1 if j < satellites % count else 0) for j in range(count)
    ]


class PlannedPlane(NamedTuple):
    """A plane of a deployment (see `Deployment.planes`)."""

    kind: Any
    """Any: The type of its satellites."""
    satellites: int
    """int: The number of satellites of the layout (its positions)."""
    orbit_spares: int
    """int: The number of spares in orbit."""


@dataclass(frozen=True)
class Launch:
    """
    A launch of satellites of a plane (see `Deployment` and
    `MissionHistory`).
    """

    time: float
    """float: The time (years from the start of development)."""
    plane: int
    """int: The plane (its index, in the order of the layouts' planes)."""
    kind: Any
    """Any: The type of its satellites."""
    first: int
    """int: The position in its plane of its first satellite (-1 for the
    launches of restorations and spares)."""
    satellites: int
    """int: The number of satellites (including spares)."""
    generation: int = 0
    """int: The generation (0 for the first, others planned replacements)."""
    purpose: str = "deployment"
    """str: `"deployment"`, `"replacement"` (planned, at the end of a
    design life), `"restoration"` (of failed satellites), or `"spare"` (a
    spare in orbit replenished)."""
    failed: bool = False
    """bool: Whether the launch failed, losing its satellites."""


@dataclass(frozen=True)
class Production:
    """
    The production of satellites of a type (see `Deployment` and
    `MissionHistory`).
    """

    kind: Any
    """Any: The type of the satellites."""
    begin: float
    """float: The start of their production (years)."""
    end: float
    """float: The end of their production (years)."""
    satellites: int
    """int: The number of satellites."""
    purpose: str = "deployment"
    """str: `"deployment"`, `"replacement"`, `"restoration"`, or `"spare"`."""


class _Manifest:
    """
    The production and launches of a deployment as they are booked, which
    plans (see `DeploymentSchedule.plan`) and simulations (of the stages they
    deploy, and of restorations) share: the production line of each type of
    satellite, the launches booked (which later ones keep the launch
    interval from, see `_next_slot`), and the production
    added.

    Args:
        deployment (Deployment): The deployment: its schedule, planes, masses,
            and values of planes, which a plan fills as it goes.
        production (list[Production]): The production, to which it adds.
        lines (dict | None): The time from which each type's line is free, or
            `None` for the end of its development.
        booked (list[float] | None): The times of the launches booked, in
            order, or `None` for none.
    """

    def __init__(
        self,
        deployment: Deployment,
        production: list[Production],
        lines: dict[Any, float] | None = None,
        booked: list[float] | None = None,
    ):
        self.deployment, self.schedule = deployment, deployment.schedule
        self.production = production
        self.lines = dict(deployment.development) if lines is None else lines
        self.booked = [] if booked is None else booked

    @classmethod
    def resume(cls, deployment: Deployment, production: list[Production]) -> _Manifest:
        """
        Gets the manifest of a planned deployment as its simulations continue
        it: each type's line free after the production of its planes (or its
        development), and the launches of the plan booked.
        """
        lines = {
            kind: max(
                (
                    batch.end
                    for batch in deployment.production
                    if batch.kind == kind and batch.purpose == "deployment"
                ),
                default=developed,
            )
            for kind, developed in deployment.development.items()
        }
        booked = sorted(launch.time for launch in deployment.launches)
        return cls(deployment, production, lines, booked)

    def plane(self, plane: int) -> tuple[int, Any]:
        """The index and values of a plane (see `DeploymentSchedule`)."""
        values = self.deployment.plane_values
        return plane, values[plane] if plane < len(values) else None

    def capacity(self, plane: int) -> float:
        """The most satellites of a launch to a plane (see `_per_launch`)."""
        kind = self.deployment.planes[plane].kind
        return _per_launch(
            self.schedule, self.deployment.masses.get(kind), self.plane(plane)
        )

    def produce(
        self, earliest: float, kind: Any, count: int, purpose: str = "deployment"
    ) -> float:
        """
        Produces satellites on their type's line, from a time (after the work
        before them), returning when they are ready.
        """
        begin = max(self.lines[kind], earliest)
        ready = begin + _production_years(self.schedule, kind, count)
        self.lines[kind] = ready
        self.production.append(Production(kind, begin, ready, count, purpose))
        return ready

    def next_slot(self, ready: float) -> float:
        """
        Gets the earliest launch at or after a time that keeps the launch
        interval (see `_next_slot`), without booking it.
        """
        return _next_slot(self.schedule, self.booked, ready)

    def slot(self, ready: float) -> float:
        """
        Books the earliest launch at or after a time that keeps the launch
        interval (see `_next_slot`), returning its time.
        """
        time = self.next_slot(ready)
        bisect.insort(self.booked, time)
        return time

    def launch(self, plane: int, ready: float) -> list[Launch]:
        """
        Launches the satellites of a plane (and its spares in orbit) when
        they are ready: in as many launches as they need (see
        `_split`), each at the next launch slot.
        """
        kind, size, spares = self.deployment.planes[plane]
        launches, first = [], 0
        for satellites in _split(
            self.schedule,
            size + spares,
            self.deployment.masses.get(kind),
            self.plane(plane),
        ):
            launches.append(Launch(self.slot(ready), plane, kind, first, satellites))
            first += satellites
        return launches

    def deploy(
        self, planes: Sequence[int], earliest: float, end: float | None = None
    ) -> list[Launch]:
        """
        Deploys planes, as a plan deploys those that are not staged and a
        simulation the planes of a stage: produces the satellites of each
        (and its spares in orbit) on its type's line from a time (after the
        work before them), in order; launches each plane when they are ready
        (see `launch`), in the order the planes are ready; and adds the
        planned replacements of their satellites (see `replace`). Returns
        the launches, then the replacements.
        """
        ready = []
        for plane in planes:
            kind, size, spares = self.deployment.planes[plane]
            if size + spares > 0:
                ready.append((self.produce(earliest, kind, size + spares), plane))
        launches = []
        for time, plane in sorted(ready, key=lambda item: item[0]):
            launches.extend(self.launch(plane, time))
        if end is None:  # a mission that begins at the first launch
            end = min(launch.time for launch in launches) + self.schedule.years
        return launches + self.replace(launches, end)

    def replace(self, launches: Sequence[Launch], end: float) -> list[Launch]:
        """
        Adds the planned replacements of the satellites of launches at the end
        of their design life (see `DeploymentSchedule`), returning their
        launches: of the positions of the layout that a launch carries (not
        of spares in orbit, which do not age), each produced just before it
        is needed (at the production rate, apart from the line) and launched
        at the next launch slot then, in the order they are needed, and
        replaced in turn a design life after that launch. A replacement that
        no slot before the end of the mission launches is not produced, but
        for a sustained constellation, whose replacements are produced after
        the end (their production begun within the development years of
        their type after it), though not launched.
        """
        schedule = self.schedule
        queue: list = []
        order = itertools.count()
        for launch in launches:
            life = _per_type(schedule.design_life, launch.kind, "design life")
            built = self.deployment.planes[launch.plane].satellites
            count = min(launch.first + launch.satellites, built) - launch.first
            if life is not None and count > 0:
                item = (launch.plane, launch.kind, launch.first, count, life)
                heapq.heappush(queue, (launch.time + life, next(order), 1, item))
        replacements = []
        while queue:
            need, _, generation, item = heapq.heappop(queue)
            plane, kind, first, count, life = item
            begin = need - _production_years(schedule, kind, count)
            time = math.inf
            if need < end - 1e-9:
                time = _next_slot(schedule, self.booked, need)
            if time < end - 1e-9:
                bisect.insort(self.booked, time)
                replacements.append(
                    Launch(time, plane, kind, first, count, generation, "replacement")
                )
                following = time + life
            elif (
                schedule.sustained
                and begin
                < end
                + float(_per_type(schedule.development_years, kind, "development"))
                - 1e-9
            ):
                following = need + life  # produced after the end, not launched
            else:
                continue
            self.production.append(Production(kind, begin, need, count, "replacement"))
            heapq.heappush(queue, (following, next(order), generation + 1, item))
        return replacements
