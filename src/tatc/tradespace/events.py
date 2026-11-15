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
    from .deployment import Deployment

PerType = float | Mapping[Any, float]
"""A value for every type of satellite, or for each (a mapping by type)."""


def per_type(value: Any, kind: Any, name: str) -> Any:
    """
    Gets the value of a type of satellite of a value for every type or a
    mapping of each (see `PerType`).

    Args:
        value (Any): The value for every type, or a mapping of each type to
            its value.
        kind (Any): The type of satellite.
        name (str): The name of the value, for errors.

    Returns:
        Any: the value of the type

    Raises:
        ValueError: if the mapping lacks the type.
    """
    if isinstance(value, Mapping):
        if kind not in value:
            raise ValueError(f"no {name} of satellites of type {kind!r}")
        return value[kind]
    return value


def per_type_values(value: Any) -> list:
    """
    Gets the values of every type of satellite of a value for every type or a
    mapping of each (see `PerType`), such as to validate them.

    Args:
        value (Any): The value for every type, or a mapping of each type to
            its value.

    Returns:
        list: the values of the mapping, or the value
    """
    return list(value.values()) if isinstance(value, Mapping) else [value]


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
    interval from, see `DeploymentSchedule.next_slot`), and the production
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
        """The most satellites of a launch to a plane (see `per_launch`)."""
        kind = self.deployment.planes[plane].kind
        return self.schedule.per_launch(
            self.deployment.masses.get(kind), self.plane(plane)
        )

    def produce(
        self, earliest: float, kind: Any, count: int, purpose: str = "deployment"
    ) -> float:
        """
        Produces satellites on their type's line, from a time (after the work
        before them), returning when they are ready.
        """
        begin = max(self.lines[kind], earliest)
        ready = begin + self.schedule.production_years(kind, count)
        self.lines[kind] = ready
        self.production.append(Production(kind, begin, ready, count, purpose))
        return ready

    def next_slot(self, ready: float) -> float:
        """
        Gets the earliest launch at or after a time that keeps the launch
        interval (see `DeploymentSchedule.next_slot`), without booking it.
        """
        return self.schedule.next_slot(self.booked, ready)

    def slot(self, ready: float) -> float:
        """
        Books the earliest launch at or after a time that keeps the launch
        interval (see `DeploymentSchedule.next_slot`), returning its time.
        """
        time = self.next_slot(ready)
        bisect.insort(self.booked, time)
        return time

    def launch(self, plane: int, ready: float) -> list[Launch]:
        """
        Launches the satellites of a plane (and its spares in orbit) when
        they are ready: in as many launches as they need (see
        `DeploymentSchedule.split`), each at the next launch slot.
        """
        kind, size, spares = self.deployment.planes[plane]
        launches, first = [], 0
        for satellites in self.schedule.split(
            size + spares, self.deployment.masses.get(kind), self.plane(plane)
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
            life = schedule.design_life_of(launch.kind)
            built = self.deployment.planes[launch.plane].satellites
            count = min(launch.first + launch.satellites, built) - launch.first
            if life is not None and count > 0:
                item = (launch.plane, launch.kind, launch.first, count, life)
                heapq.heappush(queue, (launch.time + life, next(order), 1, item))
        replacements = []
        while queue:
            need, _, generation, item = heapq.heappop(queue)
            plane, kind, first, count, life = item
            begin = need - schedule.production_years(kind, count)
            time = math.inf
            if need < end - 1e-9:
                time = schedule.next_slot(self.booked, need)
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
                + float(per_type(schedule.development_years, kind, "development"))
                - 1e-9
            ):
                following = need + life  # produced after the end, not launched
            else:
                continue
            self.production.append(Production(kind, begin, need, count, "replacement"))
            heapq.heappush(queue, (following, next(order), generation + 1, item))
        return replacements
