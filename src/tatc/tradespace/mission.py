"""
Simulated missions: the history of a deployment with failures and
restorations (see `Deployment.simulate`), the phases of a mission in which
constant sets of satellites operate, and their evaluation over the study's
representative period and combination over the mission.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import heapq
import itertools
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np
import pandas as pd

from ..schemas import Satellite
from .criteria import _OPERATORS, Requirement
from .deployment import Deployment, DeploymentSchedule, Launch, Production, _per_type
from .failures import FailureModel, RestorationPolicy


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


def _names(satellites: Sequence[Satellite | str]) -> list[str]:
    return [s.name if isinstance(s, Satellite) else str(s) for s in satellites]


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
    failures: list[tuple[float, int, int, str]] = field(default_factory=list)
    """list: The failures: their time, plane, position, and cause (`"life"`
    or `"launch"`)."""
    restorations: list[tuple[int, int, float, float, str]] = field(default_factory=list)
    """list: The restorations: their plane, position, the time of the failure
    and of the restoration, and the source (`"orbit"`, `"ground"`, or
    `"production"`)."""

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
        names = _names(satellites)
        expected = sum(built for _, built, _ in self.planes)
        if len(names) != expected:
            raise ValueError(
                f"expected {expected} satellites for the planes, not {len(names)}"
            )
        offsets = np.cumsum([0] + [built for _, built, _ in self.planes])
        intervals = {
            names[offsets[plane] + position]: spans
            for (plane, position), spans in self.operating.items()
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
        (`availability`: their years of operation, as a share of those of
        the planned deployment).

        Returns:
            dict[str, float]: the metrics
        """
        planned = self.deployment.simulate()
        operated = sum(b - a for spans in self.operating.values() for a, b in spans)
        possible = sum(b - a for spans in planned.operating.values() for a, b in spans)
        delays = [restored - failed for *_, failed, restored, _ in self.restorations]
        return {
            "failures": float(len(self.failures)),
            "launch failures": float(sum(launch.failed for launch in self.launches)),
            "restorations": float(len(self.restorations)),
            "spares used": float(
                sum(source != "production" for *_, source in self.restorations)
            ),
            "restoration years": float(np.mean(delays)) if delays else math.nan,
            "availability": operated / possible if possible > 0 else math.nan,
        }


@dataclass
class _Booking:
    """A launch of restorations or spares, open to more until it flies."""

    time: float
    plane: int
    kind: Any
    items: list = field(default_factory=list)
    order: int = 0


class _Simulation:
    """
    A discrete-event simulation of a mission (see `Deployment.simulate`).
    """

    # pylint: disable=too-many-instance-attributes

    def __init__(
        self,
        deployment: Deployment,
        failures: FailureModel,
        policy: RestorationPolicy,
        seed: int | None,
    ):
        self.deployment, self.failures, self.policy = deployment, failures, policy
        self.seed = seed
        self.schedule = deployment.schedule
        self.end = deployment.end
        self.history = MissionHistory(
            deployment, production=list(deployment.production)
        )
        self.events: list = []
        self.counter = itertools.count()
        # the state of each position: its generation, when it began to
        # operate (or None), when it failed, and whether it awaits restoration
        self.slots = {
            (plane, position): {
                "generation": -1,
                "up": None,
                "down": None,
                "pending": False,
            }
            for plane, (_, built, _) in enumerate(deployment.planes)
            for position in range(built)
        }
        self.spares = dict.fromkeys(range(len(deployment.planes)), 0)
        self.spares_since: dict[int, float] = {}
        self.ground: dict[Any, int] = {}
        self.line = {
            kind: max(
                (
                    b.end
                    for b in deployment.production
                    if b.kind == kind and b.purpose == "deployment"
                ),
                default=developed,
            )
            for kind, developed in deployment.development.items()
        }
        self.booked = sorted(launch.time for launch in deployment.launches)
        self.open: dict[int, list[_Booking]] = {}
        self.launches_of_plane: dict[int, int] = {}

    # --- randomness

    def _rng(self, *key: int) -> np.random.Generator:
        return np.random.default_rng([self.seed, *key])

    def _life(self, plane: int, position: int, generation: int) -> float:
        if self.seed is None:
            return math.inf
        life = self.failures.life_of(self.deployment.planes[plane][0])
        if life is None:
            return math.inf
        return life.sample(self._rng(0, plane, position, generation))

    def _launch_fails(self, *key: int) -> bool:
        if self.seed is None or self.failures.launch == 0:
            return False
        return bool(self._rng(1, *key).random() < self.failures.launch)

    # --- events

    def _push(self, time: float, kind: str, *payload: Any) -> None:
        heapq.heappush(self.events, (time, next(self.counter), kind, payload))

    def run(self) -> MissionHistory:
        """Runs the simulation, returning the history of the mission."""
        # pylint: disable=too-many-branches
        for kind in self.deployment.development:
            count = self.policy.ground_spares_of(kind)
            if count > 0:
                begin = self.line[kind]
                ready = begin + self._producing(kind, count)
                self.history.production.append(
                    Production(kind, begin, ready, count, "spare")
                )
                self.line[kind] = ready
                self.ground.setdefault(kind, 0)
                self._push(ready, "stock", kind, count)
        for launch in self.deployment.launches:
            self._push(launch.time, "planned", launch)
        steps = 0
        while self.events:
            time, _, kind, payload = heapq.heappop(self.events)
            if time >= self.end:
                break
            steps += 1
            if steps > 1_000_000:
                raise RuntimeError("the simulation does not progress")
            getattr(self, f"_on_{kind}")(time, *payload)
        for (plane, position), slot in self.slots.items():
            if slot["up"] is not None:
                self._close(plane, position, self.end)
        for plane, since in self.spares_since.items():
            if self.spares[plane] > 0:
                self.history.parked.append((since, self.end, self.spares[plane]))
        self.history.production.sort(key=lambda batch: batch.begin)
        self.history.launches.sort(key=lambda launch: (launch.time, launch.plane))
        return self.history

    def _on_stock(self, time: float, kind: Any, count: int) -> None:
        del time
        self.ground[kind] = self.ground.get(kind, 0) + count

    def _on_planned(self, time: float, launch: Launch) -> None:
        failed = self._launch_fails(0, launch.plane, launch.first, launch.generation)
        self.history.launches.append(replace(launch, failed=failed))
        built = self.deployment.planes[launch.plane][1]
        positions = range(launch.first, launch.first + launch.satellites)
        for position in positions:
            if position >= built:
                continue
            slot = self.slots[launch.plane, position]
            if launch.purpose == "deployment":
                if failed:
                    self._fail(time, launch.plane, position, "launch")
                else:
                    self._up(time, launch.plane, position)
            elif not failed:
                # a planned replacement takes the place of the satellite
                if slot["up"] is not None:
                    self._close(launch.plane, position, time)
                self._up(time, launch.plane, position)
        spares = sum(position >= built for position in positions)
        if spares and not failed:
            level = self.deployment.planes[launch.plane][2]
            if launch.purpose == "deployment":
                self._park(time, launch.plane, spares)
            else:
                self._park(
                    time,
                    launch.plane,
                    max(0, min(spares, level - self.spares[launch.plane])),
                )

    def _on_fail(self, time: float, plane: int, position: int, generation: int) -> None:
        slot = self.slots[plane, position]
        if slot["generation"] != generation or slot["up"] is None:
            return  # a satellite replaced before it failed
        self._close(plane, position, time)
        self._fail(time, plane, position, "life")

    def _fail(self, time: float, plane: int, position: int, cause: str) -> None:
        slot = self.slots[plane, position]
        slot["down"] = time
        self.history.failures.append((time, plane, position, cause))
        if self.policy.restore and not slot["pending"]:
            slot["pending"] = True
            self._push(time + self.policy.detection_years, "request", plane, position)

    def _on_request(self, time: float, plane: int, position: int) -> None:
        slot = self.slots[plane, position]
        if slot["up"] is not None:  # restored by a planned replacement
            slot["pending"] = False
            return
        kind = self.deployment.planes[plane][0]
        if self.spares[plane] > 0:
            self._park(time, plane, -1)
            self._push(
                time + self.policy.drift_years, "activate", plane, position, "orbit"
            )
            if self.policy.replenish:
                ready = self._produce(time, kind, 1, "spare")
                self._book(plane, kind, ready, ("spare",))
        elif self.ground.get(kind, 0) > 0:
            self.ground[kind] -= 1
            self._book(plane, kind, time, ("slot", position, "ground"))
            if self.policy.replenish:
                ready = self._produce(time, kind, 1, "spare")
                self._push(ready, "stock", kind, 1)
        else:
            ready = self._produce(time, kind, 1, "restoration")
            self._book(plane, kind, ready, ("slot", position, "production"))

    def _on_activate(self, time: float, plane: int, position: int, source: str) -> None:
        slot = self.slots[plane, position]
        if slot["up"] is not None:
            self._park(time, plane, 1)  # not needed: back to waiting
            return
        self._restore(time, plane, position, source)

    def _on_booking(self, time: float, booking: _Booking) -> None:
        self.open[booking.plane].remove(booking)
        failed = self._launch_fails(1, booking.plane, booking.order)
        spare = all(item[0] == "spare" for item in booking.items)
        self.history.launches.append(
            Launch(
                time,
                booking.plane,
                booking.kind,
                -1,
                len(booking.items),
                purpose="spare" if spare else "restoration",
                failed=failed,
            )
        )
        for item in booking.items:
            if item[0] == "spare":
                if not failed:
                    self._park(time, booking.plane, 1)
                continue
            _, position, source = item
            slot = self.slots[booking.plane, position]
            if failed:
                slot["pending"] = False
                self.history.failures.append((time, booking.plane, position, "launch"))
                if self.policy.restore:
                    slot["pending"] = True
                    self._push(
                        time + max(self.policy.detection_years, 1e-6),
                        "request",
                        booking.plane,
                        position,
                    )
            elif slot["up"] is not None:
                self._park(time, booking.plane, 1)  # not needed: a spare in orbit
            else:
                self._restore(time, booking.plane, position, source)

    # --- actions

    def _up(self, time: float, plane: int, position: int) -> None:
        slot = self.slots[plane, position]
        slot["generation"] += 1
        slot["up"], slot["pending"] = time, False
        life = self._life(plane, position, slot["generation"])
        if time + life < self.end:
            self._push(time + life, "fail", plane, position, slot["generation"])

    def _close(self, plane: int, position: int, time: float) -> None:
        slot = self.slots[plane, position]
        if time > slot["up"]:
            self.history.operating.setdefault((plane, position), []).append(
                (slot["up"], time)
            )
        slot["up"] = None

    def _restore(self, time: float, plane: int, position: int, source: str) -> None:
        slot = self.slots[plane, position]
        self.history.restorations.append((plane, position, slot["down"], time, source))
        self._up(time, plane, position)

    def _park(self, time: float, plane: int, change: int) -> None:
        if change == 0:
            return
        if self.spares[plane] > 0:
            since = self.spares_since[plane]
            if time > since:
                self.history.parked.append((since, time, self.spares[plane]))
        self.spares[plane] += change
        self.spares_since[plane] = time

    def _producing(self, kind: Any, count: int) -> float:
        rate = _per_type(self.schedule.production_rate, kind, "production rate")
        return 0.0 if rate is None else count / rate

    def _produce(self, time: float, kind: Any, count: int, purpose: str) -> float:
        begin = max(self.line[kind], time + self.policy.procurement_years)
        ready = begin + self._producing(kind, count)
        self.line[kind] = ready
        self.history.production.append(Production(kind, begin, ready, count, purpose))
        return ready

    def _slot(self, ready: float) -> float:
        """The earliest launch at or after a time, at the launch cadence."""
        interval = self.schedule.launch_interval
        time = ready
        if interval > 0:
            for booked in self.booked:
                if booked - interval < time < booked + interval:
                    time = booked + interval
        return time

    def _capacity(self, kind: Any) -> float:
        capacity = math.inf
        if self.schedule.launch_capacity is not None:
            capacity = self.schedule.launch_capacity
        if self.schedule.launch_mass_capacity is not None:
            mass = self.deployment.masses.get(kind)
            if mass is not None:
                capacity = min(
                    capacity, math.floor(self.schedule.launch_mass_capacity / mass)
                )
        return capacity

    def _book(self, plane: int, kind: Any, ready: float, item: tuple) -> None:
        """Books an item on a launch of a plane: an open one, or a new one."""
        for booking in self.open.get(plane, []):
            if booking.time >= ready and len(booking.items) < self._capacity(kind):
                booking.items.append(item)
                return
        time = self._slot(ready + self.policy.batch_years)
        order = self.launches_of_plane.get(plane, 0)
        self.launches_of_plane[plane] = order + 1
        booking = _Booking(time, plane, kind, [item], order)
        self.open.setdefault(plane, []).append(booking)
        self.booked = sorted(self.booked + [time])
        self._push(time, "booking", booking)


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
        phases (Sequence[Phase]): The phases (see `MissionHistory.phases`).
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
