"""
The discrete-event simulation of missions (see `Deployment.simulate`): the
deployment of a constellation with the failures of its satellites and
launches, the restoration of failed satellites by a policy, and the
decisions of the rules of its stages, into the history of the mission (see
`MissionHistory`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import bisect
import heapq
import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from functools import partial
from typing import TYPE_CHECKING, Any

import numpy as np

from .events import Launch, Production
from .failures import FailureModel, RestorationPolicy
from .mission import Decision, Failure, MissionHistory, Restoration

if TYPE_CHECKING:
    from .deployment import Deployment


@dataclass
class _Slot:
    """
    The state of a position of a plane: the generation of its satellite,
    when it began to operate (or `None`), when it failed, and whether it
    awaits restoration.
    """

    generation: int = -1
    up: float | None = None
    down: float | None = None
    pending: bool = False


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
        seed: int | Sequence[int] | None,
        uniforms: Callable[[tuple, int], np.ndarray] | None = None,
        stress: float = 1.0,
        launch_stress: float = 1.0,
    ):
        self.deployment, self.failures, self.policy = deployment, failures, policy
        self.seed, self.uniforms = seed, uniforms
        self.stress, self.launch_stress = stress, launch_stress
        self.log_weight = 0.0
        self.schedule = deployment.schedule
        self.end = deployment.end
        self.history = MissionHistory(
            deployment, production=list(deployment.production)
        )
        self.events: list = []
        self.counter = itertools.count()
        self.slots = {
            (plane, position): _Slot()
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
        self.decided: dict[int, float] = {}

    # --- randomness

    def _rng(self, *key: int) -> np.random.Generator:
        seed = list(self.seed) if isinstance(self.seed, Sequence) else [self.seed]
        return np.random.default_rng([*seed, *key])

    def _life(
        self, plane: int, position: int, generation: int, start: float = 0.0
    ) -> float:
        if self.seed is None:
            return math.inf
        life = self.failures.life_of(self.deployment.planes[plane][0])
        if life is None:
            return math.inf
        key = (0, plane, position, generation)
        if self.stress != 1.0:
            uniform = (
                self.uniforms(key, life.dimensions)
                if self.uniforms is not None
                else self._rng(*key).random(life.dimensions)
            )
            value, log_ratio = life.stressed(uniform, self.stress, self.end - start)
            self.log_weight += log_ratio
            return value
        if self.uniforms is not None:
            return life.from_uniform(self.uniforms(key, life.dimensions))
        return life.sample(self._rng(*key))

    def _launch_fails(self, *key: int) -> bool:
        if self.seed is None or self.failures.launch == 0:
            return False
        uniform = (
            self.uniforms((1, *key), 1)[0]
            if self.uniforms is not None
            else self._rng(1, *key).random()
        )
        p = self.failures.launch
        if self.launch_stress == 1.0:
            return bool(uniform < p)
        stressed = min(p * self.launch_stress, 0.5)
        failed = bool(uniform < stressed)
        ratio = p / stressed if failed else (1 - p) / (1 - stressed)
        self.log_weight += math.log(ratio)
        return failed

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
                ready = begin + self.schedule.production_years(kind, count)
                self.history.production.append(
                    Production(kind, begin, ready, count, "spare")
                )
                self.line[kind] = ready
                self.ground.setdefault(kind, 0)
                self._push(ready, "stock", kind, count)
        for launch in self.deployment.launches:
            self._push(launch.time, "planned", launch)
        for index, stage in enumerate(self.deployment.stages):
            for time in stage.times:
                self._push(self.deployment.begin + time, "decide", index)
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
            if slot.up is not None:
                self._close(plane, position, self.end)
        for plane, since in self.spares_since.items():
            if self.spares[plane] > 0:
                self.history.parked.append((since, self.end, self.spares[plane]))
        self.history.weight = math.exp(self.log_weight)
        self.history.stages = [
            self.decided.get(index) for index in range(len(self.deployment.stages))
        ]
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
                if slot.up is not None:
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

    def _on_decide(self, time: float, index: int) -> None:
        if index in self.decided:
            return  # deployed at an earlier decision
        stage = self.deployment.stages[index]
        if stage.rule is not None:
            operating = [
                position for position, slot in self.slots.items() if slot.up is not None
            ]
            capability = self.deployment.capability
            decision = Decision(
                time,
                index,
                len(operating),
                self.history,
                tuple(sorted(self.decided)),
                (
                    None
                    if capability is None
                    else partial(capability.at, self.deployment.planes, operating)
                ),
            )
            if not stage.rule(decision):
                return
        self.decided[index] = time
        self._deploy(time, stage)

    def _on_fail(self, time: float, plane: int, position: int, generation: int) -> None:
        slot = self.slots[plane, position]
        if slot.generation != generation or slot.up is None:
            return  # a satellite replaced before it failed
        self._close(plane, position, time)
        self._fail(time, plane, position, "life")

    def _fail(self, time: float, plane: int, position: int, cause: str) -> None:
        slot = self.slots[plane, position]
        slot.down = time
        self.history.failures.append(Failure(time, plane, position, cause))
        if self.policy.restore and not slot.pending:
            slot.pending = True
            self._push(time + self.policy.detection_years, "request", plane, position)

    def _on_request(self, time: float, plane: int, position: int) -> None:
        slot = self.slots[plane, position]
        if slot.up is not None:  # restored by a planned replacement
            slot.pending = False
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
        if slot.up is not None:
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
                slot.pending = False
                self.history.failures.append(
                    Failure(time, booking.plane, position, "launch")
                )
                if self.policy.restore:
                    slot.pending = True
                    self._push(
                        time + max(self.policy.detection_years, 1e-6),
                        "request",
                        booking.plane,
                        position,
                    )
            elif slot.up is not None:
                self._park(time, booking.plane, 1)  # not needed: a spare in orbit
            else:
                self._restore(time, booking.plane, position, source)

    # --- actions

    def _up(self, time: float, plane: int, position: int) -> None:
        slot = self.slots[plane, position]
        slot.generation += 1
        slot.up, slot.pending = time, False
        life = self._life(plane, position, slot.generation, time)
        if time + life < self.end:
            self._push(time + life, "fail", plane, position, slot.generation)

    def _close(self, plane: int, position: int, time: float) -> None:
        slot = self.slots[plane, position]
        if time > slot.up:
            self.history.operating.setdefault((plane, position), []).append(
                (slot.up, time)
            )
        slot.up = None

    def _restore(self, time: float, plane: int, position: int, source: str) -> None:
        slot = self.slots[plane, position]
        self.history.restorations.append(
            Restoration(plane, position, slot.down, time, source)
        )
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

    def _deploy(self, time: float, stage: Any) -> None:
        """
        Deploys the planes of a stage: produced after its lead time (on their
        type's line), launched at the next launch slots, and replaced at the
        end of their design life.
        """
        for plane in stage.planes:
            kind, size, spares = self.deployment.planes[plane]
            count = size + spares
            if count == 0:
                continue
            begin = max(self.line[kind], time + stage.lead_years)
            ready = begin + self.schedule.production_years(kind, count)
            self.line[kind] = ready
            self.history.production.append(Production(kind, begin, ready, count))
            position = 0
            for carried in self.schedule.split(
                count, self.deployment.masses.get(kind), self._plane(plane)
            ):
                at = self._slot(ready)
                bisect.insort(self.booked, at)
                launch = Launch(at, plane, kind, position, carried)
                self._push(at, "planned", launch)
                for production, replacement in self.schedule.replacements(
                    launch, self.end
                ):
                    self.history.production.append(production)
                    if replacement is not None:
                        bisect.insort(self.booked, replacement.time)
                        self._push(replacement.time, "planned", replacement)
                position += carried

    def _produce(self, time: float, kind: Any, count: int, purpose: str) -> float:
        begin = max(self.line[kind], time + self.policy.procurement_years)
        ready = begin + self.schedule.production_years(kind, count)
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

    def _plane(self, plane: int) -> tuple[int, Any]:
        """The index and values of a plane (see `DeploymentSchedule`)."""
        values = self.deployment.plane_values
        return plane, values[plane] if plane < len(values) else None

    def _capacity(self, plane: int, kind: Any) -> float:
        return self.schedule.per_launch(
            self.deployment.masses.get(kind), self._plane(plane)
        )

    def _book(self, plane: int, kind: Any, ready: float, item: tuple) -> None:
        """Books an item on a launch of a plane: an open one, or a new one."""
        for booking in self.open.get(plane, []):
            if booking.time >= ready and len(booking.items) < self._capacity(
                plane, kind
            ):
                booking.items.append(item)
                return
        time = self._slot(ready + self.policy.batch_years)
        order = self.launches_of_plane.get(plane, 0)
        self.launches_of_plane[plane] = order + 1
        booking = _Booking(time, plane, kind, [item], order)
        self.open.setdefault(plane, []).append(booking)
        bisect.insort(self.booked, time)
        self._push(time, "booking", booking)
