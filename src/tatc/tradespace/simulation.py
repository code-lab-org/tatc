"""
The discrete-event simulation of missions (see `Deployment.simulate`): the
deployment of a constellation with the failures of its satellites and
launches, the restoration of failed satellites by its policy, and the
decisions of the rules of its stages, into the history of the mission (see
`MissionHistory`); and the samplers of many histories of a replicate,
independent or stratified (see `Deployment.simulate_many`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import heapq
import itertools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from functools import partial
from typing import TYPE_CHECKING, Any, Literal, get_args

import numpy as np

from .events import Launch, _Manifest
from .failures import FailureModel
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


class _Draws:
    """
    The random numbers of the failures of a history: uniform numbers of each
    draw (by its key), of a sampler of the histories of a replicate (such as
    `_Sampler`) or from the seed of the history, from which the lives of
    satellites and the failures of launches are drawn, stressed (with their
    likelihood ratio) or not.
    """

    def __init__(
        self,
        seed: int | Sequence[int] | None,
        uniforms: Callable[[tuple, int], np.ndarray] | None = None,
        stress: float = 1.0,
        launch_stress: float = 1.0,
    ):
        self.seed, self.stress, self.launch_stress = seed, stress, launch_stress
        self.uniforms = self._from_seed if uniforms is None else uniforms
        self.log_weight = 0.0

    def _from_seed(self, key: tuple, dimensions: int) -> np.ndarray:
        """The uniform numbers of a draw, from the seed and its key."""
        seed = list(self.seed) if isinstance(self.seed, Sequence) else [self.seed]
        return np.random.default_rng([*seed, *key]).random(dimensions)

    @property
    def weight(self) -> float:
        """float: The likelihood ratio of the draws (1 without stress)."""
        return math.exp(self.log_weight)

    def life(self, life: Any, key: tuple, remaining: float) -> float:
        """
        Draws a life (years, infinite without a seed or a life), stressed
        over the time remaining in the mission, if the draws are.
        """
        if self.seed is None or life is None:
            return math.inf
        uniform = self.uniforms(key, life.dimensions)
        if self.stress == 1.0:
            return life.from_uniform(uniform)
        value, log_ratio = life.stressed(uniform, self.stress, remaining)
        self.log_weight += log_ratio
        return value

    def launch_fails(self, probability: float, key: tuple) -> bool:
        """Draws whether a launch fails (stressed, if the draws are)."""
        if self.seed is None or probability == 0:
            return False
        uniform = self.uniforms(key, 1)[0]
        if self.launch_stress == 1.0:
            return bool(uniform < probability)
        stressed = min(probability * self.launch_stress, 0.5)
        failed = bool(uniform < stressed)
        ratio = probability / stressed if failed else (1 - probability) / (1 - stressed)
        self.log_weight += math.log(ratio)
        return failed


class _Sampler:
    """
    The uniform numbers of the failures of a history of a replicate,
    stratified over the histories of the replicate (see
    `Deployment.simulate_many`): a Latin hypercube in which, for each draw
    (each life and launch), history `k` of `K` takes a uniform number in its
    stratum of a permutation of the `K` strata of [0, 1). Each draw is keyed
    by its key, so that the same draw of the same history of a replicate
    takes the same number in every design (common random numbers).

    The permutations of a draw are the same for every history of the
    replicate, so the samplers of its histories share them in `orders` (by
    the key and dimensions of the draw), computing each once rather than once
    for each history.
    """

    def __init__(
        self,
        seed: int,
        index: int,
        histories: int,
        orders: dict[tuple, list[np.ndarray]] | None = None,
    ):
        self.seed, self.index, self.histories = seed, index, histories
        self.orders = {} if orders is None else orders

    def __call__(self, key: tuple, dimensions: int) -> np.ndarray:
        orders = self.orders.get((key, dimensions))
        if orders is None:
            strata = np.random.default_rng([self.seed, 1, self.histories, *key])
            orders = [strata.permutation(self.histories) for _ in range(dimensions)]
            self.orders[key, dimensions] = orders
        jitter = np.random.default_rng([self.seed, 2, self.index, *key])
        offsets = jitter.random(dimensions)
        return np.array(
            [
                (order[self.index] + offset) / self.histories
                for order, offset in zip(orders, offsets)
            ]
        )


Sampling = Literal["random", "stratified"]
"""The sampling of the histories of a deployment (see
`Deployment.simulate_many`)."""

_SAMPLINGS = get_args(Sampling)
"""The samplings of the failures of histories (see `Deployment.simulate_many`)."""


def _check_sampling(
    histories: int,
    sampling: Sampling,
    stress: float,
    launch_stress: float | None = None,
) -> None:
    """
    Checks the sampling of histories (see `Deployment.simulate_many`),
    raising a `ValueError` if the histories are fewer than 1, the sampling is
    unknown, or a stress is not positive.
    """
    if histories < 1:
        raise ValueError("histories must be at least 1")
    if sampling not in _SAMPLINGS:
        raise ValueError(f"sampling must be one of {_SAMPLINGS}, not {sampling!r}")
    if stress <= 0 or (launch_stress is not None and launch_stress <= 0):
        raise ValueError("stress must be positive")


class _Simulation:
    """
    A discrete-event simulation of a mission (see `Deployment.simulate`),
    with the policy of its deployment.
    """

    # pylint: disable=too-many-instance-attributes

    def __init__(
        self,
        deployment: Deployment,
        model: FailureModel,
        seed: int | Sequence[int] | None,
        uniforms: Callable[[tuple, int], np.ndarray] | None = None,
        stress: float = 1.0,
        launch_stress: float = 1.0,
    ):
        self.deployment, self.model = deployment, model
        self.policy = deployment.policy
        self.draws = _Draws(seed, uniforms, stress, launch_stress)
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
        # the production lines and the launches booked, after the plan's
        self.manifest = _Manifest.resume(deployment, self.history.production)
        self.open: dict[int, list[_Booking]] = {}
        self.launches_of_plane: dict[int, int] = {}
        self.decided: dict[int, float] = {}

    # --- randomness

    def _life(
        self, plane: int, position: int, generation: int, start: float = 0.0
    ) -> float:
        life = self.model.life_of(self.deployment.planes[plane].kind)
        key = (0, plane, position, generation)
        return self.draws.life(life, key, self.end - start)

    def _launch_fails(self, *key: int) -> bool:
        return self.draws.launch_fails(self.model.launch, (1, *key))

    # --- events

    def _push(self, time: float, event: str, *payload: Any) -> None:
        heapq.heappush(self.events, (time, next(self.counter), event, payload))

    def run(self) -> MissionHistory:
        """Runs the simulation, returning the history of the mission."""
        # pylint: disable=too-many-branches
        for kind in self.deployment.development:
            count = self.policy.ground_spares_of(kind)
            if count > 0:
                ready = self.manifest.produce(0.0, kind, count, "spare")
                self.ground.setdefault(kind, 0)
                self._push(ready, "stock", kind, count)
        for launch in self.deployment.launches:
            self._push(launch.time, "planned", launch)
        for index, stage in enumerate(self.deployment.stages):
            for time in stage.times:
                self._push(self.deployment.begin + time, "decide", index)
        steps = 0
        while self.events:
            time, _, event, payload = heapq.heappop(self.events)
            if time >= self.end:
                break
            steps += 1
            if steps > 1_000_000:
                raise RuntimeError("the simulation does not progress")
            getattr(self, f"_on_{event}")(time, *payload)
        for (plane, position), slot in self.slots.items():
            if slot.up is not None:
                self._close(plane, position, self.end)
        for plane, since in self.spares_since.items():
            if self.spares[plane] > 0:
                self.history.parked.append((since, self.end, self.spares[plane]))
        self.history.weight = self.draws.weight
        self.history.decided = [
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
        built = self.deployment.planes[launch.plane].satellites
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
        # the spares in orbit of a deployment (planned replacements carry
        # none: spares in orbit do not age)
        spares = sum(position >= built for position in positions)
        if spares and not failed:
            self._park(time, launch.plane, spares)

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
                    # pylint: disable-next=protected-access
                    else partial(capability._at, self.deployment.planes, operating)
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
        kind = self.deployment.planes[plane].kind
        if self.spares[plane] > 0:
            self._park(time, plane, -1)
            self._push(
                time + self.policy.drift_years, "activate", plane, position, "orbit"
            )
            if self.policy.replenish:
                ready = self.manifest.produce(
                    time + self.policy.procurement_years, kind, 1, "spare"
                )
                self._book(plane, kind, ready, ("spare",))
        elif self.ground.get(kind, 0) > 0:
            self.ground[kind] -= 1
            self._book(plane, kind, time, ("slot", position, "ground"))
            if self.policy.replenish:
                ready = self.manifest.produce(
                    time + self.policy.procurement_years, kind, 1, "spare"
                )
                self._push(ready, "stock", kind, 1)
        else:
            ready = self.manifest.produce(
                time + self.policy.procurement_years, kind, 1, "restoration"
            )
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
        Deploys the planes of a stage, as a plan deploys the others (see
        `DeploymentSchedule.plan`): produced after its lead time (on their
        type's line), launched at the next launch slots, and replaced at the
        end of their design life.
        """
        start = time + stage.lead_years
        for launch in self.manifest.deploy(stage.planes, start, self.end):
            self._push(launch.time, "planned", launch)

    def _book(self, plane: int, kind: Any, ready: float, item: tuple) -> None:
        """
        Books an item that is ready at a time on a launch of a plane: an open
        one that has room and flies once the item is ready, but no later
        than a new one would (at the next launch slot after waiting the
        batch years), or else a new one.
        """
        latest = self.manifest.next_slot(ready + self.policy.batch_years)
        capacity = self.manifest.capacity(plane)
        for booking in self.open.get(plane, []):
            if ready <= booking.time <= latest and len(booking.items) < capacity:
                booking.items.append(item)
                return
        time = self.manifest.slot(ready + self.policy.batch_years)
        order = self.launches_of_plane.get(plane, 0)
        self.launches_of_plane[plane] = order + 1
        booking = _Booking(time, plane, kind, [item], order)
        self.open.setdefault(plane, []).append(booking)
        self._push(time, "booking", booking)
