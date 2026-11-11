"""
Failures of satellites and launches, and the policies that restore (or do
not restore) the satellites that fail, for simulations of missions (see
`Deployment.simulate`).

The life of a satellite (from its launch) follows a distribution whose
failure rate may vary over its lifecycle: constant (`Exponential`),
decreasing in infancy or increasing with wear (`Weibull`), both (`Bathtub`,
the earliest of competing causes), or given by intervals of age
(`Piecewise`). A `FailureModel` gives the life of each type of satellite and
the probability that a launch fails, and a `RestorationPolicy` whether and
how failed satellites are restored: by spares in orbit, by spares on the
ground, or by satellites built when needed. A policy is a decision of a
design with a fragment of a design space (`policy_space`), whose values
configure it (see `configure`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, fields
from typing import Any

import numpy as np

from .nodes import Categorical
from .space import DesignSpace
from .values import Scope


class Lifetime:
    """
    Base class of distributions of the life (years) of a satellite from its
    launch.
    """

    def cumulative_hazard(self, age: np.ndarray) -> np.ndarray:
        """
        Gets the cumulative hazard at ages (the negative logarithm of the
        survival).

        Args:
            age (numpy.ndarray): The ages (years).

        Returns:
            numpy.ndarray: the cumulative hazard
        """
        raise NotImplementedError

    def hazard(self, age: np.ndarray) -> np.ndarray:
        """
        Gets the hazard (failure rate, per year) at ages.

        Args:
            age (numpy.ndarray): The ages (years).

        Returns:
            numpy.ndarray: the hazard
        """
        raise NotImplementedError

    def survival(self, age: np.ndarray) -> np.ndarray:
        """
        Gets the probability of surviving to ages.

        Args:
            age (numpy.ndarray): The ages (years).

        Returns:
            numpy.ndarray: the survival
        """
        return np.exp(-self.cumulative_hazard(np.asarray(age, dtype=float)))

    def sample(self, rng: np.random.Generator) -> float:
        """
        Draws a life (years, infinite if the satellite never fails).

        Args:
            rng (numpy.random.Generator): The random number generator.

        Returns:
            float: the life
        """
        raise NotImplementedError


@dataclass(frozen=True)
class Exponential(Lifetime):
    """
    A life of constant failure rate (random failures).

    Args:
        rate (float): The failure rate (per year).
    """

    rate: float

    def __post_init__(self):
        if self.rate < 0:
            raise ValueError("rate must not be negative")

    def cumulative_hazard(self, age):
        return self.rate * np.asarray(age, dtype=float)

    def hazard(self, age):
        return np.full_like(np.asarray(age, dtype=float), self.rate)

    def sample(self, rng):
        return math.inf if self.rate == 0 else float(rng.exponential(1 / self.rate))


@dataclass(frozen=True)
class Weibull(Lifetime):
    """
    A life of a Weibull distribution: a failure rate that decreases with age
    (infant mortality) if its shape is less than 1, is constant if 1, and
    increases (wear-out) if greater.

    Args:
        shape (float): The shape.
        scale (float): The scale (years), the age at which 63% have failed.
    """

    shape: float
    scale: float

    def __post_init__(self):
        if self.shape <= 0 or self.scale <= 0:
            raise ValueError("shape and scale must be positive")

    def cumulative_hazard(self, age):
        return (np.asarray(age, dtype=float) / self.scale) ** self.shape

    def hazard(self, age):
        age = np.asarray(age, dtype=float)
        with np.errstate(divide="ignore"):
            return (self.shape / self.scale) * (age / self.scale) ** (self.shape - 1)

    def sample(self, rng):
        return float(self.scale * rng.weibull(self.shape))


@dataclass(frozen=True)
class Bathtub(Lifetime):
    """
    A life of competing causes of failure (the earliest of their lives), such
    as the bathtub curve of infant mortality, random failures, and wear-out:
    its hazard is the sum of theirs.

    Args:
        infant (Lifetime | None): Failures early in life.
        random (Lifetime | None): Failures at random.
        wear_out (Lifetime | None): Failures late in life.
    """

    infant: Lifetime | None = None
    random: Lifetime | None = None
    wear_out: Lifetime | None = None

    @property
    def causes(self) -> list[Lifetime]:
        """list[Lifetime]: The causes of failure."""
        return [c for c in (self.infant, self.random, self.wear_out) if c is not None]

    def cumulative_hazard(self, age):
        age = np.asarray(age, dtype=float)
        return sum((c.cumulative_hazard(age) for c in self.causes), np.zeros_like(age))

    def hazard(self, age):
        age = np.asarray(age, dtype=float)
        return sum((c.hazard(age) for c in self.causes), np.zeros_like(age))

    def sample(self, rng):
        return min((c.sample(rng) for c in self.causes), default=math.inf)


@dataclass(frozen=True)
class Piecewise(Lifetime):
    """
    A life of failure rates given by intervals of age: `rates[i]` from
    `breaks[i]` to `breaks[i + 1]` (the last beyond the last break), as from
    reliability data.

    Args:
        breaks (Sequence[float]): The ages (years) at which the intervals
            begin, from 0, increasing.
        rates (Sequence[float]): The failure rate (per year) of each.
    """

    breaks: Sequence[float]
    rates: Sequence[float]

    def __post_init__(self):
        breaks, rates = list(self.breaks), list(self.rates)
        if len(breaks) == 0 or breaks[0] != 0:
            raise ValueError("breaks must begin at 0")
        if len(breaks) != len(rates):
            raise ValueError("a rate is required for each interval")
        if any(b <= a for a, b in zip(breaks[:-1], breaks[1:])):
            raise ValueError("breaks must increase")
        if any(rate < 0 for rate in rates):
            raise ValueError("rates must not be negative")
        object.__setattr__(self, "breaks", tuple(float(b) for b in breaks))
        object.__setattr__(self, "rates", tuple(float(r) for r in rates))

    def _ends(self) -> np.ndarray:
        return np.append(np.asarray(self.breaks[1:]), math.inf)

    def cumulative_hazard(self, age):
        age = np.asarray(age, dtype=float)
        starts, ends = np.asarray(self.breaks), self._ends()
        within = np.clip(age[..., None] - starts, 0, ends - starts)
        return np.sum(within * np.asarray(self.rates), axis=-1)

    def hazard(self, age):
        age = np.asarray(age, dtype=float)
        index = np.searchsorted(np.asarray(self.breaks), age, side="right") - 1
        return np.asarray(self.rates)[np.clip(index, 0, len(self.rates) - 1)]

    def sample(self, rng):
        # the age at which the cumulative hazard reaches a unit exponential draw
        target = float(rng.exponential(1.0))
        cumulative = 0.0
        for start, end, rate in zip(self.breaks, self._ends(), self.rates):
            span = end - start
            if rate > 0 and cumulative + rate * span >= target:
                return start + (target - cumulative) / rate
            cumulative += rate * span if math.isfinite(span) else 0.0
        return math.inf


@dataclass(frozen=True)
class FailureModel:
    """
    The failures of a mission: the life of the satellites (of every type, or
    of each, by name; `None` if they do not fail) and the probability that
    a launch fails, losing its satellites.

    Args:
        life (Lifetime | Mapping[Any, Lifetime] | None): The life of the
            satellites.
        launch (float): The probability that a launch fails.
    """

    life: Lifetime | Mapping[Any, Lifetime | None] | None = None
    launch: float = 0.0

    def __post_init__(self):
        if not 0 <= self.launch <= 1:
            raise ValueError("launch must be a probability")

    def life_of(self, kind: Any) -> Lifetime | None:
        """
        Gets the life of a type of satellite.

        Args:
            kind (Any): The type.

        Returns:
            Lifetime | None: the life, or `None` if it does not fail
        """
        if isinstance(self.life, Mapping):
            if kind not in self.life:
                raise ValueError(f"no life of satellites of kind {kind!r}")
            return self.life[kind]
        return self.life


@dataclass(frozen=True)
class RestorationPolicy:
    """
    A policy for satellites that fail: to fly the others (`restore=False`,
    degrading the constellation until the satellites are replaced at the end
    of their design life, if they are), or to restore them, from the fastest
    source available when the failure is known (`detection_years` after it):

    * a spare in orbit in the failed satellite's plane (see the `spares` of
      `DeploymentSchedule.plan`), which takes its place after
      `drift_years`;
    * a spare on the ground (`ground_spares` of each type, produced after
      the first satellites and stored), launched at the next launch slot;
    * a satellite produced when needed, `procurement_years` after it is
      ordered (for long-lead items) and on its type's production line (after
      the work before it), then launched at the next launch slot.

    A restoration of a plane joins a launch already booked for the plane, if
    it has room (see the launch capacity of `DeploymentSchedule`); otherwise
    it books one at the next launch slot after waiting `batch_years` for
    others. With `replenish`, a spare that is used is produced again (and,
    in orbit, launched to its plane), procured and produced as a satellite
    produced when needed, which shares the production line with
    restorations.

    Args:
        restore (bool): `True` to restore failed satellites.
        detection_years (float): The time to know of a failure.
        drift_years (float): The time for a spare in orbit to take the place
            of a failed satellite.
        procurement_years (float): The time from ordering a satellite to the
            start of its production.
        ground_spares (int | Mapping[Any, int]): The spares on the ground of
            each type of satellite.
        replenish (bool): `True` to replace the spares that are used.
        batch_years (float): The time that restorations of a plane wait to
            share a launch.
    """

    restore: bool = True
    detection_years: float = 0.1
    drift_years: float = 0.25
    procurement_years: float = 0.0
    ground_spares: int | Mapping[Any, int] = 0
    replenish: bool = False
    batch_years: float = 0.0

    def __post_init__(self):
        for name in (
            "detection_years",
            "drift_years",
            "procurement_years",
            "batch_years",
        ):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must not be negative")
        values = (
            self.ground_spares.values()
            if isinstance(self.ground_spares, Mapping)
            else [self.ground_spares]
        )
        if any(value < 0 for value in values):
            raise ValueError("ground_spares must not be negative")

    def ground_spares_of(self, kind: Any) -> int:
        """
        Gets the spares on the ground of a type of satellite.

        Args:
            kind (Any): The type.

        Returns:
            int: the number of spares
        """
        if isinstance(self.ground_spares, Mapping):
            return int(self.ground_spares.get(kind, 0))
        return int(self.ground_spares)


def _spares_restore(d: Scope) -> bool:
    """Whether spares are only kept by a policy that restores failures."""
    return bool(d.restore) or (d.orbit_spares == 0 and d.ground_spares == 0)


def policy_space(
    *,
    restore: Any = None,
    orbit_spares: Any = 0,
    ground_spares: Any = 0,
    replenish: Any = False,
    require_restoration: bool = True,
    constraints: Sequence[Callable[[Scope], Any]] = (),
    **others: Any,
) -> DesignSpace:
    """
    Gets a fragment of a design space of policies for failed satellites (see
    `RestorationPolicy`), for trades of their cost against their
    robustness: whether to `restore` failed satellites (by default a
    variable, `False` or `True`), the spares in orbit in each plane
    (`orbit_spares`, which a deployment plans, see the `spares` of
    `DeploymentSchedule.plan`), the spares on the ground (`ground_spares`),
    whether to `replenish` used spares, and other fields of the policy (such
    as `drift_years`). Each is a node of a design space (such as a design
    variable) or a value.

    A design's policy is configured from the values of the fragment, and its
    deployment planned with their spares in orbit::

        policy = configure(BASE_POLICY, d.policy)
        deployment = SCHEDULE.plan(d.layout, spares=d.policy.orbit_spares)
        history = deployment.simulate(FAILURES, policy, d.failure_seed)

    Spares without restoration would never be used, so a constraint excludes
    them (unless `require_restoration` is `False`).

    Args:
        restore (Any): Whether to restore failed satellites, by default a
            variable.
        orbit_spares (Any): The spares in orbit in each plane.
        ground_spares (Any): The spares on the ground.
        replenish (Any): Whether to replenish used spares.
        require_restoration (bool): `True` to exclude spares without
            restoration.
        constraints (Sequence[Callable]): Other constraints of the fragment.
        **others (Any): Other fields of `RestorationPolicy`.

    Returns:
        DesignSpace: the fragment
    """
    names = {field.name for field in fields(RestorationPolicy)}
    unknown = set(others) - names
    if unknown:
        raise ValueError(f"no fields of RestorationPolicy named {sorted(unknown)}")
    children = {
        "restore": Categorical([False, True]) if restore is None else restore,
        "orbit_spares": orbit_spares,
        "ground_spares": ground_spares,
        "replenish": replenish,
        **others,
    }
    rules = ([_spares_restore] if require_restoration else []) + list(constraints)
    return DesignSpace(constraints=rules, **children)
