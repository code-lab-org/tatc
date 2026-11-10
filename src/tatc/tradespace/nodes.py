"""
Nodes of design spaces: design variables, fixed, derived, and uncertain
values, choices, and repetitions (see `tatc.tradespace.DesignSpace`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from .values import Scope


class Infeasible(Exception):
    """
    Raised while deciding a design that violates a constraint of its design
    space or has a variable without feasible values.

    Args:
        message (str): The reason.
        violation (float): The amount by which the design violates its
            constraint (positive), as for constraint-aware optimizers: the
            value of a constraint that gives one, or 1.
    """

    def __init__(self, message: str = "", violation: float = 1.0):
        super().__init__(message)
        self.violation = violation


# ---------------------------------------------------------------------------
# nodes
# ---------------------------------------------------------------------------


class Node:
    """
    Base class of the nodes of a design space.
    """


Bound = float | Callable[[Scope], float]
"""A bound of a variable: a value, or a callable of the values decided so far."""


def _resolve(bound: Any, scope: Scope) -> Any:
    """
    Resolves a bound that may be a callable of the values decided so far.
    """
    return bound(scope) if callable(bound) else bound


@dataclass(frozen=True)
class _Domain:
    """
    The values a variable may take where it is decided: a finite list of
    levels, or a continuous interval. The kind (`"real"`, `"integer"`, or
    `"categorical"`) and the step between levels inform optimizers.
    """

    levels: list | None = None
    low: float = 0.0
    high: float = 1.0
    log: bool = False
    kind: str = "categorical"
    step: float | None = None

    def from_unit(self, u: float) -> Any:
        """
        Maps a value in [0, 1] to a value of the domain: a level (each
        level spans an equal part of [0, 1]) or a point of the interval
        (uniform, or uniform in its logarithm).
        """
        u = min(max(float(u), 0.0), 1.0)
        if self.levels is not None:
            return self.levels[min(int(u * len(self.levels)), len(self.levels) - 1)]
        if self.log:
            return float(
                math.exp(
                    math.log(self.low) + u * (math.log(self.high) - math.log(self.low))
                )
            )
        return float(self.low + u * (self.high - self.low))

    def to_unit(self, value: Any) -> float:
        """
        Maps a value of the domain to [0, 1] (the center of a level's part),
        the inverse of `from_unit`.
        """
        if self.levels is not None:
            return (self.levels.index(value) + 0.5) / len(self.levels)
        if self.high == self.low:
            return 0.5
        if self.log:
            return (math.log(value) - math.log(self.low)) / (
                math.log(self.high) - math.log(self.low)
            )
        return (value - self.low) / (self.high - self.low)

    def suggest(self, trial: Any, name: str) -> Any:
        """
        Suggests a value of the domain with an Optuna trial.
        """
        if self.kind == "categorical":
            return trial.suggest_categorical(name, self.levels)
        if self.kind == "integer":
            return trial.suggest_int(
                name, self.levels[0], self.levels[-1], step=int(self.step)
            )
        if self.levels is not None:
            value = trial.suggest_float(
                name, self.levels[0], self.levels[-1], step=self.step
            )
            # the level nearest the suggestion, without rounding errors
            index = round((value - self.levels[0]) / self.step)
            return self.levels[min(max(index, 0), len(self.levels) - 1)]
        return trial.suggest_float(name, self.low, self.high, log=self.log)


class Variable(Node):
    """
    Base class of design variables: the leaves of a design space that are
    searched by enumeration, sampling, or optimization.
    """

    def domain(self, scope: Scope) -> _Domain:
        """
        Gets the values this variable may take, given the values decided so
        far.

        Args:
            scope (Scope): The values decided so far.

        Returns:
            _Domain: the domain
        """
        raise NotImplementedError


def _levels(low: float, high: float, step: float, integer: bool) -> list:
    """
    Gets the levels from `low` to `high` (inclusive, within a small
    tolerance) at intervals of `step`.
    """
    if step <= 0:
        raise ValueError("step must be positive")
    count = math.floor((high - low) / step + 1e-9) + 1
    if integer:
        return [int(low + k * step) for k in range(max(count, 0))]
    return [round(low + k * step, 12) for k in range(max(count, 0))]


class Real(Variable):
    """
    A real-valued design variable on `[low, high]`: continuous, or discrete
    at intervals of `step` from `low`. A continuous variable cannot be
    enumerated and makes every sampled design distinct, so give it a step
    unless a search needs the continuum.

    Args:
        low (float | Callable): The lower bound, or a callable of the values
            decided so far.
        high (float | Callable): The upper bound, or a callable of the values
            decided so far.
        step (float | None): The interval between levels, or `None` for a
            continuous variable.
        log (bool): `True` to sample a continuous variable uniformly in its
            logarithm (requires a positive lower bound).
    """

    def __init__(
        self, low: Bound, high: Bound, step: float | None = None, log: bool = False
    ):
        if step is not None and log:
            raise ValueError("a discrete variable cannot be logarithmic")
        self.low, self.high, self.step, self.log = low, high, step, log

    def domain(self, scope: Scope) -> _Domain:
        low, high = float(_resolve(self.low, scope)), float(_resolve(self.high, scope))
        if high < low:
            raise Infeasible(f"empty range [{low}, {high}]")
        if self.step is not None:
            levels = _levels(low, high, self.step, integer=False)
            return _Domain(levels=levels, kind="real", step=self.step)
        if self.log and low <= 0:
            raise ValueError("a logarithmic variable requires a positive lower bound")
        return _Domain(low=low, high=high, log=self.log, kind="real")

    def __repr__(self) -> str:
        return f"Real({self.low!r}, {self.high!r}, step={self.step!r})"


class Integer(Variable):
    """
    An integer design variable from `low` to `high` (inclusive) at intervals
    of `step`.

    Args:
        low (int | Callable): The lower bound, or a callable of the values
            decided so far.
        high (int | Callable): The upper bound, or a callable of the values
            decided so far.
        step (int): The interval between levels.
    """

    def __init__(self, low: Bound, high: Bound, step: int = 1):
        self.low, self.high, self.step = low, high, step

    def domain(self, scope: Scope) -> _Domain:
        low, high = int(_resolve(self.low, scope)), int(_resolve(self.high, scope))
        if high < low:
            raise Infeasible(f"empty range [{low}, {high}]")
        levels = _levels(low, high, self.step, integer=True)
        return _Domain(levels=levels, kind="integer", step=self.step)

    def __repr__(self) -> str:
        return f"Integer({self.low!r}, {self.high!r}, step={self.step!r})"


class Categorical(Variable):
    """
    A design variable that takes one of a list of unordered options (for
    example, the identifiers of instruments in the catalog). Options should
    be numbers, strings, or booleans, which identify designs (see
    `Design.key`).

    Args:
        options (Sequence | Callable): The options, or a callable of the
            values decided so far.
    """

    def __init__(self, options: Sequence | Callable[[Scope], Sequence]):
        if not callable(options) and len(options) == 0:
            raise ValueError("a categorical variable requires at least one option")
        self.options = options if callable(options) else list(options)

    def domain(self, scope: Scope) -> _Domain:
        options = list(_resolve(self.options, scope))
        if len(options) == 0:
            raise Infeasible("no options")
        return _Domain(levels=options)

    def __repr__(self) -> str:
        return f"Categorical({self.options!r})"


class Fixed(Node):
    """
    A fixed value. Plain values given as children of a design space are
    wrapped as fixed values.

    Args:
        value (Any): The value.
    """

    def __init__(self, value: Any):
        self.value = value

    def __repr__(self) -> str:
        return f"Fixed({self.value!r})"


class Derived(Node):
    """
    A value computed from the values decided so far (for example, the
    equator crossing time of an orbit plane from its index). Derived values
    are not searched, and do not identify designs (see `Design.key`). A
    derived value that depends on an uncertain parameter is computed from its
    nominal value (`NaN` or `None`) until a design is realized (see
    `Design.realize`).

    Args:
        function (Callable): The function of the values decided so far.
    """

    def __init__(self, function: Callable[[Scope], Any]):
        self.function = function

    def __repr__(self) -> str:
        return f"Derived({self.function!r})"


class Uncertain(Node):
    """
    An uncertain (nuisance) parameter: a value that is not designed but
    changes the outcome of a design, such as the phase of satellites in
    their orbit planes or the epoch of an analysis. It is drawn for each
    replicate of a design (see `Design.realize`), does not identify designs
    (see `Design.key`), and takes its nominal value until then.

    Args:
        sampler (Callable): A function that draws a value from a
            `numpy.random.Generator`.
        nominal (Any): The value before a design is realized.
    """

    def __init__(
        self, sampler: Callable[[np.random.Generator], Any], nominal: Any = math.nan
    ):
        self.sampler, self.nominal = sampler, nominal

    @classmethod
    def uniform(cls, low: float, high: float) -> Uncertain:
        """
        An uncertain parameter distributed uniformly on `[low, high)`.
        """
        return cls(lambda rng: float(rng.uniform(low, high)))

    @classmethod
    def normal(cls, mean: float, std: float) -> Uncertain:
        """
        An uncertain parameter distributed normally.
        """
        return cls(lambda rng: float(rng.normal(mean, std)))

    @classmethod
    def choice(cls, options: Sequence, p: Sequence[float] | None = None) -> Uncertain:
        """
        An uncertain parameter that takes one of `options` (with
        probabilities `p`, or equally likely).
        """
        options = list(options)
        return cls(lambda rng: options[int(rng.choice(len(options), p=p))], None)

    def __repr__(self) -> str:
        return f"Uncertain({self.sampler!r})"


def _as_node(value: Any) -> Node:
    """
    Wraps a plain value as a fixed value.
    """
    return value if isinstance(value, Node) else Fixed(value)


class Choice(Node):
    """
    A choice of exactly one of several alternatives (branches), each a node
    (typically a `DesignSpace`). The nodes of the other branches are
    inactive. The value of a choice gives the name of the chosen branch as
    `choice`, and the values of the branch as its attributes.

    Args:
        **branches (Node): The alternatives, by name.
    """

    def __init__(self, **branches: Any):
        if len(branches) == 0:
            raise ValueError("a choice requires at least one branch")
        self.branches = {name: _as_node(node) for name, node in branches.items()}

    def __repr__(self) -> str:
        return f"Choice({', '.join(self.branches)})"


class Repeat(Node):
    """
    A node repeated a (variable) number of times, such as the orbit planes of
    a constellation, each with its own design. The value of a repetition is a
    sequence of the values of its items, whose length is also given as
    `count`. Within an item, `index` gives its position.

    The items of an unordered repetition (`ordered=False`) are exchangeable:
    designs whose items are the same up to their order are the same design,
    whose items are sorted by their decisions (so their values must not
    depend on `index`).

    Args:
        item (Node): The repeated node.
        count (int | Integer): The number of items, fixed or a design variable.
        ordered (bool): `False` if designs that differ only in the order of
            their items are the same.
    """

    def __init__(self, item: Any, count: int | Integer, ordered: bool = True):
        self.item = _as_node(item)
        if isinstance(count, int) and count < 0:
            raise ValueError("count must not be negative")
        self.count, self.ordered = count, ordered

    def __repr__(self) -> str:
        return f"Repeat({self.item!r}, count={self.count!r}, ordered={self.ordered})"
