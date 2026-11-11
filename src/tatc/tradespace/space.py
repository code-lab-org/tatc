"""
Design spaces: declarative trees of design variables from which designs
(architectures) are enumerated or sampled for trade studies.

A design space is a tree of nodes. Its leaves are design variables
(`Real`, `Integer`, `Categorical`), fixed values, values derived from other
values (`Derived`), and uncertain (nuisance) parameters (`Uncertain`), which
are drawn for each replicate of a design rather than searched. Its branches
group nodes (`DesignSpace`), choose one of several alternatives (`Choice`),
or repeat a node a variable number of times (`Repeat`), which makes the
space hierarchical: for example, a constellation of a variable number of
orbit planes, each with its own design.

Nodes are decided in the order they are declared. A callable (a bound of a
variable, a derived value, or a constraint) is given a view of the values
decided so far, in which names resolve to the innermost enclosing group
that declares them (as in lexically scoped code), and `index` gives the
position of the innermost enclosing repeated item.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import hashlib
import json
import zlib
from collections.abc import Callable, Iterator, Sequence
from functools import partial
from typing import Any

import numpy as np

from .nodes import (
    Choice,
    Derived,
    Fixed,
    Infeasible,
    Integer,
    Node,
    Repeat,
    Uncertain,
    Variable,
    _as_node,
    _Domain,
)
from .values import (
    Scope,
    Values,
    _ChoiceValue,
    _RepeatValue,
)

_ABSENT = object()
"""Marker of a node that records no decision (fixed, derived, or uncertain)."""

_RESERVED_NAMES = {"index", "count", "choice"}
"""Names that cannot be given to the children of a design space."""


# ---------------------------------------------------------------------------
# deciding designs
# ---------------------------------------------------------------------------


def _decision_key(decision: Any) -> str:
    """
    Gets a canonical string of a decision.
    """
    return json.dumps(decision, sort_keys=True, default=str)


def _plain(value: Any) -> Any:
    """
    Converts a numpy scalar to a plain Python value, for decisions.
    """
    return value.item() if isinstance(value, np.generic) else value


class _Context:
    """
    How a walk of a design space decides its nodes: `choose(path, domain)`
    decides a variable, `draw(path, node)` draws an uncertain parameter (or
    `None` for nominal values), and `check` enables constraints. If given,
    `record(path, domain, value)` is called with every decision, chosen or
    prescribed. With an `odometer` (of an enumeration), the items of
    unordered repetitions are only walked in the order of their levels.
    """

    def __init__(
        self,
        choose: Callable[[str, _Domain], Any] | None = None,
        draw: Callable[[str, Uncertain], Any] | None = None,
        check: bool = True,
        record: Callable[[str, _Domain, Any], None] | None = None,
        odometer: _Odometer | None = None,
    ):
        self.choose, self.draw, self.check = choose, draw, check
        self.record, self.odometer = record, odometer

    def decide(self, path: str, domain: _Domain, prescribed: Any) -> Any:
        """
        Decides a value of a domain: the prescribed value, or a chosen one.
        """
        if prescribed is _ABSENT:
            value = _plain(_choose(self, path, domain))
        else:
            value = prescribed
        if self.record is not None:
            self.record(path, domain, value)
        return value


class _Unsorted(Infeasible):
    """
    Raised while enumerating a walk whose unordered items are not in the
    order of their levels: another walk gives the same design.
    """


def _join(path: str, name: str) -> str:
    return f"{path}.{name}" if path else name


def _walk(
    node: Node,
    path: str,
    scope: Scope,
    ctx: _Context,
    prescribed: Any = _ABSENT,
    publish: Callable[[Any], None] | None = None,
) -> tuple[Any, Any]:
    """
    Decides a node: returns its value and its decision (the values of its
    design variables, or `_ABSENT` if it has none). A prescribed decision is
    followed instead of choosing. A container publishes its (partial) value
    as soon as it is created, so that its contents are visible to the
    callables of later nodes within it.
    """
    if isinstance(node, Variable):
        if prescribed is not _ABSENT and ctx.record is None:
            return prescribed, prescribed
        value = ctx.decide(path, node.domain(scope), prescribed)
        return value, value
    if isinstance(node, Fixed):
        return node.value, _ABSENT
    if isinstance(node, Derived):
        return node.function(scope), _ABSENT
    if isinstance(node, Uncertain):
        return (node.nominal if ctx.draw is None else ctx.draw(path, node)), _ABSENT
    if isinstance(node, DesignSpace):
        return _walk_group(node, path, scope, ctx, prescribed, publish)
    if isinstance(node, Choice):
        return _walk_choice(node, path, scope, ctx, prescribed, publish)
    if isinstance(node, Repeat):
        return _walk_repeat(node, path, scope, ctx, prescribed, publish)
    raise TypeError(f"unknown node {node!r} at {path or 'the root'}")


def _choose(ctx: _Context, path: str, domain: _Domain) -> Any:
    """
    Chooses a value of a domain, if the walk can choose.
    """
    if ctx.choose is None:
        raise ValueError(f"no decision for {path}")
    return ctx.choose(path, domain)


def _walk_group(
    node: DesignSpace,
    path: str,
    scope: Scope,
    ctx: _Context,
    prescribed: Any,
    publish: Callable[[Any], None] | None,
) -> tuple[dict, dict]:
    """
    Decides the children of a design space in order, then checks its
    constraints (see `_walk`).
    """
    values: dict = {}
    if publish is not None:
        publish(values)
    inner = scope.push(values)
    decisions = {}
    for name, child in node.children.items():
        value, decision = _walk(
            child,
            _join(path, name),
            inner,
            ctx,
            _ABSENT if prescribed is _ABSENT else prescribed.get(name, _ABSENT),
            partial(values.__setitem__, name),
        )
        values[name] = value
        if decision is not _ABSENT:
            decisions[name] = decision
    if ctx.check:
        for i, constraint in enumerate(node.constraints):
            violation = _get_violation(constraint(inner))
            if violation > 0:
                label = getattr(constraint, "__name__", "<lambda>")
                raise Infeasible(
                    f"constraint {i} ({label}) of {path or 'the design space'} "
                    + "is not satisfied",
                    violation,
                )
    return values, decisions


def _get_violation(result: Any) -> float:
    """
    Gets the violation of a constraint from its result: a boolean (`True` if
    satisfied, a violation of 1 otherwise) or a number (satisfied if it is
    not positive, its value is the violation otherwise).
    """
    if isinstance(result, (bool, np.bool_)):
        return 0.0 if result else 1.0
    return float(result)


def _walk_choice(
    node: Choice,
    path: str,
    scope: Scope,
    ctx: _Context,
    prescribed: Any,
    publish: Callable[[Any], None] | None,
) -> tuple[_ChoiceValue, list]:
    """
    Decides the branch of a choice, then the branch (see `_walk`).
    """
    choice = _ABSENT if prescribed is _ABSENT else prescribed[0]
    name = ctx.decide(path, _Domain(levels=list(node.branches)), choice)
    sub = None if prescribed is _ABSENT else prescribed[1]
    result = _ChoiceValue(name)
    if publish is not None:
        publish(result)
    result.value, decision = _walk(
        node.branches[name],
        _join(path, name),
        scope,
        ctx,
        _ABSENT if sub is None else sub,
        partial(setattr, result, "value"),
    )
    return result, [name, None if decision is _ABSENT else decision]


def _walk_items(
    node: Repeat,
    path: str,
    scope: Scope,
    ctx: _Context,
    count: int,
    prescribed: Any,
    publish: Callable[[Any], None] | None,
) -> tuple[_RepeatValue, list]:
    """
    Decides the items of a repetition (see `_walk_repeat`).
    """
    result = _RepeatValue(count)
    if publish is not None:
        publish(result)
    decisions = []
    # with an odometer, permutations of unordered items give the same design:
    # enumerate only the first, whose items are in the order of their levels
    odometer = ctx.odometer if prescribed is _ABSENT and not node.ordered else None
    previous: list[int] | None = None
    for k in range(count):
        sub = _ABSENT if prescribed is _ABSENT else prescribed[k]
        start = 0 if odometer is None else odometer.position
        if odometer is not None and previous is not None:
            odometer.floors.append((start, previous))
        try:
            value, decision = _walk(
                node.item,
                f"{path}[{k}]",
                scope.push({}, k),
                ctx,
                _ABSENT if sub is None else sub,
                result.items.append,
            )
        finally:
            if odometer is not None and previous is not None:
                odometer.floors.pop()
        if len(result.items) == k:
            # a value that is not a container is not published
            result.items.append(value)
        decisions.append(None if decision is _ABSENT else decision)
        if odometer is not None:
            previous = odometer.indices[start : odometer.position]
    return result, decisions


def _walk_repeat(
    node: Repeat,
    path: str,
    scope: Scope,
    ctx: _Context,
    prescribed: Any,
    publish: Callable[[Any], None] | None,
) -> tuple[_RepeatValue, list]:
    """
    Decides the count of a repetition, then its items, which are sorted by
    their decisions and decided again in that order if they are unordered
    (see `_walk`).
    """
    if isinstance(node.count, Integer):
        count = int(
            ctx.decide(
                _join(path, "count"),
                node.count.domain(scope),
                _ABSENT if prescribed is _ABSENT else len(prescribed),
            )
        )
    else:
        count = int(node.count)
    result, decisions = _walk_items(node, path, scope, ctx, count, prescribed, publish)
    if not node.ordered:
        canonical = sorted(decisions, key=_decision_key)
        if canonical != decisions:
            result, decisions = _walk_items(
                node, path, scope, ctx, count, canonical, publish
            )
    return result, decisions


def _flatten(node: Node, value: Any, path: str, out: dict) -> None:
    """
    Flattens the values of a node into `out`, by path: the chosen branch of
    a choice at its path (and its values under the branch's name), the count
    of a repetition at `<path>.count` (and its items at `<path>[k]`), and the
    other values (except fixed values) at their paths.
    """
    if isinstance(node, Fixed):
        return
    if isinstance(node, DesignSpace):
        for name, child in node.children.items():
            _flatten(child, value[name], _join(path, name), out)
    elif isinstance(node, Choice):
        out[path] = value.choice
        _flatten(
            node.branches[value.choice], value.value, _join(path, value.choice), out
        )
    elif isinstance(node, Repeat):
        out[_join(path, "count")] = value.count
        for k, item in enumerate(value.items):
            _flatten(node.item, item, f"{path}[{k}]", out)
    else:
        out[path] = value


# ---------------------------------------------------------------------------
# designs and design spaces
# ---------------------------------------------------------------------------


class Design:
    """
    A design: one decision of every active design variable of a design
    space. Its values are given as attributes (or keys), as for the views of
    the callables of the design space, and are those of a replicate (with
    drawn uncertain parameters) if the design is realized (see `realize`).
    """

    def __init__(
        self,
        space: DesignSpace,
        decisions: dict,
        values: dict,
        replicate: int | None = None,
    ):
        self.space, self.decisions, self.replicate = space, decisions, replicate
        self._values = values
        self.key = hashlib.sha1(_decision_key(decisions).encode()).hexdigest()[:12]
        """str: An identifier of the design, from its decisions only (not
        from its derived values or uncertain parameters)."""

    @property
    def values(self) -> Values:
        """Values: the values of the design."""
        return Values(self._values)

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self.values, name)

    def __getitem__(self, name: str) -> Any:
        return self.values[name]

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Design)
            and self.key == other.key
            and self.replicate == other.replicate
        )

    def __hash__(self) -> int:
        return hash((self.key, self.replicate))

    def flat(self) -> dict[str, Any]:
        """
        Gets the values of the design by path, as for the columns of a
        table: `a.b` for the value `b` of group `a`, `c` for the chosen
        branch of choice `c` (and `c.<branch>.d` for its values), `e.count`
        for the count of repetition `e` (and `e[k].f` for the values of its
        items). Fixed values are omitted.

        Returns:
            dict: the values by path
        """
        out: dict = {}
        _flatten(self.space, self._values, "", out)
        return out

    def realize(self, replicate: int, seed: int = 0, common: bool = True) -> Design:
        """
        Realizes a replicate of the design: draws its uncertain parameters
        (and computes the values derived from them).

        Each uncertain parameter is drawn from a random number generator
        seeded by `seed`, `replicate`, and its path, so that with `common`
        random numbers, a parameter at the same path (such as the phase of
        the first orbit plane) takes the same value in the same replicate
        of every design. Comparisons between designs then measure the
        differences between the designs rather than between their draws,
        which needs fewer replicates to rank them. Otherwise, the generator
        is also seeded by the design's key, so that draws are independent
        between designs.

        Args:
            replicate (int): The replicate.
            seed (int): The seed of the trade study.
            common (bool): `True` to use common random numbers between
                designs.

        Returns:
            Design: the realized design
        """
        stream = [seed, replicate] + ([] if common else [int(self.key, 16)])

        def draw(path: str, node: Uncertain) -> Any:
            rng = np.random.default_rng(stream + [zlib.crc32(path.encode())])
            return node.sampler(rng)

        values, _ = _walk(
            self.space, "", Scope(), _Context(draw=draw, check=False), self.decisions
        )
        return Design(self.space, self.decisions, values, replicate)

    def __repr__(self) -> str:
        replicate = "" if self.replicate is None else f", replicate={self.replicate}"
        return f"Design({self.key}{replicate}: {self.flat()})"


class _Odometer:
    """
    Chooses the levels of the variables of successive walks in the order of
    an odometer whose digits are the decisions of a walk, so that walks
    enumerate every combination of levels, including those of variables
    whose existence or domain depends on earlier decisions.
    """

    def __init__(self):
        self.indices: list[int] = []
        self.sizes: list[int] = []
        self.position = 0
        self.floors: list[tuple[int, list[int]]] = []
        """The levels of the previous unordered items (and where the
        current items' levels start), which the current items must not
        precede."""

    def choose(self, path: str, domain: _Domain) -> Any:
        """
        Chooses the level of the current digit of the odometer (adding a
        digit, at the first level, if the walk is longer than before), and
        stops a walk whose unordered items precede the previous ones.
        """
        if domain.levels is None:
            raise ValueError(
                f"{path} is continuous: give it a step to enumerate the design space"
            )
        if len(domain.levels) == 0:
            raise Infeasible(f"{path} has no levels")
        if self.position == len(self.indices):
            self.indices.append(0)
            self.sizes.append(len(domain.levels))
        level = domain.levels[self.indices[self.position]]
        self.position += 1
        for start, previous in self.floors:
            current = self.indices[start : self.position]
            if current < previous[: len(current)]:
                raise _Unsorted()
        return level

    def advance(self) -> bool:
        """
        Advances to the next combination after a walk, returning `False`
        after the last.
        """
        del self.indices[self.position :], self.sizes[self.position :]
        self.position = 0
        self.floors.clear()
        while self.indices:
            if self.indices[-1] + 1 < self.sizes[-1]:
                self.indices[-1] += 1
                return True
            self.indices.pop()
            self.sizes.pop()
        return False


def _lhs_column(n: int, rng: np.random.Generator) -> np.ndarray:
    """
    Gets `n` values in [0, 1) stratified in `n` equal intervals in random
    order (a column of a Latin hypercube sample).
    """
    return (rng.permutation(n) + rng.random(n)) / n


class DesignSpace(Node):
    """
    A group of nodes: the root of a design space, or a group within it (such
    as the design of each orbit plane). Children are given by name, in the
    order they are decided; plain values are fixed values.

    Constraints are callables of the group's values (and those of enclosing
    groups), checked once the group's values are decided, that return
    whether a design satisfies them (`True` or `False`) or by how much it
    violates them: a number that is not positive if it satisfies them, and
    otherwise the amount of violation (for example, `lambda d: d.satellites -
    12` for at most 12 satellites), which guides constraint-aware optimizers
    toward feasible designs (see `optimize_optuna` and `make_pymoo_problem`).
    They must not depend on uncertain parameters.

    Args:
        constraints (Sequence[Callable]): The constraints.
        **children (Any): The nodes, by name.
    """

    def __init__(
        self, constraints: Sequence[Callable[[Scope], bool]] = (), **children: Any
    ):
        reserved = _RESERVED_NAMES.intersection(children)
        if reserved:
            raise ValueError(f"reserved names: {sorted(reserved)}")
        self.children = {name: _as_node(node) for name, node in children.items()}
        self.constraints = list(constraints)

    def __repr__(self) -> str:
        return f"DesignSpace({', '.join(self.children)})"

    def _design(self, ctx: _Context) -> Design:
        """
        Decides a design (with nominal uncertain parameters).
        """
        values, decisions = _walk(self, "", Scope(), ctx)
        return Design(self, decisions, values)

    def design(self, decisions: dict) -> Design:
        """
        Gets the design of given decisions (see `Design.decisions`), checking
        its constraints.

        Args:
            decisions (dict): The decisions.

        Returns:
            Design: the design

        Raises:
            Infeasible: if the design violates a constraint.
        """
        values, decisions = _walk(self, "", Scope(), _Context(), decisions)
        return Design(self, decisions, values)

    def enumerate(self) -> list[Design]:
        """
        Enumerates every feasible design (a full factorial design over the
        levels of the design variables, including those that exist only in
        some designs), without duplicates.

        Returns:
            list[Design]: the designs, in the order of their decisions

        Raises:
            ValueError: if a design variable is continuous.
        """
        return list(self._iterate())

    def _iterate(self) -> Iterator[Design]:
        """
        Iterates over every feasible design (see `enumerate`). The items of
        unordered repetitions are walked in only one of their orders, so
        that their permutations are not walked.
        """
        odometer = _Odometer()
        ctx = _Context(choose=odometer.choose, odometer=odometer)
        keys: set[str] = set()
        while True:
            try:
                design = self._design(ctx)
                if design.key not in keys:
                    keys.add(design.key)
                    yield design
            except Infeasible:
                pass
            if not odometer.advance():
                return

    def size(self) -> int:
        """
        Counts the feasible designs (see `enumerate`), without keeping them.

        Returns:
            int: the number of designs
        """
        return sum(1 for _ in self._iterate())

    def sample(
        self,
        n: int,
        method: str = "lhs",
        seed: int = 0,
        max_tries: int = 1000,
    ) -> list[Design]:
        """
        Samples feasible designs, which may repeat.

        With `method="lhs"`, each design variable (identified by its path,
        including those that exist only in some designs) is sampled in a
        Latin hypercube: its range is divided into `n` equal parts and each
        design samples a different part. A design that is infeasible is
        drawn again at random, which disturbs the stratification. With
        `method="random"`, every variable is sampled independently.

        Args:
            n (int): The number of designs.
            method (str): `"lhs"` (Latin hypercube) or `"random"`.
            seed (int): The seed of the random number generator.
            max_tries (int): The number of draws of a design before giving up.

        Returns:
            list[Design]: the designs

        Raises:
            Infeasible: if no feasible design is drawn in `max_tries`.
        """
        if method not in ("lhs", "random"):
            raise ValueError(f"method must be 'lhs' or 'random', not {method!r}")
        rng = np.random.default_rng(seed)
        columns: dict[str, np.ndarray] = {}
        designs = []
        for i in range(n):
            for attempt in range(max_tries):
                stratified = method == "lhs" and attempt == 0

                def unit(path: str, i: int = i, stratified: bool = stratified) -> float:
                    if not stratified:
                        return float(rng.random())
                    if path not in columns:
                        columns[path] = _lhs_column(n, rng)
                    return float(columns[path][i])

                ctx = _Context(
                    choose=lambda path, domain, unit=unit: domain.from_unit(unit(path))
                )
                try:
                    designs.append(self._design(ctx))
                    break
                except Infeasible:
                    continue
            else:
                raise Infeasible(f"no feasible design in {max_tries} draws")
        return designs

    def unit_paths(self) -> list[str]:
        """
        Gets the paths of every decision of the design space, in the order
        they are decided: its design variables, its choices, and the counts
        of its repetitions, with the items of each repetition up to its
        largest count. These are the dimensions of the unit hypercube that
        encodes its designs (see `encode` and `decode`).

        Returns:
            list[str]: the paths

        Raises:
            ValueError: if the count of a repetition has an upper bound that
                is not fixed.
        """
        return _unit_paths(self, "")

    def decode(self, unit: Sequence[float]) -> Design:
        """
        Decodes a design from a point of the unit hypercube (a value in
        [0, 1] for each of `unit_paths`): each decision maps the value of its
        path to a level (each level spans an equal part of [0, 1]) or to a
        point of its interval. Values of the paths of inactive nodes (such as
        the branches not chosen) are ignored. This encodes the design space
        for optimizers of fixed numbers of continuous variables, such as
        pymoo's.

        Args:
            unit (Sequence[float]): The values, in the order of `unit_paths`.

        Returns:
            Design: the design

        Raises:
            Infeasible: if the design violates a constraint.
        """
        index = {path: i for i, path in enumerate(self.unit_paths())}
        if len(unit) != len(index):
            raise ValueError(f"expected {len(index)} values, not {len(unit)}")
        return self._design(
            _Context(choose=lambda path, domain: domain.from_unit(unit[index[path]]))
        )

    def encode(self, design: Design) -> np.ndarray:
        """
        Encodes a design as a point of the unit hypercube (see `decode`): the
        center of the part of [0, 1] of each decided level, and 0.5 for the
        paths of inactive nodes.

        Args:
            design (Design): The design.

        Returns:
            numpy.ndarray: the values, in the order of `unit_paths`
        """
        index = {path: i for i, path in enumerate(self.unit_paths())}
        unit = np.full(len(index), 0.5)

        def record(path: str, domain: _Domain, value: Any) -> None:
            unit[index[path]] = domain.to_unit(value)

        _walk(self, "", Scope(), _Context(record=record), design.decisions)
        return unit

    def suggest(self, trial: Any) -> Design:
        """
        Decides a design with an Optuna trial (define-by-run): each decision
        is suggested by the trial, named by its path, with the kind of its
        variable (`suggest_float`, with the step of a discrete `Real`;
        `suggest_int`; or `suggest_categorical`, also for choices), so that
        Optuna's samplers search the hierarchical design space directly.

        Args:
            trial (optuna.trial.Trial): The trial.

        Returns:
            Design: the design

        Raises:
            Infeasible: if the design violates a constraint.
        """
        return self._design(
            _Context(choose=lambda path, domain: domain.suggest(trial, path))
        )


def _unit_paths(node: Node, path: str) -> list[str]:
    """
    Gets the paths of the decisions of a node (see `DesignSpace.unit_paths`).
    """
    if isinstance(node, Variable):
        return [path]
    if isinstance(node, DesignSpace):
        return [
            p
            for name, child in node.children.items()
            for p in _unit_paths(child, _join(path, name))
        ]
    if isinstance(node, Choice):
        return [path] + [
            p
            for name, branch in node.branches.items()
            for p in _unit_paths(branch, _join(path, name))
        ]
    if isinstance(node, Repeat):
        if isinstance(node.count, Integer):
            if callable(node.count.high):
                raise ValueError(
                    f"the count of {path} needs a fixed upper bound to be encoded"
                )
            paths, largest = [_join(path, "count")], int(node.count.high)
        else:
            paths, largest = [], int(node.count)
        return paths + [
            p for k in range(largest) for p in _unit_paths(node.item, f"{path}[{k}]")
        ]
    return []
