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
import math
import re
import zlib
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from functools import partial
from typing import Any

import numpy as np
import pandas as pd

from .fingerprint import _Describer
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
    `None` for nominal values), `fix(path, node)` gives a fixed value (or
    `None` for its value), and `check` enables constraints. If given,
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
        fix: Callable[[str, Fixed], Any] | None = None,
    ):
        self.choose, self.draw, self.check = choose, draw, check
        self.fix = fix
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


_ITEM_INDEX = re.compile(r"\[\d+\]")


def _generic(path: str) -> str:
    """
    Gets the path of a node with `[*]` for the index of any item of a
    repetition (as of `DesignSpace.uncertain_paths`).
    """
    return _ITEM_INDEX.sub("[*]", path)


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
        return (node.value if ctx.fix is None else ctx.fix(path, node)), _ABSENT
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
        self._scenario = ""
        """str: A description of the values of the scenario of a realized
        design (see `realize`), `""` without one."""
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

    def realize(
        self,
        replicate: int,
        seed: int = 0,
        common: bool = True,
        shift: Mapping[str, int] | None = None,
        scenario: Mapping[str, Any] | None = None,
    ) -> Design:
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

        A `shift` draws some parameters as in other replicates: a parameter
        named in it (by its path, with `[*]` for the index of any item of a
        repetition, see `DesignSpace.uncertain_paths`) is drawn as in
        replicate `replicate + shift[name]`, and the others as in
        `replicate`. Such a mixed replicate differs from `replicate` in the
        named parameters only, which measures their effects on the outcomes
        of a design (see `attribute_variance`).

        A `scenario` gives fixed values and uncertain parameters other values
        (by their paths, with `[*]` for the index of any item of a
        repetition, as of `DesignSpace.uncertain_paths`): a condition under
        which designs are compared, such as a requirement, rather than
        averaged over (see `Study`). Values derived from them are derived
        from the scenario's, and an uncertain parameter that it gives is not
        drawn.

        Args:
            replicate (int): The replicate.
            seed (int): The seed of the trade study.
            common (bool): `True` to use common random numbers between
                designs.
            shift (Mapping[str, int] | None): The offsets of the replicates
                as in which the named uncertain parameters are drawn, if any.
            scenario (Mapping[str, Any] | None): The values of fixed values
                and uncertain parameters of a scenario (by path), if any.

        Returns:
            Design: the realized design

        Raises:
            ValueError: if `shift` names a path that is not of an uncertain
                parameter of the design space, or `scenario` one that is not
                of a fixed value or an uncertain parameter.
        """
        shift = dict(shift or {})
        _check_uncertain(self.space, shift)
        scenario = dict(scenario or {})
        _check_scenario(self.space, scenario)
        tail = [] if common else [int(self.key, 16)]

        def draw(path: str, node: Uncertain) -> Any:
            name = _generic(path)
            if name in scenario:
                return scenario[name]
            drawn_in = replicate + int(shift.get(name, 0))
            rng = np.random.default_rng(
                [seed, drawn_in, *tail, zlib.crc32(path.encode())]
            )
            return node.sampler(rng)

        def fix(path: str, node: Fixed) -> Any:
            return scenario.get(_generic(path), node.value)

        ctx = _Context(draw=draw, check=False, fix=fix if scenario else None)
        values, _ = _walk(self.space, "", Scope(), ctx, self.decisions)
        realized = Design(self.space, self.decisions, values, replicate)
        # pylint: disable-next=protected-access
        realized._scenario = _describe_scenario(scenario)
        return realized

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

    def design(
        self, decisions: Mapping | Design | None = None, /, **overrides: Any
    ) -> Design:
        """
        Gets the design of given decisions (see `Design.decisions`), or of the
        decisions of another design (such as of a design space this one
        extends, see `extend`), with the decisions of some children replaced
        (by name: a value, the decisions of a group or a choice, or a design
        of their own space), checking its constraints::

            space = base.extend(policy=policy_space(orbit_spares=Integer(0, 1)))
            design = space.design(base_design, policy=policy_design)

        Args:
            decisions (Mapping | Design | None): The decisions, or a design
                whose decisions to take.
            **overrides (Any): The decisions of children, by name.

        Returns:
            Design: the design

        Raises:
            Infeasible: if the design violates a constraint.
            ValueError: if an override names no child of the design space.
        """
        if isinstance(decisions, Design):
            decisions = decisions.decisions
        merged = dict(decisions or {})
        unknown = sorted(set(overrides) - set(self.children))
        if unknown:
            raise ValueError(f"no children {unknown} of the design space")
        for name, decision in overrides.items():
            merged[name] = (
                decision.decisions if isinstance(decision, Design) else decision
            )
        values, merged = _walk(self, "", Scope(), _Context(), merged)
        return Design(self, merged, values)

    def extend(
        self, constraints: Sequence[Callable[[Scope], Any]] = (), **children: Any
    ) -> DesignSpace:
        """
        Gets a design space with other children: those named as children of
        this one replace them (in their place), and others are added after
        them, with this one's constraints and others. This one is unchanged::

            policy_trade = space.extend(policy=policy_space(orbit_spares=Integer(0, 1)))
            year = space.extend(epoch=Uncertain(draw_epoch))

        Designs of this design space are designs of the other with the
        decisions of its new children (see `design`).

        Args:
            constraints (Sequence[Callable]): Other constraints.
            **children (Any): The nodes, by name.

        Returns:
            DesignSpace: the design space
        """
        return DesignSpace(
            [*self.constraints, *constraints], **{**self.children, **children}
        )

    def designs_of(
        self, frame: pd.DataFrame, keys: Iterable[str] | None = None
    ) -> list[Design]:
        """
        Rebuilds designs from the rows of a data frame of their values by path
        (see `Design.flat`): the results of a study (see `Study.run`), or a
        summary that keeps the values of the designs (see `aggregate`). Each
        design is rebuilt once, from its first row, in the order of its first
        row (or of `keys`), and must have the key of its row (`design`, a
        column or the name of the index), if the frame has one::

            front = summary[pareto_front(summary, criteria)]
            designs = space.designs_of(front)

        Args:
            frame (pandas.DataFrame): The rows.
            keys (Iterable[str] | None): The keys of the designs to rebuild,
                or `None` for every design of the frame.

        Returns:
            list[Design]: the designs

        Raises:
            KeyError: if the frame has no rows of a key.
            ValueError: if the frame lacks the value of a decision of a
                design, or its values give a design of another key.
        """
        rows = frame
        if "design" not in rows.columns and rows.index.name == "design":
            rows = rows.reset_index()
        if "design" in rows.columns:
            first = rows.drop_duplicates("design")
            records = dict(zip(first["design"], first.to_dict("records")))
            order = list(records) if keys is None else list(dict.fromkeys(keys))
            missing = [key for key in order if key not in records]
            if missing:
                raise KeyError(f"no rows of designs {missing[:5]}")
            pairs = [(key, records[key]) for key in order]
        else:
            if keys is not None:
                raise ValueError("the frame has no keys of designs (`design`)")
            pairs = [(None, record) for record in rows.to_dict("records")]
        designs: dict[str, Design] = {}
        for key, record in pairs:
            design = self._design(_Context(choose=partial(_decided, record)))
            if key is not None and design.key != key:
                raise ValueError(
                    f"the values of design {key} give another design ({design.key}): "
                    + "are they of another design space?"
                )
            designs.setdefault(design.key, design)
        return list(designs.values())

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

    def uncertain_paths(self) -> list[str]:
        """
        Gets the paths of the uncertain parameters of the design space (see
        `Uncertain`), in the order they are drawn, with `[*]` for the index
        of any item of a repetition (so that the phases of every plane of a
        constellation are one path, such as `planes[*].phase`), and with
        the branches of choices (such as `layout.walker.phase_fraction`),
        whether a design chooses them or not.

        Returns:
            list[str]: the paths
        """
        return _uncertain_paths(self, "")

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


def _fixed_paths(node: Node, path: str) -> list[str]:
    """
    Gets the paths of the fixed values of a node, with `[*]` for the index of
    any item of a repetition (as `_uncertain_paths`).
    """
    if isinstance(node, Fixed):
        return [path]
    if isinstance(node, DesignSpace):
        return [
            p
            for name, child in node.children.items()
            for p in _fixed_paths(child, _join(path, name))
        ]
    if isinstance(node, Choice):
        return [
            p
            for name, branch in node.branches.items()
            for p in _fixed_paths(branch, _join(path, name))
        ]
    if isinstance(node, Repeat):
        return _fixed_paths(node.item, f"{path}[*]")
    return []


def _missing(value: Any) -> bool:
    """Whether a value of a data frame is missing (`None` or `NaN`)."""
    return value is None or (isinstance(value, float) and math.isnan(value))


def _decided(values: Mapping[str, Any], path: str, domain: _Domain) -> Any:
    """
    Gets the decision of a path from values by path (such as a row of the
    results of a study, see `DesignSpace.designs_of`): the level of the domain
    that equals the value (whose type a data frame may have changed, as
    integers in a column of floats), or the value of a continuous domain.
    """
    value = values.get(path)
    if _missing(value):
        raise ValueError(
            f"no value of {path}: rebuild designs from the results of a study, or "
            + "from a summary that keeps the values of designs (see aggregate)"
        )
    if domain.levels is None:
        return float(value)
    for level in domain.levels:
        try:
            if level == value:
                return level
        except (TypeError, ValueError):
            continue
    numeric = [
        level
        for level in domain.levels
        if isinstance(level, (int, float, np.number)) and not isinstance(level, bool)
    ]
    if isinstance(value, (int, float, np.number)):
        for level in numeric:
            if math.isclose(float(level), float(value), rel_tol=1e-9, abs_tol=1e-12):
                return level
    raise ValueError(f"{value!r} is not a level of {path} ({domain.levels})")


def _designs(
    space: DesignSpace,
    designs: Iterable[Design | str] | pd.DataFrame,
    *frames: pd.DataFrame,
) -> list[Design]:
    """
    Gets designs given as designs, keys, or a data frame of designs (with
    their keys, `design`), each once, in order: designs given by their keys
    are rebuilt from the values of the frame (if it has them) or of other
    frames (such as the results of a study, see `DesignSpace.designs_of`).
    """
    if isinstance(designs, pd.DataFrame):
        if "design" in designs.columns:
            items = list(designs["design"])
        elif designs.index.name == "design":
            items = list(designs.index)
        else:
            raise ValueError("the data frame has no keys of designs (`design`)")
        frames = (designs, *frames)
    else:
        items = list(designs)
    found: dict[str, Design | None] = {}
    for item in items:
        if isinstance(item, Design):
            found[item.key] = item
        elif isinstance(item, str):
            found.setdefault(item, None)
        else:
            raise TypeError(f"a design or the key of one, not {item!r}")
    for frame in frames:
        missing = [key for key, design in found.items() if design is None]
        if not missing:
            break
        rows = frame.reset_index() if "design" not in frame else frame
        if "design" not in rows:
            continue
        present = [key for key in missing if key in set(rows["design"])]
        try:
            for design in space.designs_of(rows, present):
                found[design.key] = design
        except ValueError:
            continue  # the frame lacks values of the designs: try the next
    missing = [key for key, design in found.items() if design is None]
    if missing:
        raise ValueError(
            f"the values of designs {missing[:5]} are not given: give the designs, "
            + "or their keys with the results of a study (or a summary that keeps "
            + "the values of designs, see aggregate)"
        )
    return [design for design in found.values() if design is not None]


def _check_scenario(space: DesignSpace, values: Mapping[str, Any]) -> None:
    """
    Checks that the values of a scenario are of fixed values or uncertain
    parameters of a design space (by path), raising a `ValueError` if not.
    """
    if not values:
        return
    known = [*_fixed_paths(space, ""), *space.uncertain_paths()]
    unknown = sorted(set(values) - set(known))
    if unknown:
        raise ValueError(
            "not fixed values or uncertain parameters of the design space: "
            + f"{unknown} (it has {known})"
        )


def _describe_scenario(values: Mapping[str, Any]) -> str:
    """
    Gets a description of the values of a scenario (`""` without values),
    which identifies them (by their content, see `get_definition`).
    """
    if not values:
        return ""
    described = _Describer().describe(dict(values))
    return json.dumps(described, sort_keys=True, default=str)


def _scenario_of(design: Design) -> str:
    """
    Gets the description of the values of the scenario of a realized design
    (see `Design.realize`), `""` without one.
    """
    return getattr(design, "_scenario", "")


def _check_uncertain(space: DesignSpace, paths: Iterable[str]) -> None:
    """
    Checks that paths are of uncertain parameters of a design space (see
    `DesignSpace.uncertain_paths`), raising a `ValueError` if not.
    """
    known = space.uncertain_paths()
    unknown = sorted(set(paths) - set(known))
    if unknown:
        raise ValueError(
            f"not uncertain parameters of the design space: {unknown} (it has {known})"
        )


def _uncertain_paths(node: Node, path: str) -> list[str]:
    """
    Gets the paths of the uncertain parameters of a node (see
    `DesignSpace.uncertain_paths`).
    """
    if isinstance(node, Uncertain):
        return [path]
    if isinstance(node, DesignSpace):
        return [
            p
            for name, child in node.children.items()
            for p in _uncertain_paths(child, _join(path, name))
        ]
    if isinstance(node, Choice):
        return [
            p
            for name, branch in node.branches.items()
            for p in _uncertain_paths(branch, _join(path, name))
        ]
    if isinstance(node, Repeat):
        return _uncertain_paths(node.item, f"{path}[*]")
    return []
