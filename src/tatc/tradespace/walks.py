"""
Walks of design spaces: the decisions of the nodes of a design space in the
order they are declared, which decide designs (see `DesignSpace`), and the
paths of its nodes.

A walk decides each node in turn: a design variable by a decision that is
prescribed (such as the decisions of a design) or chosen (by an enumeration,
a sample, or an optimizer), a fixed, derived, or uncertain value by its
value (or a value of a scenario, or a draw), a group by its children (then
its constraints), a choice by its branch (then the branch), and a repetition
by its count (then its items, sorted by the levels of their decisions if they
are unordered).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import copy
import math
from collections.abc import Callable, Mapping
from functools import partial
from typing import Any

import numpy as np

from .nodes import (
    Choice,
    Derived,
    Fixed,
    Group,
    Infeasible,
    Integer,
    Node,
    Repeat,
    Uncertain,
    Variable,
    _Domain,
)
from .values import Scope, _ChoiceValue, _join, _RepeatValue, is_number

_ABSENT = object()
"""Marker of a node that records no decision (fixed, derived, or uncertain),
or of a decision that is not prescribed."""

_FIRST = object()
"""Marker of a node whose decisions are their first levels (as of a branch
of a choice given by its name, see `DesignSpace.design`)."""


def _scalar(value: Any) -> Any:
    """
    Converts a numpy scalar to a plain Python value, for decisions.
    """
    return value.item() if isinstance(value, np.generic) else value


def _first_level(path: str, domain: _Domain) -> Any:
    """
    Chooses the first level of a domain (or its lower bound).
    """
    if domain.levels is not None and len(domain.levels) == 0:
        raise Infeasible(f"{path} has no levels")
    return domain.from_unit(0.0)


def _rank(domain: _Domain, value: Any) -> float:
    """
    Gets the rank of a decision in its domain, by which unordered items are
    sorted: the index of its level (or its value, if continuous).
    """
    if domain.levels is None:
        return float(value)
    try:
        return domain.levels.index(value)
    except ValueError:
        return len(domain.levels)


class _Context:
    """
    How a walk of a design space decides its nodes: `choose(path, domain)`
    decides a variable, `draw(path, node)` draws an uncertain parameter (or
    `None` for nominal values), `fix(path, node)` gives a fixed value (or
    `None` for its value), and `check` enables constraints. If given,
    `normalize(path, domain, value)` gives the value of the domain that a
    prescribed decision is (such as its level, see `_level_of`). With an
    `odometer` (of an enumeration), the items of unordered repetitions are
    only walked in the order of their levels. A walk of the decisions of a
    design (`canonical`), whose unordered items are in order, does not sort
    them again.
    """

    # pylint: disable=too-many-instance-attributes

    def __init__(
        self,
        choose: Callable[[str, _Domain], Any] | None = None,
        draw: Callable[[str, Uncertain], Any] | None = None,
        check: bool = True,
        odometer: _Odometer | None = None,
        fix: Callable[[str, Fixed], Any] | None = None,
        normalize: Callable[[str, _Domain, Any], Any] | None = None,
        canonical: bool = False,
    ):
        # pylint: disable=too-many-arguments,too-many-positional-arguments
        self.choose, self.draw, self.check = choose, draw, check
        self.fix, self.normalize = fix, normalize
        self.odometer = odometer
        self.canonical = canonical
        self.ranks: list[list[float]] = []
        """The ranks of the decisions of the unordered items being walked
        (the innermost last), by which they are sorted (see `_rank`)."""

    def decide(self, path: str, domain: _Domain, prescribed: Any) -> Any:
        """
        Decides a value of a domain: the prescribed value (normalized, if
        the walk normalizes them), or a chosen one.
        """
        if prescribed is _ABSENT:
            value = _scalar(_choose(self, path, domain))
        elif self.normalize is not None:
            value = self.normalize(path, domain, prescribed)
        else:
            value = prescribed
        if self.ranks:
            rank = _rank(domain, value)
            for ranks in self.ranks:
                ranks.append(rank)
        return value

    def first(self) -> _Context:
        """
        Gets this context, but choosing the first levels of the domains (see
        `_FIRST`).
        """
        other = copy.copy(self)
        other.choose = _first_level
        return other


class _Unsorted(Infeasible):
    """
    Raised while enumerating a walk whose unordered items are not in the
    order of their levels: another walk gives the same design.
    """


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
    followed instead of choosing (or, if `_FIRST`, the first levels are). A
    container publishes its (partial) value as soon as it is created, so
    that its contents are visible to the callables of later nodes within it.
    """
    if prescribed is _FIRST:
        ctx, prescribed = ctx.first(), _ABSENT
    if isinstance(node, Variable):
        if prescribed is not _ABSENT and ctx.normalize is None and not ctx.ranks:
            return prescribed, prescribed
        value = ctx.decide(path, node.domain(scope), prescribed)
        return value, value
    if isinstance(node, Fixed):
        return (node.value if ctx.fix is None else ctx.fix(path, node)), _ABSENT
    if isinstance(node, Derived):
        return node.function(scope), _ABSENT
    if isinstance(node, Uncertain):
        return (node.nominal if ctx.draw is None else ctx.draw(path, node)), _ABSENT
    if isinstance(node, Group):
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
    node: Group,
    path: str,
    scope: Scope,
    ctx: _Context,
    prescribed: Any,
    publish: Callable[[Any], None] | None,
) -> tuple[dict, dict]:
    """
    Decides the children of a group in order, then checks its constraints
    (see `_walk`).
    """
    if prescribed is not _ABSENT and ctx.normalize is not None:
        _check_children(node, path, prescribed)
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
            if math.isnan(violation) or violation > 0:
                label = getattr(constraint, "__name__", "<lambda>")
                raise Infeasible(
                    f"constraint {i} ({label}) of {path or 'the design space'} "
                    + ("is NaN" if math.isnan(violation) else "is not satisfied"),
                    1.0 if math.isnan(violation) else violation,
                )
    return values, decisions


def _check_children(node: Group, path: str, prescribed: Any) -> None:
    """
    Checks that the prescribed decisions of a group (such as those given to
    `DesignSpace.design`) are a mapping of the decisions of its children,
    raising a `ValueError` if not, or if they name others.
    """
    where = path or "the design space"
    if not isinstance(prescribed, Mapping):
        raise ValueError(
            f"the decisions of {where} are a mapping of its children, "
            + f"not {prescribed!r}"
        )
    unknown = [name for name in prescribed if name not in node.children]
    if unknown:
        raise ValueError(
            f"decisions of no children of {where}: {sorted(unknown, key=str)} "
            + f"(it has {list(node.children)})"
        )


def _get_violation(result: Any) -> float:
    """
    Gets the violation of a constraint from its result: a boolean (`True` if
    satisfied, a violation of 1 otherwise) or a number (satisfied if it is
    not positive, its value is the violation otherwise, and `NaN`, which is
    not satisfied, if it is `NaN`).
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
    Decides the branch of a choice, then the branch (see `_walk`). A branch
    prescribed by its name alone is decided at its first levels.
    """
    if prescribed is _ABSENT:
        choice, sub = _ABSENT, _ABSENT
    elif isinstance(prescribed, str):
        choice, sub = prescribed, _FIRST
    else:
        choice = prescribed[0]
        sub = _ABSENT if prescribed[1] is None else prescribed[1]
    name = ctx.decide(path, _Domain(levels=list(node.branches)), choice)
    result = _ChoiceValue(name)
    if publish is not None:
        publish(result)
    result.value, decision = _walk(
        node.branches[name],
        _join(path, name),
        scope,
        ctx,
        sub,
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
    ranked: bool = False,
) -> tuple[_RepeatValue, list, list[list[float]]]:
    """
    Decides the items of a repetition (see `_walk_repeat`), and the ranks of
    their decisions, if `ranked` (see `_rank`).
    """
    # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    result = _RepeatValue(count)
    if publish is not None:
        publish(result)
    decisions, ranks = [], []
    # with an odometer, permutations of unordered items give the same design:
    # enumerate only the first, whose items are in the order of their levels
    odometer = ctx.odometer if prescribed is _ABSENT and not node.ordered else None
    previous: list[int] | None = None
    for k in range(count):
        sub = _ABSENT if prescribed is _ABSENT else prescribed[k]
        start = 0 if odometer is None else odometer.position
        if odometer is not None and previous is not None:
            odometer.floors.append((start, previous))
        item_ranks: list[float] = []
        if ranked:
            ctx.ranks.append(item_ranks)
        try:
            value, decision = _walk(
                node.item,
                _join(path, k),
                scope.push({}, k),
                ctx,
                _ABSENT if sub is None else sub,
                result.items.append,
            )
        finally:
            if odometer is not None and previous is not None:
                odometer.floors.pop()
            if ranked:
                ctx.ranks.pop()
        if len(result.items) == k:
            # a value that is not a container is not published
            result.items.append(value)
        decisions.append(None if decision is _ABSENT else decision)
        ranks.append(item_ranks)
        if odometer is not None:
            previous = odometer.indices[start : odometer.position]
    return result, decisions, ranks


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
    the levels of their decisions (in the order they are decided, see
    `_rank`) and decided again in that order if they are unordered (see
    `_walk`).
    """
    if isinstance(node.count, Variable):
        count = int(
            ctx.decide(
                _join(path, "count"),
                node.count.domain(scope),
                _ABSENT if prescribed is _ABSENT else len(prescribed),
            )
        )
    elif isinstance(node.count, Fixed):
        count = int(node.count.value)
    elif isinstance(node.count, Derived):
        count = int(node.count.function(scope))
    else:
        count = int(node.count)
    if count < 0:
        raise Infeasible(f"the count of {path} is negative ({count})")
    if prescribed is not _ABSENT and len(prescribed) != count:
        raise ValueError(
            f"{len(prescribed)} items are given for {path}, whose count is {count}"
        )
    ranked = not node.ordered and not ctx.canonical
    # the ranks of the enclosing unordered items, which take those of these
    # items in their order
    marks = [len(ranks) for ranks in ctx.ranks]
    result, decisions, ranks = _walk_items(
        node, path, scope, ctx, count, prescribed, publish, ranked
    )
    if ranked:
        order = sorted(range(count), key=ranks.__getitem__)
        if order != list(range(count)):
            for enclosing, mark in zip(ctx.ranks, marks):
                del enclosing[mark:]
            canonical = [decisions[k] for k in order]
            result, decisions, _ = _walk_items(
                node, path, scope, ctx, count, canonical, publish
            )
    return result, decisions


def _flatten(node: Node, value: Any, path: str, out: dict) -> None:
    """
    Flattens the values of a node into `out`, by path: the chosen branch of
    a choice at its path (and its values under the branch's name), the count
    of a repetition at `<path>.count` (also of a fixed count) and its items
    at `<path>[k]`, and the other values (except fixed values) at their
    paths.
    """
    if isinstance(node, Fixed):
        return
    if isinstance(node, Group):
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
            _flatten(node.item, item, _join(path, k), out)
    else:
        out[path] = value


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


def _level_of(path: str, domain: _Domain, value: Any) -> Any:
    """
    Gets the value of a domain that a decision is (see `DesignSpace.design`
    and `DesignSpace.designs_of`): the level that equals it (whose type a
    data frame may have changed, as integers in a column of floats, or that
    is given as a numpy scalar), or the value (a float) in the range of a
    continuous domain, raising a `ValueError` if there is none.
    """
    if domain.levels is None:
        try:
            number = float(value)
        except (TypeError, ValueError):
            number = math.nan
        if domain.low <= number <= domain.high:
            return number
        raise ValueError(
            f"{value!r} is not in the range of {path} ({domain.low}, {domain.high})"
        )
    for level in domain.levels:
        try:
            if level == value:
                return _scalar(level)
        except (TypeError, ValueError):
            continue
    if is_number(value):
        for level in domain.levels:
            if is_number(level) and math.isclose(
                float(level), float(value), rel_tol=1e-9, abs_tol=1e-12
            ):
                return _scalar(level)
    raise ValueError(f"{value!r} is not a level of {path} ({domain.levels})")


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
            + "from a summary that keeps the values of designs (see summarize)"
        )
    return _level_of(path, domain, value)


# ---------------------------------------------------------------------------
# paths of nodes
# ---------------------------------------------------------------------------


def _unit_paths(node: Node, path: str) -> list[str]:
    """
    Gets the paths of the decisions of a node (see `DesignSpace._unit_paths`).
    """
    if isinstance(node, Variable):
        return [path]
    if isinstance(node, Group):
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
        paths, largest = [], _largest_count(node.count, path)
        if isinstance(node.count, Variable):
            paths = [_join(path, "count")]
        return paths + [
            p for k in range(largest) for p in _unit_paths(node.item, _join(path, k))
        ]
    return []


def _largest_count(count: Any, path: str) -> int:
    """
    Gets the largest count of a repetition (see `_unit_paths`), which must
    not depend on other values.
    """
    if isinstance(count, Integer) and not callable(count.high):
        return int(count.high)
    if isinstance(count, Variable) and not isinstance(count, Integer):
        options = getattr(count, "options", None)
        if options is not None and not callable(options):
            return int(max(options))
    elif isinstance(count, Fixed):
        return int(count.value)
    elif isinstance(count, int):
        return count
    raise ValueError(f"the count of {path} needs a fixed upper bound to be encoded")


def _leaf_paths(node: Node, path: str, kind: type) -> list[str]:
    """
    Gets the paths of the leaves of a kind of a node (such as its fixed
    values), with `[*]` for the index of any item of a repetition, and with
    the branches of choices, whether a design chooses them or not.
    """
    if isinstance(node, kind):
        return [path]
    if isinstance(node, Group):
        return [
            p
            for name, child in node.children.items()
            for p in _leaf_paths(child, _join(path, name), kind)
        ]
    if isinstance(node, Choice):
        return [
            p
            for name, branch in node.branches.items()
            for p in _leaf_paths(branch, _join(path, name), kind)
        ]
    if isinstance(node, Repeat):
        return _leaf_paths(node.item, f"{path}[*]", kind)
    return []


def _node_at(node: Node, path: Any) -> Node | None:
    """
    Gets the node at a path of leaves (see `_leaf_paths`): through the
    children of groups and the branches of choices by their names, and the
    items of repetitions by `[*]`, or `None` if there is none.
    """
    if not isinstance(path, str) or not path:
        return None
    for segment in path.replace("[*]", ".[*]").split("."):
        if isinstance(node, Group):
            node = node.children.get(segment)
        elif isinstance(node, Choice):
            node = node.branches.get(segment)
        elif isinstance(node, Repeat) and segment == "[*]":
            node = node.item
        else:
            return None
        if node is None:
            return None
    return node


def _missing_segment(node: Node, segments: list) -> tuple[int, list[str]] | None:
    """
    Finds the first segment of a path of values (see `Values`) that no value
    of a node's designs can have, if any: its position in the path, and the
    names that the values there have. A choice is stepped through by the
    name of a branch, or to a value of any of its branches (as of the branch
    it chooses), and a repetition by the index of an item or `count`. A leaf
    (such as a variable or a derived value) may have values of its own (such
    as a fixed group), so that a path within it is not missing.

    Args:
        node (Node): The node (such as a design space).
        segments (list): The segments of the path (see `_segments`).

    Returns:
        tuple[int, list[str]] | None: the position of the missing segment and
            the names there, or `None` if no segment is missing
    """
    # pylint: disable=too-many-return-statements
    if not segments:
        return None
    segment, rest = segments[0], segments[1:]

    def within(child: Node) -> tuple[int, list[str]] | None:
        missing = _missing_segment(child, rest)
        return None if missing is None else (missing[0] + 1, missing[1])

    if isinstance(node, Group):
        child = node.children.get(segment) if isinstance(segment, str) else None
        return (0, list(node.children)) if child is None else within(child)
    if isinstance(node, Choice):
        if segment in node.branches:
            return within(node.branches[segment])
        misses = []
        for branch in node.branches.values():
            missing = _missing_segment(branch, segments)
            if missing is None:
                return None  # a value of this branch
            misses.append(missing)
        deepest = max(misses, key=lambda missing: missing[0])
        if deepest[0] > 0:
            return deepest  # within a branch that has the first segments
        names = list(node.branches)
        for _, there in misses:
            names.extend(name for name in there if name not in names)
        return 0, names
    if isinstance(node, Repeat):
        if segment == "count":
            return None
        if isinstance(segment, int):
            return within(node.item)
        return 0, ["count"]
    return None
