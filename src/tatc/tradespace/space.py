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

import numbers
import zlib
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from functools import partial
from typing import Any, Literal

import numpy as np
import pandas as pd

from .hashing import canonical, digest
from .nodes import (
    Choice,
    Fixed,
    Group,
    Infeasible,
    Uncertain,
    _as_node,
)
from .values import Scope, Values, generic_path, missing_value
from .walks import (
    _ABSENT,
    _Context,
    _decided,
    _flatten,
    _leaf_paths,
    _level_of,
    _node_at,
    _Odometer,
    _walk,
)
from .walks import _unit_paths as _paths_of_units

STATUS_COLUMNS = ("design", "replicate", "scenario", "feasible", "error", "seconds")
"""tuple[str, ...]: The columns of the results of a study (see `Study.run`)
that are neither values of designs nor metrics, in their order: the key of
the design, the replicate, the name of the scenario (of a study with
scenarios only), whether the evaluation is feasible (and its error if not),
and the seconds it took to build and evaluate the design. Neither metrics
nor the children of design spaces may be named as them."""

RESERVED_NAMES = frozenset(
    {
        # names that the views of values give otherwise (see `Values`)
        "index",
        "count",
        "choice",
        "get",
        "items",
        "keys",
        "values",
        # attributes of designs (see `Design`)
        "decisions",
        "flat",
        "key",
        "realize",
        "replicate",
        "scenario",
        "space",
        "uncertain_values",
        # columns of the results of a study
        *STATUS_COLUMNS,
    }
)
"""frozenset[str]: The names that the children of a design space may not
have: those that designs and the views of their values give otherwise (as
attributes, such as `key`, `items`, or `index`), and the columns of the
results of a study that are not values of designs (see `STATUS_COLUMNS`)."""


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
        scenario: str | None = None,
    ):
        # pylint: disable=too-many-arguments,too-many-positional-arguments
        self.space, self.replicate = space, replicate
        self.decisions = decisions
        """dict: The decisions of the design, as the design space declares
        them: those of a group are a mapping of the decisions of its
        children, by name (omitting those without decisions: fixed,
        derived, and uncertain values); that of a design variable is its
        value; that of a choice is a list of the name of its branch and the
        decisions of the branch (`None` if it has none), such as
        `["walker", {"planes": 3}]`; and that of a repetition is a list of
        the decisions of its items (`None` for an item without decisions),
        as many as its count, whose unordered items are in the order of the
        levels of their decisions (see `Repeat`)."""
        self._values = values
        self._scenario = scenario
        self.key = digest(canonical(decisions), 12)
        """str: An identifier of the design, from its decisions only (not
        from its derived values or uncertain parameters)."""

    @property
    def scenario(self) -> str | None:
        """str | None: The name of the scenario in which the design is
        realized (see `realize`), or `None` if it is not realized in one."""
        return self._scenario

    @property
    def values(self) -> Values:
        """Values: the values of the design."""
        return Values(self._values)

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        try:
            return getattr(self.values, name)
        except AttributeError as error:
            raise AttributeError(
                f"the design has {error}", name=name, obj=self
            ) from None

    def __dir__(self) -> list[str]:
        # the values, for completion (and suggestions of misspelled names)
        return sorted({*object.__dir__(self), *self._values})

    def __getitem__(self, name: str) -> Any:
        if (
            isinstance(name, str)
            and name not in self._values
            and not ("." in name or "[" in name)
        ):
            raise KeyError(f"the design has {missing_value(name, self._values)}")
        return self.values[name]

    def __contains__(self, name: object) -> bool:
        return name in self.values

    __iter__ = None  # type: ignore[assignment]

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, Design)
            and self.key == other.key
            and self.replicate == other.replicate
            and self._scenario == other._scenario
        )

    def __hash__(self) -> int:
        return hash((self.key, self.replicate, self._scenario))

    def __fingerprint__(self) -> dict:
        """
        Gets what identifies the design in a definition (see
        `tatc.tradespace.fingerprint`), such as of a design that a function
        of a study refers to: its design space, its decisions, and its
        replicate and scenario (if it is realized).

        Returns:
            dict: the design space, decisions, replicate, and scenario
        """
        return {
            "space": self.space,
            "decisions": self.decisions,
            "replicate": self.replicate,
            "scenario": self._scenario,
        }

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

    def uncertain_values(self) -> dict[str, Any]:
        """
        Gets the values of the uncertain parameters of the design (see
        `Uncertain`) by path (as of `flat`), in the order they are drawn:
        those of the branches that it chooses and of the items of its
        repetitions, as drawn in its replicate (or given by its scenario) if
        it is realized (see `realize`), and their nominal values otherwise.
        Their paths with `[*]` for the indices of items (see
        `generic_path`) are those of `DesignSpace.uncertain_paths`::

            design.realize(3).uncertain_values()  # {"planes[0].phase": 41.2, ...}

        Returns:
            dict[str, Any]: the values, by path
        """
        paths = set(self.space.uncertain_paths())
        return {
            path: value
            for path, value in self.flat().items()
            if generic_path(path) in paths
        }

    def realize(
        self,
        replicate: int,
        *,
        seed: int = 0,
        common_random_numbers: bool = True,
        shift: Mapping[str, int] | None = None,
        overrides: Mapping[str, Any] | None = None,
        scenario: str | None = None,
    ) -> Design:
        """
        Realizes a replicate of the design: draws its uncertain parameters
        (and computes the values derived from them).

        Each uncertain parameter is drawn from a random number generator
        seeded by `seed`, `replicate`, and its path, so that with common
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

        The `overrides` of a scenario give fixed values and uncertain
        parameters other values (by their paths, with `[*]` for the index of
        any item of a repetition, as of `DesignSpace.uncertain_paths`): a
        condition under which designs are compared, such as a requirement,
        rather than averaged over (see `Study`). Values derived from them are
        derived from the overrides, an uncertain parameter that they give is
        not drawn, and the constraints of the design space are checked
        again, as the overrides may change what they depend on. The realized
        design is of the scenario named `scenario` (see `Design.scenario`).

        Args:
            replicate (int): The replicate (a non-negative integer).
            seed (int): The seed of the trade study (a non-negative integer).
            common_random_numbers (bool): `True` to use common random numbers
                between designs.
            shift (Mapping[str, int] | None): The offsets of the replicates
                as in which the named uncertain parameters are drawn (each a
                non-negative integer), if any.
            overrides (Mapping[str, Any] | None): The values of fixed values
                and uncertain parameters in a scenario (by path), if any.
            scenario (str | None): The name of the scenario of the overrides,
                if any.

        Returns:
            Design: the realized design

        Raises:
            Infeasible: if the design violates a constraint with the
                overrides.
            ValueError: if `replicate`, `seed`, or an offset of `shift` is
                not a non-negative integer, `shift` names a path that is not
                of an uncertain parameter of the design space, or
                `overrides` one that is not of a fixed value or an uncertain
                parameter.
        """
        # pylint: disable=too-many-arguments
        _check_seed(seed)
        replicate = _check_replicate(replicate)
        shift = _check_shift(shift)
        self.space._check_uncertain(shift)  # pylint: disable=protected-access
        overrides = dict(overrides or {})
        self.space._check_overrides(overrides)  # pylint: disable=protected-access
        tail = [] if common_random_numbers else [int(self.key, 16)]

        def draw(path: str, node: Uncertain) -> Any:
            name = generic_path(path)
            if name in overrides:
                return overrides[name]
            drawn_in = replicate + shift.get(name, 0)
            rng = np.random.default_rng(
                [seed, drawn_in, *tail, zlib.crc32(path.encode())]
            )
            return node.sampler(rng)

        def fix(path: str, node: Fixed) -> Any:
            return overrides.get(generic_path(path), node.value)

        ctx = _Context(
            draw=draw,
            check=bool(overrides),
            fix=fix if overrides else None,
            canonical=True,
        )
        values, _ = _walk(self.space, "", Scope(), ctx, self.decisions)
        return Design(self.space, self.decisions, values, replicate, scenario)

    def __repr__(self) -> str:
        replicate = "" if self.replicate is None else f", replicate={self.replicate}"
        scenario = "" if self._scenario is None else f", scenario={self._scenario!r}"
        return f"Design({self.key}{replicate}{scenario}: {self.flat()})"


def _check_index(value: Any, what: str) -> int:
    """
    Checks that a value that seeds uncertain parameters (a seed, a replicate,
    or the offset of a shift) is a non-negative integer, as numpy's
    generators require, raising a `ValueError` if not (naming it `what`),
    and returns it as a Python integer.
    """
    if (
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, numbers.Integral)
        or value < 0
    ):
        raise ValueError(f"{what} must be a non-negative integer, not {value!r}")
    return int(value)


def _check_seed(seed: Any) -> None:
    """
    Checks that a seed of uncertain parameters is a non-negative integer (see
    `_check_index`).
    """
    _check_index(seed, "the seed")


def _check_replicate(replicate: Any) -> int:
    """
    Checks that a replicate is a non-negative integer (see `_check_index`),
    returning it.
    """
    return _check_index(replicate, "a replicate")


def _check_shift(shift: Mapping[str, Any] | None) -> dict[str, int]:
    """
    Checks that the offsets of a shift (see `Design.realize`) are
    non-negative integers (see `_check_index`), returning them.
    """
    return {
        name: _check_index(offset, f"the offset of {name!r} in a shift")
        for name, offset in (shift or {}).items()
    }


def _as_constraints(constraints: Any) -> list[Callable[[Scope], Any]]:
    """
    Gets the constraints of a design space given as a constraint or a
    sequence of them, raising a `TypeError` if one is not callable.
    """
    found = [constraints] if callable(constraints) else list(constraints)
    for constraint in found:
        if not callable(constraint):
            raise TypeError(f"a constraint is a callable, not {constraint!r}")
    return found


def _lhs_column(n: int, rng: np.random.Generator) -> np.ndarray:
    """
    Gets `n` values in [0, 1) stratified in `n` equal intervals in random
    order (a column of a Latin hypercube sample).
    """
    return (rng.permutation(n) + rng.random(n)) / n


class DesignSpace(Group):
    """
    A group of nodes: the root of a design space, or a group within it (such
    as the design of each orbit plane). Children are given by name, in the
    order they are decided; plain values are fixed values. Children may not
    be named as the attributes of designs and of the views of their values
    (such as `key`, `items`, or `index`), nor as the columns of the results
    of a study that are not values of designs (see `STATUS_COLUMNS`): see
    `RESERVED_NAMES`.

    Constraints are callables of the group's values (and those of enclosing
    groups), checked once the group's values are decided, that return
    whether a design satisfies them (`True` or `False`) or by how much it
    violates them: a number that is not positive if it satisfies them, and
    otherwise the amount of violation (for example, `lambda d: d.satellites -
    12` for at most 12 satellites), which guides constraint-aware optimizers
    toward feasible designs (see `optimize_optuna` and `make_pymoo_problem`).
    They must not depend on uncertain parameters.

    Args:
        constraints (Callable | Sequence[Callable]): A constraint, or the
            constraints.
        **children (Any): The nodes, by name.
    """

    def __init__(
        self,
        constraints: Callable[[Scope], Any] | Sequence[Callable[[Scope], Any]] = (),
        **children: Any,
    ):
        reserved = RESERVED_NAMES.intersection(children)
        if reserved:
            raise ValueError(
                f"children may not be named {sorted(reserved)}: designs and the "
                + "views of their values give these names otherwise (as attributes, "
                + "such as `key` and `items`), as do the results of a study (as "
                + f"columns); the reserved names are {sorted(RESERVED_NAMES)}"
            )
        self.children = {name: _as_node(node) for name, node in children.items()}
        self.constraints = _as_constraints(constraints)

    def __repr__(self) -> str:
        return f"DesignSpace({', '.join(self.children)})"

    def _design(self, ctx: _Context, prescribed: Any = _ABSENT) -> Design:
        """
        Decides a design (with nominal uncertain parameters), following the
        prescribed decisions, if any.
        """
        values, decisions = _walk(self, "", Scope(), ctx, prescribed)
        return Design(self, decisions, values)

    def design(
        self, decisions: Mapping | Design | None = None, /, **overrides: Any
    ) -> Design:
        """
        Gets the design of given decisions (see `Design.decisions` for their
        shapes), or of the decisions of another design (such as of a design
        space this one extends, see `extend`), with the decisions of some
        children replaced (by name: a value, the decisions of a group or a
        choice, or a design of their own space), checking its constraints::

            space = base.extend(policy=policy_space(orbit_spares=Integer(0, 1)))
            design = space.design(base_design, policy=policy_design)
            sso = space.design(base_design, layout="sso")

        A choice may also be given the name of its branch alone, here and in
        the decisions: the branch's decisions are then those given (of a
        replaced child that chose the same branch) or, otherwise, the first
        levels of its variables (and the lower bounds of continuous ones).

        Args:
            decisions (Mapping | Design | None): The decisions, or a design
                whose decisions to take.
            **overrides (Any): The decisions of children, by name.

        Returns:
            Design: the design

        Raises:
            Infeasible: if the design violates a constraint.
            ValueError: if an override names no child of the design space,
                the decisions of a group name others than its children, or a
                decision is not a value of its variable's domain.
        """
        if isinstance(decisions, Design):
            decisions = decisions.decisions
        merged = dict(decisions or {})
        unknown = sorted(set(overrides) - set(self.children))
        if unknown:
            raise ValueError(f"no children {unknown} of the design space")
        for name, decision in overrides.items():
            if isinstance(decision, Design):
                decision = decision.decisions
            given = merged.get(name)
            if (
                isinstance(self.children[name], Choice)
                and isinstance(decision, str)
                and isinstance(given, (list, tuple))
                and len(given) == 2
                and given[0] == decision
            ):
                continue  # the branch, with the decisions given of it
            merged[name] = decision
        values, merged = _walk(self, "", Scope(), _Context(normalize=_level_of), merged)
        return Design(self, merged, values)

    def extend(
        self,
        constraints: Callable[[Scope], Any] | Sequence[Callable[[Scope], Any]] = (),
        **children: Any,
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
            constraints (Callable | Sequence[Callable]): Another constraint,
                or other constraints.
            **children (Any): The nodes, by name.

        Returns:
            DesignSpace: the design space
        """
        return DesignSpace(
            [*self.constraints, *_as_constraints(constraints)],
            **{**self.children, **children},
        )

    def designs_of(
        self, frame: pd.DataFrame, *, keys: Iterable[str] | None = None
    ) -> pd.Series:
        """
        Rebuilds designs from the rows of a data frame of their values by path
        (see `Design.flat`): the results of a study (see `Study.run`), or a
        summary that keeps the values of the designs (see `summarize`). The
        designs are those of the rows, aligned with them (with the frame's
        index), so that functions of designs map over the rows, or those of
        `keys` (indexed by them). Each design is rebuilt once, from the first
        row of its key (`design`, a column or the name of the index), whose
        key it must have, if the frame has keys::

            front = space.designs_of(summary[summary.pareto])
            realized = [study.realize(design) for design in front]

        Args:
            frame (pandas.DataFrame): The rows.
            keys (Iterable[str] | None): The keys of the designs to rebuild,
                or `None` for the design of each row of the frame.

        Returns:
            pandas.Series: the designs, by row (with the frame's index) or by
                key (of `keys`, each once)

        Raises:
            KeyError: if the frame has no rows of a key.
            ValueError: if the frame lacks the value of a decision of a
                design, or its values give a design of another key.
        """
        rows = frame
        if "design" not in rows.columns and rows.index.name == "design":
            rows = rows.reset_index()
        if "design" not in rows.columns:
            if keys is not None:
                raise ValueError("the frame has no keys of designs (`design`)")
            designs = [self._rebuilt(record) for record in rows.to_dict("records")]
            return pd.Series(designs, index=frame.index, name="design", dtype=object)
        first = rows.drop_duplicates("design")
        records = dict(zip(first["design"], first.to_dict("records")))
        wanted = list(rows["design"]) if keys is None else list(dict.fromkeys(keys))
        missing = [key for key in dict.fromkeys(wanted) if key not in records]
        if missing:
            raise KeyError(f"no rows of designs {missing[:5]}")
        rebuilt = {
            key: self._rebuilt(records[key], key) for key in dict.fromkeys(wanted)
        }
        designs = [rebuilt[key] for key in wanted]
        index = frame.index if keys is None else pd.Index(wanted, name="design")
        return pd.Series(designs, index=index, name="design", dtype=object)

    def _rebuilt(self, record: Mapping[str, Any], key: str | None = None) -> Design:
        """
        Rebuilds a design from its values by path (see `designs_of`), which
        must give its key, if given.
        """
        design = self._design(_Context(choose=partial(_decided, record)))
        if key is not None and design.key != key:
            raise ValueError(
                f"the values of design {key} give another design ({design.key}): "
                + "are they of another design space?"
            )
        return design

    def enumerate(self, designs: Iterable[Design] | None = None) -> list[Design]:
        """
        Enumerates every feasible design (a full factorial design over the
        levels of the design variables, including those that exist only in
        some designs), without duplicates, or every feasible design that
        extends given designs: designs of a design space that this one
        extends (see `extend`), each with every combination of the levels of
        the decisions it lacks, such as each policy for each of the efficient
        designs of the space without policies::

            failure_space = cost_space.extend(policy=POLICY)
            with_policies = failure_space.enumerate(cost_space.designs_of(efficient))

        Args:
            designs (Iterable[Design] | None): The designs to extend, or
                `None` for every design.

        Returns:
            list[Design]: the designs, in the order of their decisions (and of
                the given designs)

        Raises:
            ValueError: if a design variable is continuous, or the decisions
                of a given design are not those of a design of this design
                space (see `design`).
        """
        if designs is None:
            return list(self._iterate())
        found: dict[str, Design] = {}
        for given in designs:
            if not isinstance(given, Design):
                raise TypeError(f"designs to extend, not {given!r}")
            for design in self._iterate(given.decisions):
                found.setdefault(design.key, design)
        return list(found.values())

    def _iterate(self, prescribed: Any = _ABSENT) -> Iterator[Design]:
        """
        Iterates over every feasible design (see `enumerate`), or every one
        with the prescribed decisions (of some of its nodes). The items of
        unordered repetitions are walked in only one of their orders, so that
        their permutations are not walked.
        """
        odometer = _Odometer()
        ctx = _Context(
            choose=odometer.choose,
            odometer=odometer,
            normalize=None if prescribed is _ABSENT else _level_of,
        )
        keys: set[str] = set()
        while True:
            try:
                design = self._design(ctx, prescribed)
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
        *,
        method: Literal["lhs", "random"] = "lhs",
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

    def uncertain_paths(self) -> list[str]:
        """
        Gets the paths of the uncertain parameters of the design space (see
        `Uncertain`), in the order they are drawn, with `[*]` for the index
        of any item of a repetition (so that the phases of every plane of a
        constellation are one path, such as `planes[*].phase`), and with
        the branches of choices (such as `layout.walker.phase_fraction`),
        whether a design chooses them or not. These name the parameters of a
        shift (see `Design.realize`) and the overrides of a scenario (see
        `Study`).

        Returns:
            list[str]: the paths
        """
        return _leaf_paths(self, "", Uncertain)

    def _check_uncertain(self, paths: Iterable[str]) -> None:
        """
        Checks that paths are of uncertain parameters of the design space (as
        of `uncertain_paths`, with `[*]` for the index of any item of a
        repetition), such as those of a shift (see `Design.realize`), raising
        a `ValueError` if one is not.
        """
        unknown = [p for p in paths if not isinstance(_node_at(self, p), Uncertain)]
        if unknown:
            raise ValueError(
                "not uncertain parameters of the design space: "
                + f"{sorted(unknown, key=str)} (it has {self.uncertain_paths()})"
            )

    def _check_overrides(self, overrides: Iterable[str]) -> None:
        """
        Checks that the paths of the overrides of a scenario (a mapping by
        path, or the paths, see `Design.realize`) are of fixed values or
        uncertain parameters of the design space (with `[*]` for the index of
        any item of a repetition), raising a `ValueError` if one is not.
        """
        unknown = [
            p
            for p in overrides
            if not isinstance(_node_at(self, p), (Fixed, Uncertain))
        ]
        if unknown:
            known = [*_leaf_paths(self, "", Fixed), *self.uncertain_paths()]
            raise ValueError(
                "not fixed values or uncertain parameters of the design space: "
                + f"{sorted(unknown, key=str)} (it has {known})"
            )

    def _unit_paths(self) -> list[str]:
        """
        Gets the paths of every decision of the design space, in the order
        they are decided: its design variables, its choices, and the counts
        of its repetitions, with the items of each repetition up to its
        largest count. These are the dimensions of the unit hypercube that
        encodes its designs for optimizers of fixed numbers of continuous
        variables (see `_decode` and `make_pymoo_problem`), raising a
        `ValueError` if the count of a repetition has an upper bound that is
        not fixed.
        """
        return _paths_of_units(self, "")

    def _decode(self, unit: Sequence[float]) -> Design:
        """
        Decodes a design from a point of the unit hypercube (a value in
        [0, 1] for each of `_unit_paths`): each decision maps the value of
        its path to a level (each level spans an equal part of [0, 1]) or to
        a point of its interval. Values of the paths of inactive nodes (such
        as the branches not chosen) are ignored. Raises `Infeasible` if the
        design violates a constraint.
        """
        index = {path: i for i, path in enumerate(self._unit_paths())}
        if len(unit) != len(index):
            raise ValueError(f"expected {len(index)} values, not {len(unit)}")
        return self._design(
            _Context(choose=lambda path, domain: domain.from_unit(unit[index[path]]))
        )

    def _suggest(self, trial: Any) -> Design:
        """
        Decides a design with an Optuna trial (define-by-run, see
        `optimize_optuna`): each decision is suggested by the trial, named by
        its path, with the kind of its variable (`suggest_float`, with the
        step of a discrete `Real`; `suggest_int`; or `suggest_categorical`,
        also for choices), so that Optuna's samplers search the hierarchical
        design space directly. Raises `Infeasible` if the design violates a
        constraint.
        """
        return self._design(
            _Context(choose=lambda path, domain: domain.suggest(trial, path))
        )


def _design_items(designs: Iterable[Design | str] | pd.DataFrame) -> list[Design | str]:
    """
    Gets designs given as designs, keys, or a data frame of designs (with
    their keys, `design`, as a column or the name of the index), in order.
    """
    if isinstance(designs, pd.DataFrame):
        if "design" in designs.columns:
            items = list(designs["design"])
        elif designs.index.name == "design":
            items = list(designs.index)
        else:
            raise ValueError("the data frame has no keys of designs (`design`)")
    else:
        items = list(designs)
    for item in items:
        if not isinstance(item, (Design, str)):
            raise TypeError(f"a design or the key of one, not {item!r}")
    return items


def _keys_of(designs: Iterable[Design | str] | pd.DataFrame) -> list[str]:
    """
    Gets the keys of designs given as designs, keys, or a data frame of them
    (with keys `design`, a column or the name of the index), each once.
    """
    keys = [d.key if isinstance(d, Design) else d for d in _design_items(designs)]
    return list(dict.fromkeys(keys))


def _in_space(space: DesignSpace, design: Design) -> Design:
    """
    Gets the design of a design space of the decisions of a design of
    another (see `DesignSpace.design`), raising a `ValueError` if they are
    not decisions of the design space.
    """
    try:
        return space.design(design)
    except ValueError as error:
        raise ValueError(
            f"design {design.key} is of another design space, whose decisions are "
            + f"not those of a design of this one: {error}"
        ) from error


def _designs(
    space: DesignSpace,
    designs: Iterable[Design | str] | pd.DataFrame,
    *frames: pd.DataFrame,
) -> list[Design]:
    """
    Gets designs of a design space given as designs, keys, or a data frame
    of designs (with their keys, `design`), each once, in order: designs of
    other design spaces (such as of a design space that this one extends)
    are the designs of this one of their decisions (see
    `DesignSpace.design`), and designs given by their keys are rebuilt from
    the values of the frame (if it has them) or of other frames (such as the
    results of a study, see `DesignSpace.designs_of`).
    """
    if isinstance(designs, pd.DataFrame):
        frames = (designs, *frames)
    found: dict[str, Design | None] = {}
    for item in _design_items(designs):
        if isinstance(item, Design):
            design = item if item.space is space else _in_space(space, item)
            found[design.key] = design
        else:
            found.setdefault(item, None)
    error: ValueError | None = None
    for frame in frames:
        missing = [key for key, design in found.items() if design is None]
        if not missing:
            break
        rows = frame.reset_index() if "design" not in frame else frame
        if "design" not in rows:
            continue
        keys = set(rows["design"])
        present = [key for key in missing if key in keys]
        try:
            for design in space.designs_of(rows, keys=present):
                found[design.key] = design
        except ValueError as failure:
            error = failure  # the frame lacks values of the designs: try the next
    missing = [key for key, design in found.items() if design is None]
    if missing:
        raise ValueError(
            f"the values of designs {missing[:5]} are not given: give the designs, "
            + "or their keys with the results of a study (or a summary that keeps "
            + "the values of designs, see summarize)"
        ) from error
    return [design for design in found.values() if design is not None]
