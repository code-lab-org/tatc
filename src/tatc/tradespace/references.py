"""
Reading designs: references to their values, the dataclasses configured
from them, and the functions given them.

A reference names a value of a design, by which the parts of studies (such
as `ObservationStudy`, `Observations`, `Targets`, `MissionPlan`, and
`Failures`) read the designs they are given (see `value_of`). It is either
a string, which names a value by its name (such as `"points_seed"` for
`d.points_seed`) or path (as the columns of the results of a study name
values, see `Design.flat`, such as `"layout.walker.planes"`, or of a group
of values, such as `"cost_risk"`), or a function of the (realized) design
that computes the value (such as `lambda d: d.layout`). A string given where
a reference is expected is always a reference, never a value itself: give a
function that returns a string to give a string. Where a number is also
expected (such as a growth of costs), a number is the value itself, and
other values are references::

    MissionPlan(SCHEDULE, COST_MODEL, kind="field_of_regard")  # d.field_of_regard
    MissionPlan(SCHEDULE, COST_MODEL, kind=lambda d: "imager")  # always "imager"

A dataclass, such as a `RestorationPolicy` or a `DeploymentSchedule`, is
configured from the values of a design by `configure`.

A function of a study's results (such as the evaluation function of a
`Study`, or the metrics of an `ObservationStudy`) that takes a keyword named
`design` (other than its first parameter, which is given what it is of,
such as what a design builds) is given the realized design as it, as for a
function that depends on values of the design other than those it is
given.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
import inspect
import numbers
from collections.abc import Callable, Iterable, Mapping
from functools import partial
from typing import Any

from .nodes import Node
from .space import Design
from .values import ChoiceValues, Values, _join, _segments, suggestion
from .walks import _missing_segment

Reference = str | Callable[[Design], Any]
"""A reference to a value of a design: its name or path (see `Design.flat`),
or a function of the (realized) design that gives it (see `value_of`)."""


def value_of(design: Design, reference: Reference | None) -> Any:
    """
    Gets a value of a design that a reference names (see the module): a
    string names a value by its name or path (`design[reference]`, see
    `Values`: such as `"layout.walker.planes"`, the path of a column of
    results, see `Design.flat`, or of a group of values, such as
    `"cost_risk"`), a function computes it from the design, and `None` is
    `None`.

    Args:
        design (Design): The (realized) design.
        reference (Reference | None): The reference: a string, a function,
            or `None`.

    Returns:
        Any: the value

    Raises:
        KeyError: if the design has no value of a name or path.
        TypeError: if the reference is neither a string nor a function.
    """
    if reference is None:
        return None
    if callable(reference):
        return reference(design)
    if not isinstance(reference, str):
        raise TypeError(f"a reference to a value of a design, not {reference!r}")
    return design[reference]


def _check_references(space: Node, references: Iterable[tuple[str, Any]]) -> None:
    """
    Checks that the references of the parts of a study (such as the seed of
    its targets, or the type of satellite of its plan) that are names or
    paths name values that the designs of its design space have (in some
    branch of a choice), so that a misspelled name is found when the study
    is defined, rather than when it evaluates a design. References that are
    functions are not checked.

    Args:
        space (Node): The design space.
        references (Iterable[tuple[str, Any]]): The references, each with
            what it is (such as `"the seed of the targets"`).

    Raises:
        ValueError: if a reference names no value of the designs.
    """
    problems = []
    for what, reference in references:
        if not isinstance(reference, str):
            continue
        try:
            segments = _segments(reference)
        except KeyError:
            problems.append(f"{what} {reference!r} is not a name or path of a value")
            continue
        missing = _missing_segment(space, segments)
        if missing is None:
            continue
        position, names = missing
        segment, walked = segments[position], ""
        for step in segments[:position]:
            walked = _join(walked, step)
        where = f"{walked} has" if walked else "the designs have"
        problems.append(
            f"{what} {reference!r}: {where} no value {segment!r} (there are "
            + f"{names}{suggestion(segment, names)})"
        )
    if problems:
        raise ValueError(
            "references to values that the designs of the design space do not "
            + "have: "
            + "; ".join(problems)
        )


def _number(design: Design, value: Any) -> Any:
    """A number, or a value of a design that a reference names."""
    if isinstance(value, numbers.Real):
        return value
    return value_of(design, value)


def configure(base: Any, values: Any = None, **overrides: Any) -> Any:
    """
    Configures a dataclass (such as a `RestorationPolicy` or a
    `DeploymentSchedule`) from the values of a design: the fields of `base`
    that the values name (such as those of `policy_space`, or design
    variables named after fields) are replaced by them, as are those given
    as keyword arguments. Values that name no field are ignored::

        policy = configure(RestorationPolicy(), design.policy)
        schedule = configure(SCHEDULE, design, years=6)

    Args:
        base (Any): The dataclass instance, whose other fields are kept.
        values (Any): The values of a design: a `Design`, the `Values` of
            a group, a mapping, or the `ChoiceValues` of a choice (whose
            branch is a group, or a fixed `None`), or `None`.
        **overrides (Any): Other values of fields.

    Returns:
        Any: a copy of `base` with the values of its fields

    Raises:
        TypeError: if `base` is not a dataclass instance, or `values` are
            not values of a design (or a choice's branch is not a group).
        ValueError: if an override names no field of `base`.
    """
    if not dataclasses.is_dataclass(base) or isinstance(base, type):
        raise TypeError("configure requires a dataclass instance")
    if isinstance(values, ChoiceValues):
        branch = values.value
        if branch is not None and not isinstance(branch, Values):
            raise TypeError(
                "configure requires the values of a group: the choice chose "
                + f"{values.choice!r}, whose value is {branch!r}"
            )
        values = branch
    elif isinstance(values, Design):
        values = values.values
    elif values is not None and not isinstance(values, Mapping):
        raise TypeError(f"configure requires the values of a design, not {values!r}")
    names = {field.name for field in dataclasses.fields(base)}
    changes = {name: values[name] for name in (values or {}) if name in names}
    unknown = set(overrides) - names
    if unknown:
        raise ValueError(f"no fields named {sorted(unknown)}")
    return dataclasses.replace(base, **changes, **overrides)


def _takes_design(function: Callable[..., Any]) -> bool:
    """
    Whether a function takes a keyword named `design` (a parameter named
    `design` that may be given by keyword), by which the parts of studies
    give it the realized design (see `_given_design`), other than its first
    parameter, which is given what the function is of (such as what a design
    builds, or its observations), whatever its name.
    """
    try:
        parameters = list(inspect.signature(function).parameters.values())
    except (TypeError, ValueError):  # such as of some built-in functions
        return False
    return any(
        parameter.name == "design"
        and (
            parameter.kind == inspect.Parameter.KEYWORD_ONLY
            or (parameter.kind == inspect.Parameter.POSITIONAL_OR_KEYWORD and i > 0)
        )
        for i, parameter in enumerate(parameters)
    )


def _given_design(function: Callable[..., Any], design: Design) -> Callable[..., Any]:
    """
    Gets a function of its other arguments: the function given the
    (realized) design as its keyword `design`, if it takes one (see
    `_takes_design`), or the function itself::

        outcome = _given_design(evaluate, realized)(built)
    """
    return partial(function, design=design) if _takes_design(function) else function
