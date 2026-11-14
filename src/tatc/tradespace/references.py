"""
References to the values of designs, by which the parts of observation
studies (see `ObservationStudy`, `Observations`, `MissionPlan`, and
`Failures`) read them: a name (such as `"points_seed"` for `d.points_seed`,
or a path of `Design.flat`, such as `"layout.walker.planes"`), or a function
of the (realized) design (such as `lambda d: d.layout`). A string given
where a reference is expected is always a reference, never a value itself:
give a function that returns a string (or a number, where numbers are
references) to give a value::

    MissionPlan(SCHEDULE, COST_MODEL, kind="field_of_regard")  # d.field_of_regard
    MissionPlan(SCHEDULE, COST_MODEL, kind=lambda d: "imager")  # always "imager"

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numbers
from collections.abc import Callable
from typing import Any

from .space import Design

Reference = str | Callable[[Design], Any]
"""A reference to a value of a design: its name or path (see `Design.flat`),
or a function of the (realized) design that gives it (see `value_of`)."""


def value_of(design: Design, reference: Reference | None) -> Any:
    """
    Gets a value of a design that a reference names: a string names a value
    by its name or path (`design[reference]`, see `Values`: such as
    `"layout.walker.planes"`, the path of a column of results, see
    `Design.flat`, or of a group of values, such as `"cost_risk"`), a
    function computes it from the design, and `None` is `None`.

    Args:
        design (Design): The (realized) design.
        reference (Reference | None): The reference: a string, a function,
            or `None`.

    Returns:
        Any: the value

    Raises:
        KeyError: if the design has no value of a name or path.
    """
    if reference is None:
        return None
    if callable(reference):
        return reference(design)
    if not isinstance(reference, str):
        raise TypeError(f"a reference to a value of a design, not {reference!r}")
    return design[reference]


def _number(design: Design, value: Any) -> Any:
    """A number, or a value of a design that a reference names."""
    if isinstance(value, numbers.Real):
        return value
    return value_of(design, value)
