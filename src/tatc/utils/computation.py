"""
Computations run together: a computation (`TimeRequest`) yields each
Skyfield time it creates before using it, so that the Earth orientation
quantities of the times of several computations can be computed together
(see `_run_together`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Generator
from typing import Any

from skyfield.timelib import Time

from .earth_orientation import _share_earth_orientation

TimeRequest = Generator[Time | list[Time], None, Any]
"""
A computation that yields each Skyfield `Time` (or list of `Time`s) it
creates before using it, so that the costly Earth orientation quantities of
the times of several computations can be computed together (see
`_run_together`), and returns its result.
"""


def _run_together(computations: list[TimeRequest]) -> list[Any]:
    """
    Runs several computations (see `TimeRequest`) together, in steps: at each
    step, the Earth orientation quantities of all times that the computations
    yield are computed together (see `_share_earth_orientation`), and each
    computation then continues to its next time.

    Args:
        computations (list[TimeRequest]): The computations.

    Returns:
        list[Any]: the result of each computation
    """
    results: dict[int, Any] = {}
    pending: dict[int, Time | list[Time]] = {}

    def advance(i: int, first: bool = False) -> None:
        try:
            pending[i] = next(computations[i]) if first else computations[i].send(None)
        except StopIteration as stop:
            results[i] = stop.value
            pending.pop(i, None)

    for i in range(len(computations)):
        advance(i, first=True)
    while len(pending) > 0:
        _share_earth_orientation(
            [
                t
                for request in pending.values()
                for t in (request if isinstance(request, list) else [request])
            ]
        )
        for i in list(pending):
            advance(i)
    return [results[i] for i in range(len(computations))]


def _run(computation: TimeRequest) -> Any:
    """
    Runs a single computation (see `TimeRequest`).

    Args:
        computation (TimeRequest): The computation.

    Returns:
        Any: its result
    """
    return _run_together([computation])[0]
