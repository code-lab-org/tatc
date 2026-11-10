"""
Building blocks of trade studies of constellations: fragments of design
spaces for common constellation layouts (Walker constellations and planes of
sun-synchronous orbits), a builder of the satellites of a layout with a
payload (instruments, or entries of the instrument catalog), and counts of
satellites and planes for objectives and constraints.

A fragment is a `DesignSpace` to place in a larger design space, for example
as a branch of a `Choice` between layouts::

    space = DesignSpace(
        altitude=Real(500e3, 800e3, step=100e3),
        sounder=Categorical(["atms", "cris"]),
        layout=Choice(
            walker=walker_space(inclination=Real(45, 75, step=15)),
            sso=sso_planes_space(),
        ),
    )

    def build(d):
        return build_constellation(
            d.layout, d.sounder, epoch=START, altitude=d.altitude
        )

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from datetime import datetime, time
from typing import Any

import numpy as np

from ..schemas import (
    CircularOrbit,
    Instrument,
    Satellite,
    SunSynchronousOrbit,
    WalkerConstellation,
)
from .nodes import Fixed, Integer, Real, Repeat
from .space import DesignSpace
from .values import ChoiceValues, Scope, Values

WALKER = "walker"
"""The kind of the layouts of `walker_space`."""
SSO_PLANES = "sso planes"
"""The kind of the layouts of `sso_planes_space`."""


def _node(value: Any, default: Any) -> Any:
    return default if value is None else value


def walker_space(
    *,
    inclination: Any = None,
    planes: Any = None,
    per_plane: Any = None,
    relative_spacing: Any = None,
    altitude: Any = None,
    configuration: Any = "delta",
    constraints: Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of Walker constellations (see
    `tatc.schemas.WalkerConstellation`): their inclination (degrees), number
    of planes, number of satellites in each plane, relative spacing of
    satellites between planes, and (optionally, if not given to
    `build_constellation`) altitude (meters). Each is a node of a design
    space (such as a design variable) or a value.

    Args:
        inclination (Any): The inclination, by default a variable from 30 to
            90 degrees in steps of 15.
        planes (Any): The number of planes, by default a variable from 1 to 6.
        per_plane (Any): The number of satellites in each plane, by default a
            variable from 1 to 6.
        relative_spacing (Any): The relative spacing, by default a variable
            from 0 to the number of planes less 1.
        altitude (Any): The altitude, if part of the layout.
        configuration (Any): The configuration (`"delta"` or `"star"`).
        constraints (Sequence[Callable]): Constraints of the fragment.

    Returns:
        DesignSpace: the fragment
    """
    children = {
        "kind": Fixed(WALKER),
        "configuration": configuration,
        "inclination": _node(inclination, Real(30, 90, step=15)),
        "planes": _node(planes, Integer(1, 6)),
        "per_plane": _node(per_plane, Integer(1, 6)),
        "relative_spacing": _node(relative_spacing, Integer(0, lambda d: d.planes - 1)),
    }
    if altitude is not None:
        children["altitude"] = altitude
    return DesignSpace(constraints=constraints, **children)


def _distinct_crossings(d: Scope) -> bool:
    """Whether the planes of a layout have different crossing times."""
    crossings = [round(plane.crossing, 6) for plane in d.planes]
    return len(set(crossings)) == len(crossings)


def sso_planes_space(
    *,
    crossing: Any = None,
    satellites: Any = None,
    count: Any = None,
    phase: Any = 0.0,
    altitude: Any = None,
    ascending: Any = True,
    ordered: bool = False,
    distinct: bool = True,
    constraints: Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of planes of sun-synchronous orbits
    (see `tatc.schemas.SunSynchronousOrbit`), each with its own equator
    crossing time (local solar time, in hours), number of satellites (evenly
    spaced in the plane), and phase (the true anomaly of its first satellite,
    degrees), and (optionally, if not given to `build_constellation`) their
    altitude (meters). Each is a node of a design space (such as a design
    variable, a value derived from the index of the plane, or an uncertain
    parameter) or a value.

    The planes are unordered by default: designs whose planes differ only in
    their order are the same design (which requires that their values do
    not depend on their index). Planes with the same crossing time would be
    the same plane, so a constraint requires different crossing times
    (unless `distinct` is `False`).

    Args:
        crossing (Any): The crossing time of each plane, by default a variable
            from 6 to 18 hours in steps of 1.5.
        satellites (Any): The number of satellites of each plane, by default a
            variable from 1 to 4.
        count (Any): The number of planes, by default a variable from 1 to 3.
        phase (Any): The phase of each plane (an uncertain parameter, such as
            `Uncertain.uniform(0, 360)`, if it is not designed).
        altitude (Any): The altitude, if part of the layout.
        ascending (Any): Whether the crossing times are of ascending nodes.
        ordered (bool): `True` if the planes are ordered (for example, if
            their values depend on their index).
        distinct (bool): `True` to require different crossing times.
        constraints (Sequence[Callable]): Other constraints of the fragment.

    Returns:
        DesignSpace: the fragment
    """
    plane = DesignSpace(
        crossing=_node(crossing, Real(6, 18, step=1.5)),
        satellites=_node(satellites, Integer(1, 4)),
        phase=phase,
    )
    children = {
        "kind": Fixed(SSO_PLANES),
        "ascending": ascending,
        "planes": Repeat(plane, count=_node(count, Integer(1, 3)), ordered=ordered),
    }
    if altitude is not None:
        children["altitude"] = altitude
    constraints = ([_distinct_crossings] if distinct else []) + list(constraints)
    return DesignSpace(constraints=constraints, **children)


def _layout(layout: Any) -> Any:
    """Gets the values of a layout (the chosen branch of a choice)."""
    return layout.value if isinstance(layout, ChoiceValues) else layout


def _kind(layout: Any) -> str:
    try:
        return layout.kind
    except AttributeError:
        raise ValueError(
            "a layout must be the values of walker_space or sso_planes_space"
        ) from None


def count_satellites(layout: Any) -> int:
    """
    Counts the satellites of a layout (the values of `walker_space` or
    `sso_planes_space`, or of a choice of them), for objectives and
    constraints.

    Args:
        layout (Any): The values of the layout.

    Returns:
        int: the number of satellites
    """
    layout = _layout(layout)
    if _kind(layout) == WALKER:
        return int(layout.planes * layout.per_plane)
    return int(sum(plane.satellites for plane in layout.planes))


def count_planes(layout: Any) -> int:
    """
    Counts the orbit planes of a layout (the values of `walker_space` or
    `sso_planes_space`, or of a choice of them), each of which needs its own
    launch.

    Args:
        layout (Any): The values of the layout.

    Returns:
        int: the number of planes
    """
    layout = _layout(layout)
    if _kind(layout) == WALKER:
        return int(layout.planes)
    return len(layout.planes)


def hours_to_time(hours: float) -> time:
    """
    Converts a local time in hours (such as 13.5) to a time of day (13:30),
    to the nearest second.

    Args:
        hours (float): The time (hours, wrapped to [0, 24)).

    Returns:
        datetime.time: the time of day
    """
    seconds = round(hours * 3600) % 86400
    return time(seconds // 3600, seconds // 60 % 60, seconds % 60)


Payload = (
    Instrument
    | Sequence[Instrument]
    | str
    | Sequence[str]
    | Callable[[str, Any], Satellite]
)
"""A payload: instruments, catalog entries, or a function of a name and an
orbit that builds a satellite."""


def _make_satellite(
    payload: Any, name: str, orbit: Any, time_step: float | None
) -> Satellite:
    """
    Makes a satellite with a payload in an orbit (see `build_constellation`).
    """
    if callable(payload) and not isinstance(payload, Instrument):
        return payload(name, orbit)
    payloads = [payload] if isinstance(payload, (str, Instrument)) else list(payload)
    if all(isinstance(item, Instrument) for item in payloads):
        return Satellite(name=name, orbit=orbit, instruments=payloads)
    if not all(isinstance(item, str) for item in payloads):
        raise ValueError("a payload mixes instruments and catalog entries")
    # pylint: disable=import-outside-toplevel
    from ..catalog import get_entry

    instruments, attitude = [], None
    for item in payloads:
        entry_id, _, configuration = item.partition(":")
        built = get_entry(entry_id).build(
            configuration or None, orbit=orbit, time_step=time_step
        )
        instruments.extend(built.instruments)
        if attitude is None:
            attitude = built.attitude
    return Satellite(name=name, orbit=orbit, instruments=instruments, attitude=attitude)


def build_constellation(
    layout: Values | ChoiceValues,
    payload: Payload,
    *,
    epoch: datetime,
    altitude: float | None = None,
    name: str = "satellite",
    time_step: float | None = None,
) -> list[Satellite]:
    """
    Builds the satellites of a layout: the values of `walker_space` or
    `sso_planes_space` in a design (or of a choice of them), at the
    layout's altitude (or the given one), each with a payload:

    * an `Instrument` (or a list of instruments);
    * the identifier of an entry of the instrument catalog (see
      `tatc.catalog.get_entry`), optionally with a configuration after a
      colon (such as `"atms:swath"`), whose instruments are built for the
      orbit with the attitude they were validated with (and, if the
      configuration tailors a field of view, the analysis `time_step`), or
      a list of identifiers, whose instruments the satellite carries
      together (with the attitude of the first);
    * a function of a name and an orbit that builds a satellite.

    Satellites of Walker constellations are named `"<name> <k>"`, and
    those of sun-synchronous planes `"<name> <p>-<s>"` (plane `p`,
    satellite `s`).

    Args:
        layout (Values | ChoiceValues): The values of the layout.
        payload (Payload): The payload of each satellite.
        epoch (datetime.datetime): The epoch of the orbits.
        altitude (float | None): The altitude (meters), if the layout does
            not give it.
        name (str): The prefix of the names of the satellites.
        time_step (float | None): The analysis time step (seconds), for
            catalog configurations that tailor a field of view to it.

    Returns:
        list[Satellite]: the satellites
    """
    layout = _layout(layout)
    kind = _kind(layout)
    altitude = getattr(layout, "altitude", None) if altitude is None else altitude
    if altitude is None:
        raise ValueError("the layout has no altitude: give one")
    if kind == WALKER:
        members = WalkerConstellation(
            name=name,
            configuration=layout.configuration,
            orbit=CircularOrbit(
                altitude=altitude, inclination=layout.inclination, epoch=epoch
            ),
            number_planes=layout.planes,
            number_satellites=layout.planes * layout.per_plane,
            relative_spacing=layout.relative_spacing,
        ).generate_members()
        return [
            _make_satellite(payload, f"{name} {k + 1}", member.orbit, time_step)
            for k, member in enumerate(members)
        ]
    satellites = []
    for p, plane in enumerate(layout.planes):
        for s in range(plane.satellites):
            orbit = SunSynchronousOrbit(
                altitude=altitude,
                equator_crossing_time=hours_to_time(plane.crossing),
                equator_crossing_ascending=layout.ascending,
                epoch=epoch,
                true_anomaly=(plane.phase + 360 * s / plane.satellites) % 360,
            )
            satellites.append(
                _make_satellite(payload, f"{name} {p + 1}-{s + 1}", orbit, time_step)
            )
    return satellites


def remove_failures(
    satellites: Sequence[Satellite], probability: float, seed: int | None
) -> list[Satellite]:
    """
    Removes the satellites that fail, each with a probability, for the
    robustness of a design to failures: drawn for each replicate from an
    uncertain seed (see `Uncertain.seed`). Whether a satellite fails depends
    only on the seed and its position in the list, so that with common
    random numbers, the satellite at the same position fails in the same
    replicate of every design.

    Args:
        satellites (Sequence[Satellite]): The satellites, in the order they
            are built.
        probability (float): The probability that a satellite fails.
        seed (int | None): The seed of the replicate, or `None` (as before a
            design is realized) for no failures.

    Returns:
        list[Satellite]: the satellites that do not fail
    """
    if not 0 <= probability <= 1:
        raise ValueError("probability must be between 0 and 1")
    if seed is None or probability == 0:
        return list(satellites)
    return [
        satellite
        for k, satellite in enumerate(satellites)
        if np.random.default_rng([seed, k]).random() >= probability
    ]
