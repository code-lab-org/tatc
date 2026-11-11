"""
Building blocks of trade studies of constellations: fragments of design
spaces for common constellation layouts (Walker constellations, planes of
sun-synchronous orbits, and planes of inclined orbits, whose satellites are
evenly spaced or follow each other in trains), a builder of the satellites of
a layout (or of several, such as the shells of a mixed constellation) with a
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

A mixed constellation combines layouts, such as an inclined Walker
constellation and planes of sun-synchronous orbits (each layout a child of
the design space, or the items of a `Repeat` of a `Choice` of layouts), whose
satellites `build_constellation` builds together.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from datetime import datetime, time, timedelta
from typing import Any

import numpy as np

from ..schemas import (
    CircularOrbit,
    Instrument,
    Satellite,
    SunSynchronousOrbit,
    TrainConstellation,
    WalkerConstellation,
)
from .nodes import Derived, Fixed, Integer, Real, Repeat, Uncertain
from .space import DesignSpace
from .values import ChoiceValues, RepeatValues, Scope, Values

WALKER = "walker"
"""The kind of the layouts of `walker_space`."""
SSO_PLANES = "sso planes"
"""The kind of the layouts of `sso_planes_space`."""
PLANES = "planes"
"""The kind of the layouts of `planes_space`."""


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
    phase_fraction: Any = 0.0,
    node_fraction: Any = 0.0,
    constraints: Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of Walker constellations (see
    `tatc.schemas.WalkerConstellation`): their inclination (degrees), number
    of planes, number of satellites in each plane, relative spacing of
    satellites between planes, the phase and node of the lead satellite, and
    (optionally, if not given to `build_constellation`) altitude (meters).
    Each is a node of a design space (such as a design variable) or a value.

    The phase (true anomaly) and node (right ascension of the ascending
    node) of the lead satellite are given as fractions of the intervals over
    which they give distinct constellations: shifting every satellite's
    phase by `360 / per_plane` degrees moves each onto the next in its plane,
    and shifting every plane's node by `360 / planes` degrees (of a delta
    configuration) moves each onto the next plane, which, with the shift of
    phase between planes (`360 * relative_spacing / satellites`), gives the
    same constellation. Fractions in [0, 1) of both therefore give every
    distinct constellation once, without searching the same constellation
    at several phases or nodes. A star configuration has no such symmetry of
    its nodes (its planes span 180 degrees, and the next plane after the
    last is the first, in the opposite direction), so its fraction of the
    node is of 360 degrees. The fractions may be design variables (such as
    `Real(0, 0.75, step=0.25)`, without 1, which is the same as 0) or
    uncertain parameters (such as `Uncertain.uniform(0, 1)`), and the
    fragment derives the lead satellite's `true_anomaly` and `right_ascension`
    (degrees) from them.

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
        phase_fraction (Any): The phase of the lead satellite, as a fraction
            of `360 / per_plane` degrees.
        node_fraction (Any): The node of the lead satellite, as a fraction of
            `360 / planes` degrees (delta configuration) or 360 degrees (star
            configuration).
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
        "phase_fraction": phase_fraction,
        "node_fraction": node_fraction,
    }
    if altitude is not None:
        children["altitude"] = altitude
    children["true_anomaly"] = Derived(_lead_true_anomaly)
    children["right_ascension"] = Derived(_lead_right_ascension)
    return DesignSpace(constraints=constraints, **children)


def _lead_true_anomaly(d: Scope) -> float:
    """The true anomaly (degrees) of the lead satellite of a Walker layout."""
    return float(d.phase_fraction) * 360 / d.per_plane % 360


def _lead_right_ascension(d: Scope) -> float:
    """The right ascension of the ascending node (degrees) of the lead
    satellite of a Walker layout."""
    period = 360 if d.configuration == "star" else 360 / d.planes
    return float(d.node_fraction) * period % 360


def _plane_true_anomaly(d: Scope) -> float:
    """The true anomaly (degrees) of the first satellite of a plane."""
    period = 360 if d.interval is not None else 360 / d.satellites
    return float(d.phase_fraction) * period % 360


def _distinct_crossings(d: Scope) -> bool:
    """Whether the planes of a layout have different crossing times."""
    crossings = [round(plane.crossing, 6) for plane in d.planes]
    return len(set(crossings)) == len(crossings)


def sso_planes_space(
    *,
    crossing: Any = None,
    satellites: Any = None,
    count: Any = None,
    phase_fraction: Any = 0.0,
    interval: Any = None,
    altitude: Any = None,
    ascending: Any = True,
    repeat_ground_track: Any = True,
    ordered: bool = False,
    distinct: bool = True,
    constraints: Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of planes of sun-synchronous orbits
    (see `tatc.schemas.SunSynchronousOrbit`), each with its own equator
    crossing time (local solar time, in hours), number of satellites,
    spacing of its satellites, and phase (see below), and (optionally, if not
    given to
    `build_constellation`) their altitude (meters). Each is a node of a
    design space (such as a design variable, a value derived from the index
    of the plane, or an uncertain parameter) or a value.

    The satellites of a plane are evenly spaced around its orbit if its
    `interval` is `None`, or follow each other in a train (see
    `tatc.schemas.TrainConstellation`) at an `interval` (minutes): each
    crosses the equator the interval after the one ahead of it, on the same
    ground track (if `repeat_ground_track`, with its node shifted by the
    Earth's rotation, which delays its local crossing time by the interval)
    or in the same plane. Evenly spaced satellites revisit at intervals of
    the orbital period divided by their number; a train revisits within
    minutes, then not until the next orbit. The interval may be a design
    variable whose options include `None`, such as `Categorical([None, 5,
    15])`.

    The phase of each plane is a fraction of the interval over which it
    gives distinct planes (as that of a Walker constellation, see
    `walker_space`): `360 / satellites` degrees of true anomaly if its
    satellites are evenly spaced (a shift by which moves each onto the next),
    or 360 degrees if they are a train. Fractions in [0, 1) therefore give
    every distinct plane once, may be design variables or uncertain
    parameters (such as `Uncertain.uniform(0, 1)`) whatever the spacing, and
    the fragment derives the true anomaly of the plane's first satellite
    (`true_anomaly`, degrees) from them.

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
        phase_fraction (Any): The phase of each plane, as a fraction of
            `360 / satellites` degrees (evenly spaced) or 360 degrees (a
            train): an uncertain parameter, such as `Uncertain.uniform(0, 1)`,
            if it is not designed.
        interval (Any): The interval (minutes) between the satellites of a
            train in each plane, or `None` for evenly spaced satellites.
        altitude (Any): The altitude, if part of the layout.
        ascending (Any): Whether the crossing times are of ascending nodes.
        repeat_ground_track (Any): Whether the satellites of a train follow
            the same ground track (or the same plane).
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
        interval=Fixed(None) if interval is None else interval,
        phase_fraction=phase_fraction,
        true_anomaly=Derived(_plane_true_anomaly),
    )
    children = {
        "kind": Fixed(SSO_PLANES),
        "ascending": ascending,
        "repeat_ground_track": repeat_ground_track,
        "planes": Repeat(plane, count=_node(count, Integer(1, 3)), ordered=ordered),
    }
    if altitude is not None:
        children["altitude"] = altitude
    constraints = ([_distinct_crossings] if distinct else []) + list(constraints)
    return DesignSpace(constraints=constraints, **children)


def planes_space(
    *,
    inclination: Any = None,
    node: Any = None,
    satellites: Any = None,
    count: Any = None,
    phase_fraction: Any = 0.0,
    interval: Any = None,
    altitude: Any = None,
    repeat_ground_track: Any = True,
    ordered: bool = False,
    constraints: Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of planes of circular orbits (see
    `tatc.schemas.CircularOrbit`), each with its own inclination (degrees),
    node (right ascension of the ascending node at the epoch, degrees),
    number of satellites, spacing of its satellites (evenly spaced, or a
    train at an `interval`), and phase (a fraction, see `sso_planes_space`),
    and (optionally, if not given to
    `build_constellation`) their altitude (meters). Each is a node of a
    design space (such as a design variable, a value derived from the index
    of the plane, or an uncertain parameter) or a value.

    Unlike a Walker constellation, whose planes are evenly spaced in node,
    the planes may have any nodes and inclinations, as those of satellites
    launched as secondary payloads, whose nodes are chosen by the launches
    of the primary payloads: the node of each plane is therefore uncertain
    by default (uniform from 0 to 360 degrees), drawn for each replicate.
    The planes are unordered by default (see `sso_planes_space`).

    Args:
        inclination (Any): The inclination of each plane, by default a
            variable from 30 to 90 degrees in steps of 15.
        node (Any): The node of each plane, by default uncertain (uniform
            from 0 to 360 degrees).
        satellites (Any): The number of satellites of each plane, by default a
            variable from 1 to 4.
        count (Any): The number of planes, by default a variable from 1 to 3.
        phase_fraction (Any): The phase of each plane, as a fraction of
            `360 / satellites` degrees (evenly spaced) or 360 degrees (a
            train).
        interval (Any): The interval (minutes) between the satellites of a
            train in each plane, or `None` for evenly spaced satellites.
        altitude (Any): The altitude, if part of the layout.
        repeat_ground_track (Any): Whether the satellites of a train follow
            the same ground track (or the same plane).
        ordered (bool): `True` if the planes are ordered.
        constraints (Sequence[Callable]): Constraints of the fragment.

    Returns:
        DesignSpace: the fragment
    """
    plane = DesignSpace(
        inclination=_node(inclination, Real(30, 90, step=15)),
        node=_node(node, Uncertain.uniform(0, 360)),
        satellites=_node(satellites, Integer(1, 4)),
        interval=Fixed(None) if interval is None else interval,
        phase_fraction=phase_fraction,
        true_anomaly=Derived(_plane_true_anomaly),
    )
    children = {
        "kind": Fixed(PLANES),
        "repeat_ground_track": repeat_ground_track,
        "planes": Repeat(plane, count=_node(count, Integer(1, 3)), ordered=ordered),
    }
    if altitude is not None:
        children["altitude"] = altitude
    return DesignSpace(constraints=constraints, **children)


def _layouts(layout: Any) -> list[Any] | None:
    """
    Gets the layouts of several (a list or a repetition of layouts), or
    `None` for a single layout.
    """
    if isinstance(layout, (list, tuple, RepeatValues)):
        return list(layout)
    return None


def _layout(layout: Any) -> Any:
    """Gets the values of a layout (the chosen branch of a choice)."""
    return layout.value if isinstance(layout, ChoiceValues) else layout


def _kind(layout: Any) -> str:
    try:
        return layout.kind
    except AttributeError:
        raise ValueError(
            "a layout must be the values of walker_space, sso_planes_space, "
            + "or planes_space"
        ) from None


def count_satellites(layout: Any) -> int:
    """
    Counts the satellites of a layout (the values of `walker_space`,
    `sso_planes_space`, or `planes_space`, or of a choice of them), or of
    several (a list or a repetition of layouts), for objectives and
    constraints.

    Args:
        layout (Any): The values of the layout (or layouts).

    Returns:
        int: the number of satellites
    """
    several = _layouts(layout)
    if several is not None:
        return sum(count_satellites(item) for item in several)
    layout = _layout(layout)
    if _kind(layout) == WALKER:
        return int(layout.planes * layout.per_plane)
    return int(sum(plane.satellites for plane in layout.planes))


def count_planes(layout: Any) -> int:
    """
    Counts the orbit planes of a layout (the values of `walker_space`,
    `sso_planes_space`, or `planes_space`, or of a choice of them), or of
    several (a list or a repetition of layouts), each of which needs its own
    launch.

    Args:
        layout (Any): The values of the layout (or layouts).

    Returns:
        int: the number of planes
    """
    several = _layouts(layout)
    if several is not None:
        return sum(count_planes(item) for item in several)
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


def _plane_orbits(
    make: Callable[[float], Any],
    plane: Any,
    repeat_ground_track: bool,
) -> list[Any]:
    """
    Gets the orbits of the satellites of a plane: evenly spaced from the
    true anomaly of its first satellite, or a train at its interval (see
    `sso_planes_space`), with `make` the orbit of a satellite at a true
    anomaly.
    """
    interval = getattr(plane, "interval", None)
    count = int(plane.satellites)
    if interval is None:
        return [
            make((plane.true_anomaly + 360 * s / count) % 360) for s in range(count)
        ]
    interval = (
        interval
        if isinstance(interval, timedelta)
        else timedelta(minutes=float(interval))
    )
    if interval <= timedelta(0):
        raise ValueError("the interval of a train must be positive")
    train = TrainConstellation(
        name="train",
        orbit=make(plane.true_anomaly % 360),
        number_satellites=count,
        interval=interval,
        repeat_ground_track=bool(repeat_ground_track),
    )
    return [member.orbit for member in train.generate_members()]


def build_constellation(
    layout: Values | ChoiceValues | Sequence[Values | ChoiceValues] | RepeatValues,
    payload: Payload,
    *,
    epoch: datetime,
    altitude: float | None = None,
    name: str = "satellite",
    time_step: float | None = None,
) -> list[Satellite]:
    """
    Builds the satellites of a layout: the values of `walker_space`,
    `sso_planes_space`, or `planes_space` in a design (or of a choice of
    them), at the layout's altitude (or the given one), each with a payload:

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
    those of planes `"<name> <p>-<s>"` (plane `p`, satellite `s`).

    Several layouts (a list or a repetition of layouts, such as the shells
    of a mixed constellation) are built in turn, with the same payload, and
    the names of their satellites are prefixed by a letter for each layout
    (`"<name> A 1"`, `"<name> B 1-1"`). Build layouts with different
    payloads in separate calls, with different names.

    Args:
        layout (Values | ChoiceValues | Sequence): The values of the layout
            (or layouts).
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
    several = _layouts(layout)
    if several is not None:
        return [
            satellite
            for k, item in enumerate(several)
            for satellite in build_constellation(
                item,
                payload,
                epoch=epoch,
                altitude=altitude,
                name=f"{name} {chr(65 + k) if k < 26 else k + 1}",
                time_step=time_step,
            )
        ]
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
                altitude=altitude,
                inclination=layout.inclination,
                true_anomaly=getattr(layout, "true_anomaly", 0.0),
                right_ascension_ascending_node=getattr(layout, "right_ascension", 0.0),
                epoch=epoch,
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
        if kind == SSO_PLANES:

            def make(anomaly, plane=plane):
                return SunSynchronousOrbit(
                    altitude=altitude,
                    equator_crossing_time=hours_to_time(plane.crossing),
                    equator_crossing_ascending=layout.ascending,
                    epoch=epoch,
                    true_anomaly=anomaly,
                )

        else:

            def make(anomaly, plane=plane):
                return CircularOrbit(
                    altitude=altitude,
                    inclination=plane.inclination,
                    right_ascension_ascending_node=float(plane.node) % 360,
                    epoch=epoch,
                    true_anomaly=anomaly,
                )

        orbits = _plane_orbits(
            make, plane, getattr(layout, "repeat_ground_track", True)
        )
        satellites.extend(
            _make_satellite(payload, f"{name} {p + 1}-{s + 1}", orbit, time_step)
            for s, orbit in enumerate(orbits)
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
