"""
Building blocks of trade studies of constellations: fragments of design
spaces for common constellation layouts (Walker constellations, planes of
sun-synchronous orbits, and planes of inclined orbits, whose satellites are
evenly spaced or follow each other in trains, streets-of-coverage
constellations, and mutual orbiting groups), a builder of the satellites of
a layout (or of several, such as the shells of a mixed constellation) with a
payload (instruments, or entries of the instrument catalog), and counts of
satellites and planes for objectives and constraints.

The fragments share the names of their values: `planes` is the planes,
their number (of Walker and streets-of-coverage constellations) or a
repetition of them (of planes of sun-synchronous or inclined orbits), whose
`count` is their number; `per_plane` is the number of satellites in each
plane, and `satellites` the number of satellites of a layout without planes
(a mutual orbiting group); `inclination`, `node` (the right ascension of the
ascending node), and `true_anomaly` are in degrees, and `phase_fraction` and
`node_fraction` are fractions of the intervals over which the phase and
node give distinct layouts.

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

The values of a fragment in a design are those of a `Layout` (see
`layout_of`), which gives its planes and the orbits of their satellites, a
short description (such as `"Walker 60 deg: 3 x 4 (spacing 1)"`, see
`Layout.describe`), and which other layouts extend.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from datetime import datetime, time, timedelta
from functools import partial
from typing import TYPE_CHECKING, Any

import numpy as np

from ..schemas import (
    CircularOrbit,
    Instrument,
    MOGConstellation,
    Satellite,
    SOCConstellation,
    SunSynchronousOrbit,
    TrainConstellation,
    WalkerConstellation,
)
from .nodes import Derived, Fixed, Integer, Real, Repeat, Uncertain
from .space import DesignSpace, _as_constraints
from .values import ChoiceValues, RepeatValues, Scope, Values

if TYPE_CHECKING:
    from ..catalog import Grade

WALKER = "walker"
"""The pattern of the layouts of `walker_space`."""
SSO_PLANES = "sso planes"
"""The pattern of the layouts of `sso_planes_space`."""
PLANES = "planes"
"""The pattern of the layouts of `planes_space`."""
SOC = "soc"
"""The pattern of the layouts of `soc_space`."""
MOG = "mog"
"""The pattern of the layouts of `mog_space`."""


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
    constraints: Callable[[Scope], Any] | Sequence[Callable[[Scope], Any]] = (),
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
    phase between planes (`360 * relative_spacing / (planes * per_plane)`),
    gives the same constellation. Fractions in [0, 1) of both therefore give every
    distinct constellation once, without searching the same constellation
    at several phases or nodes. A star configuration has no such symmetry of
    its nodes (its planes span 180 degrees, and the next plane after the
    last is the first, in the opposite direction), so its fraction of the
    node is of 360 degrees. The fractions may be design variables (such as
    `Real(0, 0.75, step=0.25)`, without 1, which is the same as 0) or
    uncertain parameters (such as `Uncertain.uniform(0, 1)`), and the
    fragment derives the lead satellite's `true_anomaly` and `node` (degrees)
    from them.

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
        constraints (Callable | Sequence[Callable]): A constraint, or the
            constraints, of the fragment.

    Returns:
        DesignSpace: the fragment
    """
    children = {
        "pattern": Fixed(WALKER),
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
    children["node"] = Derived(_lead_node)
    return DesignSpace(constraints=constraints, **children)


def _lead_true_anomaly(d: Scope) -> float:
    """The true anomaly (degrees) of the lead satellite of a Walker layout."""
    return float(d.phase_fraction) * 360 / d.per_plane % 360


def _lead_node(d: Scope) -> float:
    """The node (right ascension of the ascending node, degrees) of the lead
    satellite of a Walker layout."""
    period = 360 if d.configuration == "star" else 360 / d.planes
    return float(d.node_fraction) * period % 360


def _plane_true_anomaly(d: Scope) -> float:
    """The true anomaly (degrees) of the first satellite of a plane."""
    period = 360 if d.interval is not None else 360 / d.per_plane
    return float(d.phase_fraction) * period % 360


def _distinct_crossings(d: Scope) -> bool:
    """Whether the planes of a layout have different crossing times."""
    crossings = [round(plane.crossing, 6) for plane in d.planes]
    return len(set(crossings)) == len(crossings)


def sso_planes_space(
    *,
    crossing: Any = None,
    planes: Any = None,
    per_plane: Any = None,
    phase_fraction: Any = 0.0,
    interval: Any = None,
    altitude: Any = None,
    ascending: Any = True,
    repeat_ground_track: Any = True,
    ordered: bool = False,
    distinct: bool = True,
    constraints: Callable[[Scope], Any] | Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of planes of sun-synchronous orbits
    (see `tatc.schemas.SunSynchronousOrbit`), each with its own equator
    crossing time (local solar time, in hours), number of satellites,
    spacing of its satellites, and phase (see below), and (optionally, if not
    given to `build_constellation`) their altitude (meters). Each is a node
    of a design space (such as a design variable, a value derived from the
    index of the plane, or an uncertain parameter) or a value.

    The planes are a `Repeat` named `planes`: the number of planes is its
    `count` (in a design, `planes.count`), and the values of each plane
    (`crossing`, `per_plane`, `interval`, `phase_fraction`, and the derived
    `true_anomaly`) are its items (`planes[0].crossing`).

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
    `walker_space`): `360 / per_plane` degrees of true anomaly if its
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
        planes (Any): The number of planes, by default a variable from 1 to 3.
        per_plane (Any): The number of satellites in each plane, by default
            a variable from 1 to 4.
        phase_fraction (Any): The phase of each plane, as a fraction of
            `360 / per_plane` degrees (evenly spaced) or 360 degrees (a
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
        constraints (Callable | Sequence[Callable]): Another constraint,
            or other constraints, of the fragment.

    Returns:
        DesignSpace: the fragment
    """
    plane = DesignSpace(
        crossing=_node(crossing, Real(6, 18, step=1.5)),
        per_plane=_node(per_plane, Integer(1, 4)),
        interval=interval,
        phase_fraction=phase_fraction,
        true_anomaly=Derived(_plane_true_anomaly),
    )
    children = {
        "pattern": Fixed(SSO_PLANES),
        "ascending": ascending,
        "repeat_ground_track": repeat_ground_track,
        "planes": Repeat(plane, count=_node(planes, Integer(1, 3)), ordered=ordered),
    }
    if altitude is not None:
        children["altitude"] = altitude
    constraints = [
        *([_distinct_crossings] if distinct else []),
        *_as_constraints(constraints),
    ]
    return DesignSpace(constraints=constraints, **children)


def planes_space(
    *,
    inclination: Any = None,
    node: Any = None,
    planes: Any = None,
    per_plane: Any = None,
    phase_fraction: Any = 0.0,
    interval: Any = None,
    altitude: Any = None,
    repeat_ground_track: Any = True,
    ordered: bool = False,
    constraints: Callable[[Scope], Any] | Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of planes of circular orbits (see
    `tatc.schemas.CircularOrbit`), each with its own inclination (degrees),
    node (right ascension of the ascending node at the epoch, degrees),
    number of satellites, spacing of its satellites (evenly spaced, or a
    train at an `interval`), and phase (a fraction, see `sso_planes_space`),
    and (optionally, if not given to `build_constellation`) their altitude
    (meters). Each is a node of a design space (such as a design variable, a
    value derived from the index of the plane, or an uncertain parameter) or
    a value. The planes are a `Repeat` named `planes`, as those of
    `sso_planes_space`.

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
        planes (Any): The number of planes, by default a variable from 1 to 3.
        per_plane (Any): The number of satellites in each plane, by default
            a variable from 1 to 4.
        phase_fraction (Any): The phase of each plane, as a fraction of
            `360 / per_plane` degrees (evenly spaced) or 360 degrees (a
            train).
        interval (Any): The interval (minutes) between the satellites of a
            train in each plane, or `None` for evenly spaced satellites.
        altitude (Any): The altitude, if part of the layout.
        repeat_ground_track (Any): Whether the satellites of a train follow
            the same ground track (or the same plane).
        ordered (bool): `True` if the planes are ordered.
        constraints (Callable | Sequence[Callable]): A constraint, or the
            constraints, of the fragment.

    Returns:
        DesignSpace: the fragment
    """
    plane = DesignSpace(
        inclination=_node(inclination, Real(30, 90, step=15)),
        node=_node(node, Uncertain.uniform(0, 360)),
        per_plane=_node(per_plane, Integer(1, 4)),
        interval=interval,
        phase_fraction=phase_fraction,
        true_anomaly=Derived(_plane_true_anomaly),
    )
    children = {
        "pattern": Fixed(PLANES),
        "repeat_ground_track": repeat_ground_track,
        "planes": Repeat(plane, count=_node(planes, Integer(1, 3)), ordered=ordered),
    }
    if altitude is not None:
        children["altitude"] = altitude
    return DesignSpace(constraints=constraints, **children)


def _soc_of(values: Any, altitude: float, epoch: datetime | None = None):
    """The streets-of-coverage constellation of the values of a layout (or of
    a scope), at an altitude and (if given) an epoch."""
    orbit = {
        "altitude": altitude,
        "inclination": values.inclination,
        "true_anomaly": getattr(values, "true_anomaly", 0.0),
        "right_ascension_ascending_node": getattr(values, "node", 0.0),
    }
    return SOCConstellation(
        name="soc",
        orbit=CircularOrbit(**orbit, **({} if epoch is None else {"epoch": epoch})),
        swath_width=values.swath_width,
        packing_distance=values.packing_distance,
        polar=values.polar,
    )


def _soc_planes(d: Scope) -> int:
    """The number of planes of a streets-of-coverage layout."""
    return _soc_of(d, d.altitude).get_number_planes()


def _soc_per_plane(d: Scope) -> int:
    """The number of satellites in each plane of a streets-of-coverage
    layout."""
    return _soc_of(d, d.altitude).get_satellites_per_plane()


def _scope_altitude(d: Scope) -> float:
    """The altitude decided before a fragment (outside it)."""
    return d.altitude


def _soc_node(d: Scope) -> float:
    """The node (right ascension of the ascending node, degrees) of the lead
    satellite of a streets-of-coverage layout."""
    period = 360 if _soc_of(d, d.altitude).is_polar() else 360 / d.planes
    return float(d.node_fraction) * period % 360


def soc_space(
    *,
    swath_width: Any,
    inclination: Any = None,
    packing_distance: Any = 1.0,
    polar: Any = None,
    altitude: Any = None,
    phase_fraction: Any = 0.0,
    node_fraction: Any = 0.0,
    constraints: Callable[[Scope], Any] | Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of streets-of-coverage constellations
    (see `tatc.schemas.SOCConstellation`): the fewest satellites whose
    footprints (of a swath width, in meters) cover the Earth continuously at
    an inclination (degrees), in a Walker delta pattern for inclined orbits
    or the polar pattern for near-polar ones (`polar`, by default chosen by
    the inclination), with a margin of overlap (`packing_distance`, 1 for
    none). Each is a node of a design space (such as a design variable, or a
    value derived from the instrument, as the swath width of the sensor of a
    design) or a value.

    The fragment derives its number of planes (`planes`) and of satellites
    in each plane (`per_plane`) from them and the altitude (meters), which
    therefore must be decided before them: a value of the fragment
    (`altitude`), or by default the altitude of the design space outside it
    (declared before it), which the fragment keeps as its own (so that
    `build_constellation` builds the satellites at the altitude its planes
    were derived at).

    The phase and node of the lead satellite are fractions, as those of a
    Walker constellation (see `walker_space`): of `360 / per_plane` degrees
    of phase, and of `360 / planes` degrees of node (inclined orbits) or 360
    degrees (polar orbits, whose planes span 180 degrees). The fragment
    derives the lead satellite's `true_anomaly` and `node` (degrees) from
    them.

    Args:
        swath_width (Any): The swath width (meters) of the footprints.
        inclination (Any): The inclination, by default a variable from 30 to
            90 degrees in steps of 15.
        packing_distance (Any): The relative distance between footprints
            (inclined orbits) or the relative radius of the footprints (polar
            orbits), in (0, 1]: less than 1 for a margin of overlap.
        polar (Any): `True` for the polar pattern, `False` for the Walker
            delta pattern, or `None` for the polar pattern of near-polar
            orbits.
        altitude (Any): The altitude, or `None` for that of the design space
            outside the fragment.
        phase_fraction (Any): The phase of the lead satellite, as a fraction
            of `360 / per_plane` degrees.
        node_fraction (Any): The node of the lead satellite, as a fraction of
            `360 / planes` degrees (inclined orbits) or 360 degrees (polar
            orbits).
        constraints (Callable | Sequence[Callable]): A constraint, or the
            constraints, of the fragment (such as a maximum of
            `planes * per_plane` satellites).

    Returns:
        DesignSpace: the fragment
    """
    children = {
        "pattern": Fixed(SOC),
        "inclination": _node(inclination, Real(30, 90, step=15)),
        "swath_width": swath_width,
        "packing_distance": packing_distance,
        "polar": polar,
        "phase_fraction": phase_fraction,
        "node_fraction": node_fraction,
        "altitude": Derived(_scope_altitude) if altitude is None else altitude,
    }
    children["planes"] = Derived(_soc_planes)
    children["per_plane"] = Derived(_soc_per_plane)
    children["true_anomaly"] = Derived(_lead_true_anomaly)
    children["node"] = Derived(_soc_node)
    return DesignSpace(constraints=constraints, **children)


def _group_true_anomaly(d: Scope) -> float:
    """The true anomaly (degrees) of the reference orbit of a mutual orbiting
    group."""
    return float(d.phase_fraction) * 360 % 360


def mog_space(
    *,
    parallel_axis: Any,
    transverse_axis: Any,
    inclination: Any = None,
    node: Any = 0.0,
    satellites: Any = None,
    clockwise: Any = True,
    altitude: Any = None,
    phase_fraction: Any = 0.0,
    constraints: Callable[[Scope], Any] | Sequence[Callable[[Scope], Any]] = (),
) -> DesignSpace:
    """
    Gets a fragment of a design space of mutual orbiting groups (see
    `tatc.schemas.MOGConstellation`): satellites that orbit each other (on
    an ellipse with axes, in meters, parallel and transverse to the
    velocity, clockwise or not) about a reference circular orbit of an
    inclination and node (degrees), and (optionally, if not given to
    `build_constellation`) altitude (meters). Each is a node of a design
    space (such as a design variable) or a value.

    The satellites of a group are launched together, as a single plane (of
    `count_planes`), although the orbit of each has its own inclination and
    node. The phase of the reference orbit is a fraction of 360 degrees,
    from which the fragment derives its `true_anomaly` (degrees).

    Args:
        parallel_axis (Any): The axis (meters) of the mutual orbit parallel
            to the velocity.
        transverse_axis (Any): The axis (meters) of the mutual orbit
            transverse to the velocity.
        inclination (Any): The inclination of the reference orbit, by default
            a variable from 30 to 90 degrees in steps of 15.
        node (Any): The node of the reference orbit.
        satellites (Any): The number of satellites, by default a variable
            from 2 to 6.
        clockwise (Any): Whether the mutual orbit is clockwise.
        altitude (Any): The altitude, if part of the layout.
        phase_fraction (Any): The phase of the reference orbit, as a fraction
            of 360 degrees.
        constraints (Callable | Sequence[Callable]): A constraint, or the
            constraints, of the fragment.

    Returns:
        DesignSpace: the fragment
    """
    children = {
        "pattern": Fixed(MOG),
        "inclination": _node(inclination, Real(30, 90, step=15)),
        "node": node,
        "satellites": _node(satellites, Integer(2, 6)),
        "parallel_axis": parallel_axis,
        "transverse_axis": transverse_axis,
        "clockwise": clockwise,
        "phase_fraction": phase_fraction,
    }
    if altitude is not None:
        children["altitude"] = altitude
    children["true_anomaly"] = Derived(_group_true_anomaly)
    return DesignSpace(constraints=constraints, **children)


@dataclass(frozen=True)
class Plane:
    """
    An orbit plane of a layout (see `Layout`): its number of satellites, and
    its values in a design.
    """

    satellites: int
    """int: The number of satellites."""
    values: Any = None
    """Any: The values of the plane in a design (such as those of each plane
    of `sso_planes_space`), or `None` (as for the planes of a Walker
    constellation), as given to the `order` and `launch_mass_capacity` of a
    `DeploymentSchedule` (functions of the index and values of a plane)."""


class Layout:
    """
    The satellites of a constellation, plane by plane: a Walker
    constellation (`WalkerLayout`), planes of sun-synchronous orbits
    (`SSOPlanesLayout`) or of inclined orbits (`PlanesLayout`), a
    streets-of-coverage constellation (`SOCLayout`), a mutual orbiting group
    (`MOGLayout`), or several layouts flown together (`Layouts`). `layout_of`
    gets the layout of the values of a design (of a fragment such as
    `walker_space`). The functions of layouts (`build_constellation`,
    `count_satellites`, `count_planes`, `plane_sizes`, and
    `DeploymentSchedule.plan`) take either.

    The satellites of a layout are those of its planes in order, and those of
    each plane in their order in it, so that the satellites that `build`
    builds are those of the planes that a deployment launches (see
    `DeploymentSchedule.plan` and `MissionHistory.phases`).

    Other layouts subclass it, implementing `planes` and `orbits` (and
    `names`, for names other than `"<name> <p>-<s>"`, and `describe`), so
    that they are built, counted, deployed, and costed as the others.

    A layout of its own altitude (such as a fragment with an `altitude`)
    is built at that altitude; the altitude given to `build_constellation`
    (or `orbits`) is that of layouts without their own.
    """

    def planes(self) -> list[Plane]:
        """
        Gets the planes, in order.

        Returns:
            list[Plane]: the planes
        """
        raise NotImplementedError

    def orbits(self, epoch: datetime, altitude: float | None = None) -> list[list[Any]]:
        """
        Gets the orbits of the satellites of each plane, in the order of the
        planes.

        Args:
            epoch (datetime.datetime): The epoch of the orbits.
            altitude (float | None): The altitude (meters) of a layout without
                its own, if any.

        Returns:
            list[list]: the orbits of the satellites of each plane

        Raises:
            ValueError: if the layout has no altitude of its own and none is
                given.
        """
        raise NotImplementedError

    def _unrealized(self) -> bool:
        """Whether the layout has values that are not realized (`NaN`)."""
        return _has_nan(getattr(self, "values", None))

    def describe(self) -> str:
        """
        Gets a short description of the layout (such as `"Walker 60 deg: 3 x
        4 (spacing 1)"` or `"SSO: 2 at 06:00, 1 at 13:30"`), as for labels of
        designs: by default, its numbers of planes and satellites.

        Returns:
            str: the description
        """
        planes = self.planes()
        return f"{len(planes)} planes, {self.satellites} satellites"

    def __str__(self) -> str:
        return self.describe()

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.describe()!r})"

    def names(self, name: str) -> list[list[str]]:
        """
        Gets the names of the satellites of each plane: `"<name> <p>-<s>"`
        (plane `p`, satellite `s`).

        Args:
            name (str): The prefix of the names.

        Returns:
            list[list[str]]: the names of the satellites of each plane
        """
        return [
            [f"{name} {p + 1}-{s + 1}" for s in range(plane.satellites)]
            for p, plane in enumerate(self.planes())
        ]

    @property
    def satellites(self) -> int:
        """int: The number of satellites."""
        return sum(plane.satellites for plane in self.planes())

    def build(
        self,
        payload: Payload,
        *,
        epoch: datetime,
        altitude: float | None = None,
        name: str = "satellite",
        time_step: float | None = None,
        start: datetime | None = None,
        end: datetime | None = None,
        grade: Grade | str | None = None,
        platform: str | None = None,
    ) -> list[Satellite]:
        """
        Builds the satellites, in the order of the planes, each with a
        payload (see `build_constellation`).

        Args:
            payload (Payload): The payload of each satellite.
            epoch (datetime.datetime): The epoch of the orbits.
            altitude (float | None): The altitude (meters), if the layout does
                not give its own.
            name (str): The prefix of the names of the satellites.
            time_step (float | None): The analysis time step (seconds), for
                catalog configurations that tailor a field of view to it.
            start (datetime.datetime | None): The start of the analysis, for
                catalog configurations that schedule their tilt by day.
            end (datetime.datetime | None): The end of the analysis, for
                catalog configurations that schedule their tilt by day.
            grade (Grade | str | None): The grade of the use, if any, which
                every catalog configuration must serve.
            platform (str | None): The identifier of the platform modeled, if
                any, which selects the catalog configurations tested on it
                (with a warning for a configuration not tested on it).

        Returns:
            list[Satellite]: the satellites

        Raises:
            ValueError: if the orbits or names of the satellites do not match
                the planes, or the layout has values that are not realized.
        """
        if self._unrealized():
            raise ValueError(
                "the layout has values that are not realized (NaN, such as of "
                + "uncertain parameters): build the satellites of a realized "
                + "design (see Design.realize and Study.realize)"
            )
        counts = [plane.satellites for plane in self.planes()]
        orbits, names = self.orbits(epoch, altitude), self.names(name)
        if [len(o) for o in orbits] != counts or [len(n) for n in names] != counts:
            raise ValueError("the orbits and names of a layout must match its planes")
        return [
            _make_satellite(
                payload,
                satellite,
                orbit,
                time_step=time_step,
                start=start,
                end=end,
                grade=grade,
                platform=platform,
            )
            for plane_names, plane_orbits in zip(names, orbits)
            for satellite, orbit in zip(plane_names, plane_orbits)
        ]


def _get_altitude(values: Any, altitude: float | None) -> float:
    """Gets the altitude of the values of a layout, or else the given one."""
    own = getattr(values, "altitude", None)
    altitude = altitude if own is None else own
    if altitude is None:
        raise ValueError("the layout has no altitude: give one")
    return altitude


def _has_nan(value: Any) -> bool:
    """Whether values (of a layout) have a missing number (`NaN`)."""
    if isinstance(value, float):
        return value != value  # pylint: disable=comparison-with-itself
    if isinstance(value, ChoiceValues):
        return _has_nan(value.value)
    if isinstance(value, RepeatValues):
        return any(_has_nan(item) for item in value)
    if isinstance(value, Values):
        return any(_has_nan(value[name]) for name in value)
    return False


def _km(altitude: Any) -> str:
    """The description of an altitude (meters) of a layout, if it has one."""
    return "" if altitude is None else f" at {float(altitude) / 1e3:g} km"


class _ValuesLayout(Layout):
    """
    A layout of the values of a fragment of a design space in a design.

    Args:
        values (Values): The values of the layout.
    """

    def __init__(self, values: Any):
        self.values = values

    def _altitude(self) -> str:
        return _km(getattr(self.values, "altitude", None))


class _MembersLayout(_ValuesLayout):
    """
    A constellation whose members a schema of TAT-C generates (see
    `_constellation`), plane by plane (the members of each plane following
    each other), its satellites named `"<name> <k>"`, in order.
    """

    _MISMATCH = "the members of the constellation do not match its planes"

    def _constellation(self, epoch: datetime, altitude: float) -> Any:
        """Gets the constellation (of `tatc.schemas`) at an epoch and altitude."""
        raise NotImplementedError

    def orbits(self, epoch: datetime, altitude: float | None = None) -> list[list[Any]]:
        altitude = _get_altitude(self.values, altitude)
        members = self._constellation(epoch, altitude).generate_members()
        sizes = [plane.satellites for plane in self.planes()]
        if len(members) != sum(sizes):
            raise ValueError(self._MISMATCH)
        starts = np.cumsum([0, *sizes])
        return [
            [member.orbit for member in members[start:stop]]
            for start, stop in zip(starts[:-1], starts[1:])
        ]

    def names(self, name: str) -> list[list[str]]:
        sizes = [plane.satellites for plane in self.planes()]
        starts = np.cumsum([0, *sizes])
        return [
            [f"{name} {k + 1}" for k in range(start, start + size)]
            for start, size in zip(starts[:-1], sizes)
        ]


class WalkerLayout(_MembersLayout):
    """
    A Walker constellation (see `tatc.schemas.WalkerConstellation`): the
    values of `walker_space` in a design. Its satellites are named
    `"<name> <k>"`, in order.

    Args:
        values (Values): The values of the layout.
    """

    def planes(self) -> list[Plane]:
        return [Plane(int(self.values.per_plane))] * int(self.values.planes)

    def _constellation(self, epoch: datetime, altitude: float) -> Any:
        values = self.values
        return WalkerConstellation(
            name="walker",
            configuration=values.configuration,
            orbit=CircularOrbit(
                altitude=altitude,
                inclination=values.inclination,
                true_anomaly=getattr(values, "true_anomaly", 0.0),
                right_ascension_ascending_node=getattr(values, "node", 0.0),
                epoch=epoch,
            ),
            number_planes=values.planes,
            number_satellites=values.planes * values.per_plane,
            relative_spacing=values.relative_spacing,
        )

    def describe(self) -> str:
        values = self.values
        star = "star " if values.configuration == "star" else ""
        return (
            f"Walker {star}{float(values.inclination):g} deg{self._altitude()}: "
            + f"{int(values.planes)} x {int(values.per_plane)} "
            + f"(spacing {int(values.relative_spacing)})"
        )


class _PlanesLayout(_ValuesLayout):
    """
    Planes of circular orbits, each with its own number of satellites,
    evenly spaced or in a train (see `sso_planes_space`).
    """

    def _describe_plane(self, plane: Any) -> str:
        """Gets a short description of a plane (its orbit)."""
        raise NotImplementedError

    def _describe_planes(self) -> str:
        """Gets a short description of the planes."""
        described = []
        for plane in self.values.planes:
            text = f"{int(plane.per_plane)} at {self._describe_plane(plane)}"
            interval = getattr(plane, "interval", None)
            if interval is not None:
                if isinstance(interval, timedelta):
                    interval = interval.total_seconds() / 60
                text += f" (train, {float(interval):g} min)"
            described.append(text)
        return ", ".join(described)

    def planes(self) -> list[Plane]:
        return [Plane(int(plane.per_plane), plane) for plane in self.values.planes]

    def _orbit(self, plane: Any, epoch: datetime, altitude: float, anomaly: float):
        """Gets the orbit of a satellite of a plane at a true anomaly."""
        raise NotImplementedError

    def orbits(self, epoch: datetime, altitude: float | None = None) -> list[list[Any]]:
        altitude = _get_altitude(self.values, altitude)
        repeat_ground_track = getattr(self.values, "repeat_ground_track", True)
        return [
            _plane_orbits(
                partial(self._orbit, plane, epoch, altitude), plane, repeat_ground_track
            )
            for plane in self.values.planes
        ]


class SSOPlanesLayout(_PlanesLayout):
    """
    Planes of sun-synchronous orbits (see `tatc.schemas.SunSynchronousOrbit`):
    the values of `sso_planes_space` in a design. Its satellites are named
    `"<name> <p>-<s>"` (plane `p`, satellite `s`).

    Args:
        values (Values): The values of the layout.
    """

    def _orbit(self, plane: Any, epoch: datetime, altitude: float, anomaly: float):
        return SunSynchronousOrbit(
            altitude=altitude,
            equator_crossing_time=hours_to_time(plane.crossing),
            equator_crossing_ascending=self.values.ascending,
            epoch=epoch,
            true_anomaly=anomaly,
        )

    def _describe_plane(self, plane: Any) -> str:
        return hours_to_time(plane.crossing).strftime("%H:%M")

    def describe(self) -> str:
        descending = "" if getattr(self.values, "ascending", True) else " descending"
        return f"SSO{descending}{self._altitude()}: {self._describe_planes()}"


class PlanesLayout(_PlanesLayout):
    """
    Planes of inclined circular orbits (see `tatc.schemas.CircularOrbit`):
    the values of `planes_space` in a design. Its satellites are named
    `"<name> <p>-<s>"` (plane `p`, satellite `s`).

    Args:
        values (Values): The values of the layout.
    """

    def _orbit(self, plane: Any, epoch: datetime, altitude: float, anomaly: float):
        return CircularOrbit(
            altitude=altitude,
            inclination=plane.inclination,
            right_ascension_ascending_node=float(plane.node) % 360,
            epoch=epoch,
            true_anomaly=anomaly,
        )

    def _describe_plane(self, plane: Any) -> str:
        node = float(plane.node)
        known = "" if node != node else f" (node {node % 360:g} deg)"
        return f"{float(plane.inclination):g} deg{known}"

    def describe(self) -> str:
        return f"Planes{self._altitude()}: {self._describe_planes()}"


class SOCLayout(_MembersLayout):
    """
    A streets-of-coverage constellation (see
    `tatc.schemas.SOCConstellation`): the values of `soc_space` in a design,
    whose planes and satellites in each plane it derived at its altitude.
    Its satellites are named `"<name> <k>"`, in order.

    Args:
        values (Values): The values of the layout.
    """

    _MISMATCH = (
        "the altitude of a streets-of-coverage layout must be the one its "
        + "planes were derived at"
    )

    def planes(self) -> list[Plane]:
        return [Plane(int(self.values.per_plane))] * int(self.values.planes)

    def _constellation(self, epoch: datetime, altitude: float) -> Any:
        return _soc_of(self.values, altitude, epoch)

    def describe(self) -> str:
        values = self.values
        polar = (
            "polar " if _soc_of(values, _get_altitude(values, None)).is_polar() else ""
        )
        return (
            f"SOC {polar}{float(values.inclination):g} deg{self._altitude()}: "
            + f"{int(values.planes)} x {int(values.per_plane)}"
        )


class MOGLayout(_MembersLayout):
    """
    A mutual orbiting group (see `tatc.schemas.MOGConstellation`): the
    values of `mog_space` in a design, a single plane (launched together).
    Its satellites are named `"<name> <k>"`, in order.

    Args:
        values (Values): The values of the layout.
    """

    def planes(self) -> list[Plane]:
        return [Plane(int(self.values.satellites))]

    def _constellation(self, epoch: datetime, altitude: float) -> Any:
        values = self.values
        return MOGConstellation(
            name="mog",
            orbit=CircularOrbit(
                altitude=altitude,
                inclination=values.inclination,
                right_ascension_ascending_node=float(values.node) % 360,
                true_anomaly=getattr(values, "true_anomaly", 0.0),
                epoch=epoch,
            ),
            parallel_axis=values.parallel_axis,
            transverse_axis=values.transverse_axis,
            clockwise=bool(values.clockwise),
            number_satellites=int(values.satellites),
        )

    def describe(self) -> str:
        values = self.values
        return (
            f"MOG {float(values.inclination):g} deg{self._altitude()}: "
            + f"{int(values.satellites)} satellites "
            + f"({float(values.parallel_axis) / 1e3:g} x "
            + f"{float(values.transverse_axis) / 1e3:g} km)"
        )


class Layouts(Layout):
    """
    Several layouts flown together (such as the shells of a mixed
    constellation): their planes in turn, whose satellites are named with a
    letter for each layout (`"<name> A 1"`, `"<name> B 1-1"`).

    Args:
        layouts (Sequence): The layouts, or the values of each (see
            `layout_of`).
    """

    def __init__(self, layouts: Sequence[Any]):
        self.layouts = [layout_of(layout) for layout in layouts]

    def planes(self) -> list[Plane]:
        return [plane for layout in self.layouts for plane in layout.planes()]

    def orbits(self, epoch: datetime, altitude: float | None = None) -> list[list[Any]]:
        return [
            orbits
            for layout in self.layouts
            for orbits in layout.orbits(epoch, altitude)
        ]

    def names(self, name: str) -> list[list[str]]:
        return [
            names
            for k, layout in enumerate(self.layouts)
            for names in layout.names(f"{name} {chr(65 + k) if k < 26 else k + 1}")
        ]

    def describe(self) -> str:
        return "; ".join(layout.describe() for layout in self.layouts)

    def _unrealized(self) -> bool:
        # pylint: disable-next=protected-access
        return any(layout._unrealized() for layout in self.layouts)


_PATTERNS = {
    WALKER: WalkerLayout,
    SSO_PLANES: SSOPlanesLayout,
    PLANES: PlanesLayout,
    SOC: SOCLayout,
    MOG: MOGLayout,
}
"""The layouts of the patterns of the fragments of design spaces (their
values `pattern`)."""


def layout_of(layout: Any) -> Layout:
    """
    Gets the layout of the values of a design: those of a fragment
    (`walker_space`, `sso_planes_space`, `planes_space`, `soc_space`, or
    `mog_space`, or of a choice of them), or of several (a list or a
    repetition of layouts, see `Layouts`). A layout is its own.

    Args:
        layout (Any): The values of the layout (or layouts), or a layout.

    Returns:
        Layout: the layout

    Raises:
        ValueError: if the values are not those of a layout.
    """
    if isinstance(layout, Layout):
        return layout
    if isinstance(layout, (list, tuple, RepeatValues)):
        return Layouts(list(layout))
    if isinstance(layout, ChoiceValues):
        layout = layout.value
    pattern = getattr(layout, "pattern", None)
    if not isinstance(pattern, str) or pattern not in _PATTERNS:
        raise ValueError(
            "a layout must be the values of walker_space, sso_planes_space, "
            + "planes_space, soc_space, or mog_space, or a Layout"
        )
    return _PATTERNS[pattern](layout)


def count_satellites(layout: Any) -> int:
    """
    Counts the satellites of a layout (the values of a fragment such as
    `walker_space`, or of a choice of them, or of several, or a `Layout`,
    see `layout_of`), for objectives and constraints.

    Args:
        layout (Any): The values of the layout (or layouts), or a layout.

    Returns:
        int: the number of satellites
    """
    return layout_of(layout).satellites


def count_planes(layout: Any) -> int:
    """
    Counts the orbit planes of a layout (the values of a fragment such as
    `walker_space`, or of a choice of them, or of several, or a `Layout`,
    see `layout_of`), each of which needs its own launch.

    Args:
        layout (Any): The values of the layout (or layouts), or a layout.

    Returns:
        int: the number of planes
    """
    return len(layout_of(layout).planes())


def plane_sizes(layout: Any) -> list[int]:
    """
    Gets the number of satellites of each orbit plane of a layout (the
    values of a fragment such as `walker_space`, or of a choice of them, or
    of several, or a `Layout`, see `layout_of`).

    Args:
        layout (Any): The values of the layout (or layouts), or a layout.

    Returns:
        list[int]: the number of satellites of each plane
    """
    return [plane.satellites for plane in layout_of(layout).planes()]


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
    payload: Any,
    name: str,
    orbit: Any,
    *,
    time_step: float | None,
    start: datetime | None,
    end: datetime | None,
    grade: Grade | str | None = None,
    platform: str | None = None,
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
    from ..catalog import build_payload

    return build_payload(
        payloads,
        orbit=orbit,
        grade=grade,
        time_step=time_step,
        platform=platform,
        start=start,
        end=end,
    ).to_satellite(orbit, name=name)


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
    count = int(plane.per_plane)
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
    layout: Any,
    payload: Payload,
    *,
    epoch: datetime,
    altitude: float | None = None,
    name: str = "satellite",
    time_step: float | None = None,
    start: datetime | None = None,
    end: datetime | None = None,
    grade: Grade | str | None = None,
    platform: str | None = None,
) -> list[Satellite]:
    """
    Builds the satellites of a layout: the values of a fragment (such as
    `walker_space`) in a design (or of a choice of them), or a `Layout` (see
    `layout_of`), at the layout's own altitude (or else the given one), each
    with a payload:

    * an `Instrument` (or a list of instruments);
    * the identifier of an entry of the instrument catalog, optionally with
      a configuration after a colon (such as `"atms:swath_noaa_20"`), or a
      list of them, whose instruments the satellite carries together with
      the attitude they share (see `tatc.catalog.build_payload`), built for
      the orbit and, if a configuration tailors a field of view, the
      analysis `time_step`, or, if it schedules its tilt by day, the
      analysis `start` and `end`; a satellite with several instruments (of
      several entries, or the parts of one view) is analyzed with
      `instrument_index=None`;
    * a function of a name and an orbit that builds a satellite.

    Satellites of Walker and streets-of-coverage constellations and of
    mutual orbiting groups are named `"<name> <k>"`, and those of planes
    `"<name> <p>-<s>"` (plane `p`, satellite `s`), in the order of the
    planes (see `Layout`).

    Several layouts (a list or a repetition of layouts, such as the shells
    of a mixed constellation, see `Layouts`) are built in turn, with the
    same payload, and the names of their satellites are prefixed by a
    letter for each layout (`"<name> A 1"`, `"<name> B 1-1"`). Build layouts
    with different payloads in separate calls, with different names.

    Args:
        layout (Any): The values of the layout (or layouts), or a layout.
        payload (Payload): The payload of each satellite.
        epoch (datetime.datetime): The epoch of the orbits.
        altitude (float | None): The altitude (meters), of layouts that do
            not give their own.
        name (str): The prefix of the names of the satellites.
        time_step (float | None): The analysis time step (seconds), for
            catalog configurations that tailor a field of view to it.
        start (datetime.datetime | None): The start of the analysis, for
            catalog configurations that schedule their tilt by day.
        end (datetime.datetime | None): The end of the analysis, for
            catalog configurations that schedule their tilt by day.
        grade (Grade | str | None): The grade of the use, if any, which
            every catalog configuration must serve (see
            `tatc.catalog.build_payload`).
        platform (str | None): The identifier of the platform modeled, if
            any, which selects the catalog configurations tested on it (with
            a warning for a configuration not tested on it).

    Returns:
        list[Satellite]: the satellites
    """
    return layout_of(layout).build(
        payload,
        epoch=epoch,
        altitude=altitude,
        name=name,
        time_step=time_step,
        start=start,
        end=end,
        grade=grade,
        platform=platform,
    )


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

    It is a snapshot, for the build function of any study: the constellation
    with some of its satellites missing, without a deployment or the time of
    a mission. The failures of an `ObservationStudy` with a `MissionPlan`
    (see `Failures`) are those over a mission instead: the lives of
    satellites, the failures of launches, and their restoration, averaged
    over the phases of the mission.

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
