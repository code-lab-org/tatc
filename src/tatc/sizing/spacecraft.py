"""
Sizing of spacecraft from their payloads and orbits: the dry mass and power
of a bus that carries a payload, the propellant of its velocity budget, and
whether its orbit lasts its life and meets its disposal rule.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from datetime import datetime
from functools import lru_cache

from pydantic import BaseModel, ConfigDict, Field, model_validator

from ..constants import EARTH_MEAN_RADIUS
from ..schemas.resources import SpacecraftResources
from ..schemas.space import Satellite
from ..schemas.space.base_constellation import BaseConstellation
from ..utils.atmosphere import SolarActivity, SolarCycle
from ..utils.decay import DEFAULT_DRAG_COEFFICIENT, compute_orbit_lifetime
from .budget import (
    DeltaVBudget,
    DisposalRule,
    compute_delta_v_budget,
    compute_propellant_mass,
)

LIFETIME_HORIZON = 200.0
"""float: The longest orbit lifetime (years) that sizing computes; longer
lifetimes are infinite."""


class SpacecraftSizing(BaseModel):
    """
    The assumptions by which to size a spacecraft (see `size_spacecraft`):
    values that a spacecraft's own resources (`SpacecraftResources`) do not
    give, such as those of a notional spacecraft in a trade study. Its dry
    mass and power may be given, or sized from those of its payload (the
    total of its instruments' resources, see `InstrumentResources`) by the
    fractions of the dry mass and power the payload takes. TAT-C gives no
    default fractions: they vary widely by mission class, so a study states
    its own, with their sources. Its drag is given by a drag area (and drag
    coefficient), or by a ballistic coefficient. Without a specific impulse
    (here or in its resources), a spacecraft has no propulsion: its orbit
    decays freely, which must last its life and meet its disposal rule.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    payload_mass_fraction: float | None = Field(
        None,
        gt=0,
        le=1,
        description="The fraction of the dry mass that is the payload's, by "
        + "which to size the dry mass, or `None` if the dry mass is given.",
    )
    payload_power_fraction: float | None = Field(
        None,
        gt=0,
        le=1,
        description="The fraction of the spacecraft's power that the payload "
        + "uses, by which to size the power, or `None`.",
    )
    dry_mass: float | None = Field(None, gt=0, description="The dry mass (kg).")
    power: float | None = Field(
        None, gt=0, description="The power (watts) of the spacecraft."
    )
    design_life: float | None = Field(
        None, gt=0, description="The design life (years)."
    )
    drag_area: float | None = Field(
        None, gt=0, description="The mean drag area (square meters)."
    )
    drag_coefficient: float | None = Field(
        None,
        gt=0,
        description="The drag coefficient, or `None` for 2.2 (see "
        + "`tatc.utils.DEFAULT_DRAG_COEFFICIENT`).",
    )
    ballistic_coefficient: float | None = Field(
        None,
        gt=0,
        description="The ballistic coefficient (kg/m^2) at the dry mass, in "
        + "place of a drag area.",
    )
    specific_impulse: float | None = Field(
        None,
        gt=0,
        description="The specific impulse (seconds) of its propulsion, or `None` "
        + "for none.",
    )
    activity: SolarActivity = Field(
        "moderate",
        description="The solar activity (see `tatc.utils.SolarActivity`).",
    )
    disposal: DisposalRule | None = Field(
        DisposalRule(),
        description="The rule of its disposal (see `DisposalRule`), or `None` "
        + "not to dispose of it.",
    )
    injection_altitude: float | None = Field(
        None,
        ge=0,
        description="The altitude (meters above the mean radius) at which a "
        + "launch injects it, from which it raises its orbit, or `None` for "
        + "its own orbit.",
    )
    other_delta_v: float = Field(
        0.0, ge=0, description="The velocity (m/s) of other maneuvers."
    )
    margin: float = Field(
        0.0, ge=0, description="The margin of the velocity budget, as a fraction."
    )

    @model_validator(mode="after")
    def _check_drag(self) -> SpacecraftSizing:
        if self.drag_area is not None and self.ballistic_coefficient is not None:
            raise ValueError("give a drag area or a ballistic coefficient, not both")
        if self.ballistic_coefficient is not None and self.drag_coefficient is not None:
            raise ValueError(
                "a ballistic coefficient includes the drag coefficient; give one"
            )
        return self


class SizedSpacecraft(BaseModel):
    """
    A spacecraft sized for its orbit and life (see `size_spacecraft`).
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    resources: SpacecraftResources = Field(
        ...,
        description="Its resources: those given, with the dry mass, power, "
        + "design life, launch mass (the dry mass and propellant), and "
        + "properties of drag and propulsion sized.",
    )
    ballistic_coefficient: float = Field(
        ..., gt=0, description="The ballistic coefficient (kg/m^2) at its dry mass."
    )
    lifetime: float = Field(
        ...,
        ge=0,
        description="The lifetime (years) of its orbit without propulsion, from "
        + "the start of its life (infinite beyond 200 years).",
    )
    budget: DeltaVBudget | None = Field(
        ..., description="Its velocity budget, or `None` without propulsion."
    )
    propellant_mass: float = Field(
        ..., ge=0, description="The propellant (kg) its budget needs."
    )
    payload_power: float | None = Field(
        ..., description="The power (watts) of its payload, if given."
    )
    violations: dict[str, float] = Field(
        default_factory=dict,
        description="The amount by which it violates each constraint that it "
        + "does, relative to the constraint: `decay` (without propulsion, its "
        + "orbit decays within its life: the shortfall over the life), "
        + "`disposal` (without propulsion, its orbit outlasts the disposal "
        + "rule's limit: the excess over the limit, or 1 for a controlled "
        + "reentry), `power` (its payload uses more than its power: the excess "
        + "over the power), and `propellant` (its budget needs more propellant "
        + "than its launch mass holds: the shortfall over the need).",
    )

    @property
    def feasible(self) -> bool:
        """bool: Whether it violates no constraint."""
        return not self.violations

    @property
    def power_margin(self) -> float:
        """float: The fraction of its power that its payload does not use, or
        NaN if either is not given."""
        if self.payload_power is None or self.resources.power is None:
            return math.nan
        return 1 - self.payload_power / self.resources.power


def _first(*values: float | None) -> float | None:
    """
    Gets the first value that is given.
    """
    return next((value for value in values if value is not None), None)


@lru_cache(maxsize=4096)
def _compute_lifetime(
    semimajor_axis: float,
    eccentricity: float,
    inclination: float,
    ballistic_coefficient: float,
    activity: SolarActivity,
    start: datetime | None,
) -> float:
    """
    Computes the lifetime of an orbit (see `compute_orbit_lifetime`), cached
    for the many spacecraft that share orbits and properties.
    """
    return compute_orbit_lifetime(
        semimajor_axis,
        eccentricity,
        inclination,
        ballistic_coefficient,
        activity=activity,
        start=start,
        horizon=LIFETIME_HORIZON,
    )


def size_spacecraft(
    spacecraft: Satellite | BaseConstellation,
    sizing: SpacecraftSizing | None = None,
    *,
    design_life: float | None = None,
    start: datetime | None = None,
) -> SizedSpacecraft:
    """
    Sizes a spacecraft (or each member of a constellation) for its orbit and
    life. Each value is the spacecraft's own (its `resources`), if given, or
    else the sizing's (see `SpacecraftSizing`): the dry mass, or the mass of
    its payload over its fraction of the dry mass; the power, or the power
    of its payload over its fraction of the power; the design life; the
    drag area and coefficient, or the ballistic coefficient; and the specific
    impulse. With a specific impulse, its velocity budget (see
    `compute_delta_v_budget`, holding its orbit against drag for its life,
    with the ballistic coefficient at its dry mass) needs propellant (see
    `compute_propellant_mass`), which a given launch mass must hold, and its
    launch mass is otherwise its dry mass and propellant. Without one, its
    orbit decays freely: its lifetime must exceed its life, and outlast it
    by no more than its disposal rule's limit (a controlled reentry is not
    possible).

    Args:
        spacecraft (Satellite | BaseConstellation): The spacecraft, with its
            orbit, instruments, and resources.
        sizing (SpacecraftSizing | None): The assumptions of values that its
            resources do not give, if any.
        design_life (float | None): The design life (years) in place of
            those of its resources and the sizing, if given.
        start (datetime.datetime | None): The start of its life, by default
            the epoch of its orbit (as of a solar cycle).

    Returns:
        SizedSpacecraft: The sized spacecraft.

    Raises:
        ValueError: If a value it needs is not given: the dry mass (or the
            payload mass and its fraction), the design life, or the drag
            area (or ballistic coefficient).
    """
    sizing = SpacecraftSizing() if sizing is None else sizing
    given = (
        SpacecraftResources() if spacecraft.resources is None else spacecraft.resources
    )
    payload = spacecraft.get_payload_resources()
    dry_mass = _first(given.dry_mass, sizing.dry_mass)
    if dry_mass is None and sizing.payload_mass_fraction is not None:
        if payload.mass is None:
            raise ValueError(
                "sizing the dry mass by the payload's fraction requires the "
                + "masses of its instruments (InstrumentResources.mass)"
            )
        dry_mass = payload.mass / sizing.payload_mass_fraction
    if dry_mass is None:
        raise ValueError(
            "the dry mass is not given: give it (in the spacecraft's resources or "
            + "the sizing), or the payload's fraction of it"
        )
    power = _first(given.power, sizing.power)
    if (
        power is None
        and sizing.payload_power_fraction is not None
        and payload.power is not None
    ):
        power = payload.power / sizing.payload_power_fraction
    life = _first(design_life, given.design_life, sizing.design_life)
    if life is None:
        raise ValueError(
            "the design life is not given: give it (in the spacecraft's resources "
            + "or the sizing)"
        )
    drag_area = _first(given.drag_area, sizing.drag_area)
    drag_coefficient = _first(given.drag_coefficient, sizing.drag_coefficient)
    if drag_area is not None:
        ballistic_coefficient = dry_mass / (
            _first(drag_coefficient, DEFAULT_DRAG_COEFFICIENT) * drag_area
        )
    elif sizing.ballistic_coefficient is not None:
        ballistic_coefficient = sizing.ballistic_coefficient
    else:
        raise ValueError(
            "the drag area is not given: give it (in the spacecraft's resources or "
            + "the sizing), or a ballistic coefficient"
        )
    specific_impulse = _first(given.specific_impulse, sizing.specific_impulse)

    orbit = spacecraft.orbit
    start = orbit.get_epoch() if start is None else start
    activity = sizing.activity
    cycle_start = start if isinstance(activity, SolarCycle) else None
    semimajor_axis = orbit.get_semimajor_axis()
    eccentricity = orbit.get_eccentricity()
    inclination = orbit.get_inclination()
    lifetime = _compute_lifetime(
        float(semimajor_axis),
        float(eccentricity),
        float(inclination),
        float(ballistic_coefficient),
        activity,
        cycle_start,
    )
    violations: dict[str, float] = {}
    if specific_impulse is None:
        budget, propellant = None, 0.0
        if lifetime < life:
            violations["decay"] = (life - lifetime) / life
        rule = sizing.disposal
        perigee = semimajor_axis * (1 - eccentricity) - EARTH_MEAN_RADIUS
        if rule is not None and perigee <= rule.region_altitude:
            if rule.perigee_altitude is not None:
                # a controlled reentry needs propulsion
                violations["disposal"] = 1.0
            elif lifetime - life > rule.limit:
                remaining = min(lifetime, LIFETIME_HORIZON) - life
                violations["disposal"] = (remaining - rule.limit) / rule.limit
    else:
        budget = compute_delta_v_budget(
            orbit,
            ballistic_coefficient,
            life,
            activity=activity,
            start=start,
            disposal=sizing.disposal,
            injection_altitude=sizing.injection_altitude,
            other=sizing.other_delta_v,
            margin=sizing.margin,
        )
        propellant = compute_propellant_mass(dry_mass, budget.total, specific_impulse)
        if given.launch_mass is not None:
            available = given.launch_mass - dry_mass
            if propellant > available:
                violations["propellant"] = (propellant - available) / propellant
    if power is not None and payload.power is not None and payload.power > power:
        violations["power"] = payload.power / power - 1
    resources = SpacecraftResources(
        **{
            **given.model_dump(),
            "dry_mass": dry_mass,
            "launch_mass": _first(given.launch_mass, dry_mass + propellant),
            "power": power,
            "design_life": life,
            "drag_area": drag_area,
            "drag_coefficient": drag_coefficient,
            "specific_impulse": specific_impulse,
        }
    )
    return SizedSpacecraft(
        resources=resources,
        ballistic_coefficient=ballistic_coefficient,
        lifetime=lifetime,
        budget=budget,
        propellant_mass=propellant,
        payload_power=payload.power,
        violations=violations,
    )
