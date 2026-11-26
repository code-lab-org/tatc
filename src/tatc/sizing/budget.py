"""
Velocity budgets of the orbits of spacecraft and the propellant they need.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from datetime import datetime, timedelta
from functools import lru_cache

from pydantic import BaseModel, ConfigDict, Field, computed_field

from ..constants import EARTH_MEAN_RADIUS, EARTH_MU, STANDARD_GRAVITY
from ..schemas.orbit import AllOrbits
from ..utils.atmosphere import YEAR_S, SolarActivity, SolarCycle
from ..utils.decay import compute_drag_delta_v, compute_orbit_lifetime
from ..utils.orbital import compute_hohmann_delta_v

MAX_ECCENTRICITY = 0.01
"""float: The largest eccentricity of the near-circular orbits whose
velocity budgets are computed (see `compute_delta_v_budget`)."""

_PERIGEE_TOLERANCE = 1e3
"""The tolerance (meters) of the perigee of a disposal orbit."""


class DisposalRule(BaseModel):
    """
    The rule by which a spacecraft in low Earth orbit is disposed of at the
    end of its life: by lowering its perigee until drag brings it down
    within a time limit (by default, 25 years, the guideline of the
    Inter-Agency Space Debris Coordination Committee, IADC-02-01, and of ISO
    24113; 5 years under the U.S. Federal Communications Commission's rule
    of 2022, FCC 22-74, and ESA's Zero Debris approach), or to a perigee
    altitude low enough to reenter at once (a controlled reentry). Orbits
    whose perigees are above the protected region of low Earth orbit (by
    default, 2000 km) are not disposed of by these rules (as geostationary
    satellites are, to graveyard orbits).
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    limit: float = Field(
        25.0,
        gt=0,
        description="The longest time (years) after the end of its life that "
        + "a spacecraft may remain in orbit.",
    )
    perigee_altitude: float | None = Field(
        None,
        ge=0,
        description="The perigee altitude (meters above the mean radius) to "
        + "which to lower the orbit for a controlled reentry, or `None` to "
        + "lower it only as far as the time limit requires.",
        examples=[50e3],
    )
    region_altitude: float = Field(
        2000e3,
        gt=0,
        description="The highest perigee altitude (meters above the mean "
        + "radius) of the orbits to dispose of.",
    )


class DeltaVBudget(BaseModel):
    """
    The velocity budget (meters per second) of a spacecraft's orbit over its
    life (see `compute_delta_v_budget`).
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    raising: float = Field(
        0.0, ge=0, description="To raise its orbit from its injection by a launch."
    )
    drag: float = Field(
        0.0, ge=0, description="To hold its orbit against drag over its life."
    )
    disposal: float = Field(0.0, ge=0, description="To dispose of it.")
    other: float = Field(
        0.0, ge=0, description="For other maneuvers, such as collision avoidance."
    )
    margin: float = Field(0.0, ge=0, description="The margin on the others.")

    @computed_field  # type: ignore[prop-decorator]
    @property
    def total(self) -> float:
        """float: The total velocity (meters per second)."""
        return self.raising + self.drag + self.disposal + self.other + self.margin


def compute_propellant_mass(
    dry_mass: float, delta_v: float, specific_impulse: float
) -> float:
    """
    Computes the mass of propellant that gives a spacecraft a velocity, by
    the rocket equation: `dry_mass * (exp(delta_v / (specific_impulse * g0))
    - 1)`, for standard gravity `g0`.

    Args:
        dry_mass (float): The dry mass (kg).
        delta_v (float): The velocity (meters per second).
        specific_impulse (float): The specific impulse (seconds).

    Returns:
        float: The mass of propellant (kg).

    Raises:
        ValueError: If a mass or velocity is negative or the specific
            impulse is not positive.
    """
    if dry_mass < 0 or delta_v < 0:
        raise ValueError("the dry mass and velocity must not be negative")
    if specific_impulse <= 0:
        raise ValueError("the specific impulse must be positive")
    return dry_mass * math.expm1(delta_v / (specific_impulse * STANDARD_GRAVITY))


def _get_lowering_delta_v(radius: float, perigee_radius: float) -> float:
    """
    Gets the velocity of a burn that lowers the perigee of a circular orbit
    of a radius to another radius.
    """
    circular = math.sqrt(EARTH_MU / radius)
    apogee = math.sqrt(EARTH_MU * (2 / radius - 2 / (radius + perigee_radius)))
    return circular - apogee


def compute_disposal_delta_v(
    altitude: float,
    inclination: float,
    ballistic_coefficient: float,
    *,
    rule: DisposalRule | None = None,
    activity: SolarActivity = "moderate",
    start: datetime | None = None,
) -> float:
    """
    Computes the velocity to dispose of a spacecraft in a circular orbit by a
    rule (see `DisposalRule`): to lower its perigee to the rule's perigee
    altitude, if given; or else none if drag brings it down within the
    rule's time limit, and otherwise to lower its perigee to the highest
    altitude (within 1 km) from which drag brings it down within the limit
    (see `tatc.utils.compute_orbit_lifetime`).

    Args:
        altitude (float): The altitude (meters above the mean radius).
        inclination (float): The inclination (degrees).
        ballistic_coefficient (float): The ballistic coefficient (kg/m^2) of
            the spacecraft at its disposal.
        rule (DisposalRule | None): The rule, or `None` for the default rule.
        activity (SolarActivity): The solar activity after its disposal.
        start (datetime.datetime | None): The time of its disposal, required
            for a solar cycle.

    Returns:
        float: The velocity (meters per second).
    """
    rule = DisposalRule() if rule is None else rule
    if altitude > rule.region_altitude:
        return 0.0
    if not isinstance(activity, SolarCycle):
        start = None  # the same for any start
    return _compute_disposal_delta_v(
        float(altitude),
        float(inclination),
        float(ballistic_coefficient),
        rule,
        activity,
        start,
    )


@lru_cache(maxsize=4096)
def _compute_disposal_delta_v(
    altitude: float,
    inclination: float,
    ballistic_coefficient: float,
    rule: DisposalRule,
    activity: SolarActivity,
    start: datetime | None,
) -> float:
    """
    Computes the velocity to dispose of a spacecraft (see
    `compute_disposal_delta_v`), whose searches of perigees are cached for
    the many spacecraft that share orbits and properties.
    """
    radius = EARTH_MEAN_RADIUS + altitude
    if rule.perigee_altitude is not None:
        if rule.perigee_altitude >= altitude:
            return 0.0
        return _get_lowering_delta_v(radius, EARTH_MEAN_RADIUS + rule.perigee_altitude)

    def lifetime(perigee: float) -> float:
        perigee_radius = EARTH_MEAN_RADIUS + perigee
        return compute_orbit_lifetime(
            (radius + perigee_radius) / 2,
            (radius - perigee_radius) / (radius + perigee_radius),
            inclination,
            ballistic_coefficient,
            activity=activity,
            start=start,
            horizon=rule.limit,
        )

    if lifetime(altitude) <= rule.limit:
        return 0.0
    # the highest perigee (lower) that comes down within the limit, and the
    # lowest (upper) that does not
    lower, upper = 100e3, altitude
    while upper - lower > _PERIGEE_TOLERANCE:
        middle = (lower + upper) / 2
        if lifetime(middle) <= rule.limit:
            lower = middle
        else:
            upper = middle
    return _get_lowering_delta_v(radius, EARTH_MEAN_RADIUS + lower)


def compute_delta_v_budget(
    orbit: AllOrbits,
    ballistic_coefficient: float,
    design_life: float,
    *,
    activity: SolarActivity = "moderate",
    start: datetime | None = None,
    disposal: DisposalRule | None = DisposalRule(),
    injection_altitude: float | None = None,
    other: float = 0.0,
    margin: float = 0.0,
) -> DeltaVBudget:
    """
    Computes the velocity budget of a spacecraft in a near-circular orbit
    over its life: to raise its orbit from the altitude at which a launch
    injects it (a Hohmann transfer), to hold its altitude against drag (see
    `tatc.utils.compute_drag_delta_v`), to dispose of it at the end of its
    life (see `compute_disposal_delta_v`), for other maneuvers, and a margin
    on these. Its altitude is the orbit's mean altitude, held for its life,
    so that the solar activity of its disposal is that at the end of its
    life.

    Args:
        orbit (AllOrbits): The orbit, with an eccentricity of at most 0.01
            (see `MAX_ECCENTRICITY`).
        ballistic_coefficient (float): The ballistic coefficient (kg/m^2).
        design_life (float): The life (years) for which the orbit is held.
        activity (SolarActivity): The solar activity.
        start (datetime.datetime | None): The start of its life, by default
            the epoch of the orbit (as of a solar cycle).
        disposal (DisposalRule | None): The rule of its disposal, or `None`
            not to dispose of it.
        injection_altitude (float | None): The altitude (meters above the mean
            radius) of the circular orbit into which a launch injects it, or
            `None` for its own orbit.
        other (float): The velocity (meters per second) of other maneuvers.
        margin (float): The margin, as a fraction of the other velocities.

    Returns:
        DeltaVBudget: The budget.

    Raises:
        ValueError: If the orbit is not near-circular, the design life or a
            velocity is negative, or the altitude is below 100 km.
    """
    eccentricity = orbit.get_eccentricity()
    if eccentricity > MAX_ECCENTRICITY:
        raise ValueError(
            f"velocity budgets are of near-circular orbits (eccentricity at most "
            + f"{MAX_ECCENTRICITY}), not {eccentricity:.4f}"
        )
    if design_life < 0 or other < 0 or margin < 0:
        raise ValueError(
            "the design life, other velocity, and margin must not be negative"
        )
    radius = orbit.get_semimajor_axis()
    altitude = radius - EARTH_MEAN_RADIUS
    inclination = orbit.get_inclination()
    start = orbit.get_epoch() if start is None else start
    raising = (
        0.0
        if injection_altitude is None
        else compute_hohmann_delta_v(EARTH_MEAN_RADIUS + injection_altitude, radius)
    )
    drag = compute_drag_delta_v(
        altitude,
        inclination,
        ballistic_coefficient,
        design_life,
        activity=activity,
        start=start,
    )
    disposed = (
        0.0
        if disposal is None
        else compute_disposal_delta_v(
            altitude,
            inclination,
            ballistic_coefficient,
            rule=disposal,
            activity=activity,
            start=start + timedelta(seconds=design_life * YEAR_S),
        )
    )
    return DeltaVBudget(
        raising=raising,
        drag=drag,
        disposal=disposed,
        other=other,
        margin=margin * (raising + drag + disposed + other),
    )
