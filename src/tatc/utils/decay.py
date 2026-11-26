"""
Methods for the decay of orbits by atmospheric drag: the rates at which drag
shrinks and circularizes an orbit, its decay over time and lifetime, and the
velocity a satellite must add to hold its altitude against drag.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from datetime import datetime

import numpy as np
import numpy.typing as npt
import pandas as pd

from ..constants import EARTH_MEAN_RADIUS, EARTH_MU, EARTH_ROTATION_RATE
from .atmosphere import (
    MAX_DENSITY_ALTITUDE,
    MIN_DENSITY_ALTITUDE,
    YEAR_S,
    SolarActivity,
    SolarCycle,
    get_solar_activity,
    get_years_between,
    interpolate_density,
)

DEFAULT_DRAG_COEFFICIENT = 2.2
"""float: The drag coefficient of a satellite whose coefficient is not given:
2.2, the value commonly assumed for compact satellites in free molecular flow
(see, for example, Vallado, "Fundamentals of Astrodynamics and Applications",
4th ed., 2013, section 8.6.2)."""

_PEAK_NODES = np.polynomial.legendre.leggauss(48)
"""Gauss-Legendre nodes and weights on [-1, 1] near the perigee of an
eccentric orbit, where most of its drag acts."""

_TAIL_NODES = np.polynomial.legendre.leggauss(16)
"""Gauss-Legendre nodes and weights on [-1, 1] over the rest of its orbit."""

_STEP_SCALE_HEIGHTS = 0.05
"""The largest fraction of the density scale height at its perigee by which
an orbit may decay in one step of `compute_orbit_decay`."""

_STEPS_PER_CYCLE = 50
"""The fewest steps of `compute_orbit_decay` per period of a solar cycle."""

_SCALE_HEIGHT_DELTA = 5e3
"""The difference of altitudes (meters) by which to compute scale heights."""


def _get_nodes(
    a: float, e: float, scale_height: float
) -> tuple[npt.NDArray, npt.NDArray]:
    """
    Gets the eccentric anomalies (radians) and weights of the quadrature of
    an integral over a half revolution (from perigee, 0, to apogee, pi),
    concentrated within the anomalies near perigee where the density falls by
    a few scale heights, `sqrt(2 H / (a e))`.
    """
    width = math.sqrt(2 * scale_height / (a * e)) if e > 0 else math.inf
    split = min(math.pi, 10 * width)
    x, w = _PEAK_NODES
    anomalies = [split / 2 * (x + 1)]
    weights = [split / 2 * w]
    if split < math.pi:
        x, w = _TAIL_NODES
        anomalies.append(split + (math.pi - split) / 2 * (x + 1))
        weights.append((math.pi - split) / 2 * w)
    return np.concatenate(anomalies), np.concatenate(weights)


def _get_scale_height(altitude: float, activity: float) -> float:
    """
    Gets the scale height (meters) of the density at an altitude (and above
    it), by which it falls by a factor of e, at an index of activity.
    """
    lower, upper = interpolate_density(
        [altitude, altitude + _SCALE_HEIGHT_DELTA], activity
    )
    if not upper > 0 or not lower > upper:
        return math.inf
    return _SCALE_HEIGHT_DELTA / math.log(lower / upper)


def _get_rotation_factor(
    radius: npt.NDArray, speed: npt.NDArray, inclination: float
) -> npt.NDArray:
    """
    Gets the factor `(1 - r w cos(i) / v) ** 2` of the drag of an atmosphere
    that rotates with the Earth (at rate w) on a satellite at radius `r` and
    speed `v` in an orbit inclined at `i`: the square of its speed relative to
    the air in the direction of its motion, relative to its speed (King-Hele,
    "Satellite Orbits in an Atmosphere", 1987).
    """
    return (
        1 - radius * EARTH_ROTATION_RATE * math.cos(math.radians(inclination)) / speed
    ) ** 2


def compute_decay_rates(
    semimajor_axis: float,
    eccentricity: float,
    inclination: float,
    ballistic_coefficient: float,
    activity: float = 0.5,
) -> tuple[float, float]:
    """
    Computes the rates at which drag decays an orbit, averaged over a
    revolution: the changes per revolution of its semimajor axis,
    `-(a ** 2 / B) * integral(rho * F * (1 + e cos E) ** 1.5 / (1 - e cos E) ** 0.5 dE)`,
    and of its eccentricity,
    `-(a (1 - e ** 2) / B) * integral(rho * F * ((1 + e cos E) / (1 - e cos E)) ** 0.5 * cos E dE)`,
    over its eccentric anomaly `E` (King-Hele, "Satellite Orbits in an
    Atmosphere", 1987), divided by its period, for a ballistic coefficient
    `B` (mass over the drag coefficient times area), the density `rho` of
    the thermosphere (see `compute_atmospheric_density`), and the factor `F`
    of the rotation of the atmosphere with the Earth. A circular orbit
    shrinks at `-rho * F * sqrt(mu * a) / B`; an eccentric orbit shrinks and
    circularizes, mostly by drag near its perigee.

    Args:
        semimajor_axis (float): The semimajor axis (meters).
        eccentricity (float): The eccentricity.
        inclination (float): The inclination (degrees).
        ballistic_coefficient (float): The ballistic coefficient (kg/m^2):
            the mass over the drag coefficient times the drag area.
        activity (float): The index of solar activity, from 0 (low) to 1
            (high, see `tatc.utils.get_solar_activity`).

    Returns:
        tuple[float, float]: The rates of change of the semimajor axis
            (meters per second) and the eccentricity (per second), not
            positive, and NaN for a perigee below 100 km.
    """
    a, e = float(semimajor_axis), max(float(eccentricity), 0.0)
    perigee = a * (1 - e) - EARTH_MEAN_RADIUS
    if perigee < MIN_DENSITY_ALTITUDE:
        return math.nan, math.nan
    if perigee > MAX_DENSITY_ALTITUDE:
        return 0.0, 0.0
    period = 2 * math.pi * math.sqrt(a**3 / EARTH_MU)
    if e < 1e-9:
        radius = np.array([a])
        speed = np.sqrt(EARTH_MU / radius)
        density = interpolate_density(radius - EARTH_MEAN_RADIUS, activity)
        factor = _get_rotation_factor(radius, speed, inclination)
        integral_a = float(2 * math.pi * density[0] * factor[0])
        return -(a**2) * integral_a / ballistic_coefficient / period, 0.0
    anomalies, weights = _get_nodes(a, e, _get_scale_height(perigee, activity))
    cosine = np.cos(anomalies)
    radius = a * (1 - e * cosine)
    speed = np.sqrt(EARTH_MU * (2 / radius - 1 / a))
    density = interpolate_density(radius - EARTH_MEAN_RADIUS, activity)
    drag = density * _get_rotation_factor(radius, speed, inclination)
    ratio = (1 + e * cosine) / (1 - e * cosine)
    # integrands are symmetric about perigee: twice the half revolution
    integral_a = 2 * np.sum(weights * drag * ratio**1.5 * (1 - e * cosine))
    integral_e = 2 * np.sum(weights * drag * np.sqrt(ratio) * cosine)
    rate_a = -(a**2) * integral_a / ballistic_coefficient / period
    rate_e = -a * (1 - e**2) * integral_e / ballistic_coefficient / period
    return float(rate_a), float(min(rate_e, 0.0))


def _get_activity_function(activity: SolarActivity, start: datetime | None):
    """
    Gets a function of the index of activity at a number of years from the
    start, and the longest step (years) that resolves its changes.
    """
    if isinstance(activity, SolarCycle):
        if start is None:
            raise ValueError("the decay under a solar cycle requires a start")
        offset = get_years_between(activity.minimum, start)
        return (
            lambda years: activity.get_activity_after(offset + years),
            activity.period / _STEPS_PER_CYCLE,
        )
    index = get_solar_activity(activity)
    return lambda years: index, math.inf


def compute_orbit_decay(
    semimajor_axis: float,
    eccentricity: float,
    inclination: float,
    ballistic_coefficient: float,
    *,
    activity: SolarActivity = "moderate",
    start: datetime | None = None,
    horizon: float = 200.0,
) -> pd.DataFrame:
    """
    Computes the decay of an orbit by drag (see `compute_decay_rates`) until
    its perigee falls to 100 km (see `MIN_DENSITY_ALTITUDE`), where it
    reenters within days, or until a horizon. It integrates the rates
    averaged over revolutions with steps (second-order Runge-Kutta) in which
    the perigee falls by at most a twentieth of the density's scale height
    there, and that resolve a solar cycle. Its accuracy is that of the
    density: the global mean of the thermosphere at a level of activity, or
    under a mean solar cycle, about which real lifetimes vary by a factor of
    about 2 with the timing and strength of solar cycles and storms.

    Args:
        semimajor_axis (float): The semimajor axis (meters).
        eccentricity (float): The eccentricity.
        inclination (float): The inclination (degrees).
        ballistic_coefficient (float): The ballistic coefficient (kg/m^2):
            the mass over the drag coefficient times the drag area.
        activity (SolarActivity): The solar activity: a reference level, an
            index from 0 (low) to 1 (high), or a solar cycle.
        start (datetime.datetime | None): The start of the decay, required
            for a solar cycle.
        horizon (float): The longest time (years) to compute.

    Returns:
        pandas.DataFrame: a row per step: its `time` (years from the start),
            `semimajor_axis` (meters), `eccentricity`, `perigee_altitude`, and
            `apogee_altitude` (meters above the mean radius); the last row is
            at the reentry (a perigee of 100 km) or the horizon.

    Raises:
        ValueError: If the perigee is below 100 km, the ballistic coefficient
            is not positive, or a solar cycle has no start.
    """
    if ballistic_coefficient <= 0:
        raise ValueError("the ballistic coefficient must be positive")
    a, e = float(semimajor_axis), float(eccentricity)
    if a * (1 - e) - EARTH_MEAN_RADIUS < MIN_DENSITY_ALTITUDE:
        raise ValueError("the perigee must be at least 100 km above the mean radius")
    activity_at, longest = _get_activity_function(activity, start)
    longest = min(longest, horizon)

    def rates(years: float, a: float, e: float) -> tuple[float, float]:
        rate_a, rate_e = compute_decay_rates(
            a, e, inclination, ballistic_coefficient, activity_at(years)
        )
        return rate_a * YEAR_S, rate_e * YEAR_S

    def perigee_of(a: float, e: float) -> float:
        return a * (1 - e) - EARTH_MEAN_RADIUS

    time = 0.0
    rows = [(time, a, e)]
    while time < horizon:
        perigee = perigee_of(a, e)
        if perigee > MAX_DENSITY_ALTITUDE:
            break  # no drag: no decay
        rate_a, rate_e = rates(time, a, e)
        falling = -(rate_a * (1 - e) - a * rate_e)  # meters per year
        height = _get_scale_height(perigee, activity_at(time))
        step = min(longest, horizon - time)
        if falling > 0:
            step = min(step, _STEP_SCALE_HEIGHTS * height / falling)
        # second-order Runge-Kutta (midpoint), shortened to stay above 100 km
        while True:
            a_mid = a + step / 2 * rate_a
            e_mid = max(e + step / 2 * rate_e, 0.0)
            mid_a, mid_e = rates(time + step / 2, a_mid, e_mid)
            if not math.isnan(mid_a):
                break
            step /= 2
        a_next = a + step * mid_a
        e_next = max(e + step * mid_e, 0.0)
        next_perigee = perigee_of(a_next, e_next)
        if next_perigee <= MIN_DENSITY_ALTITUDE:
            # the time of reentry, interpolating the perigee linearly
            fraction = (perigee - MIN_DENSITY_ALTITUDE) / (perigee - next_perigee)
            a_end = a + fraction * (a_next - a)
            e_end = 1 - (EARTH_MEAN_RADIUS + MIN_DENSITY_ALTITUDE) / a_end
            rows.append((time + fraction * step, a_end, max(e_end, 0.0)))
            break
        time, a, e = time + step, a_next, e_next
        rows.append((time, a, e))
    if rows[-1][0] < horizon and perigee_of(*rows[-1][1:]) > MIN_DENSITY_ALTITUDE:
        rows.append((horizon, rows[-1][1], rows[-1][2]))
    result = pd.DataFrame(rows, columns=["time", "semimajor_axis", "eccentricity"])
    result["perigee_altitude"] = (
        result.semimajor_axis * (1 - result.eccentricity) - EARTH_MEAN_RADIUS
    )
    result["apogee_altitude"] = (
        result.semimajor_axis * (1 + result.eccentricity) - EARTH_MEAN_RADIUS
    )
    return result


def compute_orbit_lifetime(
    semimajor_axis: float,
    eccentricity: float,
    inclination: float,
    ballistic_coefficient: float,
    *,
    activity: SolarActivity = "moderate",
    start: datetime | None = None,
    horizon: float = 200.0,
) -> float:
    """
    Computes the lifetime of an orbit decayed by drag: the time until its
    perigee falls to 100 km (see `compute_orbit_decay`).

    Args:
        semimajor_axis (float): The semimajor axis (meters).
        eccentricity (float): The eccentricity.
        inclination (float): The inclination (degrees).
        ballistic_coefficient (float): The ballistic coefficient (kg/m^2):
            the mass over the drag coefficient times the drag area.
        activity (SolarActivity): The solar activity: a reference level, an
            index from 0 (low) to 1 (high), or a solar cycle.
        start (datetime.datetime | None): The start of the decay, required
            for a solar cycle.
        horizon (float): The longest lifetime (years) to compute.

    Returns:
        float: The lifetime (years), infinite if longer than the horizon.
    """
    decay = compute_orbit_decay(
        semimajor_axis,
        eccentricity,
        inclination,
        ballistic_coefficient,
        activity=activity,
        start=start,
        horizon=horizon,
    )
    last = decay.iloc[-1]
    if last.perigee_altitude > MIN_DENSITY_ALTITUDE + 1:
        return math.inf
    return float(last.time)


def compute_drag_delta_v(
    altitude: float,
    inclination: float,
    ballistic_coefficient: float,
    duration: float,
    *,
    activity: SolarActivity = "moderate",
    start: datetime | None = None,
) -> float:
    """
    Computes the velocity a satellite must add to hold a circular orbit at an
    altitude against drag for a duration: the integral of its drag
    acceleration, `rho * F * v ** 2 / (2 * B)`, over the duration, for its
    speed `v`, the density `rho` of the thermosphere (see
    `compute_atmospheric_density`), the factor `F` of the rotation of the
    atmosphere with the Earth, and the ballistic coefficient `B`.

    Args:
        altitude (float): The altitude (meters) above the mean radius.
        inclination (float): The inclination (degrees).
        ballistic_coefficient (float): The ballistic coefficient (kg/m^2):
            the mass over the drag coefficient times the drag area.
        duration (float): The duration (years).
        activity (SolarActivity): The solar activity: a reference level, an
            index from 0 (low) to 1 (high), or a solar cycle.
        start (datetime.datetime | None): The start, required for a solar
            cycle.

    Returns:
        float: The velocity (meters per second).

    Raises:
        ValueError: If the altitude is below 100 km, the ballistic
            coefficient is not positive, the duration is negative, or a solar
            cycle has no start.
    """
    if ballistic_coefficient <= 0:
        raise ValueError("the ballistic coefficient must be positive")
    if duration < 0:
        raise ValueError("the duration must not be negative")
    if altitude < MIN_DENSITY_ALTITUDE:
        raise ValueError("the altitude must be at least 100 km above the mean radius")
    if duration == 0 or altitude > MAX_DENSITY_ALTITUDE:
        return 0.0
    activity_at, longest = _get_activity_function(activity, start)
    count = 2 if math.isinf(longest) else max(2, math.ceil(duration / longest * 4) + 1)
    years = np.linspace(0, duration, count)
    density = interpolate_density(
        altitude, np.array([activity_at(year) for year in years])
    )
    # the trapezoidal rule
    mean_density = float(
        np.sum((density[1:] + density[:-1]) / 2 * np.diff(years)) / duration
    )
    radius = EARTH_MEAN_RADIUS + altitude
    speed = math.sqrt(EARTH_MU / radius)
    factor = float(
        _get_rotation_factor(np.array([radius]), np.array([speed]), inclination)[0]
    )
    acceleration = mean_density * factor * speed**2 / (2 * ballistic_coefficient)
    return acceleration * duration * YEAR_S
