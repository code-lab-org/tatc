"""
Atmospheric utility functions: the pressure of the lower atmosphere (the U.S.
Standard Atmosphere, 1976) and the density of the thermosphere (NRLMSIS 2.1),
which drags low orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import os
from datetime import datetime, timezone
from functools import cache
from typing import Literal, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
from numba import njit
from pydantic import BaseModel, ConfigDict, Field

from ..constants import resources_dir

# U.S. Standard Atmosphere (1976) layer definitions, valid from 0 to 86 km
# geopotential altitude: for each of the seven layers, the base altitude
# (m), base temperature (K), lapse rate (K/m), and base pressure (Pa).
_LAYER_BASE_ALTITUDE = (0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0)
_LAYER_BASE_TEMPERATURE = (288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65)
_LAYER_LAPSE_RATE = (-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002)
_LAYER_BASE_PRESSURE = (
    101325.0,
    22632.634393022883,
    5475.157308245898,
    868.0881576113082,
    110.91901929163447,
    66.9471149988061,
    3.957095749749069,
)

_STANDARD_GRAVITY = 9.80665  # m/s^2
_MOLAR_MASS_AIR = 0.0289644  # kg/mol
_GAS_CONSTANT = 8.3144598  # J/(mol*K)
_MAX_ALTITUDE = 86000.0  # m, upper bound of the standard atmosphere model
_MIN_PRESSURE = 0.3023942287612965  # Pa, pressure at _MAX_ALTITUDE (layer 6 formula)


@njit(cache=True)
def altitude_to_pressure(altitude: float) -> float:
    """
    Fast computation of atmospheric pressure at a specified altitude, using
    the U.S. Standard Atmosphere (1976) model. Inputs outside the model's
    valid range (0 to 86 km) saturate to the pressure at the nearest bound.

    Args:
        altitude (float): Altitude (meters) above mean sea level.

    Returns:
        float: Atmospheric pressure (Pa). Divide by 100 to convert to the
            hPa (millibar) units common in atmospheric science.
    """
    # as floats, since integers can overflow when compiled
    altitude = altitude * 1.0
    h = min(max(altitude, 0.0), _MAX_ALTITUDE)
    # select the highest layer whose base altitude is at or below h
    layer = 0
    for i in range(1, len(_LAYER_BASE_ALTITUDE)):
        if h >= _LAYER_BASE_ALTITUDE[i]:
            layer = i
    h_b = _LAYER_BASE_ALTITUDE[layer]
    t_b = _LAYER_BASE_TEMPERATURE[layer]
    l_b = _LAYER_LAPSE_RATE[layer]
    p_b = _LAYER_BASE_PRESSURE[layer]
    if l_b == 0.0:
        return p_b * np.exp(
            -_STANDARD_GRAVITY * _MOLAR_MASS_AIR * (h - h_b) / (_GAS_CONSTANT * t_b)
        )
    return p_b * (t_b / (t_b + l_b * (h - h_b))) ** (
        _STANDARD_GRAVITY * _MOLAR_MASS_AIR / (_GAS_CONSTANT * l_b)
    )


@njit(cache=True)
def pressure_to_altitude(pressure: float) -> float:
    """
    Fast computation of altitude at a specified atmospheric pressure, using
    the U.S. Standard Atmosphere (1976) model. This is the inverse of
    `altitude_to_pressure`; inputs outside the model's valid range (down to
    the pressure at 86 km) saturate to the altitude at the nearest bound.

    Args:
        pressure (float): Atmospheric pressure (Pa). Multiply hPa
            (millibar) values, common in atmospheric science, by 100 to
            convert to Pa.

    Returns:
        float: Altitude (meters) above mean sea level.
    """
    # as floats, since integers can overflow when compiled
    pressure = pressure * 1.0
    p = min(max(pressure, _MIN_PRESSURE), _LAYER_BASE_PRESSURE[0])
    # select the highest layer whose base pressure is at or above p
    # (pressure decreases monotonically with altitude/layer)
    layer = 0
    for i in range(1, len(_LAYER_BASE_PRESSURE)):
        if p <= _LAYER_BASE_PRESSURE[i]:
            layer = i
    h_b = _LAYER_BASE_ALTITUDE[layer]
    t_b = _LAYER_BASE_TEMPERATURE[layer]
    l_b = _LAYER_LAPSE_RATE[layer]
    p_b = _LAYER_BASE_PRESSURE[layer]
    if l_b == 0.0:
        return h_b - (_GAS_CONSTANT * t_b) / (
            _STANDARD_GRAVITY * _MOLAR_MASS_AIR
        ) * np.log(p / p_b)
    return h_b + (t_b / l_b) * (
        (p / p_b) ** (-_GAS_CONSTANT * l_b / (_STANDARD_GRAVITY * _MOLAR_MASS_AIR)) - 1
    )


_DENSITY_FILE = os.path.join(resources_dir, "atmosphere", "nrlmsis_density.csv")
"""
The global-mean mass density (kg/m^3) of the thermosphere of NRLMSIS 2.1
(Emmert et al., 2022, https://doi.org/10.1029/2022SW003219) at heights from
100 to 1500 km above the Earth's mean radius, every 5 km, averaged over the
sphere, local times, and seasons, for the long-term reference levels of
solar and geomagnetic activity of ECSS-E-ST-10-04C (Table 6-3); see
`validation/atmosphere/build_density_table.py` of the repository.
"""

MIN_DENSITY_ALTITUDE = 100e3
"""float: The lowest altitude (meters) of the density of the thermosphere
(see `compute_atmospheric_density`), at which an orbit is taken to reenter."""

MAX_DENSITY_ALTITUDE = 1500e3
"""float: The highest altitude (meters) of the density of the thermosphere
(see `compute_atmospheric_density`), above which it is taken to be zero."""

SOLAR_ACTIVITY_LEVELS = {"low": 0.0, "moderate": 0.5, "high": 1.0}
"""dict[str, float]: The index of each reference level of solar and
geomagnetic activity (see `compute_atmospheric_density`): the long-term
values of ECSS-E-ST-10-04C, Table 6-3, low (F10.7 = 65, Ap = 0), moderate
(F10.7 = 140, Ap = 15), and high (F10.7 = 250, Ap = 45)."""

YEAR_S = 365.25 * 86400
"""float: The seconds of a (Julian) year, the unit of the times of solar
cycles and orbit lifetimes."""


class SolarCycle(BaseModel):
    """
    A solar cycle, by which the activity of the Sun (and the density of the
    thermosphere) varies over time as a cosine between the low level of
    activity at its minima and the high level halfway between them: the index
    of activity (see `SolarActivity`) is
    `(1 - cos(2 * pi * (time - minimum) / period)) / 2`. Real cycles differ
    in length (9 to 14 years) and strength, and rise faster than they fall;
    the cosine is a mean cycle for predictions over years, such as orbit
    lifetimes. By default, its minimum is that of December 2019, which began
    solar cycle 25 (NASA, "Solar Cycle 25 Is Here", 15 September 2020,
    https://www.nasa.gov/news-release/solar-cycle-25-is-here-nasa-noaa-scientists-explain-what-that-means/).
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    minimum: datetime = Field(
        datetime(2019, 12, 15, tzinfo=timezone.utc),
        description="The time of a minimum of the cycle (UTC if without time zone).",
    )
    period: float = Field(
        11.0, gt=0, description="The period (years) of the cycle, minimum to minimum."
    )

    def get_activity(self, time: datetime | npt.ArrayLike) -> float | npt.NDArray:
        """
        Gets the index of activity at times: 0 at a minimum, 1 at a maximum.

        Args:
            time (datetime.datetime | numpy.typing.ArrayLike): The time(s)
                (UTC if without time zone).

        Returns:
            float | numpy.typing.NDArray: The index of activity.
        """
        return self.get_activity_after(get_years_between(self.minimum, time))

    def get_activity_after(self, years: float | npt.ArrayLike) -> float | npt.NDArray:
        """
        Gets the index of activity a number of years after a minimum.

        Args:
            years (float | numpy.typing.ArrayLike): The years after a minimum.

        Returns:
            float | numpy.typing.NDArray: The index of activity.
        """
        activity = (1 - np.cos(2 * np.pi * np.asarray(years) / self.period)) / 2
        return float(activity) if np.ndim(activity) == 0 else activity


SolarActivity = Union[Literal["low", "moderate", "high"], float, SolarCycle]
"""
The activity of the Sun that sets the density of the thermosphere: a
reference level (`"low"`, `"moderate"`, or `"high"`, see
`SOLAR_ACTIVITY_LEVELS`), an index of activity from 0 (low) to 1 (high), or a
solar cycle (`SolarCycle`) by which it varies over time.
"""


def get_years_between(
    start: datetime, time: datetime | npt.ArrayLike
) -> float | npt.NDArray:
    """
    Gets the (Julian) years from a time to other times.

    Args:
        start (datetime.datetime): The time from which to count (UTC if
            without time zone).
        time (datetime.datetime | numpy.typing.ArrayLike): The time(s) (UTC if
            without time zone).

    Returns:
        float | numpy.typing.NDArray: The years (negative before the start).
    """
    origin = pd.Timestamp(start)
    origin = origin.tz_localize("UTC") if origin.tzinfo is None else origin
    times = pd.to_datetime(time, utc=True)
    years = np.asarray((times - origin) / pd.Timedelta(seconds=1)) / YEAR_S
    return float(years) if np.ndim(years) == 0 else years


def get_solar_activity(
    activity: SolarActivity, time: datetime | npt.ArrayLike | None = None
) -> float | npt.NDArray:
    """
    Gets the index of solar activity, from 0 (low) to 1 (high), of a level, an
    index, or a solar cycle at times.

    Args:
        activity (SolarActivity): The activity (see `SolarActivity`).
        time (datetime.datetime | numpy.typing.ArrayLike | None): The time(s)
            of a solar cycle, ignored for a level or an index.

    Returns:
        float | numpy.typing.NDArray: The index of activity.

    Raises:
        ValueError: If a level is not known, an index is not between 0 and 1,
            or a solar cycle has no time.
    """
    if isinstance(activity, SolarCycle):
        if time is None:
            raise ValueError("the activity of a solar cycle requires a time")
        return activity.get_activity(time)
    if isinstance(activity, str):
        if activity not in SOLAR_ACTIVITY_LEVELS:
            raise ValueError(
                f"unknown level of solar activity {activity!r}; "
                + f"expected one of {list(SOLAR_ACTIVITY_LEVELS)}"
            )
        return SOLAR_ACTIVITY_LEVELS[activity]
    index = float(activity)
    if not 0 <= index <= 1:
        raise ValueError(f"the index of solar activity must be in [0, 1], not {index}")
    return index


@cache
def _load_densities() -> tuple[npt.NDArray, npt.NDArray]:
    """
    Loads the table of densities.

    Returns:
        tuple[numpy.typing.NDArray, numpy.typing.NDArray]: the altitudes
            (meters), and the natural logarithms of the densities (kg/m^3) at
            each altitude (rows) and the low, moderate, and high levels of
            activity (columns)
    """
    table = pd.read_csv(_DENSITY_FILE, comment="#")
    altitudes = table["altitude_km"].to_numpy(dtype=float) * 1e3
    log_densities = np.log(table[list(SOLAR_ACTIVITY_LEVELS)].to_numpy(dtype=float))
    return altitudes, log_densities


def interpolate_density(
    altitude: npt.ArrayLike, activity: float | npt.ArrayLike
) -> npt.NDArray:
    """
    Interpolates the density (kg/m^3) of the thermosphere at altitudes
    (meters) and indices of activity (broadcast together, see
    `compute_atmospheric_density`), log-linearly in both: NaN below the
    table, and zero above it. It is the vectorized core of
    `compute_atmospheric_density`, for indices already computed.

    Args:
        altitude (numpy.typing.ArrayLike): The altitude(s) (meters) above the
            mean radius.
        activity (float | numpy.typing.ArrayLike): The index (or indices) of
            activity, from 0 (low) to 1 (high).

    Returns:
        numpy.typing.NDArray: The density (kg/m^3).
    """
    altitudes, log_densities = _load_densities()
    altitude, activity = np.broadcast_arrays(
        np.asarray(altitude, dtype=float), np.asarray(activity, dtype=float)
    )
    levels = [
        np.interp(altitude, altitudes, log_densities[:, k])
        for k in range(log_densities.shape[1])
    ]
    # piecewise linear in the index of activity: low (0), moderate (0.5),
    # high (1)
    upper_half = activity > 0.5
    lower = np.where(upper_half, levels[1], levels[0])
    upper = np.where(upper_half, levels[2], levels[1])
    fraction = np.where(upper_half, activity - 0.5, activity) / 0.5
    density = np.exp(lower + fraction * (upper - lower))
    density = np.where(altitude > MAX_DENSITY_ALTITUDE, 0.0, density)
    return np.where(altitude < MIN_DENSITY_ALTITUDE, np.nan, density)


def compute_atmospheric_density(
    altitude: float | npt.ArrayLike,
    activity: SolarActivity = "moderate",
    time: datetime | npt.ArrayLike | None = None,
) -> float | npt.NDArray:
    """
    Computes the mass density of the thermosphere at altitudes above the
    Earth's mean radius (the mean altitude of an orbit, see
    `tatc.schemas.CircularOrbit`), by which drag decays low orbits: the mean
    of NRLMSIS 2.1 (Emmert et al., 2022,
    https://doi.org/10.1029/2022SW003219) over the sphere, local times, and
    seasons, which an orbit averages over its revolutions, interpolated
    log-linearly from a table every 5 km (errors at most 2%, at 125 km, and
    0.04% on average) and between the reference levels of activity (see
    `SolarActivity`). The local density varies about the mean by a factor
    of about 2 (from about 0.4 to 2 times the mean, 5th to 95th percentiles,
    at 400 km), and the mean by the solar activity (by about 24 times from
    low to high at 400 km).

    Args:
        altitude (float | numpy.typing.ArrayLike): The altitude(s) (meters)
            above the mean radius.
        activity (SolarActivity): The solar activity: a reference level, an
            index from 0 (low) to 1 (high), or a solar cycle.
        time (datetime.datetime | numpy.typing.ArrayLike | None): The time(s)
            of a solar cycle (broadcast with the altitudes), ignored for a
            level or an index.

    Returns:
        float | numpy.typing.NDArray: The density (kg/m^3): NaN below 100 km
            (see `MIN_DENSITY_ALTITUDE`), and zero above 1500 km (see
            `MAX_DENSITY_ALTITUDE`).
    """
    index = get_solar_activity(activity, time)
    density = interpolate_density(altitude, index)
    return float(density) if np.ndim(density) == 0 else density
