"""
Atmospheric utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
from numba import njit

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


@njit
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


@njit
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
        return h_b - (_GAS_CONSTANT * t_b) / (_STANDARD_GRAVITY * _MOLAR_MASS_AIR) * np.log(
            p / p_b
        )
    return h_b + (t_b / l_b) * (
        (p / p_b) ** (-_GAS_CONSTANT * l_b / (_STANDARD_GRAVITY * _MOLAR_MASS_AIR)) - 1
    )
