"""
Coarse radio link budgets: the rate at which a transmitter can send data to a
receiver over a distance, from free-space path loss, a clear-sky atmospheric
loss scaled from zenith by elevation, and the energy per bit to noise
density that the link's modulation and coding require (or the capacity of
its channel less a gap). Suited to concept studies with few known
parameters: rain fades, interference, and polarization are not modeled
(beyond the losses and margin given).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import functools
import os
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
import numpy.typing as npt
import yaml

from ..constants import resources_dir

if TYPE_CHECKING:  # pragma: no cover
    from ..schemas.communication import LinkMode, Receiver, ShannonModem, Transmitter

SPEED_OF_LIGHT = 299792458.0
"""The speed of light (meters per second) in vacuum."""

BOLTZMANN_CONSTANT = 10 * np.log10(1.380649e-23)
"""Boltzmann's constant (dBW per kelvin per hertz), about -228.6."""

_COMMUNICATIONS_DIR = os.path.join(resources_dir, "communications")


@dataclass(frozen=True)
class BandDefaults:
    """
    The nominal values of a frequency band of space communications.

    Attributes:
        band (str): The band.
        frequency (float): The nominal carrier frequency (hertz).
        zenith_attenuation (float): The nominal clear-sky attenuation (dB) of
            the atmosphere at zenith at a ground station.
        sources (dict[str, str]): The source of each value.
    """

    band: str
    frequency: float
    zenith_attenuation: float
    sources: dict[str, str]


@functools.cache
def _load_table(name: str) -> dict:
    """Loads a table of the communications resources."""
    with open(os.path.join(_COMMUNICATIONS_DIR, name), encoding="utf-8") as file:
        return yaml.safe_load(file)


def get_band(band: str) -> BandDefaults:
    """
    Gets the nominal values of a frequency band (`S`, `X`, `Ku`, or `Ka`),
    with their sources.

    Args:
        band (str): The band.

    Returns:
        BandDefaults: The nominal values.

    Raises:
        ValueError: If the band is not known.
    """
    bands = _load_table("bands.yml")["bands"]
    if band not in bands:
        raise ValueError(f"unknown band {band!r}; known bands are {sorted(bands)}")
    entry = bands[band]
    return BandDefaults(
        band=band,
        frequency=float(entry["frequency"]["value"]),
        zenith_attenuation=float(entry["zenith_attenuation"]["value"]),
        sources={
            name: entry[name]["source"] for name in ("frequency", "zenith_attenuation")
        },
    )


def get_link_modes(
    standard: Literal["dvb-s2", "ccsds"],
    symbol_rate: float,
    names: Sequence[str] | None = None,
    implementation_loss: float = 0,
) -> list[LinkMode]:
    """
    Gets the modes (modulations and codes) of a standard at a symbol rate:
    the rate of information of each is the symbol rate times its spectral
    efficiency (information bits per symbol), and the energy per bit to noise
    density it requires is its published threshold plus an implementation
    loss. The thresholds of DVB-S2 (ETSI EN 302 307-1) are ideal, for normal
    frames without pilots; those of CCSDS codes are at the error rates their
    sources state.

    Args:
        standard (Literal["dvb-s2", "ccsds"]): The standard.
        symbol_rate (float): The symbol rate (symbols per second).
        names (Sequence[str] | None): The names of the modes to get, or `None`
            for all of them.
        implementation_loss (float): The loss (dB) of a modem from the
            published thresholds.

    Returns:
        list[LinkMode]: The modes, by increasing rate.

    Raises:
        ValueError: If the standard or a name is not known.
    """
    from ..schemas.communication import (  # pylint: disable=import-outside-toplevel
        LinkMode,
    )

    standards = _load_table("modcods.yml")["standards"]
    if standard not in standards:
        raise ValueError(
            f"unknown standard {standard!r}; known standards are {sorted(standards)}"
        )
    table = {mode["name"]: mode for mode in standards[standard]["modes"]}
    if names is None:
        names = list(table)
    unknown = [name for name in names if name not in table]
    if unknown:
        raise ValueError(
            f"unknown modes {unknown} of {standard!r}; known modes are {list(table)}"
        )
    modes = []
    for name in names:
        mode = table[name]
        efficiency = float(mode["spectral_efficiency"])
        if "required_ebn0" in mode:
            ebn0 = float(mode["required_ebn0"])
        else:
            ebn0 = float(mode["required_esn0"]) - 10 * np.log10(efficiency)
        modes.append(
            LinkMode(
                name=name,
                data_rate=symbol_rate * efficiency,
                required_ebn0=ebn0 + implementation_loss,
            )
        )
    return sorted(modes, key=lambda mode: mode.data_rate)


def compute_free_space_path_loss(
    distance: npt.ArrayLike, frequency: float
) -> float | npt.NDArray:
    """
    Computes the free-space path loss, `20 log10(4 pi d f / c)`.

    Args:
        distance (numpy.typing.ArrayLike): The distance (meters).
        frequency (float): The frequency (hertz).

    Returns:
        float | numpy.typing.NDArray: The loss (dB).
    """
    loss = 20 * np.log10(
        4 * np.pi * np.asarray(distance, dtype=float) * frequency / SPEED_OF_LIGHT
    )
    return float(loss) if np.ndim(loss) == 0 else loss


def compute_antenna_gain(
    diameter: float, frequency: float, efficiency: float = 0.55
) -> float:
    """
    Computes the gain of a parabolic antenna, `efficiency (pi D f / c) ** 2`.

    Args:
        diameter (float): The diameter (meters).
        frequency (float): The frequency (hertz).
        efficiency (float): The aperture efficiency.

    Returns:
        float: The gain (dBi).
    """
    return float(
        10 * np.log10(efficiency * (np.pi * diameter * frequency / SPEED_OF_LIGHT) ** 2)
    )


def compute_atmospheric_loss(
    elevation: npt.ArrayLike, zenith_attenuation: float, min_elevation: float = 5
) -> float | npt.NDArray:
    """
    Computes the loss through the atmosphere at an elevation angle from its
    attenuation at zenith, scaled by the length of the path through a flat
    layer (`zenith_attenuation / sin(elevation)`), which holds above low
    elevations: lower ones take the loss at `min_elevation`.

    Args:
        elevation (numpy.typing.ArrayLike): The elevation angle (degrees).
        zenith_attenuation (float): The attenuation (dB) at zenith.
        min_elevation (float): The lowest elevation angle (degrees) of the
            scaling.

    Returns:
        float | numpy.typing.NDArray: The loss (dB).
    """
    angle = np.maximum(np.asarray(elevation, dtype=float), min_elevation)
    loss = zenith_attenuation / np.sin(np.radians(angle))
    return float(loss) if np.ndim(loss) == 0 else loss


def compute_carrier_to_noise_density(
    eirp: float,
    path_loss: npt.ArrayLike,
    g_over_t: float,
    losses: npt.ArrayLike = 0,
) -> float | npt.NDArray:
    """
    Computes the carrier to noise density of a link,
    `EIRP - path loss - losses + G/T - 10 log10(k)`.

    Args:
        eirp (float): The transmitter's effective isotropic radiated power (dBW).
        path_loss (numpy.typing.ArrayLike): The path loss (dB).
        g_over_t (float): The receiver's figure of merit (dB/K).
        losses (numpy.typing.ArrayLike): Other losses (dB).

    Returns:
        float | numpy.typing.NDArray: The carrier to noise density (dB-Hz).
    """
    cn0 = (
        eirp
        - np.asarray(path_loss, dtype=float)
        - np.asarray(losses, dtype=float)
        + g_over_t
        - BOLTZMANN_CONSTANT
    )
    return float(cn0) if np.ndim(cn0) == 0 else cn0


def compute_mode_rate(
    carrier_to_noise_density: npt.ArrayLike,
    modes: Sequence[LinkMode],
    margin: float = 0,
) -> float | npt.NDArray:
    """
    Computes the rate of the fastest mode that closes a link: whose required
    energy per bit to noise density plus the margin is at most the carrier to
    noise density less the rate (`Eb/N0 = C/N0 - 10 log10(R)`).

    Args:
        carrier_to_noise_density (numpy.typing.ArrayLike): The carrier to noise
            density (dB-Hz).
        modes (Sequence[LinkMode]): The modes.
        margin (float): The margin (dB).

    Returns:
        float | numpy.typing.NDArray: The rate (bits per second), 0 where no
            mode closes.
    """
    cn0 = np.asarray(carrier_to_noise_density, dtype=float)
    rate = np.zeros_like(cn0)
    for mode in modes:
        closes = cn0 - 10 * np.log10(mode.data_rate) >= mode.required_ebn0 + margin
        rate = np.where(closes, np.maximum(rate, mode.data_rate), rate)
    return float(rate) if np.ndim(rate) == 0 else rate


def compute_shannon_rate(
    carrier_to_noise_density: npt.ArrayLike,
    modem: ShannonModem,
    margin: float = 0,
) -> float | npt.NDArray:
    """
    Computes the rate of a modem that follows the capacity of its channel
    less a gap, `B log2(1 + S/N / gap)` with the signal to noise ratio
    `S/N = C/N0 - margin - 10 log10(B)`, at most its maximum spectral
    efficiency times its bandwidth `B`.

    Args:
        carrier_to_noise_density (numpy.typing.ArrayLike): The carrier to noise
            density (dB-Hz).
        modem (ShannonModem): The modem.
        margin (float): The margin (dB).

    Returns:
        float | numpy.typing.NDArray: The rate (bits per second).
    """
    snr = (
        np.asarray(carrier_to_noise_density, dtype=float)
        - margin
        - 10 * np.log10(modem.bandwidth)
        - modem.gap
    )
    rate = modem.bandwidth * np.log2(1 + 10 ** (snr / 10))
    if modem.max_spectral_efficiency is not None:
        rate = np.minimum(rate, modem.bandwidth * modem.max_spectral_efficiency)
    return float(rate) if np.ndim(rate) == 0 else rate


def get_frequency(transmitter: Transmitter) -> float:
    """
    Gets the carrier frequency of a transmitter: its `frequency`, or the
    nominal frequency of its band.

    Args:
        transmitter (Transmitter): The transmitter.

    Returns:
        float: The frequency (hertz).
    """
    if transmitter.frequency is not None:
        return transmitter.frequency
    return get_band(transmitter.band).frequency


def get_eirp(transmitter: Transmitter) -> float:
    """
    Gets the effective isotropic radiated power of a transmitter: its `eirp`,
    or its power into its antenna's gain (given, or of a parabolic antenna of
    its diameter at its frequency).

    Args:
        transmitter (Transmitter): The transmitter.

    Returns:
        float: The EIRP (dBW).

    Raises:
        ValueError: If neither the EIRP nor the power is given.
    """
    if transmitter.eirp is not None:
        return transmitter.eirp
    if transmitter.power is None:
        raise ValueError("the transmitter gives neither its eirp nor its power")
    gain = transmitter.antenna_gain
    if gain is None:
        gain = compute_antenna_gain(
            transmitter.antenna_diameter,
            get_frequency(transmitter),
            transmitter.antenna_efficiency,
        )
    return float(10 * np.log10(transmitter.power) + gain)


def get_g_over_t(receiver: Receiver, frequency: float) -> float:
    """
    Gets the figure of merit of a receiver: its `g_over_t`, or the gain of a
    parabolic antenna of its diameter at a frequency over its system noise
    temperature.

    Args:
        receiver (Receiver): The receiver.
        frequency (float): The frequency (hertz).

    Returns:
        float: The G/T (dB/K).

    Raises:
        ValueError: If neither the G/T nor the antenna is given.
    """
    if receiver.g_over_t is not None:
        return receiver.g_over_t
    if receiver.antenna_diameter is None:
        raise ValueError("the receiver gives neither its g_over_t nor its antenna")
    gain = compute_antenna_gain(
        receiver.antenna_diameter, frequency, receiver.antenna_efficiency
    )
    return float(gain - 10 * np.log10(receiver.system_noise_temperature))


def compute_link_rate(
    transmitter: Transmitter,
    receiver: Receiver | None,
    distance: npt.ArrayLike,
    elevation: npt.ArrayLike | None = None,
) -> float | npt.NDArray:
    """
    Computes the rate at which a transmitter sends data to a receiver. If
    both give what a link budget needs (see `Transmitter.has_budget` and
    `Receiver.has_budget`), the rate is that of the fastest mode that closes
    the link with the transmitter's margin (or of its modem), from the
    carrier to noise density over the distance, with the atmospheric loss at
    the elevation angle for a ground station's receiver. Otherwise, it is the
    highest rate the transmitter sends (see `Transmitter.get_max_data_rate`).
    Either way, it is at most the receiver's `max_data_rate`.

    Args:
        transmitter (Transmitter): The transmitter.
        receiver (Receiver | None): The receiver, or `None` for one that is
            not known.
        distance (numpy.typing.ArrayLike): The distance (meters).
        elevation (numpy.typing.ArrayLike | None): The elevation angle
            (degrees) of the transmitter seen from a ground station, or `None`
            for a link between satellites (without an atmosphere).

    Returns:
        float | numpy.typing.NDArray: The rate (bits per second), 0 where the
            link does not close.
    """
    distance = np.asarray(distance, dtype=float)
    if transmitter.has_budget() and receiver is not None and receiver.has_budget():
        frequency = get_frequency(transmitter)
        losses = np.asarray(
            compute_free_space_path_loss(distance, frequency), dtype=float
        ) + (transmitter.losses + receiver.losses)
        if elevation is not None:
            zenith = receiver.zenith_attenuation
            if zenith is None:
                zenith = get_band(receiver.band).zenith_attenuation
            losses = losses + compute_atmospheric_loss(elevation, zenith)
        cn0 = compute_carrier_to_noise_density(
            get_eirp(transmitter), losses, get_g_over_t(receiver, frequency)
        )
        if transmitter.modes:
            rate = compute_mode_rate(cn0, transmitter.modes, transmitter.margin)
        else:
            rate = compute_shannon_rate(cn0, transmitter.modem, transmitter.margin)
        rate = np.minimum(rate, transmitter.get_max_data_rate())
    else:
        rate = np.full(distance.shape, transmitter.get_max_data_rate())
    if receiver is not None and receiver.max_data_rate is not None:
        rate = np.minimum(rate, receiver.max_data_rate)
    rate = np.asarray(rate, dtype=float)
    return float(rate) if np.ndim(rate) == 0 else rate
