"""
Object schemas for the communications of spacecraft, ground stations, and
relays: coarse radio transmitters and receivers from which link budgets
estimate the rates at which data flow, the generation of data onboard, and
the options of data flow analyses (see `tatc.analysis.compute_data_flow`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, PositiveInt, model_validator

Band = Literal["S", "X", "Ku", "Ka"]
"""A radio frequency band of space communications (see
`tatc.utils.get_band` for its nominal frequency and atmospheric
attenuation)."""


class LinkMode(BaseModel):
    """
    One modulation and coding of a link: the rate of information it carries
    and the energy per bit to noise density (Eb/N0) it requires at the
    target error rate. A transmitter with several modes adapts its rate to
    the link (adaptive coding and modulation); one with one mode has a fixed
    rate.
    """

    model_config = ConfigDict(extra="forbid")

    name: str | None = Field(
        None, description="The name of the mode, if any.", examples=["QPSK 3/4"]
    )
    data_rate: float = Field(
        ...,
        gt=0,
        description="The rate (bits per second) of information.",
        examples=[150e6],
    )
    required_ebn0: float = Field(
        ...,
        description="The energy per bit to noise density (dB) required at "
        + "the target error rate.",
        examples=[4.5],
    )


class ShannonModem(BaseModel):
    """
    A modem whose rate follows the capacity of its channel less a gap, for
    when the modes of a link are unknown: the rate is
    `bandwidth * log2(1 + snr / gap)` (with the signal-to-noise ratio `snr`
    and `gap` as ratios), at most `max_spectral_efficiency * bandwidth`.
    """

    model_config = ConfigDict(extra="forbid")

    bandwidth: float = Field(
        ...,
        gt=0,
        description="The bandwidth (hertz) of the channel.",
        examples=[375e6],
    )
    gap: float = Field(
        3,
        ge=0,
        description="The gap (dB) between the rate and the capacity of the "
        + "channel, for the modulations and codes of the modem.",
    )
    max_spectral_efficiency: float | None = Field(
        None,
        gt=0,
        description="The maximum spectral efficiency (bits per second per "
        + "hertz) of the modem's modes, or `None` for none.",
        examples=[4.45],
    )


class Transmitter(BaseModel):
    """
    A radio transmitter of a spacecraft, which sends data down to ground
    stations (a `downlink`) or across to the satellites of a relay (a
    `crosslink`). Its effective isotropic radiated power (EIRP) is given
    directly (`eirp`) or as the `power` into an antenna of `antenna_gain` (or
    of `antenna_diameter`, a parabolic antenna). Without either, or without a
    receiver to close the link, its rate is the highest it can send (the
    `max_data_rate`, or the rate of its fastest mode) whenever it is in
    contact.
    """

    model_config = ConfigDict(extra="forbid")

    type: Literal["radio"] = Field("radio", description="The type of transmitter.")
    name: str | None = Field(
        None,
        description="The name of the transmitter, or `None` to name it by its "
        + "band and link.",
        examples=["X-band downlink"],
    )
    link: Literal["downlink", "crosslink"] = Field(
        "downlink",
        description="The link of the transmitter: down to ground stations or "
        + "across to the satellites of a relay.",
    )
    band: Band = Field(..., description="The frequency band.", examples=["X"])
    frequency: float | None = Field(
        None,
        gt=0,
        description="The carrier frequency (hertz), or `None` for the nominal "
        + "frequency of the band.",
        examples=[8.2e9],
    )
    eirp: float | None = Field(
        None,
        description="The effective isotropic radiated power (dBW), if given "
        + "rather than the power and antenna.",
        examples=[20],
    )
    power: float | None = Field(
        None,
        gt=0,
        description="The radio frequency power (watts) into the antenna.",
        examples=[5],
    )
    antenna_gain: float | None = Field(
        None,
        description="The gain (dBi) of the antenna.",
        examples=[15],
    )
    antenna_diameter: float | None = Field(
        None,
        gt=0,
        description="The diameter (meters) of a parabolic antenna.",
        examples=[0.3],
    )
    antenna_efficiency: float = Field(
        0.55,
        gt=0,
        le=1,
        description="The aperture efficiency of a parabolic antenna.",
    )
    losses: float = Field(
        0,
        ge=0,
        description="The losses (dB) of the link not counted in the EIRP or the "
        + "receiver's G/T, such as pointing and polarization losses.",
        examples=[1],
    )
    margin: float = Field(
        3,
        ge=0,
        description="The margin (dB) that the link must keep above the "
        + "energy per bit to noise density its mode requires.",
    )
    modes: list[LinkMode] = Field(
        [],
        description="The modulations and codes of the link (one for a fixed "
        + "rate; the fastest that closes is used).",
    )
    modem: ShannonModem | None = Field(
        None,
        description="A modem whose rate follows the capacity of the channel, "
        + "for when the modes are unknown.",
    )
    max_data_rate: float | None = Field(
        None,
        gt=0,
        description="The maximum rate (bits per second) of information.",
        examples=[150e6],
    )

    @model_validator(mode="after")
    def _check_power(self) -> Transmitter:
        antenna = (self.antenna_gain is not None, self.antenna_diameter is not None)
        if self.eirp is not None:
            if self.power is not None or any(antenna):
                raise ValueError(
                    "give either the eirp or the power and antenna, not both"
                )
        elif self.power is not None or any(antenna):
            if self.power is None or sum(antenna) != 1:
                raise ValueError(
                    "the power requires exactly one of antenna_gain or "
                    + "antenna_diameter"
                )
        if not self.modes and self.modem is None and self.max_data_rate is None:
            raise ValueError("give at least one of modes, modem, or max_data_rate")
        return self

    def get_name(self) -> str:
        """
        Gets the name of the transmitter: its `name`, or its band and link.

        Returns:
            str: The name.
        """
        return self.name if self.name is not None else f"{self.band}-band {self.link}"

    def get_max_data_rate(self) -> float:
        """
        Gets the highest rate (bits per second) the transmitter can send: the
        least of its `max_data_rate`, its fastest mode, and its modem's
        maximum, of those given.

        Returns:
            float: The rate (bits per second), infinite if only a modem without
                a maximum spectral efficiency limits it.
        """
        limits = []
        if self.max_data_rate is not None:
            limits.append(self.max_data_rate)
        if self.modes:
            limits.append(max(mode.data_rate for mode in self.modes))
        if self.modem is not None and self.modem.max_spectral_efficiency is not None:
            limits.append(self.modem.bandwidth * self.modem.max_spectral_efficiency)
        return min(limits) if limits else float("inf")

    def has_budget(self) -> bool:
        """
        Checks whether the transmitter gives what a link budget needs: its
        EIRP (or power and antenna) and its modes or modem.

        Returns:
            bool: `True` if a link budget can estimate its rate
        """
        return (self.eirp is not None or self.power is not None) and bool(
            self.modes or self.modem is not None
        )


class Receiver(BaseModel):
    """
    A radio receiver of a ground station or relay satellite in one band. Its
    figure of merit, the ratio of its antenna's gain to its system noise
    temperature (G/T), is given directly (`g_over_t`) or as the
    `antenna_diameter` of a parabolic antenna and its
    `system_noise_temperature`. Without either, a link to it is not budgeted
    (its rate is the transmitter's highest).
    """

    model_config = ConfigDict(extra="forbid")

    type: Literal["radio"] = Field("radio", description="The type of receiver.")
    band: Band = Field(..., description="The frequency band.", examples=["X"])
    g_over_t: float | None = Field(
        None,
        description="The figure of merit (dB/K): the antenna's gain over the "
        + "system noise temperature, if given rather than the antenna and "
        + "temperature.",
        examples=[31.0],
    )
    antenna_diameter: float | None = Field(
        None,
        gt=0,
        description="The diameter (meters) of a parabolic antenna.",
        examples=[11.3],
    )
    antenna_efficiency: float = Field(
        0.55,
        gt=0,
        le=1,
        description="The aperture efficiency of a parabolic antenna.",
    )
    system_noise_temperature: float | None = Field(
        None,
        gt=0,
        description="The system noise temperature (kelvins).",
        examples=[150],
    )
    losses: float = Field(
        0,
        ge=0,
        description="The losses (dB) of the receiver not counted in its G/T.",
    )
    zenith_attenuation: float | None = Field(
        None,
        ge=0,
        description="The attenuation (dB) of the atmosphere at zenith at a "
        + "ground station, or `None` for the nominal attenuation of the band "
        + "(not applied between satellites).",
        examples=[0.2],
    )
    max_data_rate: float | None = Field(
        None,
        gt=0,
        description="The maximum rate (bits per second) the receiver (or the "
        + "service it provides) accepts.",
        examples=[300e6],
    )

    @model_validator(mode="after")
    def _check_merit(self) -> Receiver:
        antenna = (
            self.antenna_diameter is not None,
            self.system_noise_temperature is not None,
        )
        if self.g_over_t is not None and any(antenna):
            raise ValueError(
                "give either g_over_t or the antenna and noise temperature, not both"
            )
        if self.g_over_t is None and any(antenna) and not all(antenna):
            raise ValueError(
                "the antenna_diameter requires the system_noise_temperature "
                + "(and vice versa)"
            )
        return self

    def has_budget(self) -> bool:
        """
        Checks whether the receiver gives what a link budget needs: its G/T
        (or antenna and noise temperature).

        Returns:
            bool: `True` if a link budget can use it
        """
        return self.g_over_t is not None or self.antenna_diameter is not None


AllTransmitters = Annotated[Transmitter, Field(discriminator="type")]
"""Transmitters of a spacecraft (only radio ones for now)."""

AllReceivers = Annotated[Receiver, Field(discriminator="type")]
"""Receivers of a ground station or relay (only radio ones for now)."""


class DataGeneration(BaseModel):
    """
    How a spacecraft generates data: at a `rate` (by default, the total data
    rate of its instruments; see `SpaceSystem.get_payload_resources`) for a
    fraction `duty_cycle` of the time, while it is active (`always`, while
    the spacecraft is `sunlit`, or while the point below it is in daylight,
    `day`).
    """

    model_config = ConfigDict(extra="forbid")

    rate: float | None = Field(
        None,
        ge=0,
        description="The rate (bits per second) of generation while active, or "
        + "`None` for the total data rate of the instruments.",
        examples=[10e6],
    )
    duty_cycle: float = Field(
        1,
        gt=0,
        le=1,
        description="The fraction of the active time in which data is "
        + "generated (spread evenly).",
    )
    activity: Literal["always", "sunlit", "day"] = Field(
        "always",
        description="When the spacecraft is active: always, while it is "
        + "sunlit, or while the point below it is in daylight.",
    )


class SpacecraftCommunications(BaseModel):
    """
    The communications of a spacecraft: its transmitters, and how it
    generates the data they send. Its storage is its resources' `storage`
    (see `SpacecraftResources`), unlimited if not given.
    """

    model_config = ConfigDict(extra="forbid")

    transmitters: list[AllTransmitters] = Field(
        [], description="The transmitters (each sends over one link at a time)."
    )
    generation: DataGeneration = Field(
        DataGeneration(), description="How the spacecraft generates data."
    )


class StationCommunications(BaseModel):
    """
    The communications of a ground station: the number of antennas (each
    receiving from one satellite at a time) and the receivers in each band.
    Without receivers, a station receives in every band without a link
    budget.
    """

    model_config = ConfigDict(extra="forbid")

    antennas: PositiveInt = Field(
        1,
        description="The number of antennas, each receiving from one satellite "
        + "at a time.",
    )
    receivers: list[AllReceivers] = Field(
        [],
        description="The receivers (at most one per band), or none to receive "
        + "in every band without a link budget.",
    )

    @model_validator(mode="after")
    def _check_bands(self) -> StationCommunications:
        bands = [receiver.band for receiver in self.receivers]
        if len(bands) != len(set(bands)):
            raise ValueError("give at most one receiver per band")
        return self


class RelayCommunications(BaseModel):
    """
    The communications of a relay: the number of single-access antennas of
    each relay satellite (each receiving from one satellite at a time), the
    receivers of the return link from users in each band, and the rate at
    which each relay satellite forwards data to the ground (its feeder link),
    shared by its users. The relay stores nothing: data flows through it as
    it is received.
    """

    model_config = ConfigDict(extra="forbid")

    antennas: PositiveInt = Field(
        1,
        description="The number of antennas of each relay satellite, each "
        + "receiving from one satellite at a time.",
    )
    receivers: list[AllReceivers] = Field(
        [],
        description="The receivers of each relay satellite (at most one per "
        + "band), or none to receive in every band without a link budget.",
    )
    feeder_rate: float | None = Field(
        None,
        gt=0,
        description="The rate (bits per second) at which each relay satellite "
        + "forwards data to the ground, shared by its users, or `None` for no "
        + "limit.",
        examples=[600e6],
    )

    @model_validator(mode="after")
    def _check_bands(self) -> RelayCommunications:
        bands = [receiver.band for receiver in self.receivers]
        if len(bands) != len(set(bands)):
            raise ValueError("give at most one receiver per band")
        return self


class DataFlow(BaseModel):
    """
    Options of a data flow analysis (see `tatc.analysis.compute_data_flow`):
    the step of the grid of time on which links are scheduled and storage is
    tracked (the accuracy of times), the backlog of data stored on each
    spacecraft at the start, and when links are released.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    step: timedelta = Field(
        timedelta(seconds=30),
        gt=timedelta(0),
        description="The step of the grid of time.",
    )
    initial_backlog: float = Field(
        0,
        ge=0,
        description="The volume (bits) of data stored on each spacecraft at "
        + "the start.",
    )
    release: Literal["empty", "end"] = Field(
        "empty",
        description="When a link is released: once the spacecraft's storage "
        + "empties (then requested again once more than `min_backlog` is "
        + "stored), or at the end of its contact (as for a scheduled pass).",
    )
    min_backlog: float = Field(
        0,
        ge=0,
        description="The volume (bits) stored above which a spacecraft "
        + 'requests a link, with `release="empty"`.',
    )
    strict: bool = Field(
        False,
        description="Whether to raise an error for a satellite whose data "
        + "flow cannot be modeled (without a rate of generation or of any "
        + "link), rather than to warn and find its latencies by timing alone.",
    )
