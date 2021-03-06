# -*- coding: utf-8 -*-
"""
Object schemas for satellites.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""
from __future__ import annotations

import copy
from datetime import datetime, timedelta, timezone
from enum import Enum
import math
import numpy as np
from pydantic import BaseModel, Field, root_validator
from sgp4.api import Satrec, WGS72
from sgp4.conveniences import sat_epoch_datetime
from sgp4 import exporter
from typing import Optional, List, Union
from typing_extensions import Literal

from .instrument import Instrument
from .orbit import TwoLineElements, CircularOrbit, SunSynchronousOrbit, KeplerianOrbit


class SpaceSystem(BaseModel):
    """
    Base class for space systems.
    """

    name: str = Field(
        ...,
        description="Space system name.",
        example="International Space Station",
    )
    orbit: Union[
        TwoLineElements, CircularOrbit, SunSynchronousOrbit, KeplerianOrbit
    ] = Field(..., description="Orbit specification.")
    instruments: List[Instrument] = Field(
        [], description="List of assigned instruments."
    )


class Satellite(SpaceSystem):
    """
    Single satellite.
    """

    type: Literal["satellite"] = Field(
        "satellite", description="Space system type discriminator."
    )

    def generate_members(self) -> List[Satellite]:
        """
        Generate space system member satellites (returns a list containing this satellite).

        Returns:
            List[Satellite]: the member satellites
        """
        return [self]


class TrainConstellation(Satellite):
    """
    A constellation that arranges member satellites in sequence.
    """

    type: Literal["train"] = Field(
        "train", description="Space system type discriminator."
    )
    orbit: Union[
        TwoLineElements, SunSynchronousOrbit, CircularOrbit, KeplerianOrbit
    ] = Field(..., description="Lead orbit for this constellation.")
    number_satellites: int = Field(
        1, description="The count of the number of satellites.", ge=1
    )
    interval: timedelta = Field(
        ...,
        description="The local time interval between satellites in a train constellation.",
    )
    repeat_ground_track: bool = Field(
        True,
        description="True, if the train satellites should repeat the same ground track.",
    )

    def get_delta_mean_anomaly(self) -> float:
        """
        Gets the difference in mean anomaly (decimal degrees) for adjacent
        member satellites.

        Returns:
            float: the difference in mean anomaly
        """
        return 360 * self.orbit.get_mean_motion() * (self.interval / timedelta(days=1))

    def get_delta_raan(self) -> float:
        """
        Gets the difference in right ascension of ascending node (decimal
        degrees) for adjacent member satellites.

        Returns:
            float: the difference in right ascension of ascending node
        """
        if self.repeat_ground_track:
            return -1 * 360 * (self.interval / timedelta(days=1))
        else:
            return 0

    def generate_members(self) -> List[Satellite]:
        """
        Generate space system member satellites.

        Returns:
            List[Satellite]: the member satellites
        """
        return [
            Satellite(
                name=f"{self.name} #{i+1:02d}",
                orbit=self.orbit.get_derived_orbit(
                    i * self.get_delta_mean_anomaly(), i * self.get_delta_raan()
                ),
                instruments=copy.deepcopy(self.instruments),
            )
            for i in range(self.number_satellites)
        ]


class WalkerConfiguration(str, Enum):
    delta = "delta"
    star = "star"


class WalkerConstellation(Satellite):
    """
    A constellation that arranges member satellites following the Walker pattern.
    """

    type: Literal["walker"] = Field(
        "walker", description="Space system type discriminator."
    )
    configuration: WalkerConfiguration = Field(
        WalkerConfiguration.delta, description="Walker configuration."
    )
    orbit: Union[
        TwoLineElements, SunSynchronousOrbit, CircularOrbit, KeplerianOrbit
    ] = Field(..., description="Lead orbit for this constellation.")
    number_satellites: int = Field(
        1, description="Number of satellites in the constellation.", ge=1
    )
    number_planes: int = Field(
        1,
        description="The number of equally-spaced planes in a Walker Delta constellation. Ranges from 1 to (number of satellites).",
        ge=1,
    )
    relative_spacing: int = Field(
        0,
        description="Relative spacing of satellites between plans for a Walker Delta constellation. Ranges from 0 for equal true anomaly to (number of planes) - 1. For example, `relative_spacing=1` means the true anomaly is shifted by `360/number_satellites` between adjacent planes.",
        ge=0,
    )

    @root_validator
    def number_planes_le_number_satellites(cls, values):
        """
        Validates the number of planes given the number of satellites.
        """
        p, t = values.get("number_planes"), values.get("number_satellites")
        if p is not None and t is not None and p > t:
            raise ValueError("number planes exceeds number satellites")
        return values

    @root_validator
    def relative_spacing_lt_number_planes(cls, values):
        """
        Validates the relative spacing given the number of planes.
        """
        p, f = values.get("number_planes"), values.get("relative_spacing")
        if p is not None and f is not None and f >= p:
            raise ValueError("relative spacing exceeds number planes - 1")
        return values

    def get_satellites_per_plane(self) -> int:
        """
        Gets the (max) number of satellites per plane.

        Returns:
            int: number of satellites per plane
        """
        return math.ceil(self.number_satellites / self.number_planes)

    def get_delta_mean_anomaly_within_planes(self) -> float:
        """
        Gets the difference in mean anomaly (decimal degrees) for adjacent
        member satellites within a single plane.

        Returns:
            float: difference in mean anomaly
        """
        return 360 / self.get_satellites_per_plane()

    def get_delta_mean_anomaly_between_planes(self) -> float:
        """
        Gets the difference in mean anomaly (decimal degrees) for adjacent
        member satellites between adjacent planes.

        Returns:
            float: difference in mean anomaly
        """
        return 360 * self.relative_spacing / self.number_satellites

    def get_delta_raan_between_planes(self) -> float:
        """
        Gets the difference in right ascension of ascending node (decimal
        degrees) for adjacent planes of member satellites.

        Returns:
            float: difference in right ascension of ascending node
        """
        if self.configuration == WalkerConfiguration.delta:
            return 360 / self.number_planes
        else:
            return 180 / self.number_planes

    def generate_members(self) -> List[Satellite]:
        """
        Generate space system member satellites.

        Returns:
            List[Satellite]: the member satellites
        """
        return [
            Satellite(
                name=f"{self.name} #{i+1}",
                orbit=self.orbit.get_derived_orbit(
                    np.mod(i, self.get_satellites_per_plane())
                    * self.get_delta_mean_anomaly_within_planes()
                    + (i // self.get_satellites_per_plane())
                    * self.get_delta_mean_anomaly_between_planes(),
                    (i // self.get_satellites_per_plane())
                    * self.get_delta_raan_between_planes(),
                ),
                instruments=copy.deepcopy(self.instruments),
            )
            for i in range(self.number_satellites)
        ]
