# -*- coding: utf-8 -*-
"""
Object schemas for satellites.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""
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
    Representation of a space system.
    """

    type: str
    name: str = Field(
        ...,
        description="Name of this satellite.",
        example="International Space Station",
    )


class Satellite(SpaceSystem):
    """
    Representation of a satellite in orbit.
    """

    type: Literal["satellite"] = Field("satellite")
    orbit: Union[
        TwoLineElements, CircularOrbit, SunSynchronousOrbit, KeplerianOrbit
    ] = Field(..., description="Orbit specification.")
    instruments: List[Instrument] = Field(
        [], description="List of assigned instruments."
    )

    def generate_members(self) -> List[SpaceSystem]:
        return [self]


class TrainConstellation(Satellite):
    """
    A constellation that arranges member satellites sequentially in one orbit plane.
    """

    type: Literal["train"] = Field("train")
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

    def generate_members(self) -> List[Satellite]:
        members = []
        lead_orbit = self.orbit.to_tle()
        lead_tle = Satrec.twoline2rv(lead_orbit.tle[0], lead_orbit.tle[1])
        epoch = sat_epoch_datetime(lead_tle)
        for satellite in range(self.number_satellites):
            delta_mean_anomaly = (
                lead_tle.no_kozai * satellite * (self.interval / timedelta(minutes=1))
            )
            if self.repeat_ground_track:
                delta_raan = 2 * np.pi * satellite * (self.interval / timedelta(days=1))
            else:
                delta_raan = 0
            satrec = Satrec()
            satrec.sgp4init(
                WGS72,
                "i",
                0,
                (epoch - datetime(1949, 12, 31, tzinfo=timezone.utc))
                / timedelta(days=1),
                0,
                0.0,
                0.0,
                lead_tle.ecco,
                lead_tle.argpo,
                lead_tle.inclo,
                np.mod(lead_tle.mo + delta_mean_anomaly, 2 * np.pi),
                lead_tle.no_kozai,
                np.mod(lead_tle.nodeo - delta_raan, 2 * np.pi),
            )
            tle1, tle2 = exporter.export_tle(satrec)
            members.append(
                Satellite(
                    name=f"{self.name} #{satellite+1:02d}",
                    orbit=TwoLineElements(tle=[tle1.replace("\x00", "U"), tle2]),
                    instruments=copy.deepcopy(self.instruments),
                )
            )
        return members


class WalkerConfiguration(str, Enum):
    delta = "delta"
    star = "star"


class WalkerConstellation(Satellite):
    """
    A constellation that arranges member satellites following the Walker Delta pattern.
    """

    type: Literal["walker"] = Field("walker")
    configuration: WalkerConfiguration = Field(
        WalkerConfiguration.delta, description="Walker constellation configuration."
    )
    orbit: Union[
        TwoLineElements, SunSynchronousOrbit, CircularOrbit, KeplerianOrbit
    ] = Field(..., description="Lead orbit for this constellation.")
    number_satellites: int = Field(
        1, description="The count of the number of satellites.", ge=1
    )
    number_planes: int = Field(
        1,
        description="The number of equally-spaced planes in a Walker Delta constellation. Ranges from 1 to (number of satellites).",
        ge=1,
    )
    relative_spacing: int = Field(
        0,
        description="Relative spacing of satellites between plans for a Walker Delta constellation. Ranges from 0 for equal true anomaly to (number of planes) - 1.",
        ge=0,
    )

    @root_validator
    def number_planes_le_number_satellites(cls, values):
        p, t = values.get("number_planes"), values.get("number_satellites")
        if p is not None and t is not None and p > t:
            raise ValueError("number planes exceeds number satellites")
        return values

    @root_validator
    def relative_spacing_lt_number_planes(cls, values):
        p, f = values.get("number_planes"), values.get("relative_spacing")
        if p is not None and f is not None and f >= p:
            raise ValueError("relative spacing exceeds number planes - 1")
        return values

    def generate_members(self) -> List[Satellite]:
        members = []
        lead_orbit = self.orbit.to_tle()
        lead_tle = Satrec.twoline2rv(lead_orbit.tle[0], lead_orbit.tle[1])
        epoch = sat_epoch_datetime(lead_tle)
        satellites_per_plane = math.ceil(self.number_satellites / self.number_planes)
        for satellite in range(self.number_satellites):
            plane = satellite // satellites_per_plane
            delta_mean_anomaly = (
                (
                    np.mod(satellite, satellites_per_plane) * self.number_planes
                    + self.relative_spacing * plane
                )
                * 2
                * np.pi
                / (satellites_per_plane * self.number_planes)
            )
            delta_raan = (
                plane
                * (
                    2 * np.pi
                    if self.configuration == WalkerConfiguration.delta
                    else np.pi
                )
                / self.number_planes
            )
            satrec = Satrec()
            satrec.sgp4init(
                WGS72,
                "i",
                0,
                (epoch - datetime(1949, 12, 31, tzinfo=timezone.utc))
                / timedelta(days=1),
                0,
                0.0,
                0.0,
                lead_tle.ecco,
                lead_tle.argpo,
                lead_tle.inclo,
                np.mod(lead_tle.mo + delta_mean_anomaly, 2 * np.pi),
                lead_tle.no_kozai,
                np.mod(lead_tle.nodeo + delta_raan, 2 * np.pi),
            )
            tle1, tle2 = exporter.export_tle(satrec)
            members.append(
                Satellite(
                    name=f"{self.name} #{satellite+1}",
                    orbit=TwoLineElements(tle=[tle1.replace("\x00", "U"), tle2]),
                    instruments=copy.deepcopy(self.instruments),
                )
            )
        return members
