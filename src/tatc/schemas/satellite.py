# -*- coding: utf-8 -*-
"""
Object schemas for satellites.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""
from __future__ import annotations

import copy
from datetime import timedelta
from enum import Enum
import math
from typing import List, Union

import numpy as np
from pydantic import BaseModel, Field, root_validator
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
        return -360 * self.orbit.get_mean_motion() * (self.interval / timedelta(days=1))

    def get_delta_raan(self) -> float:
        """
        Gets the difference in right ascension of ascending node (decimal
        degrees) for adjacent member satellites.

        Returns:
            float: the difference in right ascension of ascending node
        """
        if self.repeat_ground_track:
            return 360 * (self.interval / timedelta(days=1))
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
    """
    Enumeration of different Walker constellation configurations.
    """

    DELTA = "delta"
    STAR = "star"


class WalkerConstellation(Satellite):
    """
    A constellation that arranges member satellites following the Walker pattern.
    """

    type: Literal["walker"] = Field(
        "walker", description="Space system type discriminator."
    )
    configuration: WalkerConfiguration = Field(
        WalkerConfiguration.DELTA, description="Walker configuration."
    )
    orbit: Union[
        TwoLineElements, SunSynchronousOrbit, CircularOrbit, KeplerianOrbit
    ] = Field(..., description="Lead orbit for this constellation.")
    number_satellites: int = Field(
        1, description="Number of satellites in the constellation.", ge=1
    )
    number_planes: int = Field(
        1,
        description="The number of equally-spaced planes in a Walker Delta "
        + "constellation. Ranges from 1 to (number of satellites).",
        ge=1,
    )
    relative_spacing: int = Field(
        0,
        description="Relative spacing of satellites between plans for a Walker Delta "
        + "constellation. Ranges from 0 for equal true anomaly to "
        + "(number of planes) - 1. For example, `relative_spacing=1` "
        + "means the true anomaly is shifted by `360/number_satellites` "
        + "between adjacent planes.",
        ge=0,
    )

    @root_validator
    def number_planes_le_number_satellites(cls, values):
        """
        Validates the number of planes given the number of satellites.
        """
        planes = values.get("number_planes")
        count = values.get("number_satellites")
        if planes is not None and count is not None and planes > count:
            raise ValueError("number planes exceeds number satellites")
        return values

    @root_validator
    def relative_spacing_lt_number_planes(cls, values):
        """
        Validates the relative spacing given the number of planes.
        """
        planes = values.get("number_planes")
        spacing = values.get("relative_spacing")
        if planes is not None and spacing is not None and spacing >= planes:
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
        if self.configuration == WalkerConfiguration.DELTA:
            return 360 / self.number_planes
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

class MOGConstellation(Satellite):
    """
    A constellation that arranges member satellites following the mutual orbiting group pattern.
    """

    type: Literal["mog"] = Field(
        "mog", description="Space system type discriminator."
    )
    orbit: Union[
        TwoLineElements, SunSynchronousOrbit, CircularOrbit, KeplerianOrbit
    ] = Field(..., description="Lead orbit for this constellation.")
    number_satellites: int = Field(
        2, description="The count of the number of satellites.", ge=2, le=2
    )
    delta: float = Field(
        float(np.radians(5)), description="The separation of the angular momentum vectors of the reference orbiter and the mutually orbiting satellite", ge=0
    )
    theta: float = Field(
        float(np.radians(5)), description="", ge=0
    )
    delta_anomaly: int = Field(
        50, description="delta mean anomaly in degrees for the satellites with respect to the reference orbiter", le=360, ge=0
    )


    x_vector = [1, 0, 0]
    y_vector = [0, 1, 0]
    z_vector = [0, 0, 1]

    def get_angular_momentum_direction_of_orbiter(self) -> List[float]:
        """
        Gets the angular momentum direction of the reference circular orbiter
        
        Returns:
            float: the angular momentum direction of orbiter
        """
        inclination = (
            self.orbit.inclination 
            if isinstance(self.orbit, CircularOrbit) 
            else self.orbit.get_inclination()
        )
        cos_inclination = float(np.cos(np.radians(inclination)))
        sin_inclination = float(np.sin(np.radians(inclination)))
        return [(cos_inclination * z) - (sin_inclination * y) for z, y in zip(self.z_vector, self.y_vector)]

    def get_angular_momentum_direction_of_satellite(self) -> List[float]:
        """
        Gets the angular momentum direction of the mutual orbiting satellite
        
        Returns:
            float: the angular momentum direction of the mutual orbiting satellite
        """
        cos_delta = float(np.cos(self.delta))
        sin_delta = float(np.cos(self.delta))
        cos_theta = float(np.cos(self.theta))
        sin_theta = float(np.cos(self.theta))

        p1 = [cos_delta * l for l in self.get_angular_momentum_direction_of_orbiter()]
        p2 = [(sin_delta * cos_theta) * m for m in (np.cross(self.get_angular_momentum_direction_of_orbiter(), self.x_vector))]
        p3 = [(sin_delta * sin_theta) * x for x in self.x_vector]

        return ([p1[i] + p2[i] + p3[i] for i in range(3)])

    def get_satellite_inclination(self) -> float:
        """
        Gets the inclination of the mutual orbiting satellite
        
        Returns:
            float: the inclination of the mutual orbiting satellite
        """
        return np.arccos(np.dot(self.get_angular_momentum_direction_of_satellite(), self.z_vector))

    def get_delta_raan(self) -> float:
        """
        Gets the longitude of the right ascending node of the mutual orbiter
        Returns:
            float: longitude of the mutual orbiter's right ascension of the ascending node
        """
        raan_mo = np.arctan2(np.dot(self.get_angular_momentum_direction_of_satellite(), self.x_vector), 
                             np.dot([-1 * x for x in self.get_angular_momentum_direction_of_satellite()], self.y_vector)
                            )
        return raan_mo

    def get_delta_mean_anomaly(self) -> float:
        return np.radians(self.delta_anomaly)

    def generate_members(self) -> List[Satellite]:
        """
        Generate space system member satellites.
        Returns:
            List[Satellite]: the member satellites
        """
        return [
            Satellite(
                name=f"{self.name} #{0+1:02d}",
                orbit=self.orbit.get_derived_orbit(
                    -0.5 * (self.get_delta_mean_anomaly()), -1 * self.get_delta_raan()
                ),
                instruments=copy.deepcopy(self.instruments),
            ),
            Satellite(
                name=f"{self.name} #{1+1:02d}",
                orbit=self.orbit.get_derived_orbit(
                    (0.5 * self.get_delta_mean_anomaly()), 1 * self.get_delta_raan()
                ),
                instruments=copy.deepcopy(self.instruments),
            )
        ]