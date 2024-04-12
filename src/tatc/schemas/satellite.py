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

from ..constants import EARTH_MEAN_RADIUS
from .instrument import Instrument
from .orbit import TwoLineElements, CircularOrbit, SunSynchronousOrbit, KeplerianOrbit
from tatc.utils import zero_pad


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
                name=zero_pad(self.name, self.number_satellites, i + 1),
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
                name=zero_pad(self.name, self.number_satellites, i + 1),
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

    Based on Stephen Leroy, Riley Fitzgerald, Kerri Cahoy, James Abel, and James Clark (2020).
    "Orbital Maintenance of a Constellation of CubeSats for Internal Gravity Wave Tomography,"
    IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 13,
    pp. 307-317. doi: 10.1109/JSTARS.2019.2961084
    """

    type: Literal["mog"] = Field("mog", description="Space system type discriminator.")
    orbit: CircularOrbit = Field(
        ..., description="Reference circular orbit for this constellation."
    )
    parallel_axis: float = Field(
        ...,
        description="Mutual orbit axis length (m) parallel to velocity vector.",
        gt=0,
    )
    transverse_axis: float = Field(
        ...,
        description="Mutual orbit axis length (m) transverse to velocity vector.",
        gt=0,
    )
    clockwise: bool = Field(True, description="True, if the mutual orbit is clockwise.")
    number_satellites: int = Field(
        2, description="Number of equally-spaced mutually orbiting satellites.", gt=0
    )

    def generate_members(self) -> List[Satellite]:
        """
        Generate space system member satellites.
        Returns:
            List[Satellite]: the member satellites
        """

        orbits = []

        # semimajor axis (m)
        a = self.orbit.altitude + EARTH_MEAN_RADIUS

        # angle of separation (radians) of angular momentum vectors for reference and mutual orbiter
        delta = self.transverse_axis / (2 * a)

        # eccentricity of the mutual orbit
        e = self.parallel_axis / (4 * a)

        # direction of the mutual orbit (1: clockwise, -1: counter-clockwise)
        s = 1 if self.clockwise else -1

        # inclination (radians) of reference orbit
        i_0 = np.radians(self.orbit.inclination)

        # right ascension of ascending node (radians) of reference orbit
        omega_0 = np.radians(self.orbit.right_ascension_ascending_node)

        # direction of angular momentum for reference orbit [Eq. (21) in Leroy et al. (2020)]
        l_0 = np.cos(i_0) * np.array([0, 0, 1]) - np.sin(i_0) * np.array([0, 1, 0])

        # angle describing position of the mutual orbiter w.r.t. reference orbiter at ascending node
        for theta in np.linspace(0, 2 * np.pi, self.number_satellites, endpoint=False):

            # direction of angular momentum of mutual orbit [Eq. (22) in Leroy et al. (2020)]
            l = (
                np.cos(delta) * l_0
                + np.sin(delta) * np.cos(theta) * np.cross(l_0, np.array([1, 0, 0]))
                + np.sin(delta) * np.sin(theta) * np.array([1, 0, 0])
            )

            # inclination (radians) of mutual orbiting satellite [Eq. (23) in Leroy et al. (2020)]
            i = np.arccos(np.dot(l, np.array([0, 0, 1])))

            # raan (radians) of mutual orbiting satellite [Eq. (24) in Leroy et al. (2020)]
            omega = omega_0 + np.arctan2(
                np.dot(l, np.array([1, 0, 0])), np.dot(-l, np.array([0, 1, 0]))
            )

            # direction of mutual orbit ascending node w.r.t. Earth's center of mass [Eq. (25) in Leroy et al. (2020)]
            p_node = np.cos(omega - omega_0) * np.array([1, 0, 0]) + np.sin(
                omega - omega_0
            ) * np.array([0, 1, 0])

            # direction of mutual and reference orbit intersection ([Eq. (26) in Leroy et al. (2020)])
            t = np.cross(l, l_0) / np.sin(delta)

            # perigee direction [Eq. (27) in Leroy et al. (2020)]
            p_peri = s * np.cross(t, l)

            # argument of perigee [Eq. (28) in Leroy et al. (2020)]
            w = np.arctan2(np.dot(p_peri, np.cross(l, p_node)), np.dot(p_peri, p_node))

            # time from mutual orbiter passing through its perigee to time when circular orbiter
            # passes through its ascending node [Eq. (31) in Leroy et al. (2020)]
            n = np.arctan2(
                s * np.dot(t, np.array([1, 0, 0])),
                s * np.dot(t, np.cross(l_0, np.array([1, 0, 0]))),
            )

            # eccentric anomaly (radians) implicit equation [Eq. (32) in Leroy et al. (2020)]
            psi_ = 0
            psi = 0.1  # initial guess
            while np.abs(psi - psi_) > 1e-6:  # convergence criterion
                psi_ = psi
                psi = n + e * np.sin(psi_)

            # true anomaly (radians) [Eq. (33a) in Leroy et al. (2020)]
            nu = np.arctan2(np.sin(psi) * np.sqrt(1 - e**2), np.cos(psi) - e)

            orbits.append(
                KeplerianOrbit(
                    altitude=self.orbit.altitude,
                    inclination=(360 + np.degrees(i)) % 360,
                    eccentricity=e,
                    right_ascension_ascending_node=(360 + np.degrees(omega)) % 360,
                    perigee_argument=(360 + np.degrees(w)) % 360,
                    true_anomaly=(360 + np.degrees(nu)) % 360,
                    epoch=self.orbit.epoch,
                )
            )

        return [
            Satellite(
                name=zero_pad(self.name, self.number_satellites, i + 1),
                orbit=orbit,
                instruments=copy.deepcopy(self.instruments),
            )
            for i, orbit in enumerate(orbits)
        ]
