"""
Object schema for mutual orbiting group (MOG) constellations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import copy
from typing import Literal

import numpy as np
from pydantic import Field

from tatc.utils.formatting import zero_pad

from ..orbit import CircularOrbit, KeplerianOrbit
from .base_constellation import BaseConstellation
from .satellite import Satellite


class MOGConstellation(BaseConstellation):
    """
    A constellation that arranges member satellites following the mutual orbiting group pattern.

    Based on Stephen Leroy, Riley Fitzgerald, Kerri Cahoy, James Abel, and James Clark (2020).
    "Orbital Maintenance of a Constellation of CubeSats for Internal Gravity Wave Tomography,"
    IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 13,
    pp. 307-317. doi: 10.1109/JSTARS.2019.2961084
    """

    type: Literal["mog"] = Field(default="mog", description="Space system type discriminator.")
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

    def generate_members(self) -> list[Satellite]:
        """
        Generate space system member satellites.

        Returns:
            list[Satellite]: the member satellites
        """

        orbits = []

        # semimajor axis (m)
        a = self.orbit.get_semimajor_axis()

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

            # direction of mutual orbit ascending node w.r.t. 
            # Earth's center of mass [Eq. (25) in Leroy et al. (2020)]
            p_node = np.cos(omega - omega_0) * np.array([1, 0, 0]) + np.sin(
                omega - omega_0
            ) * np.array([0, 1, 0])

            # direction of mutual and reference orbit intersection 
            # ([Eq. (26) in Leroy et al. (2020)])
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
                    semimajor_axis=self.orbit.get_semimajor_axis(),
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
