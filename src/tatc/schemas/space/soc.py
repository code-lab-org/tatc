"""
Object schema for streets-of-coverage (SOC) constellations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import math
from typing import Literal

from pydantic import Field

from tatc.utils.observation import (
    compute_min_elevation_angle,
    swath_width_to_field_of_regard,
)

from ...constants import EARTH_MEAN_RADIUS
from ..orbit import CircularOrbit
from .base_constellation import BaseConstellation
from .satellite import Satellite
from .walker import WalkerConstellation


class SOCConstellation(BaseConstellation):
    """
    A constellation that arranges member satellites following the streets of coverage pattern.

    Based on Joshua F. Anderson, Michel-Alexandre Cardin, and Paul T. Grogan (2022).
    "Design and analysis of flexible multi-layer staged deployment for satellite
    mega-constellations under demand uncertainty"
    Acta Astronautica, vol. 198,
    pp. 179-193. doi: 10.1016/j.actaastro.2022.05.022
    """

    type: Literal["soc"] = Field(
        default="soc", description="Space system type discriminator."
    )
    orbit: CircularOrbit = Field(
        ..., description="Reference circular orbit for this constellation."
    )
    swath_width: float = Field(
        ..., description="Observation diameter (meters) at specified elevation.", gt=0
    )
    packing_distance: float = Field(
        ..., description="Relative distance between footprint centers", gt=0, le=1
    )

    def generate_walker(self) -> WalkerConstellation:
        """
        Generate a WalkerConstellation fitting the Streets of Coverage description.

        Returns:
            WalkerConstellation: the member satellites following the Walker pattern.
        """
        # compute min elevation angle
        e = compute_min_elevation_angle(
            altitude=self.orbit.mean_altitude,
            field_of_regard=swath_width_to_field_of_regard(
                altitude=self.orbit.mean_altitude, swath_width=self.swath_width
            ),
        )

        # nadir angle (degrees) [Eq. (19) in Anderson et al. (2022)]
        eta = math.degrees(
            math.asin(
                (EARTH_MEAN_RADIUS / (EARTH_MEAN_RADIUS + self.orbit.mean_altitude))
                * math.cos(math.radians(e))
            )
        )

        # compute gamma (earth central angle) [Eq. (20) in Anderson et al. (2022)]
        gamma = 90 - e - eta

        # satellite footprint radius (km) [Eq. (21) in Anderson et al. (2022)]
        r_foot = EARTH_MEAN_RADIUS * math.sin(math.radians(gamma))

        # distance between adjacent footprint centers (km) [Eq. (23) in Anderson et al. (2022)]
        d_f = 2 * r_foot * self.packing_distance

        # distance between adjacent planes (km) [Eq. (24) in Anderson et al. (2022)]
        d_p = math.sqrt(3) * r_foot * self.packing_distance

        # angle (degrees) between footprint centers [Eq. (25) in Anderson et al. (2022)]
        gamma_f = 2 * math.asin((0.5 * d_f) / (EARTH_MEAN_RADIUS))

        # number of satellites per plane [Eq. (26) in Anderson et al. (2022)]
        satellites_per_plane = math.ceil((2 * math.pi) / gamma_f)

        # angle (degrees) between adjacent planes [Eq. (27) in Anderson et al. (2022)]
        gamma_p = 2 * math.asin((0.5 * d_p) / (EARTH_MEAN_RADIUS))

        # number of planes [Eq. (28) in Anderson et al. (2022)]
        number_planes = math.ceil((2 * math.pi) / gamma_p)

        number_satellites = satellites_per_plane * number_planes

        return WalkerConstellation(
            name=self.name,
            orbit=self.orbit,
            instruments=self.instruments,
            number_satellites=number_satellites,
            number_planes=number_planes,
        )

    def generate_members(self) -> list[Satellite]:
        """
        Generate member satellites for this streets-of-coverage constellation.

        Returns:
            list[Satellite]: the member satellites
        """
        return self.generate_walker().generate_members()
