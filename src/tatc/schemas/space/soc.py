"""
Object schema for streets-of-coverage (SOC) constellations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import copy
import math
from typing import Literal

from pydantic import Field

from tatc.utils.formatting import zero_pad
from tatc.utils.observation import (
    compute_min_elevation_angle,
    swath_width_to_field_of_regard,
)

from ...constants import EARTH_MEAN_RADIUS
from ..orbit import CircularOrbit
from .base_constellation import BaseConstellation
from .satellite import Satellite
from .walker import WalkerConstellation

POLAR_INCLINATION_TOLERANCE = 10
"""Maximum difference (degrees) between a near-polar inclination and 90 degrees."""


class SOCConstellation(BaseConstellation):
    """
    A constellation that arranges member satellites following the streets of coverage pattern.

    For inclined orbits, the satellites follow a Walker delta pattern based on Joshua F.
    Anderson, Michel-Alexandre Cardin, and Paul T. Grogan (2022). "Design and analysis
    of flexible multi-layer staged deployment for satellite mega-constellations under
    demand uncertainty" Acta Astronautica, vol. 198, pp. 179-193.
    doi: 10.1016/j.actaastro.2022.05.022, with footprint centers on the hexagonal
    lattice that covers a plane with the fewest circles (rather than the hexagonal
    lattice of touching circles, which leaves gaps between them).

    For near-polar orbits, the satellites follow the polar streets-of-coverage pattern
    of Lloyd Rider (1985). "Optimized polar orbit constellations for redundant earth
    coverage" Journal of the Astronautical Sciences, vol. 33, pp. 147-161, and John
    G. Adams and Lloyd Rider (1987). "Circular polar constellations providing continuous
    single or multiple coverage above a specified latitude" Journal of the Astronautical
    Sciences, vol. 35, pp. 155-192: planes spanning 180 degrees of right ascension,
    in which the ascending (co-rotating) sides of adjacent planes are farther apart than
    the ascending and descending (counter-rotating) sides of the first and last planes
    across the "seam".
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
        ...,
        description="Relative distance between footprint centers (inclined orbits), "
        + "where 1 places the footprints on the hexagonal lattice that just "
        + "covers continuously, or relative footprint radius used as a coverage "
        + "margin (polar orbits). Smaller values add overlap as a coverage margin.",
        gt=0,
        le=1,
    )
    polar: bool | None = Field(
        default=None,
        description="True, to use the polar streets-of-coverage pattern; False, to "
        + "use the Walker delta pattern; None (default), to use the polar pattern "
        + "for near-polar orbits (inclination within "
        + f"{POLAR_INCLINATION_TOLERANCE} degrees of 90 degrees).",
    )

    def is_polar(self) -> bool:
        """
        Checks whether this constellation follows the polar streets-of-coverage pattern.

        Returns:
            bool: True, if the polar pattern is used.
        """
        if self.polar is not None:
            return self.polar
        return abs(self.orbit.get_inclination() - 90) <= POLAR_INCLINATION_TOLERANCE

    def get_footprint_angle(self) -> float:
        """
        Gets the Earth central angle (degrees) of the footprint radius
        [Eqs. (19)-(20) in Anderson et al. (2022)].

        Returns:
            float: the footprint Earth central angle
        """
        # minimum elevation angle (degrees)
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
        # earth central angle [Eq. (20) in Anderson et al. (2022)]
        return 90 - e - eta

    def get_polar_design(self) -> tuple[int, int, float, float]:
        """
        Gets the polar streets-of-coverage design for continuous single coverage
        with the fewest satellites (Adams and Rider, 1987). With S satellites per
        plane, footprints of Earth central angle gamma (reduced by the packing
        distance as a coverage margin) overlap along each plane in a street of
        half width c, where cos(c) = cos(gamma) / cos(180 / S). Adjacent
        co-rotating planes may be up to gamma + c apart and the counter-rotating
        planes across the seam up to 2c apart, so that P planes span 180 degrees
        if (P - 1) (gamma + c) + 2 c >= 180. Both spacings are scaled down
        equally to span exactly 180 degrees.

        Returns:
            tuple[int, int, float, float]: the number of satellites per plane,
                the number of planes, the spacing (degrees) between co-rotating
                planes, and the spacing (degrees) across the seam.
        """
        gamma = self.get_footprint_angle() * self.packing_distance
        if gamma <= 0:
            raise ValueError("Footprint too small for a streets-of-coverage design.")
        best, best_satellites = None, math.inf
        # the fewest satellites per plane whose footprints overlap
        min_satellites_per_plane = math.floor(180 / gamma) + 1
        for satellites_per_plane in range(
            min_satellites_per_plane, 3 * min_satellites_per_plane + 1
        ):
            # street half width (degrees)
            c = math.degrees(
                math.acos(
                    math.cos(math.radians(gamma))
                    / math.cos(math.pi / satellites_per_plane)
                )
            )
            number_planes = 1 + math.ceil((180 - 2 * c) / (gamma + c) - 1e-9)
            if satellites_per_plane * number_planes < best_satellites:
                best_satellites = satellites_per_plane * number_planes
                scale = 180 / ((number_planes - 1) * (gamma + c) + 2 * c)
                best = (
                    satellites_per_plane,
                    number_planes,
                    (gamma + c) * scale,
                    2 * c * scale,
                )
        return best  # type: ignore

    def get_satellites_per_plane(self) -> int:
        """
        Gets the number of satellites per plane.

        Returns:
            int: the number of satellites per plane
        """
        if self.is_polar():
            return self.get_polar_design()[0]
        walker = self.generate_walker()
        return walker.number_satellites // walker.number_planes

    def get_number_planes(self) -> int:
        """
        Gets the number of planes.

        Returns:
            int: the number of planes
        """
        if self.is_polar():
            return self.get_polar_design()[1]
        return self.generate_walker().number_planes

    def generate_walker(self) -> WalkerConstellation:
        """
        Generate a WalkerConstellation fitting the Streets of Coverage description
        for inclined orbits (the Walker delta pattern). Polar streets-of-coverage
        designs, whose planes are unequally spaced, are not Walker constellations.

        Returns:
            WalkerConstellation: the member satellites following the Walker pattern.
        """
        if self.is_polar():
            raise ValueError(
                "A polar streets-of-coverage design is not a Walker constellation; "
                + "use generate_members() or set polar=False."
            )
        # satellite footprint radius (m) [Eq. (21) in Anderson et al. (2022)]
        r_foot = EARTH_MEAN_RADIUS * math.sin(math.radians(self.get_footprint_angle()))

        # distance between adjacent footprint centers (m) [cf. Eq. (23) in Anderson
        # et al. (2022)], on the hexagonal covering lattice: footprints overlap
        # along each plane, leaving no gaps with the footprints of adjacent planes
        # (Eq. (23) spaced footprints by their diameter, so that they only touched)
        d_f = math.sqrt(3) * r_foot * self.packing_distance

        # distance between adjacent planes (m) [cf. Eq. (24) in Anderson et al. (2022)]
        d_p = 1.5 * r_foot * self.packing_distance

        # angle (radians) between footprint centers [Eq. (25) in Anderson et al. (2022)]
        gamma_f = 2 * math.asin((0.5 * d_f) / (EARTH_MEAN_RADIUS))

        # number of satellites per plane [Eq. (26) in Anderson et al. (2022)]
        satellites_per_plane = math.ceil((2 * math.pi) / gamma_f)

        # angle (radians) between adjacent planes [Eq. (27) in Anderson et al. (2022)]
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
            # offset adjacent planes by half a within-plane satellite
            # spacing, so the d_p row spacing (derived above for the
            # hexagonal covering lattice) actually yields a staggered
            # hexagonal layout rather than a plain rectangular grid of planes
            relative_spacing=number_planes // 2,
        )

    def generate_members(self) -> list[Satellite]:
        """
        Generate member satellites for this streets-of-coverage constellation.
        In the polar pattern, adjacent co-rotating planes are offset by half the
        in-plane spacing, so that each plane's satellites fill the gaps between
        the footprints of its neighbors.

        Returns:
            list[Satellite]: the member satellites
        """
        if not self.is_polar():
            return self.generate_walker().generate_members()
        satellites_per_plane, number_planes, spacing, _ = self.get_polar_design()
        number_satellites = satellites_per_plane * number_planes
        return [
            Satellite(
                name=zero_pad(self.name, number_satellites, i + 1),
                orbit=self.orbit.get_derived_orbit(
                    (i % satellites_per_plane) * 360 / satellites_per_plane
                    + (i // satellites_per_plane) * 180 / satellites_per_plane,
                    (i // satellites_per_plane) * spacing,
                ),
                instruments=copy.deepcopy(self.instruments),
            )
            for i in range(number_satellites)
        ]
