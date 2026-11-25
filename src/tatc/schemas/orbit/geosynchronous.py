"""
Object schemas for geosynchronous orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import timedelta
from typing import Literal

import numpy as np
from pydantic import AliasChoices, Field
from skyfield.api import wgs84

from ... import constants, utils
from ...utils.cache import get_cached
from .base_circular import CircularOrbitBase
from .gp import GeneralPerturbationsOrbit

# mean altitude (meters) yielding an orbit period of exactly one sidereal
# day, via Kepler's third law
_GEOSYNCHRONOUS_MEAN_ALTITUDE = (
    utils.orbital.mean_motion_to_semimajor_axis(360 / constants.EARTH_SIDEREAL_DAY_S)
    - constants.EARTH_MEAN_RADIUS
)
# the corrections of the node that place the sub-satellite point at epoch, at most,
# and the error (degrees) at which they stop
_NODE_ITERATIONS, _NODE_TOLERANCE = 5, 1e-9


class GeosynchronousOrbit(CircularOrbitBase):
    """
    Orbit defined by geosynchronous parameters: a fixed longitude rather
    than a right ascension of ascending node.

    The orbit is propagated as maintained at its longitude (as by
    station-keeping), on a repeat ground track with a repeat cycle of one
    nodal day (see `GeneralPerturbationsOrbit.repeat_cycle`): its mean motion
    is adjusted so that it completes one orbit per nodal day, and every time
    is propagated within a day of its epoch. Otherwise, the Earth's
    oblateness and the Moon and Sun would move its longitude by tenths of a
    degree within a month.
    """

    type: Literal["geosynchronous"] = Field(
        default="geosynchronous", description="Orbit type discriminator."
    )
    mean_altitude: float = Field(
        default=_GEOSYNCHRONOUS_MEAN_ALTITUDE,
        description="Mean altitude (meters).",
        ge=0,
        validation_alias=AliasChoices("mean_altitude", "altitude"),
    )
    inclination: float = Field(
        default=0, description="Inclination (degrees).", ge=0, lt=180
    )
    longitude: float = Field(
        ...,
        description="Longitude (decimal degrees) of the sub-satellite point "
        + "in the WGS 84 coordinate system.",
        ge=-180,
        le=180,
    )

    def get_inclination(self) -> float:
        """
        Gets the inclination.

        Returns:
            float: the inclination (degrees)
        """
        return self.inclination

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node, derived from this orbit's
        fixed longitude, since a geosynchronous satellite remains above a fixed
        Earth-relative longitude rather than a fixed inertial-frame RAAN: the node
        (in the TEME frame of SGP4) that places the satellite's sub-satellite
        point at epoch at the orbit's longitude, rotated by its argument of
        latitude (its true anomaly and argument of perigee), as SGP4 propagates
        it (see `to_gp_orbit`).

        The node starts from the Earth's rotation angle that orients TEME
        (Greenwich Mean Sidereal Time) at epoch plus the longitude, and is
        corrected for the periodic lunar and solar terms of SGP4's deep-space
        propagation, which displace the satellite at its epoch from the position
        of its mean elements (by up to about 0.04 degrees, by the date). As with
        SunSynchronousOrbit's equator-crossing-time convention, the longitude is
        that of the node, independent of the true anomaly.

        Returns:
            float: the right ascension of ascending node (degrees)
        """
        return get_cached(
            self,
            "right_ascension_ascending_node",
            self.model_dump_json(),
            self._compute_right_ascension_ascending_node,
        )

    def _compute_right_ascension_ascending_node(self) -> float:
        """
        Computes the right ascension of ascending node (see
        `get_right_ascension_ascending_node`), without caching.

        Returns:
            float: the right ascension of ascending node (degrees)
        """
        epoch_time = constants.timescale.from_datetime(self.epoch)
        # the longitude of the sub-satellite point, from the node's longitude and
        # the argument of latitude
        latitude_argument = np.radians(
            self.get_perigee_argument() + self.get_true_anomaly()
        )
        inclination = np.radians(self.inclination)
        expected = self.longitude + np.degrees(
            np.arctan2(
                np.cos(inclination) * np.sin(latitude_argument),
                np.cos(latitude_argument),
            )
        )
        node = epoch_time.gmst * 15 + self.longitude  # hours -> degrees
        for _ in range(_NODE_ITERATIONS):
            track = self._build_gp_orbit(node).get_orbit_track(self.epoch)
            error = (expected - wgs84.subpoint_of(track).longitude.degrees + 180) % 360
            error -= 180
            node += error
            if abs(error) < _NODE_TOLERANCE:
                break
        return float(node % 360)

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> GeosynchronousOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and
        longitude (equivalent, for a geosynchronous orbit, to a right
        ascension of ascending node perturbation, since both epoch and
        Earth's rotation angle stay fixed).

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            GeosynchronousOrbit: the derived orbit
        """
        true_anomaly = utils.orbital.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360)
        )
        longitude = ((self.longitude + delta_raan + 180) % 360) - 180
        return GeosynchronousOrbit(
            mean_altitude=self.mean_altitude,
            inclination=self.inclination,
            longitude=longitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
        )

    def _compute_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Computes a general perturbations orbit representation of this orbit,
        maintained at its longitude on a repeat ground track with a repeat
        cycle of one nodal day.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        return self._build_gp_orbit(self.get_right_ascension_ascending_node())

    def _build_gp_orbit(self, node: float) -> GeneralPerturbationsOrbit:
        """
        Builds the general perturbations orbit of this orbit with a right
        ascension of ascending node, maintained at its longitude on a repeat
        ground track with a repeat cycle of one nodal day.

        Args:
            node (float): The right ascension of ascending node (degrees).

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        # deferred import to avoid a circular import (keplerian.py imports the
        # orbit base module)
        # pylint: disable-next=import-outside-toplevel,cyclic-import
        from .keplerian import KeplerianOrbit

        elements = KeplerianOrbit(
            semimajor_axis=self.get_semimajor_axis(),
            inclination=self.get_inclination(),
            right_ascension_ascending_node=node % 360,
            eccentricity=self.get_eccentricity(),
            perigee_argument=self.get_perigee_argument(),
            true_anomaly=self.get_true_anomaly(),
            epoch=self.get_epoch(),
        ).to_gp_orbit()
        return GeneralPerturbationsOrbit(
            elements=elements.elements,
            remove_drag=True,
            repeat_cycle=timedelta(days=1),
        )
