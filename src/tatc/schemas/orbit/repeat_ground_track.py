"""
Object schema for circular orbits maintained on a repeat ground track.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timedelta
from typing import Literal

import numpy as np
from pydantic import AwareDatetime, BaseModel, ConfigDict, Field, model_validator

from ... import constants, utils
from ...utils.cache import get_cached
from .gp import GeneralPerturbationsOrbit
from .gp_elements import GeneralPerturbationsElements


class AscendingNode(BaseModel):
    """
    An ascending node of an orbit: a time at which the satellite crosses the
    Earth's equatorial plane from south to north in the Earth-fixed (ITRS)
    frame, and the longitude of the crossing.
    """

    model_config = ConfigDict(extra="forbid")

    time: AwareDatetime = Field(..., description="Time of the crossing.")
    longitude: float = Field(
        ...,
        description="Longitude (decimal degrees) of the crossing.",
        ge=-180,
        le=180,
    )


class RepeatGroundTrackOrbit(BaseModel):
    """
    Circular orbit maintained on a repeat ground track, as an operator
    maintains it against drag (and, for a sun-synchronous orbit, at a
    constant local time of its ascending node): the satellite completes
    `revolutions` orbits (between ascending nodes) in `days` nodal days (for
    a sun-synchronous orbit, mean solar days), after which it returns over
    the same ground track, pinned by one of its ascending nodes
    (`ascending_node`).

    Unlike a two-line element set, which describes a satellite near its
    epoch and diverges from it as the satellite is maneuvered, this orbit
    describes the ground track the satellite is kept on, so that it can be
    built for any time, including years from its ascending node. Its
    altitude, inclination, and the local time of its ascending node are
    derived from its repeat, its inclination (if not sun-synchronous), and
    its node, rather than given.

    The orbit is propagated as a `GeneralPerturbationsOrbit` with one
    element (circular, without drag, at the time of the node), maintained on
    its repeat ground track (`repeat_cycle`), so that every time is
    propagated within one repeat cycle of the node. Its repeat cycle is a
    whole number of days of the International System of Units (in
    Terrestrial Time), whereas the Earth's rotation keeps Universal Time
    (UT1): the difference between the two accumulates (a few tenths of a
    second per year recently, and up to about a second per year
    historically) and delays or advances the crossings far from the node by
    that much, on the same ground track.
    """

    model_config = ConfigDict(extra="forbid")

    type: Literal["repeat"] = Field(
        default="repeat", description="Orbit type discriminator."
    )
    revolutions: int = Field(
        ...,
        description="Number of orbits (between ascending nodes) in the repeat cycle.",
        gt=0,
        examples=[227],
    )
    days: int = Field(
        ...,
        description="Number of nodal days (for a sun-synchronous orbit, mean "
        + "solar days) in the repeat cycle.",
        gt=0,
        examples=[16],
    )
    inclination: float | None = Field(
        default=None,
        description="Inclination (decimal degrees), or None (the default) for "
        + "a sun-synchronous orbit, whose inclination is derived so that its "
        + "nodal day is a mean solar day.",
        gt=0,
        lt=180,
    )
    ascending_node: AscendingNode = Field(
        ...,
        description="An ascending node of the orbit, which pins its ground "
        + "track and the times at which the satellite follows it.",
    )

    @model_validator(mode="after")
    def _validate_repeat(self) -> RepeatGroundTrackOrbit:
        """
        Validates that the repeat gives an orbit above the Earth's surface
        and in low or medium Earth orbit, and builds (and caches) its
        general perturbations orbit, which checks that its revolutions are
        consistent with its days.
        """
        if self.revolutions < self.days:
            raise ValueError(
                f"{self.revolutions} revolutions in {self.days} days is "
                "beyond geosynchronous orbit (use a GeosynchronousOrbit)"
            )
        period = self.days * constants.EARTH_SOLAR_DAY_S / self.revolutions
        if period < 5000:
            raise ValueError(
                f"{self.revolutions} revolutions in {self.days} days gives an "
                f"orbit period of {period:.0f} s, below the Earth's surface"
            )
        self.to_gp_orbit()
        return self

    def is_sun_synchronous(self) -> bool:
        """
        Checks whether this orbit is sun-synchronous (has no `inclination`).

        Returns:
            bool: True, if this orbit is sun-synchronous
        """
        return self.inclination is None

    def to_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Converts this orbit to a general perturbations orbit representation:
        one circular element without drag at the time of the ascending node,
        with a mean motion (and, for a sun-synchronous orbit, inclination)
        for which `revolutions` nodal periods span `days` nodal (or mean
        solar) days, and a mean anomaly and right ascension of ascending node
        for which the satellite crosses the ascending node, as propagated,
        at its time and longitude; maintained on its repeat ground track
        with a repeat cycle of `days`. The conversion is cached.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        return get_cached(
            self, "gp_orbit", self.model_dump_json(), self._compute_gp_orbit
        )

    def _compute_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Computes the general perturbations orbit, without caching (see
        `to_gp_orbit`).

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        element = self._compute_repeat_element()
        repeat_cycle = timedelta(days=self.days)
        nodal_period, _ = element.get_nodal_period_and_day()
        refined = element.refine_repeat_cycle(repeat_cycle).total_seconds()
        if round(refined / nodal_period) != self.revolutions:
            raise ValueError(
                f"the orbit maintained on a repeat cycle of {self.days} days "
                f"completes {round(refined / nodal_period)} revolutions, not "
                f"{self.revolutions}"
            )
        node_time = self.ascending_node.time
        for _ in range(10):
            orbit = GeneralPerturbationsOrbit(
                elements=[element], remove_drag=True, repeat_cycle=repeat_cycle
            )
            # the crossing nearest the node, within half an orbit
            nodes = orbit.get_ascending_nodes(
                node_time - timedelta(seconds=nodal_period / 2),
                node_time + timedelta(seconds=nodal_period / 2),
            )
            if len(nodes) == 0:
                raise ValueError("the orbit does not cross its ascending node")
            nearest = nodes.iloc[
                int(np.argmin(np.abs((nodes.time - node_time).dt.total_seconds())))
            ]
            delta_time = (nearest.time - node_time).total_seconds()
            delta_longitude = (
                self.ascending_node.longitude - nearest.longitude + 180
            ) % 360 - 180
            if abs(delta_time) < 1e-4 and abs(delta_longitude) < 1e-8:
                break
            # a later mean anomaly advances the crossing (east, by the Earth's
            # rotation in the time), and a larger right ascension of
            # ascending node moves it east
            element = element.model_copy(
                update={
                    "mean_anomaly": float(
                        np.mod(
                            element.mean_anomaly + 360 * delta_time / nodal_period,
                            360,
                        )
                    ),
                    "ra_of_asc_node": float(
                        np.mod(
                            element.ra_of_asc_node
                            + delta_longitude
                            - np.degrees(constants.EARTH_ROTATION_RATE) * delta_time,
                            360,
                        )
                    ),
                }
            )
        return orbit

    def _compute_repeat_element(self) -> GeneralPerturbationsElements:
        """
        Computes the element at the time of the ascending node with the mean
        motion (and, for a sun-synchronous orbit, inclination) of the repeat,
        before its phase is adjusted to cross the node.

        Returns:
            GeneralPerturbationsElements: the element
        """
        node_time = self.ascending_node.time
        period = self.days * constants.EARTH_SOLAR_DAY_S / self.revolutions
        mean_motion = 360 / period
        semimajor_axis = utils.orbital.mean_motion_to_semimajor_axis(mean_motion)
        if self.inclination is None:
            # classical J2 estimate of the sun-synchronous inclination
            inclination = float(
                np.degrees(np.arccos(-np.power(semimajor_axis / 12352000, 7 / 2)))
            )
        else:
            inclination = self.inclination
        epoch_time = constants.timescale.from_datetime(node_time)
        element = GeneralPerturbationsElements(
            epoch=node_time,
            mean_motion=mean_motion,
            eccentricity=0,
            inclination=inclination,
            ra_of_asc_node=float(
                np.mod(epoch_time.gmst * 15 + self.ascending_node.longitude, 360)
            ),
            arg_of_pericenter=0,
            mean_anomaly=0,
        )
        for _ in range(20):
            nodal_period, nodal_day = element.get_nodal_period_and_day()
            day = constants.EARTH_SOLAR_DAY_S if self.inclination is None else nodal_day
            target = self.days * day / self.revolutions
            update = {}
            if abs(nodal_period - target) >= 1e-7:
                # the nodal period varies (nearly) inversely with mean motion
                update["mean_motion"] = element.mean_motion * nodal_period / target
            if (
                self.inclination is None
                and abs(nodal_day - constants.EARTH_SOLAR_DAY_S) >= 1e-7
            ):
                # the node precession rate is (nearly) proportional to the
                # cosine of inclination
                ratio = (
                    constants.EARTH_ROTATION_RATE
                    - 2 * np.pi / constants.EARTH_SOLAR_DAY_S
                ) / (constants.EARTH_ROTATION_RATE - 2 * np.pi / nodal_day)
                update["inclination"] = float(
                    np.degrees(
                        np.arccos(
                            np.clip(np.cos(np.radians(element.inclination)) * ratio, -1, 1)
                        )
                    )
                )
            if not update:
                break
            element = element.model_copy(update=update)
        return element

    def get_epoch(self) -> datetime:
        """
        Gets the epoch: the time of the ascending node.

        Returns:
            datetime: the epoch
        """
        return self.ascending_node.time

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters).

        Returns:
            float: the semimajor axis
        """
        return self.to_gp_orbit().get_semimajor_axis()

    def get_mean_altitude(self) -> float:
        """
        Gets the mean altitude (meters) above the WGS 84 mean radius.

        Returns:
            float: the mean altitude
        """
        return self.to_gp_orbit().get_mean_altitude()

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion (degrees per second).

        Returns:
            float: the mean motion
        """
        return self.to_gp_orbit().get_mean_motion()

    def get_orbit_period(self) -> timedelta:
        """
        Gets the orbit period (of the mean motion).

        Returns:
            timedelta: the orbit period
        """
        return self.to_gp_orbit().get_orbit_period()

    def get_inclination(self) -> float:
        """
        Gets the inclination (decimal degrees), derived for a sun-synchronous
        orbit.

        Returns:
            float: the inclination
        """
        return self.to_gp_orbit().get_inclination()

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity, which is zero.

        Returns:
            float: the eccentricity
        """
        return 0.0

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node (decimal degrees) of the
        element at the epoch.

        Returns:
            float: the right ascension of ascending node
        """
        return self.to_gp_orbit().get_right_ascension_ascending_node()

    def get_perigee_argument(self) -> float:
        """
        Gets the argument of perigee (decimal degrees) of the element at the
        epoch.

        Returns:
            float: the argument of perigee
        """
        return self.to_gp_orbit().get_perigee_argument()

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees) of the element at the epoch.

        Returns:
            float: the mean anomaly
        """
        return self.to_gp_orbit().get_mean_anomaly()

    def get_true_anomaly(self) -> float:
        """
        Gets the true anomaly (decimal degrees) of the element at the epoch.

        Returns:
            float: the true anomaly
        """
        return self.to_gp_orbit().get_true_anomaly()

    def get_local_time_ascending_node(self) -> float:
        """
        Gets the local mean solar time (hours) of the ascending node: the
        universal time (UT1) of the node plus its longitude (15 degrees per
        hour). For a sun-synchronous orbit, it is the same at every node.

        Returns:
            float: the local mean solar time (hours, from 0 to 24)
        """
        node_time = constants.timescale.from_datetime(self.ascending_node.time)
        universal_hours = np.mod(node_time.ut1 + 0.5, 1) * 24
        return float(np.mod(universal_hours + self.ascending_node.longitude / 15, 24))

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> RepeatGroundTrackOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node: its ascending node is earlier by the
        fraction of a nodal period of the mean anomaly, and east by the right
        ascension of ascending node and by the Earth's rotation relative to
        the orbit plane in that time.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            RepeatGroundTrackOrbit: the derived orbit
        """
        element = self.to_gp_orbit().elements[0]
        nodal_period, nodal_day = element.get_nodal_period_and_day()
        if self.is_sun_synchronous():
            nodal_day = constants.EARTH_SOLAR_DAY_S
        delta_time = nodal_period * delta_mean_anomaly / 360
        longitude = (
            self.ascending_node.longitude + delta_raan + 360 * delta_time / nodal_day
        )
        return RepeatGroundTrackOrbit(
            revolutions=self.revolutions,
            days=self.days,
            inclination=self.inclination,
            ascending_node=AscendingNode(
                time=self.ascending_node.time - timedelta(seconds=delta_time),
                longitude=(longitude + 180) % 360 - 180,
            ),
        )
