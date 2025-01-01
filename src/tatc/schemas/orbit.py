# -*- coding: utf-8 -*-
"""
Object schemas for satellite orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, time, timedelta, timezone
from typing import List, Optional, Union
import re

import numpy as np
from pydantic import BaseModel, Field, field_validator
from sgp4.api import Satrec, WGS72
from sgp4 import exporter
from sgp4.conveniences import sat_epoch_datetime
from skyfield.api import EarthSatellite, Time, wgs84
from skyfield.positionlib import Geocentric
from skyfield.framelib import itrs
from typing_extensions import Literal

from .point import Point
from .. import config, constants, utils


class TwoLineElements(BaseModel):
    """
    Orbit defined with standard two line elements.
    """

    type: Literal["tle"] = Field("tle", description="Orbit type discriminator.")
    tle: List[str] = Field(
        ...,
        description="Two line elements.",
        min_length=2,
        max_length=2,
        examples=[
            [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        ],
    )

    def get_catalog_number(self) -> int:
        """
        Gets the TLE catalog number.

        Returns:
            int: the catalog number
        """
        # pylint: disable=E1136
        return int(self.tle[0][2:7])

    def get_classification(self) -> str:
        """
        Gets the TLE classification type (U: unclassified; C: classified).

        Returns:
            str: the classification type
        """
        # pylint: disable=E1136
        return self.tle[0][7]

    def get_international_designator(self) -> str:
        """
        Gets the TLE international designator.

        Returns:
            str: the international designator
        """
        # pylint: disable=E1136
        year = self.tle[0][9:11]
        launch = self.tle[0][11:14]
        piece = self.tle[0][14:17]
        return str(
            ("" if not year.isdigit() else "20" if int(year) < 57 else "19")
            + year
            + "-"
            + launch
            + piece
        ).strip()

    def get_epoch(self) -> datetime:
        """
        Gets the TLE epoch time.

        Returns:
            datetime: the epoch datetime
        """
        # pylint: disable=E1136
        return sat_epoch_datetime(Satrec.twoline2rv(self.tle[0], self.tle[1]))

    def get_first_derivative_mean_motion(self) -> float:
        """
        Gets the first derivative of mean motion (ballistic coefficient).

        Returns:
            float: the first derivative of mean motion
        """
        # pylint: disable=E1136
        return float(self.tle[0][33:43])

    def get_second_derivative_mean_motion(self) -> float:
        """
        Gets the second derivative of mean motion.

        Returns:
            float: the second derivative of mean motion
        """
        # pylint: disable=E1136
        return float("0." + self.tle[0][44:50].strip()) * 10 ** (
            int(self.tle[0][50:52])
        )

    def get_b_star(self) -> float:
        """
        Gets the b-star term (drag or radiation pressure coefficient).

        Returns:
            float: the b-star term
        """
        # pylint: disable=E1136
        return float("0." + self.tle[0][53:59].strip()) * 10 ** (
            int(self.tle[0][59:61])
        )

    def get_ephemeris_type(self) -> int:
        """
        Gets the TLE ephemeris type.

        Returns:
            int: the ephemeris type
        """
        # pylint: disable=E1136
        return int(self.tle[0][62])

    def get_element_set_number(self) -> int:
        """
        Gets the TLE element set number.

        Returns:
            int: the element set number
        """
        # pylint: disable=E1136
        return int(self.tle[0][64:68])

    def get_inclination(self) -> float:
        """
        Gets the orbit inclination (decimal degrees).

        Returns:
            float: the inclination
        """
        # pylint: disable=E1136
        return float(self.tle[1][8:16])

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node (decimal degrees).

        Returns:
            float: the right ascension of ascending node
        """
        # pylint: disable=E1136
        return float(self.tle[1][17:25])

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity.

        Returns:
            float: the eccentricity
        """
        # pylint: disable=E1136
        return float("0." + self.tle[1][26:33].strip())

    def get_perigee_argument(self) -> float:
        """
        Gets the argument of perigee (decimal degrees).

        Returns:
            float: the argument of perigee
        """
        # pylint: disable=E1136
        return float(self.tle[1][34:42])

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        # pylint: disable=E1136
        return float(self.tle[1][43:51])

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion (revolutions per day).

        Returns:
            float: the mean motion
        """
        # pylint: disable=E1136
        return float(self.tle[1][52:63])

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period.

        Returns:
            timedelta: the orbit period
        """
        return timedelta(days=1 / self.get_mean_motion())

    def get_revolution_number_at_epoch(self) -> int:
        """
        Gets the revolution number at epoch.

        Returns:
            timedelta: the revolution number
        """
        # pylint: disable=E1136
        return int(self.tle[1][63:68])

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters).

        Returns:
            float: the semimajor axis
        """
        mean_motion_rad_s = self.get_mean_motion() * 2 * np.pi / 86400
        return np.power(
            constants.EARTH_MU / mean_motion_rad_s**2,
            1 / 3,
        )

    def get_altitude(self) -> float:
        """
        Gets the altitude (meters).

        Returns:
            float: the altitude
        """
        return self.get_semimajor_axis() - constants.EARTH_MEAN_RADIUS

    def get_true_anomaly(self) -> float:
        """
        Gets the true anomaly (decimal degrees).

        Returns:
            float: the true anomaly
        """
        return utils.mean_anomaly_to_true_anomaly(
            self.get_mean_anomaly(), self.get_eccentricity()
        )

    @field_validator("tle")
    @classmethod
    def valid_tle(cls, v):
        """
        Validate the two line element set.
        """
        # based on orekit's TLE.isFormatOK function
        if len(v[0]) != 69:
            raise ValueError("Invalid tle: line 1 incorrect length.")
        if len(v[1]) != 69:
            raise ValueError("Invalid tle: line 2 incorrect length.")

        line_1_pattern = (
            r"1 [ 0-9A-HJ-NP-Z][ 0-9]{4}[A-Z] [ 0-9]{5}[ A-Z]{3} "
            + r"[ 0-9]{5}[.][ 0-9]{8} (?:(?:[ 0+-][.][ 0-9]{8})|(?: "
            + r"[ +-][.][ 0-9]{7})) [ +-][ 0-9]{5}[+-][ 0-9] "
            + r"[ +-][ 0-9]{5}[+-][ 0-9] [ 0-9] [ 0-9]{4}[ 0-9]"
        )
        if re.match(line_1_pattern, v[0]) is None:
            raise ValueError("Invalid tle: line 1 does not match pattern.")
        line_2_pattern = (
            r"2 [ 0-9A-HJ-NP-Z][ 0-9]{4} [ 0-9]{3}[.][ 0-9]{4} "
            + r"[ 0-9]{3}[.][ 0-9]{4} [ 0-9]{7} [ 0-9]{3}[.][ 0-9]{4} "
            + r"[ 0-9]{3}[.][ 0-9]{4} [ 0-9]{2}[.][ 0-9]{13}[ 0-9]"
        )
        if re.match(line_2_pattern, v[1]) is None:
            raise ValueError("Invalid tle: line 2 does not match pattern.")

        def checksum(line):
            the_sum = 0
            for i in range(68):
                if line[i].isdigit():
                    the_sum += int(line[i])
                elif line[i] == "-":
                    the_sum += 1
            return the_sum % 10

        if int(v[0][68]) != checksum(v[0]):
            raise ValueError("Invalid tle: line 1 checksum failed.")
        if int(v[1][68]) != checksum(v[1]):
            raise ValueError("Invalid tle: line 2 checksum failed.")
        return v

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> TwoLineElements:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            TwoLineElements: the derived orbit
        """
        # pylint: disable=E1136
        lead_tle = Satrec.twoline2rv(self.tle[0], self.tle[1])
        epoch = sat_epoch_datetime(lead_tle)
        satrec = Satrec()
        satrec.sgp4init(
            WGS72,
            "i",
            0,
            (epoch - datetime(1949, 12, 31, tzinfo=timezone.utc)) / timedelta(days=1),
            lead_tle.bstar,
            lead_tle.ndot,
            lead_tle.nddot,
            lead_tle.ecco,
            lead_tle.argpo,
            lead_tle.inclo,
            np.mod(lead_tle.mo + np.radians(delta_mean_anomaly), 2 * np.pi),
            lead_tle.no_kozai,
            np.mod(lead_tle.nodeo + np.radians(delta_raan), 2 * np.pi),
        )
        tle1, tle2 = exporter.export_tle(satrec)
        return TwoLineElements(tle=[tle1.replace("\x00", "U"), tle2])

    def as_skyfield(self):
        """
        Converts this orbit to a Skyfield `EarthSatellite`.

        Returns:
            skyfield.api.EarthSatellite: the Skyfield EarthSatellite
        """
        # pylint: disable=E1136
        return EarthSatellite(self.tle[0], self.tle[1])

    def get_repeat_cycle(
        self,
        max_delta_position: float = None,
        max_delta_velocity: float = None,
        min_elevation_angle: float = None,
        max_search_duration: timedelta = None,
        lazy_load: bool = None,
    ) -> timedelta:
        """
        Compute the orbit repeat cycle. Lazy-loads a previously-computed repeat cycle if available.

        Args:
            max_delta_position (float): the maximum difference in position (m) allowed for a repeat.
            max_delta_velocity (float): the maximum difference in velocity (m/s) allowed for a repeat.
            min_elevation_angle (float): the minimum elevation angle (deg) for screening repeats.
            max_search_duration (timedelta): the maximum period of time to search for repeats.
            lazy_load (bool): True, if the previously-computed repeat cycle should be loaded.

        Returns:
            timedelta: the repeat cycle duration (if it exists)
        """
        # load defaults
        if max_delta_position is None:
            max_delta_position = config.rc.repeat_cycle_delta_position_m
        if max_delta_velocity is None:
            max_delta_velocity = config.rc.repeat_cycle_delta_velocity_m_per_s
        if min_elevation_angle is None:
            min_elevation_angle = config.rc.repeat_cycle_search_elevation_deg
        if max_search_duration is None:
            max_search_duration = timedelta(
                days=config.rc.repeat_cycle_search_duration_days
            )
        if lazy_load is None:
            lazy_load = config.rc.repeat_cycle_lazy_load

        if lazy_load:
            repeat_cycle = self.__dict__.get("repeat_cycle")
        else:
            repeat_cycle = None
        if repeat_cycle is None:
            # extract the orbit epoch time
            epoch = self.get_epoch()
            # record the initial position and velocity in Earth-centered Earth-fixed frame
            datum = wgs84.subpoint_of(
                self.as_skyfield().at(constants.timescale.from_datetime(epoch))
            )
            position_0, velocity_0 = (
                self.as_skyfield()
                .at(constants.timescale.from_datetime(epoch))
                .frame_xyz_and_velocity(itrs)
            )
            # find candidate repeat events
            ts, es = self.as_skyfield().find_events(
                datum,
                constants.timescale.from_datetime(epoch + timedelta(minutes=10)),
                constants.timescale.from_datetime(epoch + max_search_duration),
                min_elevation_angle,
            )
            # compute position and velocity at culmination in Earth-centered Earth-fixed frame
            position, velocity = (
                self.as_skyfield().at(ts[es == 1]).frame_xyz_and_velocity(itrs)
            )
            # apply validity conditions on position and velocity error norms
            is_valid = np.logical_and(
                np.linalg.norm((position.m.T - position_0.m.T).T, axis=0)
                < max_delta_position,
                np.linalg.norm((velocity.m_per_s.T - velocity_0.m_per_s.T).T, axis=0)
                < max_delta_velocity,
            )
            if np.any(is_valid):
                # assign repeat cycle
                repeat_cycle = ts[es == 1][is_valid][0].utc_datetime() - epoch
            else:
                # assign zero repeat cycle value to avoid recalculation
                repeat_cycle = timedelta(0)
            self.__dict__["repeat_cycle"] = repeat_cycle
        if repeat_cycle > timedelta(0):
            return repeat_cycle
        return None

    def get_orbit_track(
        self, times: Union[datetime, List[datetime]], try_repeat: bool = None
    ) -> Geocentric:
        """
        Gets the orbit track of this orbit using Skyfield.

        Args:
            times (Union[datetime, List[datetime]]): time(s) at which to compute position/velocity.
            try_repeat (bool): True, if a repeat orbit should be used to improve long-term accuracy.

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        # load defaults
        if try_repeat is None:
            try_repeat = config.rc.repeat_cycle_for_orbit_track
        # create skyfield Time
        if isinstance(times, datetime):
            ts_times = constants.timescale.from_datetime(times)
        else:
            ts_times = constants.timescale.from_datetimes(times)
        if try_repeat:
            # try to compute repeat cycle positions
            repeat_cycle = self.get_repeat_cycle()
            if repeat_cycle is not None:
                epoch = self.get_epoch()
                if isinstance(times, datetime):
                    repeat_times = constants.timescale.from_datetime(
                        epoch + np.mod(times - epoch, repeat_cycle)
                    )
                else:
                    repeat_times = constants.timescale.from_datetimes(
                        epoch + np.mod(np.array(times) - epoch, repeat_cycle)
                    )
                repeat_track = self.as_skyfield().at(repeat_times)
                return Geocentric(
                    repeat_track.position.au, repeat_track.velocity.au_per_d, ts_times
                )
        # compute satellite positions
        return self.as_skyfield().at(ts_times)

    def get_observation_events(
        self,
        point: Point,
        start: datetime,
        end: datetime,
        min_elevation_angle: float,
        try_repeat: bool = None,
    ) -> tuple:
        """
        Gets the observation events of this orbit using Skyfield.

        Args:
            point (Point): Target location to observe.
            start (datetime): Start time of the observation period.
            end (datetime): End time of the observation period.
            min_elevation_angle (float): Minimum elevation angle (deg) to constrain observation.
            try_repeat (bool): True, if a repeat orbit should be used to improve long-term accuracy.

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        # load defaults
        if try_repeat is None:
            try_repeat = config.rc.repeat_cycle_for_observation_events
        # create skyfield Time
        t_0 = constants.timescale.from_datetime(start)
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        if try_repeat:
            # try to compute repeat cycle events
            repeat_cycle = self.get_repeat_cycle()
            if repeat_cycle is not None and repeat_cycle < end - start:
                repeat_t_1 = constants.timescale.from_datetime(start + repeat_cycle)
                times, events = self.as_skyfield().find_events(
                    topos, t_0, repeat_t_1, min_elevation_angle
                )
                number_cycles = int(np.ceil((end - start) / repeat_cycle))
                if len(times) == 0:
                    return (Time([], []), np.array([], dtype=int))
                times_py = np.concatenate(
                    [
                        times.utc_datetime() + i * repeat_cycle
                        for i in range(number_cycles)
                    ]
                )
                events_py = np.concatenate([events for _ in range(number_cycles)])
                if len(times_py) == 0:
                    return Time([], []), np.array([], dtype=int)
                return (
                    constants.timescale.from_datetimes(times_py[times_py <= end]),
                    events_py[times_py <= end],
                )
        # compute observation events
        t_1 = constants.timescale.from_datetime(end)
        # pylint: disable=E1101
        return self.as_skyfield().find_events(topos, t_0, t_1, min_elevation_angle)

    def to_tle(self) -> TwoLineElements:
        """
        Converts this orbit to a two line elements representation.

        Returns:
            TwoLineElements: the two line elements orbit
        """
        return self


class OrbitBase(BaseModel):
    """
    Base class for orbits.
    """

    altitude: float = Field(..., description="Mean altitude (meters).")
    true_anomaly: float = Field(0, description="True anomaly (degrees).", ge=0, lt=360)
    epoch: Optional[datetime] = Field(
        datetime.now(tz=timezone.utc),
        description="Timestamp (epoch) of the initial orbital state.",
    )

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.true_anomaly_to_mean_anomaly(self.true_anomaly)

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis (meters).

        Returns:
            float: the semimajor axis
        """
        return constants.EARTH_MEAN_RADIUS + self.altitude

    def get_mean_motion(self) -> float:
        """
        Gets the mean motion (revolutions per day).

        Returns:
            float: the mean motion
        """
        return 1 / (self.get_orbit_period() / timedelta(days=1))

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period.

        Returns:
            timedelta: the orbit period
        """
        return timedelta(seconds=utils.compute_orbit_period(self.altitude))


class CircularOrbit(OrbitBase):
    """
    Orbit specification using Keplerian elements for elliptical motion -- circular motion case.
    """

    type: Literal["circular"] = Field(
        "circular", description="Orbit type discriminator."
    )
    inclination: float = Field(0, description="Inclination (degrees).", ge=0, lt=180)
    right_ascension_ascending_node: float = Field(
        0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> CircularOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            CircularOrbit: the derived orbit
        """
        true_anomaly = utils.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360)
        )
        raan = np.mod(self.right_ascension_ascending_node + delta_raan, 360)
        return CircularOrbit(
            altitude=self.altitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
        )

    def to_tle(self, lazy_load: bool = None) -> TwoLineElements:
        """
        Converts this orbit to a two line elements representation.

        Args:
            lazy_load (bool): True, if this tle should be lazy-loaded.

        Returns:
            TwoLineElements: the two line elements orbit
        """
        if lazy_load is None:
            lazy_load = config.rc.orbit_tle_lazy_load
        if lazy_load:
            tle = self.__dict__.get("tle")
        else:
            tle = None
        if tle is None:
            tle = KeplerianOrbit(
                altitude=self.altitude,
                true_anomaly=self.true_anomaly,
                epoch=self.epoch,
                inclination=self.inclination,
                right_ascension_ascending_node=self.right_ascension_ascending_node,
            ).to_tle()
            self.__dict__["tle"] = tle
        return tle


class SunSynchronousOrbit(OrbitBase):
    """
    Orbit defined by sun synchronous parameters.
    """

    type: Literal["sso"] = Field("sso", description="Orbit type discriminator.")
    altitude: float = Field(
        ...,
        description="Mean altitude (meters).",
        ge=0,
        lt=12352000 - constants.EARTH_MEAN_RADIUS,
    )
    equator_crossing_time: time = Field(
        ..., description="Equator crossing time (local solar time)."
    )
    equator_crossing_ascending: bool = Field(
        True,
        description="True, if the equator crossing time is ascending (south-to-north).",
    )

    def get_inclination(self) -> float:
        """
        Gets the inclination (decimal degrees).

        Returns:
            float: the inclination
        """
        return np.degrees(
            np.arccos(-np.power(self.get_semimajor_axis() / 12352000, 7 / 2))
        )

    def get_right_ascension_ascending_node(self) -> float:
        """
        Gets the right ascension of ascending node (decimal degrees).

        Returns:
            float: the right ascension of ascending node
        """
        # pylint: disable=E1101
        ect_day = timedelta(
            hours=self.equator_crossing_time.hour,
            minutes=self.equator_crossing_time.minute,
            seconds=self.equator_crossing_time.second,
            microseconds=self.equator_crossing_time.microsecond,
        ) / timedelta(days=1)
        epoch_time = constants.timescale.from_datetime(self.epoch)
        sun = constants.de421["sun"]
        earth = constants.de421["earth"]
        right_ascension, _, _ = earth.at(epoch_time).observe(sun).radec()
        # pylint: disable=W0212
        return (
            right_ascension._degrees
            + 360 * ect_day
            + 180 * self.equator_crossing_ascending
        ) % 360

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> CircularOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            CircularOrbit: the derived orbit
        """
        true_anomaly = utils.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360)
        )
        raan = np.mod(self.get_right_ascension_ascending_node() + delta_raan, 360)
        # TODO the resulting orbit *is* still sun-synchronous; however, with a
        # different equator-crossing time. Need to do the math to determine.
        # For now, simply return a circular orbit with the correct parameters.
        return CircularOrbit(
            altitude=self.altitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.get_inclination(),
            right_ascension_ascending_node=raan,
        )

    def to_tle(self, lazy_load: bool = None) -> TwoLineElements:
        """
        Converts this orbit to a two line elements representation.

        Args:
            lazy_load (bool): True, if this tle should be lazy-loaded.

        Returns:
            TwoLineElements: the two line elements orbit
        """
        if lazy_load is None:
            lazy_load = config.rc.orbit_tle_lazy_load
        if lazy_load:
            tle = self.__dict__.get("tle")
        else:
            tle = None
        if tle is None:
            tle = KeplerianOrbit(
                altitude=self.altitude,
                inclination=self.get_inclination(),
                right_ascension_ascending_node=self.get_right_ascension_ascending_node(),
                true_anomaly=self.true_anomaly,
                epoch=self.epoch,
            ).to_tle()
            self.__dict__["tle"] = tle
        return tle


class KeplerianOrbit(CircularOrbit):
    """
    Orbit specification using Keplerian elements for elliptical motion.
    """

    type: Literal["keplerian"] = Field(
        "keplerian", description="Orbit type discriminator."
    )
    eccentricity: float = Field(0, description="Eccentricity.", ge=0)
    perigee_argument: float = Field(
        0, description="Perigee argument (degrees).", ge=0, lt=360
    )

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.true_anomaly_to_mean_anomaly(self.true_anomaly, self.eccentricity)

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> KeplerianOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            KeplerianOrbit: the derived orbit
        """
        true_anomaly = utils.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360),
            eccentricity=self.eccentricity,
        )
        raan = np.mod(self.right_ascension_ascending_node + delta_raan, 360)
        return KeplerianOrbit(
            altitude=self.altitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
            eccentricity=self.eccentricity,
            perigee_argument=self.perigee_argument,
        )

    def to_tle(self, lazy_load: bool = None) -> TwoLineElements:
        """
        Converts this orbit to a two line elements representation.

        Args:
            lazy_load (bool): True, if this tle should be lazy-loaded.

        Returns:
            TwoLineElements: the two line elements orbit
        """
        if lazy_load is None:
            lazy_load = config.rc.orbit_tle_lazy_load
        if lazy_load:
            tle = self.__dict__.get("tle")
        else:
            tle = None
        if tle is None:
            satrec = Satrec()
            satrec.sgp4init(
                WGS72,
                "i",
                0,
                (self.epoch - datetime(1949, 12, 31, tzinfo=timezone.utc))
                / timedelta(days=1),
                0,
                0.0,
                0.0,
                self.eccentricity,
                np.radians(self.perigee_argument),
                np.radians(self.inclination),
                np.radians(self.get_mean_anomaly()),
                self.get_mean_motion() * 2 * np.pi / 1440,
                np.radians(self.right_ascension_ascending_node),
            )
            tle1, tle2 = exporter.export_tle(satrec)
            tle = TwoLineElements(tle=[tle1.replace("\x00", "U"), tle2])
            self.__dict__["tle"] = tle
        return tle


class MolniyaOrbit(OrbitBase):
    """
    Orbit defined by Molniya parameters. The altitude parameter is considered the altitude at perigee.
    """

    type: Literal["molniya"] = Field("molniya", description="Orbit type discriminator.")
    inclination: float = Field(0, description="Inclination (degrees).", ge=0, lt=180)
    right_ascension_ascending_node: float = Field(
        0, description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )
    orbital_period: float = Field(43080, description="Orbital period (seconds).", gt=0)
    perigee_argument: float = Field(
        0, description="Perigee argument (degrees).", ge=0, lt=360
    )

    def get_semimajor_axis(self) -> float:
        gravitational_constant = 6.6743e-11
        earth_mass = 5.9736e24

        return np.cbrt(
            gravitational_constant
            * earth_mass
            * (self.orbital_period**2)
            / (4 * (np.pi**2))
        )

    def get_eccentricity(self) -> float:
        """
        Gets the eccentricity (float between 0 and 1).

        Returns:
            float: the eccentricity
        """
        semimajor_axis = self.get_semimajor_axis()
        return (
            semimajor_axis - (self.altitude + constants.EARTH_EQUATORIAL_RADIUS)
        ) / semimajor_axis

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.true_anomaly_to_mean_anomaly(
            self.true_anomaly, self.get_eccentricity()
        )

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> MolniyaOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            Molniya Orbit: the derived orbit
        """
        true_anomaly = utils.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360),
            eccentricity=self.get_eccentricity(),
        )
        raan = np.mod(self.right_ascension_ascending_node + delta_raan, 360)
        return MolniyaOrbit(
            altitude=self.altitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
            eccentricity=self.get_eccentricity(),
            perigee_argument=self.perigee_argument,
        )

    def to_tle(self, lazy_load: bool = None) -> TwoLineElements:
        """
        Converts this orbit to a two line elements representation.

        Args:
            lazy_load (bool): True, if this tle should be lazy-loaded.

        Returns:
            TwoLineElements: the two line elements orbit
        """
        if lazy_load is None:
            lazy_load = config.rc.orbit_tle_lazy_load
        if lazy_load:
            tle = self.__dict__.get("tle")
        else:
            tle = None
        if tle is None:
            satrec = Satrec()
            satrec.sgp4init(
                WGS72,
                "i",
                0,
                (self.epoch - datetime(1949, 12, 31, tzinfo=timezone.utc))
                / timedelta(days=1),
                0,
                0.0,
                0.0,
                self.get_eccentricity(),
                np.radians(self.perigee_argument),
                np.radians(self.inclination),
                np.radians(self.get_mean_anomaly()),
                2 * np.pi / (self.orbital_period / 60),
                np.radians(self.right_ascension_ascending_node),
            )
            tle1, tle2 = exporter.export_tle(satrec)
            tle = TwoLineElements(tle=[tle1.replace("\x00", "U"), tle2])
            self.__dict__["tle"] = tle
        return tle


class TundraOrbit(MolniyaOrbit):
    """
    Orbit defined by Tundra parameters, inherits the Molniya class.
    """

    type: Literal["tundra"] = Field("tundra", description="Orbit type discriminator.")
    inclination: float = Field(63.4, description="Inclination (degrees).", ge=0, lt=180)
    orbital_period: float = Field(86400, description="Orbital period (seconds).", gt=0)
    perigee_argument: float = Field(
        270, description="Perigee argument (degrees).", ge=0, lt=360
    )
    eccentricity: float = Field(0.2, description="Eccentricity.", ge=0)

    def get_eccentricity(self):
        return self.eccentricity

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> TundraOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            Tundra Orbit: the derived orbit
        """
        true_anomaly = utils.mean_anomaly_to_true_anomaly(
            np.mod(self.get_mean_anomaly() + delta_mean_anomaly, 360),
            eccentricity=self.get_eccentricity(),
        )
        raan = np.mod(self.right_ascension_ascending_node + delta_raan, 360)
        return TundraOrbit(
            altitude=self.altitude,
            true_anomaly=true_anomaly,
            epoch=self.epoch,
            inclination=self.inclination,
            right_ascension_ascending_node=raan,
            eccentricity=self.eccentricity,
            perigee_argument=self.perigee_argument,
        )
