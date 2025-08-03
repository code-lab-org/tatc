# -*- coding: utf-8 -*-
"""
Object schemas for satellite orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, time, timedelta, timezone
from typing import List, Optional, Tuple, Union

import numpy as np
from astropy.time import Time
from pydantic import AfterValidator, BaseModel, Field
from sgp4.api import Satrec, WGS72
from sgp4 import exporter
from sgp4.conveniences import sat_epoch_datetime
from skyfield.api import EarthSatellite, Time, wgs84
from skyfield.positionlib import Geocentric
from skyfield.framelib import itrs
from typing_extensions import Annotated, Literal

from .point import Point
from .. import config, constants, utils


class TwoLineElements(BaseModel):
    """
    Orbit defined with standard two line elements.
    """

    type: Literal["tle"] = Field("tle", description="Orbit type discriminator.")
    tle: Annotated[
        List[str],
        AfterValidator(utils.is_chronological_tle),
        AfterValidator(utils.is_valid_tle),
        AfterValidator(utils.is_even_length_list),
    ] = Field(
        ...,
        description="Two line elements. Multiple TLEs must be in chronological order.",
        min_length=2,
        examples=[
            [
                "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
                "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
            ]
        ],
    )

    def get_tle_count(self) -> int:
        """
        Gets the number of TLEs specified for this orbit.

        Returns:
            int: the number of TLEs
        """
        return len(self.tle) // 2

    def get_catalog_number(self, tle_i=0) -> int:
        """
        Gets the TLE catalog number.

        Args:
            tle_i (int): the TLE index.

        Returns:
            int: the catalog number
        """
        # pylint: disable=E1136
        return int(self.tle[2 * tle_i][2:7])

    def get_classification(self, tle_i=0) -> str:
        """
        Gets the TLE classification type (U: unclassified; C: classified).

        Args:
            tle_i (int): the TLE index.

        Returns:
            str: the classification type
        """
        # pylint: disable=E1136
        return self.tle[2 * tle_i][7]

    def get_international_designator(self, tle_i=0) -> str:
        """
        Gets the TLE international designator.

        Args:
            tle_i (int): the TLE index.

        Returns:
            str: the international designator
        """
        # pylint: disable=E1136
        year = self.tle[2 * tle_i][9:11]
        launch = self.tle[2 * tle_i][11:14]
        piece = self.tle[2 * tle_i][14:17]
        return str(
            ("" if not year.isdigit() else "20" if int(year) < 57 else "19")
            + year
            + "-"
            + launch
            + piece
        ).strip()

    def get_epoch(self, tle_i=0) -> datetime:
        """
        Gets the TLE epoch time.

        Args:
            tle_i (int): the TLE index.

        Returns:
            datetime: the epoch datetime
        """
        # pylint: disable=E1136
        return sat_epoch_datetime(
            Satrec.twoline2rv(self.tle[2 * tle_i], self.tle[2 * tle_i + 1])
        )

    def get_first_derivative_mean_motion(self, tle_i=0) -> float:
        """
        Gets the first derivative of mean motion (ballistic coefficient).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the first derivative of mean motion
        """
        # pylint: disable=E1136
        return float(self.tle[2 * tle_i][33:43])

    def get_second_derivative_mean_motion(self, tle_i=0) -> float:
        """
        Gets the second derivative of mean motion.

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the second derivative of mean motion
        """
        # pylint: disable=E1136
        return float("0." + self.tle[2 * tle_i][44:50].strip()) * 10 ** (
            int(self.tle[2 * tle_i][50:52])
        )

    def get_b_star(self, tle_i=0) -> float:
        """
        Gets the b-star term (drag or radiation pressure coefficient).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the b-star term
        """
        # pylint: disable=E1136
        return float("0." + self.tle[2 * tle_i][53:59].strip()) * 10 ** (
            int(self.tle[2 * tle_i][59:61])
        )

    def get_ephemeris_type(self, tle_i=0) -> int:
        """
        Gets the TLE ephemeris type.

        Args:
            tle_i (int): the TLE index.

        Returns:
            int: the ephemeris type
        """
        # pylint: disable=E1136
        return int(self.tle[2 * tle_i][62])

    def get_element_set_number(self, tle_i=0) -> int:
        """
        Gets the TLE element set number.

        Args:
            tle_i (int): the TLE index.

        Returns:
            int: the element set number
        """
        # pylint: disable=E1136
        return int(self.tle[2 * tle_i][64:68])

    def get_inclination(self, tle_i=0) -> float:
        """
        Gets the orbit inclination (decimal degrees).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the inclination
        """
        # pylint: disable=E1136
        return float(self.tle[2 * tle_i + 1][8:16])

    def get_right_ascension_ascending_node(self, tle_i=0) -> float:
        """
        Gets the right ascension of ascending node (decimal degrees).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the right ascension of ascending node
        """
        # pylint: disable=E1136
        return float(self.tle[2 * tle_i + 1][17:25])

    def get_eccentricity(self, tle_i=0) -> float:
        """
        Gets the eccentricity.

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the eccentricity
        """
        # pylint: disable=E1136
        return float("0." + self.tle[2 * tle_i + 1][26:33].strip())

    def get_perigee_argument(self, tle_i=0) -> float:
        """
        Gets the argument of perigee (decimal degrees).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the argument of perigee
        """
        # pylint: disable=E1136
        return float(self.tle[2 * tle_i + 1][34:42])

    def get_mean_anomaly(self, tle_i=0) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the mean anomaly
        """
        # pylint: disable=E1136
        return float(self.tle[2 * tle_i + 1][43:51])

    def get_mean_motion(self, tle_i=0) -> float:
        """
        Gets the mean motion (revolutions per day).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the mean motion
        """
        # pylint: disable=E1136
        return float(self.tle[2 * tle_i + 1][52:63])

    def get_orbit_period(self, tle_i=0) -> timedelta:
        """
        Gets the approximate orbit period.

        Args:
            tle_i (int): the TLE index.

        Returns:
            timedelta: the orbit period
        """
        return timedelta(days=1 / self.get_mean_motion(tle_i))

    def get_revolution_number_at_epoch(self, tle_i=0) -> int:
        """
        Gets the revolution number at epoch.

        Args:
            tle_i (int): the TLE index.

        Returns:
            timedelta: the revolution number
        """
        # pylint: disable=E1136
        return int(self.tle[2 * tle_i + 1][63:68])

    def get_semimajor_axis(self, tle_i=0) -> float:
        """
        Gets the semimajor axis (meters).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the semimajor axis
        """
        mean_motion_rad_s = self.get_mean_motion(tle_i) * 2 * np.pi / 86400
        return np.power(
            constants.EARTH_MU / mean_motion_rad_s**2,
            1 / 3,
        )

    def get_altitude(self, tle_i=0) -> float:
        """
        Gets the altitude (meters).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the altitude
        """
        return self.get_semimajor_axis(tle_i) - constants.EARTH_MEAN_RADIUS

    def get_true_anomaly(self, tle_i=0) -> float:
        """
        Gets the true anomaly (decimal degrees).

        Args:
            tle_i (int): the TLE index.

        Returns:
            float: the true anomaly
        """
        return utils.mean_anomaly_to_true_anomaly(
            self.get_mean_anomaly(tle_i), self.get_eccentricity(tle_i)
        )

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
        # pylint: disable=E1101
        tles = self.tle.copy()
        for i in range(self.get_tle_count()):
            # pylint: disable=E1136
            lead_tle = Satrec.twoline2rv(self.tle[2 * i], self.tle[2 * i + 1])
            epoch = sat_epoch_datetime(lead_tle)
            satrec = Satrec()
            satrec.sgp4init(
                WGS72,
                "i",
                0,
                (epoch - datetime(1949, 12, 31, tzinfo=timezone.utc))
                / timedelta(days=1),
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
            tles[2 * i] = tle1.replace("\x00", "U")
            tles[2 * i + 1] = tle2
        return TwoLineElements(tle=tles)

    def as_skyfield(self, tle_i=0):
        """
        Converts this orbit to a Skyfield `EarthSatellite`.

        Args:
            tle_i (int): the TLE index.

        Returns:
            skyfield.api.EarthSatellite: the Skyfield EarthSatellite
        """
        # pylint: disable=E1136
        return EarthSatellite(self.tle[2 * tle_i], self.tle[2 * tle_i + 1])

    def get_closest_tle_index(
        self, at_times: Union[datetime, List[datetime]]
    ) -> Union[int, List[int]]:
        """
        Gets the closest TLE index to specified time(s).

        Args:
            at_times (datetime or List[datetime]): specified times

        Returns:
            int or List[int]: closest TLE index or indices
        """

        if at_times is None:
            return 0
        # lazy-load epochs
        tle_epochs = self.__dict__.get("tle_epochs")
        if tle_epochs is None:
            # extract the orbit epoch time
            tle_epochs = np.array(
                [self.get_epoch(i) for i in range(self.get_tle_count())]
            )
            self.__dict__["tle_epochs"] = tle_epochs
        # handle scalar
        if isinstance(at_times, datetime):
            idx = np.searchsorted(tle_epochs, at_times, side="left")
            return (
                int(idx - 1)
                if idx > 0
                and (
                    idx == len(tle_epochs)
                    or abs(at_times - tle_epochs[idx - 1])
                    < abs(at_times - tle_epochs[idx])
                )
                else int(idx)
            )
        # handle vector
        indices = np.searchsorted(tle_epochs, at_times, side="left")
        return [
            (
                int(idx - 1)
                if idx > 0
                and (
                    idx == len(tle_epochs)
                    or abs(at_times[i] - tle_epochs[idx - 1])
                    < abs(at_times[i] - tle_epochs[idx])
                )
                else int(idx)
            )
            for i, idx in enumerate(indices)
        ]

    def partition_by_tle_index(
        self, start: datetime, end: datetime
    ) -> Tuple[List[datetime], List[int]]:
        """
        Partition a timeline based on closest TLE index.

        Args:
            start (datetime): Start time.
            end (datetime): End time.

        Returns:
            Tuple[List[datetime], List[int]]: list of partitioned times and assigned TLE indices
        """
        if self.get_tle_count() <= 1:
            return [start, end], [0, 0]
        epochs = [self.get_epoch(i) for i in range(self.get_tle_count())]
        sorted_epochs = np.sort(epochs)
        tle_i = np.argsort(epochs)
        epoch_midpoints = (
            sorted_epochs[1:] + (sorted_epochs[:-1] - sorted_epochs[1:]) / 2
        )
        midpoint_indices = [i for i, t in enumerate(epoch_midpoints) if start < t < end]
        return (
            [start] + list(epoch_midpoints[midpoint_indices]) + [end],
            (
                [self.get_closest_tle_index(start)]
                + list(map(int, tle_i[midpoint_indices]))
                + [self.get_closest_tle_index(end)]
            ),
        )

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

        if self.get_tle_count() > 1:
            # try to use use multiple TLEs
            if isinstance(times, datetime):
                return self.as_skyfield(self.get_closest_tle_index(times)).at(
                    constants.timescale.from_datetime(times)
                )
            tle_i = self.get_closest_tle_index(times)
            tracks = [
                self.as_skyfield(i).at(constants.timescale.from_datetime(t))
                for i, t in zip(tle_i, times)
            ]
            return Geocentric(
                np.array([track.position.au for track in tracks]).T,
                np.array([track.velocity.au_per_d for track in tracks]).T,
                constants.timescale.from_datetimes(times),
            )
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
                    offset = times - epoch
                    repeat_times = constants.timescale.from_datetime(
                        epoch
                        + np.sign(offset / timedelta(1))
                        * np.mod(np.abs(offset), repeat_cycle)
                    )
                else:
                    offset = np.array(times) - epoch
                    repeat_times = constants.timescale.from_datetimes(
                        epoch
                        + np.sign(offset / timedelta(1))
                        * np.mod(np.abs(offset), repeat_cycle)
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
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        if self.get_tle_count() > 1:
            # try to use use multiple TLEs
            part_ts, tle_is = self.partition_by_tle_index(start, end)
            events = [
                self.as_skyfield(tle_is[i]).find_events(
                    topos,
                    constants.timescale.from_datetime(part_ts[i]),
                    constants.timescale.from_datetime(part_ts[i + 1]),
                    min_elevation_angle,
                )
                for i in range(len(part_ts) - 1)
            ]
            return (
                constants.timescale.from_datetimes(
                    [t.utc_datetime() for e in events for t in e[0]]
                ),
                np.array([v for e in events for v in e[1]]),
            )
        # create skyfield Time
        t_0 = constants.timescale.from_datetime(start)
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

class GeosynchronousOrbit(CircularOrbit):
    """
    Geosynchronous orbit defined by longitude to be observed.
    Geostationary orbit can be specified with inclination equal to 0.
    """

    type: Literal["geosynchronous"] = Field(
        "geosynchronous", description="Orbit type discriminator."
    )
    altitude: float = Field(35786000, description="Mean altitude (meters).")
    longitude: float = Field(0, description="Longitude (degrees).", ge=0, lt=360)

    def get_true_anomaly(self) -> float:
        """
        Gets the true anomaly for the satellite over a given longitude at a given epoch.

        Returns:
            float: True anomaly in degrees.
        """
        # convert to astropy time
        t = Time(self.epoch, scale="utc")
        # get greenwich sidereal time in degrees
        gst = t.sidereal_time('mean', 'greenwich').deg
        # return earth-centered inertial angle
        return (self.longitude + gst) % 360

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
            tle = CircularOrbit(
                altitude=self.altitude,
                inclination=self.inclination,
                right_ascension_ascending_node=self.right_ascension_ascending_node,
                true_anomaly=self.compute_true_anomaly(),
                epoch=self.epoch,
            ).to_tle()
            self.__dict__["tle"] = tle
        return tle