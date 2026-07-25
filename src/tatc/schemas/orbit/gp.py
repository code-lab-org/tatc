"""
Object schemas for general perturbations orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import BaseModel, Field
from sgp4 import exporter, omm
from sgp4.api import WGS72, Satrec
from sgp4.conveniences import sat_epoch_datetime
from sgp4.io import twoline2rv
from skyfield.api import EarthSatellite, Time, wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from typing_extensions import Literal

from ... import config, constants, utils
from ..point import Point


class GeneralPerturbationsElements(BaseModel):
    """General perturbations orbital elements for a satellite."""

    object_name: str | None = Field(None, description="Object name.")
    epoch: datetime = Field(..., description="Epoch.")
    mean_motion: float = Field(..., description="Mean motion (degrees/second).", gt=0)
    eccentricity: float = Field(..., description="Eccentricity.", ge=0, le=1)
    inclination: float = Field(..., description="Inclination (degrees).", ge=-90, le=90)
    ra_of_asc_node: float = Field(..., description="Right ascension of ascending node (degrees).", ge=0, lt=360)
    arg_of_pericenter: float = Field(..., description="Argument of pericenter (degrees).", ge=0, lt=360)
    mean_anomaly: float = Field(..., description="Mean anomaly (degrees).", ge=0, lt=360)
    norad_cat_id: int = Field(0, description="NORAD catalog identifier.", ge=0)
    bstar: float = Field(0, description="Starred ballistic coefficient.")
    mean_motion_dot: float = Field(0, description="First derivative of mean motion (degrees/second^2).")
    mean_motion_ddot: float = Field(0, description="Second derivative of mean motion (degrees/second^3).")

    @classmethod
    def from_satrec(cls, satrec: Satrec) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from a Satrec object.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        return GeneralPerturbationsElements(
            object_name=satrec.satname,
            epoch=sat_epoch_datetime(satrec),
            mean_motion=np.degrees(satrec.no_kozai) / 60,
            eccentricity=satrec.ecco,
            inclination=np.degrees(satrec.inclo),
            ra_of_asc_node=np.degrees(satrec.nodeo),
            arg_of_pericenter=np.degrees(satrec.argpo),
            mean_anomaly=np.degrees(satrec.mo),
            norad_cat_id=satrec.satnum,
            bstar=satrec.bstar,
            mean_motion_dot=np.degrees(satrec.ndot) / 60**2,
            mean_motion_ddot=np.degrees(satrec.nddot) / 60**3
        )
    
    def to_satrec(self) -> Satrec:
        """
        Converts this GP elements object to a Satrec object.

        Returns:
            Satrec: the Satrec object
        """
        satrec = Satrec()
        satrec.sgp4init(
            WGS72,
            "i",
            self.norad_cat_id,
            (self.epoch - datetime(1949, 12, 31, tzinfo=timezone.utc)) / timedelta(days=1),
            self.bstar,
            np.radians(self.mean_motion_dot) * 60**2,
            np.radians(self.mean_motion_ddot) * 60**3,
            self.eccentricity,
            np.radians(self.arg_of_pericenter),
            np.radians(self.inclination),
            np.radians(self.mean_anomaly),
            np.radians(self.mean_motion) * 60,
            np.radians(self.ra_of_asc_node),
        )
        return satrec

    def get_orbit_period(self) -> timedelta:
        """
        Gets the approximate orbit period.

        Returns:
            timedelta: the orbit period
        """
        return timedelta(days=1 / self.mean_motion)

    def get_semimajor_axis(self) -> float:
        """
        Gets the semimajor axis.

        Returns:
            float: the semimajor axis (meters)
        """

        return utils.orbital.mean_motion_to_semimajor_axis(self.mean_motion)

    def get_mean_altitude(self) -> float:
        """
        Gets the mean altitude.

        Returns:
            float: the mean altitude (meters)
        """
        return self.get_semimajor_axis() - constants.EARTH_MEAN_RADIUS

    def get_true_anomaly(self) -> float:
        """
        Gets the true anomaly.

        Returns:
            float: the true anomaly (degrees)
        """
        return utils.orbital.mean_anomaly_to_true_anomaly(
            self.mean_anomaly, self.eccentricity
        )

    @classmethod
    def from_tle(cls, tle_lines: list[str]) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from two line element (TLE) lines.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        return GeneralPerturbationsElements.from_satrec(
            twoline2rv(tle_lines[0], tle_lines[1])
        )

    def to_tle(self) -> list[str]:
        """
        Converts this GP elements object to a two line element (TLE) representation.

        Returns:
            list[str]: the two line elements
        """
        return exporter.export_tle(self.to_satrec())
    
    @classmethod
    def from_omm_dict(cls, omm_dict: dict) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from an OMM dictionary.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        satrec = Satrec()
        omm.initialize(satrec, omm_dict)
        return GeneralPerturbationsElements.from_satrec(satrec)

    def to_omm_dict(self) -> dict:
        """
        Converts this GP elements object to an OMM dictionary.

        Returns:
            dict: the OMM dictionary
        """
        return exporter.export_omm(self.to_satrec(), self.object_name)
    
    @classmethod
    def from_omm_csv(cls, omm_csv: list[str]) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from OMM CSV lines.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        for fields in csv.DictReader(omm_csv):
            return GeneralPerturbationsElements.from_omm_dict(fields)
    
    @classmethod
    def from_omm_json(cls, omm_json: str) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from OMM JSON lines.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        for fields in json.loads(omm_json):
            return GeneralPerturbationsElements.from_omm_dict(fields)

    def to_skyfield(self):
        """
        Converts this GP elements object to a Skyfield `EarthSatellite`.

        Returns:
            skyfield.api.EarthSatellite: the Skyfield EarthSatellite
        """
        return EarthSatellite.from_omm(constants.timescale, self.to_omm_dict())

class GeneralPerturbationsOrbit(BaseModel):
    """
    Orbit defined with general perturbations (GP) elements.
    """

    type: Literal["gp"] = Field("gp", description="Orbit type discriminator.")
    elements: list[GeneralPerturbationsElements] = Field(..., description="General perturbations elements.")

    def get_mean_altitude(self, index: int = 0) -> float:
        """
        Gets the mean altitude of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the mean altitude (meters)
        """
        return self.elements[index].get_mean_altitude()

    def get_inclination(self, index: int = 0) -> float:
        """
        Gets the inclination of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the inclination (degrees)
        """
        return self.elements[index].inclination
    
    def get_true_anomaly(self, index: int = 0) -> float:
        """
        Gets the true anomaly of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the true anomaly (degrees)
        """
        return self.elements[index].get_true_anomaly()
    
    def get_right_ascension_ascending_node(self, index: int = 0) -> float:
        """
        Gets the right ascension of ascending node of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the right ascension of ascending node (degrees)
        """
        return self.elements[index].ra_of_asc_node
    
    def get_perigee_argument(self, index: int = 0) -> float:
        """
        Gets the argument of perigee of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the argument of perigee (degrees)
        """
        return self.elements[index].arg_of_pericenter

    @classmethod
    def from_tle(cls, tle_lines: list[str]) -> GeneralPerturbationsOrbit:
        """
        Creates a GP orbit from two line element (TLE) lines.

        Returns:
            GeneralPerturbationsOrbit: the GP orbit
        """
        return GeneralPerturbationsOrbit(
            elements=[
                GeneralPerturbationsElements.from_tle(tle_lines[i:i+2])
                for i in range(0, len(tle_lines), 2)
            ]
        )

    @classmethod
    def from_omm_csv(cls, omm_csv: list[str]) -> GeneralPerturbationsOrbit:
        """
        Creates a GP orbit from a OMM CSV lines.

        Returns:
            GeneralPerturbationsOrbit: the GP orbit
        """
        return GeneralPerturbationsOrbit(
            elements=[
                GeneralPerturbationsElements.from_omm_dict(fields)
                for fields in csv.DictReader(omm_csv)
            ]
        )

    @classmethod
    def from_omm_json(cls, omm_json: str) -> GeneralPerturbationsOrbit:
        """
        Creates a GP orbit from OMM JSON lines.

        Returns:
            GeneralPerturbationsOrbit: the GP orbit
        """
        return GeneralPerturbationsOrbit(
            elements=[
                GeneralPerturbationsElements.from_omm_dict(fields)
                for fields in json.loads(omm_json)
            ]
        )

    def get_element_epochs(self) -> list[datetime]:
        """
        Lazy-loads the epoch times for all elements in this orbit.

        Returns:
            list[datetime]: the epoch times
        """
        # lazy-load epochs
        element_epochs = self.__dict__.get("element_epochs")
        if element_epochs is None:
            # extract the element epoch times
            element_epochs = np.array(
                [el.epoch for el in self.elements]
            )
            self.__dict__["element_epochs"] = element_epochs
        return element_epochs

    def get_derived_orbit(
        self, delta_mean_anomaly: float, delta_raan: float
    ) -> GeneralPerturbationsOrbit:
        """
        Gets a derived orbit with perturbations to the mean anomaly and right
        ascension of ascending node.

        Args:
            delta_mean_anomaly (float):  Delta mean anomaly (degrees).
            delta_raan (float):  Delta right ascension of ascending node (degrees).

        Returns:
            GeneralPerturbationsOrbit: the derived orbit
        """
        derived_elements = []
        for original_el in self.elements:
            derived_el = original_el.model_copy(deep=True)
            derived_el.mean_anomaly = np.mod(derived_el.mean_anomaly + delta_mean_anomaly, 360)
            derived_el.ra_of_asc_node = np.mod(derived_el.ra_of_asc_node + delta_raan, 360)
            derived_elements.append(derived_el)
        return GeneralPerturbationsOrbit(elements=derived_elements)

    def get_closest_element_index(
        self, at_times: datetime | list[datetime]
    ) -> int | list[int]:
        """
        Gets the closest element index to specified time(s).

        Args:
            at_times (datetime | list[datetime]): specified times

        Returns:
            int | list[int]: closest element index or indices
        """

        if at_times is None:
            return 0
        # lazy-load element epochs
        element_epochs = self.get_element_epochs()
        # handle scalar
        if isinstance(at_times, datetime):
            idx = np.searchsorted(element_epochs, at_times, side="left")
            return (
                int(idx - 1)
                if idx > 0
                and (
                    idx == len(element_epochs)
                    or abs(at_times - element_epochs[idx - 1])
                    < abs(at_times - element_epochs[idx])
                )
                else int(idx)
            )
        # handle vector
        indices = np.searchsorted(element_epochs, at_times, side="left")
        return [
            (
                int(idx - 1)
                if idx > 0
                and (
                    idx == len(element_epochs)
                    or abs(at_times[i] - element_epochs[idx - 1])
                    < abs(at_times[i] - element_epochs[idx])
                )
                else int(idx)
            )
            for i, idx in enumerate(indices)
        ]

    def partition_by_element_index(
        self, start: datetime, end: datetime
    ) -> tuple[list[datetime], list[int]]:
        """
        Partition a timeline based on closest element index.

        Args:
            start (datetime): Start time.
            end (datetime): End time.

        Returns:
            tuple[list[datetime], list[int]]: list of partitioned times and assigned element indices
        """
        if len(self.elements) <= 1:
            return [start, end], [0, 0]
        element_epochs = self.get_element_epochs()
        sorted_epochs = np.sort(element_epochs)
        element_indices = np.argsort(element_epochs)
        epoch_midpoints = (
            sorted_epochs[1:] + (sorted_epochs[:-1] - sorted_epochs[1:]) / 2
        )
        midpoint_indices = [i for i, t in enumerate(epoch_midpoints) if start < t < end]
        return (
            [start] + list(epoch_midpoints[midpoint_indices]) + [end],
            (
                [self.get_closest_element_index(start)]
                + list(map(int, element_indices[midpoint_indices]))
                + [self.get_closest_element_index(end)]
            ),
        )

    def get_repeat_cycle(
        self,
        max_delta_position: float | None = None,
        max_delta_velocity: float | None = None,
        min_elevation_angle: float | None = None,
        max_search_duration: timedelta | None = None,
        lazy_load: bool | None = None,
    ) -> timedelta:
        """
        Compute the orbit repeat cycle. Lazy-loads a previously-computed repeat cycle if available.

        Args:
            max_delta_position (float | None): the maximum difference in position (m) allowed for a repeat.
            max_delta_velocity (float | None): the maximum difference in velocity (m/s) allowed for a repeat.
            min_elevation_angle (float | None): the minimum elevation angle (deg) for screening repeats.
            max_search_duration (timedelta | None): the maximum period of time to search for repeats.
            lazy_load (bool | None): True, if the previously-computed repeat cycle should be loaded.

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
        if repeat_cycle is None and len(self.elements) > 0:
            # extract the orbit epoch time from the first element
            epoch = self.elements[0].epoch
            satellite = self.elements[0].to_skyfield()
            # record the initial position and velocity in Earth-centered Earth-fixed frame
            datum = wgs84.subpoint_of(
                satellite.at(constants.timescale.from_datetime(epoch))
            )
            position_0, velocity_0 = (
                satellite
                .at(constants.timescale.from_datetime(epoch))
                .frame_xyz_and_velocity(itrs)
            )
            # find candidate repeat events
            ts, es = satellite.find_events(
                datum,
                constants.timescale.from_datetime(epoch + timedelta(minutes=10)),
                constants.timescale.from_datetime(epoch + max_search_duration),
                min_elevation_angle,
            )
            # compute position and velocity at culmination in Earth-centered Earth-fixed frame
            position, velocity = (
                satellite.at(ts[es == 1]).frame_xyz_and_velocity(itrs)
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
        self, times: datetime | list[datetime], try_repeat: bool | None = None
    ) -> Geocentric:
        """
        Gets the orbit track of this orbit using Skyfield.

        Args:
            times (datetime | list[datetime]): time(s) at which to compute position/velocity.
            try_repeat (bool | None): True, if a repeat orbit should be used to improve long-term accuracy.

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        # load defaults
        if try_repeat is None:
            try_repeat = config.rc.repeat_cycle_for_orbit_track

        if len(self.elements) > 1:
            # try to use use multiple TLEs
            if isinstance(times, datetime):
                nearest_index = self.get_closest_element_index(times)
                return self.elements[nearest_index].to_skyfield().at(
                    constants.timescale.from_datetime(times)
                )
            nearest_indices = self.get_closest_element_index(times)
            tracks = [
                self.elements[i].to_skyfield().at(constants.timescale.from_datetime(t))
                for i, t in zip(nearest_indices, times)
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
                repeat_track = self.elements[0].to_skyfield().at(repeat_times)
                return Geocentric(
                    repeat_track.position.au, repeat_track.velocity.au_per_d, ts_times
                )
        # compute satellite positions
        return self.elements[0].to_skyfield().at(ts_times)

    def get_observation_events(
        self,
        point: Point,
        start: datetime,
        end: datetime,
        min_elevation_angle: float,
        try_repeat: bool | None = None,
    ) -> tuple:
        """
        Gets the observation events of this orbit using Skyfield.

        Args:
            point (Point): Target location to observe.
            start (datetime): Start time of the observation period.
            end (datetime): End time of the observation period.
            min_elevation_angle (float): Minimum elevation angle (deg) to constrain observation.
            try_repeat (bool | None): True, if a repeat orbit should be used to improve long-term accuracy.

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        # load defaults
        if try_repeat is None:
            try_repeat = config.rc.repeat_cycle_for_observation_events
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        if len(self.elements) > 1:
            # try to use use multiple TLEs
            part_ts, element_is = self.partition_by_element_index(start, end)
            events = [
                element_is[i].to_skyfield().find_events(
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
                times, events = self.elements[0].to_skyfield().find_events(
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
        return self.elements[0].to_skyfield().find_events(topos, t_0, t_1, min_elevation_angle)

    def to_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Converts this orbit to a general perturbations orbit representation.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        return self
