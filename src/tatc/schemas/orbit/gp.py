"""
Object schemas for general perturbations orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta, timezone
from typing import Literal, overload

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, model_validator
from sgp4 import exporter, omm
from sgp4.api import WGS72, Satrec
from sgp4.conveniences import sat_epoch_datetime
from skyfield.api import EarthSatellite, Time, wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ... import config, constants, utils
from ..surface import Point


class GeneralPerturbationsElements(BaseModel):
    """General perturbations orbital elements for a satellite."""

    object_name: str | None = Field(default=None, description="Object name.")
    epoch: datetime = Field(..., description="Epoch.")
    mean_motion: float = Field(..., description="Mean motion (degrees/second).", gt=0)
    eccentricity: float = Field(..., description="Eccentricity.", ge=0, le=1)
    inclination: float = Field(..., description="Inclination (degrees).", ge=0, le=180)
    ra_of_asc_node: float = Field(
        ..., description="Right ascension of ascending node (degrees).", ge=0, lt=360
    )
    arg_of_pericenter: float = Field(
        ..., description="Argument of pericenter (degrees).", ge=0, lt=360
    )
    mean_anomaly: float = Field(
        ..., description="Mean anomaly (degrees).", ge=0, lt=360
    )
    norad_cat_id: int = Field(default=0, description="NORAD catalog identifier.", ge=0)
    bstar: float = Field(default=0, description="Starred ballistic coefficient.")
    mean_motion_dot: float = Field(
        default=0, description="First derivative of mean motion (degrees/second^2)."
    )
    mean_motion_ddot: float = Field(
        default=0, description="Second derivative of mean motion (degrees/second^3)."
    )
    classification: str = Field(default="U", description="Classification type.")
    international_designator: str = Field(
        default="00000A", description="International designator."
    )
    ephemeris_type: int = Field(default=0, description="Ephemeris type.")
    element_set_num: int = Field(default=0, description="Element set number.")
    revolution_num: int = Field(default=0, description="Revolution number at epoch.")

    @classmethod
    def from_satrec(cls, satrec: Satrec) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from a Satrec object.

        Args:
            satrec (Satrec): The Satrec object.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        return GeneralPerturbationsElements(
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
            mean_motion_ddot=np.degrees(satrec.nddot) / 60**3,
            classification=satrec.classification,
            international_designator=satrec.intldesg,
            ephemeris_type=satrec.ephtype,
            element_set_num=satrec.elnum,
            revolution_num=satrec.revnum,
        )

    def to_satrec(self) -> Satrec:
        """
        Converts this GP elements object to a Satrec object.

        Returns:
            Satrec: the Satrec object
        """
        satrec = Satrec()
        satrec.classification = self.classification
        satrec.intldesg = self.international_designator
        satrec.ephtype = self.ephemeris_type
        satrec.elnum = self.element_set_num
        satrec.revnum = self.revolution_num
        satrec.sgp4init(
            WGS72,
            "i",
            self.norad_cat_id,
            (self.epoch - datetime(1949, 12, 31, tzinfo=timezone.utc))
            / timedelta(days=1),
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
        return timedelta(
            seconds=utils.orbital.mean_motion_to_orbit_period(self.mean_motion)
        )

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
    def from_tle(cls, tle_lines: tuple[str, str]) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from two line element (TLE) lines.

        Args:
            tle_lines (tuple[str, str]): The two TLE lines.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        return GeneralPerturbationsElements.from_satrec(
            Satrec.twoline2rv(tle_lines[0], tle_lines[1])
        )

    def to_tle(self) -> tuple[str, str]:
        """
        Converts this GP elements object to a two line element (TLE) representation.

        Returns:
            tuple[str, str]: the two line elements
        """
        return exporter.export_tle(self.to_satrec())

    @classmethod
    def from_omm_dict(cls, omm_dict: dict) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from an OMM dictionary.

        Args:
            omm_dict (dict): The OMM dictionary.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        satrec = Satrec()
        omm.initialize(satrec, omm_dict)
        elements = GeneralPerturbationsElements.from_satrec(satrec)
        # object_name has no equivalent on Satrec, so from_satrec can never
        # recover it; restore it directly from the OMM dictionary
        return elements.model_copy(update={"object_name": omm_dict.get("OBJECT_NAME")})

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
        Creates a GP elements object from OMM CSV lines. Only the first
        data row is used; all subsequent rows are ignored.

        Args:
            omm_csv (list[str]): The OMM CSV lines, including a header row.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        for fields in csv.DictReader(omm_csv):
            return GeneralPerturbationsElements.from_omm_dict(fields)
        raise ValueError("No OMM CSV lines found.")

    @classmethod
    def from_omm_json(cls, omm_json: str) -> GeneralPerturbationsElements:
        """
        Creates a GP elements object from an OMM JSON string. Only the
        first entry in the JSON array is used; all subsequent entries
        are ignored.

        Args:
            omm_json (str): The OMM JSON string, encoding a list of OMM
                records.

        Returns:
            GeneralPerturbationsElements: the GP elements
        """
        for fields in json.loads(omm_json):
            return GeneralPerturbationsElements.from_omm_dict(fields)
        raise ValueError("No OMM JSON lines found.")

    def to_skyfield(self) -> EarthSatellite:
        """
        Converts this GP elements object to a Skyfield `EarthSatellite`,
        which can be used to propagate this orbital state via SGP4.

        Returns:
            skyfield.api.EarthSatellite: the Skyfield EarthSatellite
        """
        return EarthSatellite.from_omm(constants.timescale, self.to_omm_dict())


class GeneralPerturbationsOrbit(BaseModel):
    """
    Orbit defined with general perturbations (GP) elements.
    """

    type: Literal["gp"] = Field(default="gp", description="Orbit type discriminator.")
    elements: list[GeneralPerturbationsElements] = Field(
        ..., description="General perturbations elements.", min_length=1
    )

    @model_validator(mode="after")
    def _sort_elements_by_epoch(self) -> GeneralPerturbationsOrbit:
        """
        Sorts elements by epoch (ascending) after construction.
        get_closest_element_index relies on np.searchsorted, which
        silently returns incorrect results if its input is not sorted, so
        this guarantees that precondition holds regardless of the order
        elements were provided in. This only covers elements provided at
        construction time; directly mutating self.elements in place
        afterward (e.g. via .append()) bypasses this validator and can
        reintroduce unsorted order.
        """
        self.elements.sort(key=lambda el: el.epoch)
        return self

    def get_semimajor_axis(self, index: int = 0) -> float:
        """
        Gets the semimajor axis of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the semimajor axis (meters)
        """
        return self.elements[index].get_semimajor_axis()

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

    def get_eccentricity(self, index: int = 0) -> float:
        """
        Gets the eccentricity of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the eccentricity
        """
        return self.elements[index].eccentricity

    def get_epoch(self, index: int = 0) -> datetime:
        """
        Gets the epoch of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            datetime: the epoch
        """
        return self.elements[index].epoch

    def get_mean_motion(self, index: int = 0) -> float:
        """
        Gets the mean motion of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the mean motion (degrees/second)
        """
        return self.elements[index].mean_motion

    def get_mean_anomaly(self, index: int = 0) -> float:
        """
        Gets the mean anomaly of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the mean anomaly (degrees)
        """
        return self.elements[index].mean_anomaly

    def get_orbit_period(self, index: int = 0) -> timedelta:
        """
        Gets the approximate orbit period of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            timedelta: the orbit period
        """
        return self.elements[index].get_orbit_period()

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

    def get_catalog_number(self, index: int = 0) -> int:
        """
        Gets the NORAD catalog number of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            int: the NORAD catalog number
        """
        return self.elements[index].norad_cat_id

    def get_bstar(self, index: int = 0) -> float:
        """
        Gets the starred ballistic coefficient of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the starred ballistic coefficient
        """
        return self.elements[index].bstar

    def get_mean_motion_dot(self, index: int = 0) -> float:
        """
        Gets the first derivative of mean motion of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the first derivative of mean motion (degrees/second^2)
        """
        return self.elements[index].mean_motion_dot

    def get_mean_motion_ddot(self, index: int = 0) -> float:
        """
        Gets the second derivative of mean motion of the specified element.

        Args:
            index (int): the index of the element

        Returns:
            float: the second derivative of mean motion (degrees/second^3)
        """
        return self.elements[index].mean_motion_ddot

    @classmethod
    def from_tle(cls, tle_lines: list[str]) -> GeneralPerturbationsOrbit:
        """
        Creates a GP orbit from two line element (TLE) lines. Multiple
        TLEs (e.g. a history of element sets for one satellite) may be
        concatenated into a single flat list of lines, two per element
        set, to construct an orbit with multiple elements.

        Args:
            tle_lines (list[str]): the two line element lines

        Returns:
            GeneralPerturbationsOrbit: the GP orbit
        """
        if len(tle_lines) % 2 != 0:
            raise ValueError(
                f"Expected an even number of TLE lines (two per element set), "
                f"got {len(tle_lines)}."
            )
        return GeneralPerturbationsOrbit(
            elements=[
                GeneralPerturbationsElements.from_tle((tle_lines[i], tle_lines[i + 1]))
                for i in range(0, len(tle_lines), 2)
            ]
        )

    @classmethod
    def from_omm_csv(cls, omm_csv: list[str]) -> GeneralPerturbationsOrbit:
        """
        Creates a GP orbit from OMM CSV lines, using every data row
        (unlike GeneralPerturbationsElements.from_omm_csv, which only
        uses the first) to build one element per row.

        Args:
            omm_csv (list[str]): The OMM CSV lines, including a header row.

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
        Creates a GP orbit from an OMM JSON string, using every entry
        (unlike GeneralPerturbationsElements.from_omm_json, which only
        uses the first) to build one element per entry.

        Args:
            omm_json (str): The OMM JSON string, encoding a list of OMM
                records.

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
        Lazy-loads the epoch times for all elements in this orbit. The
        cache is invalidated (and recomputed) if the number of elements
        has changed since it was last computed (e.g. after appending or
        removing an element), but not if an element is replaced in place
        at the same list position with a different epoch -- such a
        same-length in-place replacement is unusual/unsupported usage
        that this lightweight invalidation check cannot detect without
        recomputing on every call, which would defeat the purpose of
        caching.

        Returns:
            list[datetime]: the epoch times
        """
        # lazy-load epochs, invalidating the cache if the element count changed
        element_epochs = self.__dict__.get("element_epochs")
        if element_epochs is None or len(element_epochs) != len(self.elements):
            # extract the element epoch times
            element_epochs = [el.epoch for el in self.elements]
            self.__dict__["element_epochs"] = element_epochs  # type: ignore
        return element_epochs  # type: ignore

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
            derived_el.mean_anomaly = np.mod(
                derived_el.mean_anomaly + delta_mean_anomaly, 360
            )
            derived_el.ra_of_asc_node = np.mod(
                derived_el.ra_of_asc_node + delta_raan, 360
            )
            derived_elements.append(derived_el)
        return GeneralPerturbationsOrbit(elements=derived_elements)

    @overload
    def get_closest_element_index(self, at_times: None) -> int: ...

    @overload
    def get_closest_element_index(self, at_times: datetime) -> int: ...

    @overload
    def get_closest_element_index(self, at_times: list[datetime]) -> list[int]: ...

    @overload
    def get_closest_element_index(
        self, at_times: npt.NDArray[np.datetime64]
    ) -> list[int]: ...

    def get_closest_element_index(
        self, at_times: datetime | list[datetime] | npt.NDArray[np.datetime64] | None
    ) -> int | list[int]:
        """
        Gets the closest element index to specified time(s), assuming
        elements are sorted by epoch (guaranteed for any orbit built
        through the constructor, since elements are sorted at
        construction time; see _sort_elements_by_epoch).

        Args:
            at_times (datetime | list[datetime] | npt.NDArray[np.datetime64] | None):
                specified times, or None to always select the first element (index 0)

        Returns:
            int | list[int]: closest element index or indices
        """

        if at_times is None:
            return 0
        # lazy-load element epochs
        element_epochs = utils.to_datetime64_ns(self.get_element_epochs())
        # handle scalar
        if isinstance(at_times, datetime):
            at_time = utils.to_datetime64_ns(at_times)
            idx = np.searchsorted(element_epochs, at_time, side="left")
            return (
                int(idx - 1)
                if idx > 0
                and (
                    idx == len(element_epochs)
                    or abs(at_time - element_epochs[idx - 1])
                    < abs(at_time - element_epochs[idx])
                )
                else int(idx)
            )
        # handle vector
        at_time_array = utils.to_datetime64_ns(at_times)
        indices = np.searchsorted(element_epochs, at_time_array, side="left")
        return [
            (
                int(idx - 1)
                if idx > 0
                and (
                    idx == len(element_epochs)
                    or abs(at_time_array[i] - element_epochs[idx - 1])
                    < abs(at_time_array[i] - element_epochs[idx])
                )
                else int(idx)
            )
            for i, idx in enumerate(indices)
        ]

    @overload
    def get_closest_element(self, at_times: None) -> GeneralPerturbationsElements: ...

    @overload
    def get_closest_element(
        self, at_times: datetime
    ) -> GeneralPerturbationsElements: ...

    @overload
    def get_closest_element(
        self, at_times: list[datetime]
    ) -> list[GeneralPerturbationsElements]: ...

    @overload
    def get_closest_element(
        self, at_times: npt.NDArray[np.datetime64]
    ) -> list[GeneralPerturbationsElements]: ...

    def get_closest_element(
        self, at_times: datetime | list[datetime] | npt.NDArray[np.datetime64] | None
    ) -> GeneralPerturbationsElements | list[GeneralPerturbationsElements]:
        """
        Gets the closest element to specified time(s).

        Args:
            at_times (datetime | list[datetime] | npt.NDArray[np.datetime64] | None):
                specified times, or None to always select the first element (index 0)

        Returns:
            GeneralPerturbationsElements | list[GeneralPerturbationsElements]: closest element or elements
        """
        indices = self.get_closest_element_index(at_times)
        if isinstance(indices, int):
            return self.elements[indices]
        return [self.elements[i] for i in indices]

    def partition_by_element_index(
        self, start: datetime, end: datetime
    ) -> tuple[list[datetime], list[int]]:
        """
        Partitions the time range [start, end] into consecutive segments,
        each assigned the index of whichever element is closest
        throughout that segment. Uses get_element_epochs() (benefiting
        from its lazy-load cache) rather than re-reading each element's
        epoch directly. The midpoint between each pair of consecutive
        elements' epochs is where the closest element switches from one
        to the next (elements are guaranteed sorted by epoch by the
        constructor's _sort_elements_by_epoch validator), so only
        midpoints strictly inside (start, end) become segment boundaries.
        Each segment's element index is determined by querying
        get_closest_element_index at that segment's own midpoint, so it
        is correct by construction rather than tracked separately.

        Args:
            start (datetime): Start time.
            end (datetime): End time.

        Returns:
            tuple[list[datetime], list[int]]: segment boundary times
                (length N+1, including start and end) and the element
                index for each of the N segments between consecutive
                boundaries (length N).
        """
        epochs = self.get_element_epochs()
        epoch_midpoints = [
            epochs[i] + (epochs[i + 1] - epochs[i]) / 2 for i in range(len(epochs) - 1)
        ]
        boundary_times = (
            [start] + [t for t in epoch_midpoints if start < t < end] + [end]
        )
        segment_midpoints = [
            boundary_times[i] + (boundary_times[i + 1] - boundary_times[i]) / 2
            for i in range(len(boundary_times) - 1)
        ]
        return boundary_times, self.get_closest_element_index(segment_midpoints)

    def get_repeat_cycle(
        self,
        max_delta_position: float | None = None,
        max_delta_velocity: float | None = None,
        min_elevation_angle: float | None = None,
        max_search_duration: timedelta | None = None,
        lazy_load: bool | None = None,
    ) -> timedelta | None:
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
            position_0, velocity_0 = satellite.at(
                constants.timescale.from_datetime(epoch)
            ).frame_xyz_and_velocity(itrs)
            # find candidate repeat events
            ts, es = satellite.find_events(
                datum,
                constants.timescale.from_datetime(epoch + timedelta(minutes=10)),
                constants.timescale.from_datetime(epoch + max_search_duration),
                min_elevation_angle,
            )
            # compute position and velocity at culmination in Earth-centered Earth-fixed frame
            position, velocity = satellite.at(ts[es == 1]).frame_xyz_and_velocity(itrs)
            p_m = np.array(position.m)
            v_m_per_s = np.array(velocity.m_per_s)
            p_0_m = np.array(position_0.m)
            v_0_m_per_s = np.array(velocity_0.m_per_s)
            # apply validity conditions on position and velocity error norms
            is_valid = np.logical_and(
                np.linalg.norm((p_m.T - p_0_m.T).T, axis=0) < max_delta_position,
                np.linalg.norm((v_m_per_s.T - v_0_m_per_s.T).T, axis=0)
                < max_delta_velocity,
            )
            if np.any(is_valid):
                # assign repeat cycle
                repeat_cycle = ts[es == 1][is_valid][0].utc_datetime() - epoch
            else:
                # assign zero repeat cycle value to avoid recalculation
                repeat_cycle = timedelta(0)
            self.__dict__["repeat_cycle"] = repeat_cycle  # type: ignore
        if repeat_cycle is not None and repeat_cycle > timedelta(0):
            return repeat_cycle
        return None

    def get_orbit_track_at_time(self, t: Time) -> Geocentric:
        """
        Gets the true (directly propagated) orbit track of this orbit at given
        Skyfield time(s), in the inertial (GCRS) frame.

        Prefer this method over `get_orbit_track` when a Skyfield `Time` is
        already in hand (e.g. while iterating a Skyfield search such as
        `skyfield.searchlib.find_discrete`). Converting a `Time` to Python
        `datetime` objects and back (as `get_orbit_track` must, since it only
        accepts `datetime`) builds a new `Time` instance that starts without any
        of the per-instant quantities Skyfield caches on a `Time` object (such as
        nutation angles), forcing Skyfield to recompute them from scratch.

        Args:
            t (skyfield.timelib.Time): time(s) at which to compute position/velocity.

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        if len(self.elements) > 1:
            # try to use multiple TLEs
            nearest_indices = self.get_closest_element_index(t.utc_datetime())
            if isinstance(nearest_indices, int):
                return self.elements[nearest_indices].to_skyfield().at(t)  # type: ignore
            nearest_indices = np.asarray(nearest_indices)
            position_au = np.empty((3,) + t.shape)
            velocity_au_per_d = np.empty((3,) + t.shape)
            for element_index in np.unique(nearest_indices):
                # propagate each distinct nearest TLE across all its assigned
                # times in one vectorized call, rather than one time at a time
                mask = nearest_indices == element_index
                track = self.elements[element_index].to_skyfield().at(t[mask])
                position_au[:, mask] = track.position.au
                velocity_au_per_d[:, mask] = track.velocity.au_per_d
            return Geocentric(position_au, velocity_au_per_d, t)
        # compute satellite positions directly at the given time(s)
        return self.elements[0].to_skyfield().at(t)  # type: ignore

    def get_orbit_track(self, times: datetime | list[datetime]) -> Geocentric:
        """
        Gets the true (directly propagated) orbit track of this orbit using
        Skyfield, in the inertial (GCRS) frame.

        Args:
            times (datetime | list[datetime]): time(s) at which to compute position/velocity.

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        t = (
            constants.timescale.from_datetime(times)
            if isinstance(times, datetime)
            else constants.timescale.from_datetimes(times)
        )
        return self.get_orbit_track_at_time(t)

    def get_geographic_position_at_time(
        self, t: Time, try_repeat: bool | None = None
    ) -> GeographicPosition:
        """
        Gets the geodetic (WGS84) position of this orbit at given Skyfield
        time(s), in an Earth-fixed frame.

        Unlike `get_orbit_track_at_time`, this method may substitute a detected
        repeat cycle to improve long-term accuracy: rather than directly
        propagating to a possibly-distant `t`, it propagates near this orbit's
        epoch (reducing `t`'s offset from epoch modulo the repeat cycle) and
        relies on the orbit's ground track repeating with that period. Because
        the result is a `GeographicPosition` -- a location descriptor, not a
        frozen inertial state vector -- it can be freely reused afterward (e.g.
        `.at(some_time)` for a look angle or Sun angle at any moment) without
        carrying forward any inaccuracy from the substitution.

        Args:
            t (skyfield.timelib.Time): time(s) at which to compute geodetic position.
            try_repeat (bool | None): True, if a repeat orbit should be used to improve long-term accuracy.

        Returns:
            skyfield.toposlib.GeographicPosition: the geodetic position
        """
        if try_repeat is None:
            try_repeat = config.rc.repeat_cycle_for_orbit_track
        if try_repeat and len(self.elements) == 1:
            repeat_cycle = self.get_repeat_cycle()
            if repeat_cycle is not None:
                epoch = self.get_epoch()
                offset = t.utc_datetime() - epoch
                repeat_offset = np.multiply(
                    np.sign(offset / timedelta(1)),
                    np.mod(np.abs(offset / timedelta(1)), repeat_cycle / timedelta(1)),
                )
                repeat_times = (
                    constants.timescale.from_datetime(epoch + repeat_offset)
                    if t.shape == ()
                    else constants.timescale.from_datetimes(epoch + repeat_offset)
                )
                return wgs84.geographic_position_of(
                    self.elements[0].to_skyfield().at(repeat_times)
                )
        # compute geodetic position from a true, directly propagated orbit track
        return wgs84.geographic_position_of(self.get_orbit_track_at_time(t))

    def get_geographic_position(
        self, times: datetime | list[datetime], try_repeat: bool | None = None
    ) -> GeographicPosition:
        """
        Gets the geodetic (WGS84) position of this orbit at given time(s), in
        an Earth-fixed frame.

        Args:
            times (datetime | list[datetime]): time(s) at which to compute geodetic position.
            try_repeat (bool | None): True, if a repeat orbit should be used to improve long-term accuracy.

        Returns:
            skyfield.toposlib.GeographicPosition: the geodetic position
        """
        t = (
            constants.timescale.from_datetime(times)
            if isinstance(times, datetime)
            else constants.timescale.from_datetimes(times)
        )
        return self.get_geographic_position_at_time(t, try_repeat)

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
                self.elements[element_is[i]]
                .to_skyfield()
                .find_events(
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
                times, events = (
                    self.elements[0]
                    .to_skyfield()
                    .find_events(topos, t_0, repeat_t_1, min_elevation_angle)
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
        return (
            self.elements[0]
            .to_skyfield()
            .find_events(topos, t_0, t_1, min_elevation_angle)
        )

    def to_gp_orbit(self) -> GeneralPerturbationsOrbit:
        """
        Converts this orbit to a general perturbations orbit representation.

        Returns:
            GeneralPerturbationsOrbit: the general perturbations orbit
        """
        return self
