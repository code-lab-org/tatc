"""
Object schema for general perturbations orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta
from typing import Literal, overload

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, model_validator
from skyfield.api import Time, wgs84
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ... import config, constants, utils
from ..surface import Point
from .gp_elements import GeneralPerturbationsElements


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
        max_search_duration: timedelta | None = None,
        lazy_load: bool | None = None,
        consistency_threshold: timedelta | None = None,
    ) -> timedelta | None:
        """
        Compute the orbit's repeat cycle, if every element agrees on one.
        Lazy-loads a previously-computed repeat cycle if available.

        Each element's own repeat cycle is computed independently (see
        `GeneralPerturbationsElements.get_repeat_cycle` for how), since a
        `GeneralPerturbationsOrbit` with multiple elements may span a
        significant maneuver (altitude change, plane change, etc.)
        partway through its history -- in which case there may be no
        single repeat cycle that legitimately describes the whole orbit.
        This method reports a repeat cycle for the orbit only if every
        element has one and they all agree within `consistency_threshold`
        of each other; otherwise it returns None. For a single-element
        orbit (the common case) this is equivalent to just asking that
        one element, since there is nothing to compare against.

        Args:
            max_delta_position (float | None): the maximum difference in position (m) allowed for a repeat.
            max_delta_velocity (float | None): the maximum difference in velocity (m/s) allowed for a repeat.
            max_search_duration (timedelta | None): the maximum period of time to search for repeats.
            lazy_load (bool | None): True, if the previously-computed repeat cycle should be loaded.
            consistency_threshold (timedelta | None): the maximum allowed spread between elements' repeat cycles.

        Returns:
            timedelta: the repeat cycle duration (if every element agrees on one)
        """
        if lazy_load is None:
            lazy_load = config.get_rc().repeat_cycle_lazy_load
        if consistency_threshold is None:
            consistency_threshold = timedelta(
                seconds=config.get_rc().repeat_cycle_consistency_threshold_s
            )

        if lazy_load:
            repeat_cycle = self.__dict__.get("repeat_cycle")
        else:
            repeat_cycle = None
        if repeat_cycle is None:
            repeat_cycle = timedelta(0)
            min_cycle = max_cycle = None
            for element in self.elements:
                cycle = element.get_repeat_cycle(
                    max_delta_position,
                    max_delta_velocity,
                    max_search_duration,
                    lazy_load,
                )
                if cycle is None:
                    repeat_cycle = timedelta(0)
                    break
                min_cycle = cycle if min_cycle is None else min(min_cycle, cycle)
                max_cycle = cycle if max_cycle is None else max(max_cycle, cycle)
                if max_cycle - min_cycle > consistency_threshold:
                    repeat_cycle = timedelta(0)
                    break
                # keep the most recent element's cycle as the orbit's
                # representative value, once every element seen so far agrees
                repeat_cycle = cycle
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

        If this orbit has multiple elements (e.g. built from a historical
        archive of TLEs via `from_tle` with multiple pairs), each query time
        is independently propagated using whichever element's epoch is
        closest to it (see `get_closest_element_index`), not always the
        first or most recent element. A vectorized `t` may therefore draw
        from different elements for different entries.

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

        Accepts plain Python `datetime`(s) for convenience; builds a Skyfield
        `Time` and delegates to `get_orbit_track_at_time`, which documents
        the multi-element selection behavior that also applies here. Prefer
        calling `get_orbit_track_at_time` directly when a Skyfield `Time` is
        already in hand, to avoid rebuilding one (see that method's
        docstring for why that matters).

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
            try_repeat = config.get_rc().repeat_cycle_for_orbit_track
        if try_repeat and len(self.elements) == 1:
            repeat_cycle = self.get_repeat_cycle()
            if repeat_cycle is not None:
                epoch = self.get_epoch()
                offset = t.utc_datetime() - epoch
                repeat_offset = np.multiply(
                    np.sign(offset / timedelta(1)),
                    np.mod(np.abs(offset / timedelta(1)), repeat_cycle / timedelta(1)),
                ) * timedelta(1)
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
        Gets the observation events (rise/culminate/set) of this orbit
        with respect to a ground point, between `start` and `end`, using
        Skyfield's `find_events`.

        Tries three strategies, in order, and uses the first that applies:

        1. If `try_repeat` and this orbit has a repeat cycle shorter than
           `end - start` (see `get_repeat_cycle`, which validates that
           *every* element agrees on the same cycle if there are several),
           events are computed once over a single repeat cycle starting at
           `start` and then copy-pasted forward for as many cycles as
           needed to cover the full period, rather than propagating the
           whole span directly. Because a validated repeat cycle means the
           whole orbit -- not just one element -- repeats identically,
           this does not need to consider which element is closest to
           each time the way strategy 2 does; the single element closest
           to `start` is enough to compute the one cycle's worth of events
           that every subsequent cycle repeats.
        2. Otherwise, if this orbit has multiple elements, the requested
           period is partitioned by whichever element's epoch is closest
           at each point in time (`partition_by_element_index`), and
           events are computed separately over each segment using its
           assigned element.
        3. Otherwise (a single element with no usable repeat cycle),
           events are computed directly over the whole period.

        Args:
            point (Point): Target location to observe.
            start (datetime): Start time of the observation period.
            end (datetime): End time of the observation period.
            min_elevation_angle (float): Minimum elevation angle (deg) to constrain observation.
            try_repeat (bool | None): True, if a repeat orbit should be used to improve long-term accuracy.

        Returns:
            tuple[skyfield.timelib.Time, numpy.ndarray]: event times and their rise (0) / culminate (1) / set (2) codes
        """
        # load defaults
        if try_repeat is None:
            try_repeat = config.get_rc().repeat_cycle_for_observation_events
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        t_0 = constants.timescale.from_datetime(start)
        if try_repeat:
            # try to compute repeat cycle events
            repeat_cycle = self.get_repeat_cycle()
            if repeat_cycle is not None and repeat_cycle < end - start:
                repeat_t_1 = constants.timescale.from_datetime(start + repeat_cycle)
                times, events = (
                    self.get_closest_element(start)
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
        if len(self.elements) > 1:
            # partition the period by whichever element is closest at each time
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
        # compute observation events directly, over the whole period
        t_1 = constants.timescale.from_datetime(end)
        return (
            self.elements[0]
            .to_skyfield()
            .find_events(topos, t_0, t_1, min_elevation_angle)
        )

    def to_gp_orbit(self, lazy_load: bool | None = None) -> GeneralPerturbationsOrbit:
        """
        Converts this orbit to a general perturbations orbit representation.
        Since this orbit already is one, this is always just `self` -- no
        computation or caching is needed, unlike `OrbitBase.to_gp_orbit`,
        which fits a `GeneralPerturbationsOrbit` from other orbit
        representations (e.g. via SGP4 fitting) and so benefits from
        lazy-loading a previously-computed result.

        Args:
            lazy_load (bool | None): accepted, but has no effect, for
                interface parity with `OrbitBase.to_gp_orbit`: callers
                that only know an orbit as `AllOrbits` (e.g.
                `Satellite.orbit`) can call `to_gp_orbit(lazy_load=...)`
                uniformly without checking which concrete type it is.

        Returns:
            GeneralPerturbationsOrbit: this orbit, unchanged
        """
        return self
