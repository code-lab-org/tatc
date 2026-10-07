"""
Object schema for general perturbations orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import csv
import json
import warnings
from collections.abc import Callable
from datetime import datetime, timedelta
from typing import Literal, overload

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, model_validator
from skyfield.api import Time, wgs84
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from skyfield.sgp4lib import EarthSatellite
from skyfield.toposlib import GeographicPosition

from ... import config, constants, utils
from ...utils.cache import get_cached
from ..surface import Point
from .gp_elements import GeneralPerturbationsElements

_FIRST_REPEAT = -1
"""Source code of times propagated with the first element's repeat track."""
_LAST_REPEAT = -2
"""Source code of times propagated with the last element's repeat track."""


def _to_time(times: datetime | list[datetime]) -> Time:
    """
    Converts time(s) to a Skyfield `Time`.

    Args:
        times (datetime | list[datetime]): The time(s).

    Returns:
        skyfield.timelib.Time: the Skyfield time(s)
    """
    if isinstance(times, datetime):
        return constants.timescale.from_datetime(times)
    return constants.timescale.from_datetimes(times)


def _bisect(
    excess: Callable[[npt.NDArray], npt.NDArray],
    lower: npt.NDArray,
    upper: npt.NDArray,
) -> npt.NDArray:
    """
    Bisects brackets of a change of sign of a function to a millisecond.

    Args:
        excess (Callable[[numpy.typing.NDArray], numpy.typing.NDArray]): The
            function of time (TT Julian date).
        lower (numpy.typing.NDArray): The lower ends of the brackets.
        upper (numpy.typing.NDArray): The upper ends of the brackets.

    Returns:
        numpy.typing.NDArray: the times of the changes of sign
    """
    if len(lower) == 0:
        return lower
    f_lower = excess(lower)
    for _ in range(64):
        if np.max(upper - lower) < 1e-3 / 86400:
            break
        middle = (lower + upper) / 2
        f_middle = excess(middle)
        same = np.sign(f_middle) == np.sign(f_lower)
        lower = np.where(same, middle, lower)
        f_lower = np.where(same, f_middle, f_lower)
        upper = np.where(same, upper, middle)
    return (lower + upper) / 2


def _find_crossing(
    excess: Callable[[npt.NDArray], npt.NDArray], start: float, stop: float
) -> float | None:
    """
    Finds the first time, from `start` (where a function is not negative)
    toward `stop`, when the function becomes negative, if it does so once
    between them.

    Args:
        excess (Callable[[numpy.typing.NDArray], numpy.typing.NDArray]): The
            function of time (TT Julian date).
        start (float): The time from which to search.
        stop (float): The time toward which to search (excluded).

    Returns:
        float | None: the time of the crossing, if found
    """
    samples = np.linspace(start, stop, 26)[:-1]
    negative = np.flatnonzero(excess(samples) < 0)
    if len(negative) == 0:
        return None
    bracket = np.sort(samples[negative[0] - 1 : negative[0] + 1])
    return float(_bisect(excess, bracket[:1], bracket[1:])[0])


def _find_events(
    satellite: EarthSatellite,
    topos: GeographicPosition,
    t_0: Time,
    t_1: Time,
    min_elevation_angle: float,
) -> tuple[Time, npt.NDArray]:
    """
    Find the rise, culminate, and set events of a satellite with respect to a
    ground position using Skyfield's `find_events`, with the rise and set
    times refined by bisection and the rise or set events of passes that
    culminate outside the period added.

    Skyfield's search stops refining all rise and set brackets once the first
    one converges, which assumes they start with equal widths; over long
    periods they do not, so some rise and set times can be several seconds
    early or late. Each rise or set event lies between the preceding event
    (or `t_0`) and the reported time, where the elevation angle crosses the
    minimum, so that bracket is bisected to a millisecond.

    Skyfield finds rise and set events around culminations, so it misses the
    set of a pass that culminates before `t_0` (and the rise of one that
    culminates after `t_1`): if the satellite is above the minimum elevation
    angle at `t_0` but the first event found is a rise (or none is found and
    it is below at `t_1`), its set in between is added (and, conversely, a
    rise before `t_1`).

    Args:
        satellite (skyfield.sgp4lib.EarthSatellite): The satellite.
        topos (skyfield.toposlib.GeographicPosition): The ground position.
        t_0 (skyfield.timelib.Time): The start time.
        t_1 (skyfield.timelib.Time): The end time.
        min_elevation_angle (float): The minimum elevation angle (degrees).

    Returns:
        tuple[skyfield.timelib.Time, numpy.ndarray]: event times and their
            rise (0) / culminate (1) / set (2) codes
    """
    times, events = satellite.find_events(topos, t_0, t_1, min_elevation_angle)
    jd = np.array(times.tt, dtype=float, ndmin=1)
    events = np.array(events, dtype=int, ndmin=1)
    relative = satellite - topos

    def excess(x: npt.NDArray) -> npt.NDArray:
        t = constants.timescale.tt_jd(x)
        return relative.at(t).altaz()[0].degrees - min_elevation_angle

    refine = np.flatnonzero(events != 1)
    if len(refine) > 0:
        lower = np.concatenate(([t_0.tt], jd))[refine]
        upper = jd[refine]
        bracketed = np.sign(excess(lower)) != np.sign(excess(upper))
        jd[refine[bracketed]] = _bisect(excess, lower[bracketed], upper[bracketed])
    crossing_jd, crossing_events = jd[events != 1], events[events != 1]
    f_0, f_1 = excess(np.array([t_0.tt, t_1.tt]))
    added = []
    if f_0 >= 0 and (len(crossing_events) == 0 or crossing_events[0] == 0):
        if len(crossing_events) > 0 or f_1 < 0:
            stop = crossing_jd[0] if len(crossing_events) > 0 else t_1.tt
            added.append((_find_crossing(excess, t_0.tt, stop), 2))
    if f_1 >= 0 and (len(crossing_events) == 0 or crossing_events[-1] == 2):
        if len(crossing_events) > 0 or f_0 < 0:
            stop = crossing_jd[-1] if len(crossing_events) > 0 else t_0.tt
            added.append((_find_crossing(excess, t_1.tt, stop), 0))
    added = [(time, code) for time, code in added if time is not None]
    if len(added) == 0:
        return constants.timescale.tt_jd(jd), events
    jd = np.concatenate((jd, [time for time, _ in added]))
    events = np.concatenate((events, [code for _, code in added]))
    order = np.argsort(jd, kind="stable")
    return constants.timescale.tt_jd(jd[order]), events[order]


class _RepeatTrack:
    """
    The orbit track of an element maintained on its repeat ground track (see
    `GeneralPerturbationsElements.get_repeat_element`), repeated with its
    repeat cycle: each time is shifted by a whole number of repeat cycles to
    the repeat cycle just after the element's epoch (or, for a time before
    the epoch, just before it), so that the element is propagated no more
    than one repeat cycle from its epoch, however far the time is from it.
    The maintained element returns to its initial Earth-fixed position after
    each repeat cycle, so that the repeated cycles join.
    """

    def __init__(self, element: GeneralPerturbationsElements, repeat_cycle: timedelta):
        """
        Initializes a repeat track.

        Args:
            element (GeneralPerturbationsElements): The element maintained on its repeat ground track.
            repeat_cycle (timedelta): The (exact) repeat cycle of the maintained element.
        """
        self.satellite = element.to_skyfield()
        self.epoch = element.epoch
        self.repeat_cycle = repeat_cycle
        self.epoch_time = constants.timescale.from_datetime(element.epoch)

    def shift(self, t: Time) -> Time:
        """
        Shifts times by whole repeat cycles to the repeat cycle just after
        (or, for times before the epoch, just before) the epoch.

        Args:
            t (skyfield.timelib.Time): The time(s).

        Returns:
            skyfield.timelib.Time: the shifted time(s)
        """
        cycle = self.repeat_cycle.total_seconds() / 86400
        offset = (t.whole - self.epoch_time.whole) + (
            t.tt_fraction - self.epoch_time.tt_fraction
        )
        shift = np.trunc(offset / cycle) * cycle
        return constants.timescale.tt_jd(t.whole, t.tt_fraction - shift)

    def at(self, t: Time) -> Geocentric:
        """
        Gets the orbit track at given time(s): the satellite's Earth-fixed
        position and velocity are those propagated to the shifted time(s)
        (see `shift`), expressed in the inertial (GCRS) frame at the time(s)
        themselves, so that quantities that depend on time (such as the
        Sun's position) are unaffected.

        Args:
            t (skyfield.timelib.Time): The time(s).

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        shifted = self.shift(t)
        if np.all(shifted.tt_fraction == t.tt_fraction):
            return self.satellite.at(t)  # type: ignore
        position, velocity = self.satellite.at(shifted).frame_xyz_and_velocity(itrs)
        # invert the Earth-fixed transformation (r_itrs = R r, v_itrs = R v + V r_itrs)
        rotation = itrs.rotation_at(t)
        rate = itrs._dRdt_times_RT_at(t)  # pylint: disable=protected-access
        v_rotating = velocity.au_per_d - np.einsum(
            "ij...,j...->i...", rate, position.au
        )
        return Geocentric(
            np.einsum("ji...,j...->i...", rotation, position.au),
            np.einsum("ji...,j...->i...", rotation, v_rotating),
            t,
            center=399,
        )

    def get_boundaries(self, start: datetime, end: datetime) -> list[datetime]:
        """
        Gets the ends of repeat cycles (other than the epoch) strictly
        between `start` and `end`, at which the repeated orbit track joins
        the next cycle.

        Args:
            start (datetime): Start time.
            end (datetime): End time.

        Returns:
            list[datetime]: the ends of repeat cycles
        """
        first = int(np.ceil((start - self.epoch) / self.repeat_cycle))
        last = int(np.floor((end - self.epoch) / self.repeat_cycle))
        return [
            self.epoch + cycles * self.repeat_cycle
            for cycles in range(first, last + 1)
            if cycles != 0 and start < self.epoch + cycles * self.repeat_cycle < end
        ]

    def find_events(
        self,
        topos: GeographicPosition,
        start: datetime,
        end: datetime,
        min_elevation_angle: float,
    ) -> list[tuple[datetime, int]]:
        """
        Finds the observation events between `start` and `end`, which must
        both be on the same side of the epoch, by computing those of the
        repeat cycle just after (or before) the epoch once, over the parts of
        it that the period (shifted by whole repeat cycles) covers, and
        shifting them back to the period.

        Args:
            topos (skyfield.toposlib.GeographicPosition): Target location to observe.
            start (datetime): Start time of the observation period.
            end (datetime): End time of the observation period.
            min_elevation_angle (float): Minimum elevation angle (deg) to constrain observation.

        Returns:
            list[tuple[datetime, int]]: event times and their rise (0) / culminate (1) / set (2) codes
        """
        after = start >= self.epoch
        cycle = self.repeat_cycle
        # pieces of the period within each repeat cycle, as their shift and
        # their bounds once shifted
        pieces = []
        first = int(np.trunc((start - self.epoch) / cycle))
        last = int(np.trunc((end - self.epoch) / cycle))
        for cycles in range(first, last + 1):
            shift = cycles * cycle
            cycle_start = self.epoch + shift - (timedelta(0) if after else cycle)
            lower = max(start, cycle_start)
            upper = min(end, cycle_start + cycle)
            if lower < upper:
                pieces.append((shift, lower - shift, upper - shift))
        if len(pieces) == 0:
            return []
        times, codes = _find_events(
            self.satellite,
            topos,
            constants.timescale.from_datetime(min(p[1] for p in pieces)),
            constants.timescale.from_datetime(max(p[2] for p in pieces)),
            min_elevation_angle,
        )
        if len(codes) == 0:
            return []
        times_py = np.atleast_1d(times.utc_datetime())
        codes = np.atleast_1d(codes)
        events = []
        for shift, lower, upper in pieces:
            selected = (times_py >= lower) & (times_py <= upper)
            events.extend(zip(times_py[selected] + shift, codes[selected]))
        return events


class GeneralPerturbationsOrbit(BaseModel):
    """
    Orbit defined with general perturbations (GP) elements.

    The orbit is propagated with SGP4, using whichever element's epoch is
    closest to each time. If a repeat cycle is used (by default, see the
    `repeat_cycle_for_orbit_track` and `repeat_cycle_for_observation_events`
    settings) and either the orbit declares one (`repeat_cycle`) or one is
    found (see `GeneralPerturbationsElements.get_repeat_cycle`), times before
    the first element's epoch and after the last element's epoch are instead
    propagated with that element maintained on its repeat ground track and
    repeated with its repeat cycle (see `get_repeat_element`), as for an
    orbit maintained on its repeat ground track, which SGP4 does not model.
    For an orbit with a single element, this applies to all times.
    """

    type: Literal["gp"] = Field(default="gp", description="Orbit type discriminator.")
    elements: list[GeneralPerturbationsElements] = Field(
        ..., description="General perturbations elements.", min_length=1
    )
    remove_drag: bool = Field(
        default=False,
        description="True, to propagate the orbit without drag, as for an "
        + "orbit maintained against drag: the elements' drag terms (B* and "
        + "the derivatives of mean motion) are ignored, while the "
        + "Earth's oblateness (and, for deep-space orbits, the Moon and Sun) "
        + "still perturb the orbit.",
    )
    repeat_cycle: timedelta | None = Field(
        default=None,
        description="Repeat cycle of an orbit maintained on a repeat ground "
        + "track (for example, 91 days), which overrides the computed repeat "
        + "cycle; it is refined to the nearest whole number of nodal days "
        + "(or, for a sun-synchronous orbit, mean solar days) of the elements "
        + "maintained on the repeat ground track (see "
        + "`GeneralPerturbationsElements.get_repeat_element`). Requires `remove_drag`.",
        gt=timedelta(0),
    )

    @model_validator(mode="after")
    def _sort_elements_by_epoch(self) -> GeneralPerturbationsOrbit:
        """
        Sorts elements by epoch (ascending) after construction, which
        `get_closest_element_index` requires. This only covers elements
        provided at construction time; mutating `elements` in place
        afterward (e.g. via `append()`) can reintroduce unsorted order.
        """
        self.elements.sort(key=lambda el: el.epoch)
        return self

    @model_validator(mode="after")
    def _validate_repeat_cycle(self) -> GeneralPerturbationsOrbit:
        """
        Validates that a repeat cycle is only declared without drag (with
        drag, the orbit would not repeat).
        """
        if self.repeat_cycle is not None and not self.remove_drag:
            raise ValueError("repeat_cycle requires remove_drag=True.")
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
    def from_tle(
        cls,
        tle_lines: list[str],
        remove_drag: bool = False,
        repeat_cycle: timedelta | None = None,
    ) -> GeneralPerturbationsOrbit:
        """
        Creates a GP orbit from two line element (TLE) lines. Multiple
        TLEs (e.g. a history of element sets for one satellite) may be
        concatenated into a single flat list of lines, two per element
        set, to construct an orbit with multiple elements.

        Args:
            tle_lines (list[str]): the two line element lines
            remove_drag (bool): True, to propagate the orbit without drag.
            repeat_cycle (timedelta | None): The repeat cycle of a maintained
                orbit (requires `remove_drag`).

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
            ],
            remove_drag=remove_drag,
            repeat_cycle=repeat_cycle,
        )

    @classmethod
    def from_omm_csv(
        cls,
        omm_csv: list[str],
        remove_drag: bool = False,
        repeat_cycle: timedelta | None = None,
    ) -> GeneralPerturbationsOrbit:
        """
        Creates a GP orbit from OMM CSV lines, using every data row
        (unlike GeneralPerturbationsElements.from_omm_csv, which only
        uses the first) to build one element per row.

        Args:
            omm_csv (list[str]): The OMM CSV lines, including a header row.
            remove_drag (bool): True, to propagate the orbit without drag.
            repeat_cycle (timedelta | None): The repeat cycle of a maintained
                orbit (requires `remove_drag`).

        Returns:
            GeneralPerturbationsOrbit: the GP orbit
        """
        return GeneralPerturbationsOrbit(
            elements=[
                GeneralPerturbationsElements.from_omm_dict(fields)
                for fields in csv.DictReader(omm_csv)
            ],
            remove_drag=remove_drag,
            repeat_cycle=repeat_cycle,
        )

    @classmethod
    def from_omm_json(
        cls,
        omm_json: str,
        remove_drag: bool = False,
        repeat_cycle: timedelta | None = None,
    ) -> GeneralPerturbationsOrbit:
        """
        Creates a GP orbit from an OMM JSON string, using every entry
        (unlike GeneralPerturbationsElements.from_omm_json, which only
        uses the first) to build one element per entry.

        Args:
            omm_json (str): The OMM JSON string, encoding a list of OMM
                records.
            remove_drag (bool): True, to propagate the orbit without drag.
            repeat_cycle (timedelta | None): The repeat cycle of a maintained
                orbit (requires `remove_drag`).

        Returns:
            GeneralPerturbationsOrbit: the GP orbit
        """
        return GeneralPerturbationsOrbit(
            elements=[
                GeneralPerturbationsElements.from_omm_dict(fields)
                for fields in json.loads(omm_json)
            ],
            remove_drag=remove_drag,
            repeat_cycle=repeat_cycle,
        )

    def get_element_epochs(self) -> list[datetime]:
        """
        Gets the epochs of all elements in this orbit. The epochs are cached
        until the elements change (e.g. after appending or removing an
        element), but not if an element's epoch is changed in place.

        Returns:
            list[datetime]: the epoch times
        """
        return get_cached(
            self,
            "element_epochs",
            tuple(id(el) for el in self.elements),
            lambda: [el.epoch for el in self.elements],
        )

    def _get_element_epochs_ns(self) -> npt.NDArray[np.datetime64]:
        """
        Gets the epochs of all elements in this orbit as a cached array (see
        `get_element_epochs`).

        Returns:
            numpy.typing.NDArray[numpy.datetime64]: the epoch times
        """
        return get_cached(
            self,
            "element_epochs_ns",
            tuple(id(el) for el in self.elements),
            lambda: utils.to_datetime64_ns(self.get_element_epochs()),
        )

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
        return GeneralPerturbationsOrbit(
            elements=derived_elements,
            remove_drag=self.remove_drag,
            repeat_cycle=self.repeat_cycle,
        )

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

    @overload
    def get_closest_element_index(self, at_times: Time) -> int | list[int]: ...

    def get_closest_element_index(
        self,
        at_times: datetime | list[datetime] | npt.NDArray[np.datetime64] | Time | None,
    ) -> int | list[int]:
        """
        Gets the closest element index to specified time(s), assuming
        elements are sorted by epoch (guaranteed for any orbit built
        through the constructor, since elements are sorted at
        construction time; see _sort_elements_by_epoch). A time exactly
        halfway between two epochs selects the later element.

        Args:
            at_times (datetime | list[datetime] | npt.NDArray[np.datetime64] | Time | None):
                specified times (a Skyfield `Time` may be a scalar or an array),
                or None to always select the first element (index 0)

        Returns:
            int | list[int]: closest element index (for a scalar time) or indices
        """
        if at_times is None:
            return 0
        element_epochs = self._get_element_epochs_ns()
        at_time_array = utils.to_datetime64_ns(at_times)
        is_scalar = np.ndim(at_time_array) == 0
        at_time_array = np.atleast_1d(at_time_array)
        indices = np.searchsorted(element_epochs, at_time_array, side="left")
        # select the preceding epoch when it is strictly closer than the
        # following one (or when there is no following epoch)
        preceding = element_epochs[np.maximum(indices - 1, 0)]
        following = element_epochs[np.minimum(indices, len(element_epochs) - 1)]
        use_preceding = (indices > 0) & (
            (indices == len(element_epochs))
            | (np.abs(at_time_array - preceding) < np.abs(at_time_array - following))
        )
        closest = np.where(use_preceding, indices - 1, indices)
        return int(closest[0]) if is_scalar else closest.tolist()

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
        each assigned the index of whichever element is closest throughout
        that segment. The closest element switches at the midpoint between
        consecutive elements' epochs, so only midpoints strictly inside
        (start, end) become segment boundaries.

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
        max_delta_semimajor_axis: float | None = None,
    ) -> timedelta | None:
        """
        Gets the orbit's repeat cycle, if every element agrees on one.

        If the orbit declares a repeat cycle (`repeat_cycle`, for an orbit
        maintained on a repeat ground track without drag), it is refined to
        the whole number of nodal (or mean solar) days of the first element
        maintained on its repeat ground track (see `get_repeat_element`).

        Otherwise, each element's own repeat cycle is computed independently
        (see `GeneralPerturbationsElements.get_repeat_cycle`), since an
        orbit with multiple elements may span a significant maneuver
        (altitude change, plane change, etc.), in which case there may be no
        single repeat cycle that describes the whole orbit. A repeat cycle is
        reported only if every element has one and they all agree within
        `consistency_threshold` of each other (the most recent element's is
        reported). Propagation does not depend on this agreement: times
        before the first and after the last element's epoch use those
        elements' own repeat cycles.

        Args:
            max_delta_position (float | None): the maximum difference in position (m) allowed for a repeat.
            max_delta_velocity (float | None): the maximum difference in velocity (m/s) allowed for a repeat.
            max_search_duration (timedelta | None): the maximum period of time to search for repeats.
            lazy_load (bool | None): True, if the previously-computed repeat cycle should be loaded.
            consistency_threshold (timedelta | None): the maximum allowed spread between elements' repeat cycles.
            max_delta_semimajor_axis (float | None): the maximum difference (m) between the
                semimajor axis and that of an exact repeat for a candidate repeat.

        Returns:
            timedelta: the repeat cycle duration (if every element agrees on one)
        """
        if self.repeat_cycle is not None:
            return self._get_repeat_track(0).repeat_cycle  # type: ignore
        if consistency_threshold is None:
            consistency_threshold = timedelta(
                seconds=config.get_rc().repeat_cycle_consistency_threshold_s
            )
        cycles = []
        for element in self.elements:
            cycle = element.get_repeat_cycle(
                max_delta_position,
                max_delta_velocity,
                max_search_duration,
                lazy_load,
                max_delta_semimajor_axis,
            )
            if cycle is None:
                return None
            cycles.append(cycle)
            if max(cycles) - min(cycles) > consistency_threshold:
                return None
        return cycles[-1]

    def _get_approximate_repeat_cycle(self, index: int) -> timedelta | None:
        """
        Gets the approximate repeat cycle of an element: the orbit's declared
        repeat cycle, if any, or otherwise the element's computed one.

        Args:
            index (int): the index of the element

        Returns:
            timedelta | None: the approximate repeat cycle, if any
        """
        if self.repeat_cycle is not None:
            return self.repeat_cycle
        return self.elements[index].get_repeat_cycle()

    def get_repeat_element(self, index: int = 0) -> GeneralPerturbationsElements | None:
        """
        Gets the specified element maintained on its repeat ground track (see
        `GeneralPerturbationsElements.get_repeat_element`): without drag,
        with its mean motion adjusted so that a whole number of orbits spans
        the repeat cycle exactly, and, for a sun-synchronous orbit, with its
        inclination adjusted so that its nodal day is a mean solar day. The
        first and last elements' are used to propagate times before the
        first and after the last element's epoch, repeated with the repeat
        cycle.

        Args:
            index (int): the index of the element

        Returns:
            GeneralPerturbationsElements | None: the maintained element, if
                the orbit declares a repeat cycle or the element has one
        """
        repeat_cycle = self._get_approximate_repeat_cycle(index)
        if repeat_cycle is None:
            return None
        return self.elements[index].get_repeat_element(repeat_cycle)

    def _get_repeat_track(self, index: int) -> _RepeatTrack | None:
        """
        Gets the repeat track of the specified element maintained on its
        repeat ground track (see `get_repeat_element`), which is cached.

        Args:
            index (int): the index of the element

        Returns:
            _RepeatTrack | None: the repeat track, if the element has a repeat cycle
        """
        repeat_cycle = self._get_approximate_repeat_cycle(index)
        if repeat_cycle is None:
            return None
        element = self.elements[index].get_repeat_element(repeat_cycle)
        return get_cached(
            element,
            "repeat_track",
            repeat_cycle,
            lambda: _RepeatTrack(element, element.refine_repeat_cycle(repeat_cycle)),
        )

    def _get_repeat_tracks(
        self, try_repeat: bool
    ) -> tuple[_RepeatTrack | None, _RepeatTrack | None]:
        """
        Gets the repeat tracks used before the first element's epoch and
        after the last element's epoch (the same for a single element).

        Args:
            try_repeat (bool): True, if repeat tracks should be used.

        Returns:
            tuple[_RepeatTrack | None, _RepeatTrack | None]: the first and
                last elements' repeat tracks, if used
        """
        if not try_repeat:
            return None, None
        first = self._get_repeat_track(0)
        if len(self.elements) == 1:
            return first, first
        return first, self._get_repeat_track(-1)

    def _propagate(
        self,
        source: int,
        t: Time,
        first: _RepeatTrack | None,
        last: _RepeatTrack | None,
    ) -> Geocentric:
        """
        Propagates times from one source: an element index (propagated
        directly) or a repeat track code (`_FIRST_REPEAT` or `_LAST_REPEAT`).

        Args:
            source (int): The source.
            t (skyfield.timelib.Time): time(s) at which to compute position/velocity.
            first (_RepeatTrack | None): The first element's repeat track.
            last (_RepeatTrack | None): The last element's repeat track.

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        if source == _FIRST_REPEAT:
            return first.at(t)  # type: ignore
        if source == _LAST_REPEAT:
            return last.at(t)  # type: ignore
        return self.elements[source].to_skyfield(self.remove_drag).at(t)  # type: ignore

    def get_orbit_track_at_time(
        self, t: Time, try_repeat: bool | None = None
    ) -> Geocentric:
        """
        Gets the orbit track of this orbit at given Skyfield time(s), in the
        inertial (GCRS) frame.

        Each time is propagated with whichever element's epoch is closest to
        it (see `get_closest_element_index`), so a vectorized `t` may draw
        from different elements for different entries. If `try_repeat`,
        times before the first element's epoch and after the last element's
        epoch are instead propagated with that element's repeat track, if it
        has a repeat cycle: the satellite's Earth-fixed position and velocity
        are those of the element maintained on its repeat ground track at
        the time shifted by whole repeat cycles to within one repeat cycle of
        its epoch (see `get_repeat_element`), expressed in the inertial frame
        at the time itself. For an orbit with a single element, this applies
        to all times.

        Prefer this method over `get_orbit_track` when a Skyfield `Time` is
        already in hand (e.g. while iterating a Skyfield search such as
        `skyfield.searchlib.find_discrete`): converting a `Time` to Python
        `datetime` objects and back discards the per-instant quantities
        Skyfield caches on it (such as nutation angles).

        Args:
            t (skyfield.timelib.Time): time(s) at which to compute position/velocity.
            try_repeat (bool | None): True, if repeat tracks should be used
                (by default, the `repeat_cycle_for_orbit_track` setting).

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        if try_repeat is None:
            try_repeat = config.get_rc().repeat_cycle_for_orbit_track
        first, last = self._get_repeat_tracks(try_repeat)
        # the source of each time: the closest element's index, or a repeat track
        sources = np.asarray(self.get_closest_element_index(t))
        if last is not None:
            sources = np.where(t.tt >= last.epoch_time.tt, _LAST_REPEAT, sources)
        if first is not None:
            sources = np.where(
                t.tt < first.epoch_time.tt,
                _LAST_REPEAT if first is last else _FIRST_REPEAT,
                sources,
            )
        if try_repeat:
            self._warn_if_propagating_with_drag(t, sources)
        unique = np.unique(sources)
        if len(unique) == 1:
            return self._propagate(int(unique[0]), t, first, last)
        position_au = np.empty((3,) + t.shape)
        velocity_au_per_d = np.empty((3,) + t.shape)
        # Skyfield caches per-instant quantities on a `Time` object, but
        # slicing a `Time` does not carry them over. Compute the costly ones
        # used to rotate SGP4 (TEME) results into GCRS (sidereal time and the
        # precession-nutation matrix) once for all of `t` and share them with
        # each slice, which also leaves them cached on `t` for later frame
        # conversions (e.g. to ITRS).
        gast, precession_nutation = t.gast, t.M
        for source in unique:
            # propagate each source across all its assigned times in one
            # vectorized call, rather than one time at a time
            mask = sources == source
            t_mask = t[mask]
            t_mask.gast = gast[mask]
            t_mask.M = precession_nutation[:, :, mask]
            track = self._propagate(int(source), t_mask, first, last)
            position_au[:, mask] = track.position.au
            velocity_au_per_d[:, mask] = track.velocity.au_per_d
        return Geocentric(position_au, velocity_au_per_d, t)

    def get_orbit_track(
        self, times: datetime | list[datetime], try_repeat: bool | None = None
    ) -> Geocentric:
        """
        Gets the orbit track of this orbit at given time(s), in the inertial
        (GCRS) frame (see `get_orbit_track_at_time`).

        Args:
            times (datetime | list[datetime]): time(s) at which to compute position/velocity.
            try_repeat (bool | None): True, if repeat tracks should be used
                (by default, the `repeat_cycle_for_orbit_track` setting).

        Returns:
            skyfield.positionlib.Geocentric: the orbit track position/velocity
        """
        return self.get_orbit_track_at_time(_to_time(times), try_repeat)

    def get_geographic_position_at_time(
        self, t: Time, try_repeat: bool | None = None
    ) -> GeographicPosition:
        """
        Gets the geodetic (WGS84) position of this orbit at given Skyfield
        time(s), in an Earth-fixed frame (see `get_orbit_track_at_time`).

        Args:
            t (skyfield.timelib.Time): time(s) at which to compute geodetic position.
            try_repeat (bool | None): True, if repeat tracks should be used
                (by default, the `repeat_cycle_for_orbit_track` setting).

        Returns:
            skyfield.toposlib.GeographicPosition: the geodetic position
        """
        return wgs84.geographic_position_of(self.get_orbit_track_at_time(t, try_repeat))

    def get_geographic_position(
        self, times: datetime | list[datetime], try_repeat: bool | None = None
    ) -> GeographicPosition:
        """
        Gets the geodetic (WGS84) position of this orbit at given time(s), in
        an Earth-fixed frame (see `get_orbit_track_at_time`).

        Args:
            times (datetime | list[datetime]): time(s) at which to compute geodetic position.
            try_repeat (bool | None): True, if repeat tracks should be used
                (by default, the `repeat_cycle_for_orbit_track` setting).

        Returns:
            skyfield.toposlib.GeographicPosition: the geodetic position
        """
        return self.get_geographic_position_at_time(_to_time(times), try_repeat)

    def _warn_if_propagating_with_drag(self, t: Time, sources: npt.NDArray) -> None:
        """
        Warns if elements are propagated directly with drag to times farther
        from their epoch than the `repeat_cycle_search_duration_days`
        setting: drag lowers the propagated orbit, which then diverges from a
        maintained orbit.

        Args:
            t (skyfield.timelib.Time): The propagated times.
            sources (numpy.typing.NDArray): The source of each time: an
                element index if propagated directly (see `get_orbit_track_at_time`).
        """
        if self.remove_drag or not any(el.has_drag for el in self.elements):
            return
        sources = np.atleast_1d(sources)
        direct = sources >= 0
        if not np.any(direct):
            return
        offsets = (
            np.atleast_1d(utils.to_datetime64_ns(t))[direct]
            - self._get_element_epochs_ns()[sources[direct]]
        )
        limit = timedelta(days=config.get_rc().repeat_cycle_search_duration_days)
        if np.max(np.abs(offsets)) > np.timedelta64(limit):
            warnings.warn(
                "Propagating general perturbations elements with drag, "
                + "without a repeat cycle, more than "
                + f"{limit.total_seconds() / 86400:g} days from their "
                + "epoch: the propagated orbit decays and can diverge "
                + "from an orbit maintained against drag. To model a "
                + "maintained orbit, set `remove_drag=True` and declare "
                + "its `repeat_cycle`.",
                stacklevel=3,
            )

    def _get_boundary_events(
        self,
        topos: GeographicPosition,
        boundaries: list[datetime],
        min_elevation_angle: float,
        try_repeat: bool,
    ) -> list[tuple[datetime, int]]:
        """
        Gets the rise and set events at boundaries where the propagation
        switches between sources (elements, repeat tracks, or repeat cycles),
        whose orbit tracks differ slightly. A satellite above the minimum
        elevation angle just before a boundary but not just after it sets at
        the boundary (or, conversely, rises), although `find_events`, which
        searches each source separately, reports neither.

        Args:
            topos (skyfield.toposlib.GeographicPosition): Target location to observe.
            boundaries (list[datetime]): The boundaries.
            min_elevation_angle (float): Minimum elevation angle (deg) to constrain observation.
            try_repeat (bool): True, if repeat tracks are used.

        Returns:
            list[tuple[datetime, int]]: event times and their rise (0) / set (2) codes
        """
        if len(boundaries) == 0:
            return []
        epsilon = timedelta(milliseconds=1)
        t = _to_time([b + d for b in boundaries for d in (-epsilon, epsilon)])
        altitude = (self.get_orbit_track_at_time(t, try_repeat) - topos.at(t)).altaz()[
            0
        ]
        up = np.reshape(altitude.degrees >= min_elevation_angle, (-1, 2))
        return [
            (boundary, 2 if before else 0)
            for boundary, (before, after) in zip(boundaries, up)
            if before != after
        ]

    def get_observation_events(
        self,
        point: Point,
        start: datetime,
        end: datetime,
        min_elevation_angle: float,
        try_repeat: bool | None = None,
    ) -> tuple[Time, npt.NDArray]:
        """
        Gets the observation events (rise/culminate/set) of this orbit
        with respect to a ground point, between `start` and `end`, using
        Skyfield's `find_events` with refined rise and set times.

        Times are propagated as by `get_orbit_track_at_time`: the parts of
        the period before the first element's epoch and after the last
        element's epoch, if `try_repeat` and the element has a repeat cycle,
        with its repeat track (the events of the repeat cycle just before or
        after the epoch are computed once and repeated, shifted by whole
        repeat cycles, to cover the period), and the rest by partitioning
        the period by whichever element's epoch is closest at each time (see
        `partition_by_element_index`). Where a pass spans a switch between
        these sources, whose orbit tracks differ slightly, the satellite may
        rise or set at the switch (see `_get_boundary_events`).

        Args:
            point (Point): Target location to observe.
            start (datetime): Start time of the observation period.
            end (datetime): End time of the observation period.
            min_elevation_angle (float): Minimum elevation angle (deg) to constrain observation.
            try_repeat (bool | None): True, if repeat tracks should be used (by
                default, the `repeat_cycle_for_observation_events` setting).

        Returns:
            tuple[skyfield.timelib.Time, numpy.ndarray]: event times and their rise (0) / culminate (1) / set (2) codes
        """
        if try_repeat is None:
            try_repeat = config.get_rc().repeat_cycle_for_observation_events
        topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
        first, last = self._get_repeat_tracks(try_repeat)
        # events, and boundaries between the sources that propagate them
        events, boundaries = [], []
        if first is not None and start < first.epoch:
            events.extend(
                first.find_events(
                    topos, start, min(end, first.epoch), min_elevation_angle
                )
            )
            boundaries.extend(first.get_boundaries(start, min(end, first.epoch)))
        direct_start = start if first is None else max(start, first.epoch)
        direct_end = end if last is None else min(end, last.epoch)
        if first is not None and first is not last and start < first.epoch < end:
            boundaries.append(first.epoch)
        if last is not None and first is not last and start < last.epoch < end:
            boundaries.append(last.epoch)
        if direct_start < direct_end:
            if try_repeat:
                self._warn_if_propagating_with_drag(
                    _to_time([direct_start, direct_end]),
                    np.array(
                        self.get_closest_element_index([direct_start, direct_end])
                    ),
                )
            parts, indices = self.partition_by_element_index(direct_start, direct_end)
            boundaries.extend(parts[1:-1])
            for i, index in enumerate(indices):
                times, codes = _find_events(
                    self.elements[index].to_skyfield(self.remove_drag),
                    topos,
                    constants.timescale.from_datetime(parts[i]),
                    constants.timescale.from_datetime(parts[i + 1]),
                    min_elevation_angle,
                )
                if len(codes) > 0:
                    events.extend(
                        zip(np.atleast_1d(times.utc_datetime()), np.atleast_1d(codes))
                    )
        if last is not None and end > last.epoch:
            events.extend(
                last.find_events(
                    topos, max(start, last.epoch), end, min_elevation_angle
                )
            )
            boundaries.extend(last.get_boundaries(max(start, last.epoch), end))
        events.extend(
            self._get_boundary_events(
                topos, boundaries, min_elevation_angle, try_repeat
            )
        )
        # sort by time, removing duplicates at the ends of repeated cycles
        events = sorted(set((t, int(code)) for t, code in events))
        if len(events) == 0:
            return constants.timescale.tt_jd(np.array([])), np.array([], dtype=int)
        return (
            constants.timescale.from_datetimes([t for t, _ in events]),
            np.array([code for _, code in events], dtype=int),
        )

    def to_gp_orbit(self, lazy_load: bool | None = None) -> GeneralPerturbationsOrbit:
        """
        Converts this orbit to a general perturbations orbit representation,
        which it already is.

        Args:
            lazy_load (bool | None): accepted, but has no effect, for
                interface parity with `OrbitBase.to_gp_orbit`.

        Returns:
            GeneralPerturbationsOrbit: this orbit, unchanged
        """
        return self
