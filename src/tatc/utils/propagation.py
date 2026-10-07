"""
Orbit propagation utilities for general perturbations (SGP4) orbits:
observation event search, propagation on a repeat ground track, and repeat
cycle search.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timedelta
from typing import TYPE_CHECKING, NamedTuple

import numpy as np
import numpy.typing as npt
from skyfield.api import Time
from skyfield.framelib import itrs
from skyfield.positionlib import Geocentric
from skyfield.searchlib import find_minima
from skyfield.sgp4lib import EarthSatellite
from skyfield.toposlib import GeographicPosition

from .. import config, constants

if TYPE_CHECKING:
    from ..schemas.orbit.gp_elements import GeneralPerturbationsElements


class RepeatCycleSearch(NamedTuple):
    """
    Options to search for a repeat cycle, which default to the runtime
    configuration (see `resolve`).
    """

    max_delta_position: float
    """Maximum difference in position (m) allowed for a repeat."""
    max_delta_velocity: float
    """Maximum difference in velocity (m/s) allowed for a repeat."""
    max_search_duration: timedelta
    """Maximum period of time to search for repeats."""
    max_delta_semimajor_axis: float
    """Maximum difference (m) between the semimajor axis and that of an exact repeat."""

    @classmethod
    def resolve(
        cls,
        max_delta_position: float | None = None,
        max_delta_velocity: float | None = None,
        max_search_duration: timedelta | None = None,
        max_delta_semimajor_axis: float | None = None,
    ) -> RepeatCycleSearch:
        """
        Resolves search options, taking those not specified from the runtime
        configuration, so that a cached repeat cycle (keyed by the options)
        is recomputed if the configuration changes.

        Args:
            max_delta_position (float | None): the maximum difference in position (m) allowed for a repeat.
            max_delta_velocity (float | None): the maximum difference in velocity (m/s) allowed for a repeat.
            max_search_duration (timedelta | None): the maximum period of time to search for repeats.
            max_delta_semimajor_axis (float | None): the maximum difference (m) between the
                semimajor axis and that of an exact repeat for a candidate repeat.

        Returns:
            RepeatCycleSearch: the search options
        """
        rc = config.get_rc()
        return cls(
            (
                rc.repeat_cycle_delta_position_m
                if max_delta_position is None
                else max_delta_position
            ),
            (
                rc.repeat_cycle_delta_velocity_m_per_s
                if max_delta_velocity is None
                else max_delta_velocity
            ),
            (
                timedelta(days=rc.repeat_cycle_search_duration_days)
                if max_search_duration is None
                else max_search_duration
            ),
            (
                rc.repeat_cycle_delta_semimajor_axis_m
                if max_delta_semimajor_axis is None
                else max_delta_semimajor_axis
            ),
        )


def _get_closest_return(
    satellite: EarthSatellite, epoch: datetime, center: datetime, period: float
) -> tuple[float, float] | None:
    """
    Gets the differences of Earth-fixed position and velocity between a
    satellite's closest return, within half an orbit of a time, and its
    initial state at an epoch.

    Args:
        satellite (skyfield.api.EarthSatellite): The satellite.
        epoch (datetime): The epoch of the initial state.
        center (datetime): The time near which to search for the closest return.
        period (float): The orbit period (seconds).

    Returns:
        tuple[float, float] | None: the differences of position (m) and
            velocity (m/s), or None if no closest return is found
    """
    # initial position and velocity in the Earth-fixed frame
    position_0, velocity_0 = satellite.at(
        constants.timescale.from_datetime(epoch)
    ).frame_xyz_and_velocity(itrs)
    p_0_m = np.array(position_0.m)
    v_0_m_per_s = np.array(velocity_0.m_per_s)

    def position_error(t: Time) -> npt.NDArray[np.float64]:
        position, _ = satellite.at(t).frame_xyz_and_velocity(itrs)
        return np.linalg.norm((np.array(position.m).T - p_0_m.T).T, axis=0)

    position_error.rough_period = period / 86400  # type: ignore
    window = timedelta(seconds=period / 2)
    times, errors = find_minima(
        constants.timescale.from_datetime(center - window),
        constants.timescale.from_datetime(center + window),
        position_error,
    )
    if len(times) == 0:
        return None
    t_min = times[np.argmin(errors)]
    position, velocity = satellite.at(t_min).frame_xyz_and_velocity(itrs)
    return (
        float(np.linalg.norm(np.array(position.m) - p_0_m)),
        float(np.linalg.norm(np.array(velocity.m_per_s) - v_0_m_per_s)),
    )


def _compute_repeat_element(
    element: GeneralPerturbationsElements, repeat_cycle: timedelta
) -> GeneralPerturbationsElements:
    """
    Computes an element maintained on the repeat ground track with the
    approximate repeat cycle, without caching (see
    `GeneralPerturbationsElements.get_repeat_element`).

    Args:
        element (GeneralPerturbationsElements): The element.
        repeat_cycle (timedelta): The approximate repeat cycle.

    Returns:
        GeneralPerturbationsElements: the maintained element
    """
    maintained = element.without_drag()
    sun_synchronous = maintained.is_sun_synchronous()
    nodal_period, _ = maintained.get_nodal_period_and_day()
    orbits = max(
        1,
        round(
            maintained.refine_repeat_cycle(repeat_cycle).total_seconds() / nodal_period
        ),
    )
    for _ in range(10):
        nodal_period, nodal_day = maintained.get_nodal_period_and_day()
        refined = maintained.refine_repeat_cycle(repeat_cycle).total_seconds()
        update = {}
        residual = orbits * nodal_period - refined
        if abs(residual) >= 1e-6:
            # the nodal period varies (nearly) inversely with mean motion
            update["mean_motion"] = maintained.mean_motion * (1 + residual / refined)
        if sun_synchronous and abs(nodal_day - constants.EARTH_SOLAR_DAY_S) >= 1e-6:
            # the node precession rate is (nearly) proportional to the
            # cosine of inclination
            ratio = (
                constants.EARTH_ROTATION_RATE - 2 * np.pi / constants.EARTH_SOLAR_DAY_S
            ) / (constants.EARTH_ROTATION_RATE - 2 * np.pi / nodal_day)
            update["inclination"] = float(
                np.degrees(
                    np.arccos(
                        np.clip(
                            np.cos(np.radians(maintained.inclination)) * ratio, -1, 1
                        )
                    )
                )
            )
        if not update:
            break
        maintained = maintained.model_copy(update=update)
    return maintained


def _search_repeat_cycle(
    element: GeneralPerturbationsElements, search: RepeatCycleSearch
) -> timedelta | None:
    """
    Searches for an element's repeat cycle, without caching (see
    `GeneralPerturbationsElements.get_repeat_cycle`).

    Args:
        element (GeneralPerturbationsElements): The element.
        search (RepeatCycleSearch): The search options.

    Returns:
        timedelta | None: the repeat cycle duration (if it exists)
    """
    # analytic repeat ground track candidates: how many nodal days (D)
    # are needed for a whole number of orbits (C) to elapse
    drag_free = element.without_drag()
    nodal_period, nodal_day = drag_free.get_nodal_period_and_day()
    orbits_per_day = nodal_day / nodal_period
    max_days = int(search.max_search_duration.total_seconds() / nodal_day)
    days_range = np.arange(1, max_days + 1)
    orbit_counts = np.round(orbits_per_day * days_range)
    residual_orbits = orbits_per_day * days_range - orbit_counts
    # ground-track drift (m) at the equator implied by missing a whole
    # orbit count by residual_orbits
    ground_track_spacing = 2 * np.pi * constants.EARTH_MEAN_RADIUS / orbits_per_day
    drift = np.abs(residual_orbits) * ground_track_spacing
    # difference (m) between the semimajor axis and that of the exact
    # repeat, from the sensitivity of orbits per nodal day to the
    # semimajor axis (by a finite difference of mean motion)
    perturbed = drag_free.model_copy(
        update={"mean_motion": drag_free.mean_motion * (1 + 1e-6)}
    )
    perturbed_period, perturbed_day = perturbed.get_nodal_period_and_day()
    sensitivity = (perturbed_day / perturbed_period - orbits_per_day) / (
        perturbed.get_semimajor_axis() - drag_free.get_semimajor_axis()
    )
    delta_semimajor_axis = (orbit_counts / days_range - orbits_per_day) / sensitivity
    # exact repeats of up to D nodal days are spaced by about 1/D^2
    # orbits per day: apply the semimajor axis tolerance only where they
    # are spaced by at least three times the tolerance, beyond which it
    # would admit chance near-repeats
    short = 3 * search.max_delta_semimajor_axis * abs(sensitivity) * days_range**2 <= 1
    # generous margin: the drift is estimated from secular rates, while
    # the verification below also includes periodic terms
    near_drift = drift < 3 * search.max_delta_position
    near_semimajor_axis = short & (
        np.abs(delta_semimajor_axis) < search.max_delta_semimajor_axis
    )

    def is_confirmed(result: tuple[float, float] | None) -> bool:
        return (
            result is not None
            and result[0] < search.max_delta_position
            and result[1] < search.max_delta_velocity
        )

    for i in np.flatnonzero(near_drift | near_semimajor_axis):
        candidate = timedelta(seconds=int(days_range[i]) * nodal_day)
        maintained = element.get_repeat_element(candidate)
        # verify a candidate near the element's mean motion with the
        # element itself (without drag, as for an orbit maintained against
        # drag), and one near the semimajor axis of the exact repeat with
        # the element maintained on its repeat ground track
        confirmed = near_drift[i] and is_confirmed(
            _get_closest_return(
                drag_free.to_skyfield(),
                element.epoch,
                element.epoch + candidate,
                nodal_period,
            )
        )
        if not confirmed and near_semimajor_axis[i]:
            confirmed = is_confirmed(
                _get_closest_return(
                    maintained.to_skyfield(),
                    element.epoch,
                    element.epoch + maintained.refine_repeat_cycle(candidate),
                    nodal_period,
                )
            )
        if confirmed:
            return maintained.refine_repeat_cycle(candidate)
    return None


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
