"""
Orbit propagation utilities for general perturbations (SGP4) orbits:
observation event search, propagation on a repeat ground track, and repeat
cycle search.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Generator
from datetime import datetime, timedelta, timezone
from typing import TYPE_CHECKING, Any, NamedTuple

import numpy as np
import numpy.typing as npt
from skyfield.api import Time
from skyfield.framelib import itrs
from skyfield.nutationlib import iau2000a_radians
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


def _to_time_from_offsets(reference: datetime, seconds: npt.ArrayLike) -> Time:
    """
    Converts offsets from a reference time to a Skyfield `Time`, without
    building a `datetime` for each offset: equivalent to `_to_time` of
    `reference + timedelta(seconds=x)` for each offset `x` (as UTC calendar
    arithmetic, so each time's leap second offset is that of its UTC day).

    Args:
        reference (datetime): The (timezone-aware) reference time.
        seconds (numpy.typing.ArrayLike): The offsets (seconds).

    Returns:
        skyfield.timelib.Time: the Skyfield time(s)
    """
    utc = reference.astimezone(timezone.utc)
    # split each time into whole days after the reference's UTC day and
    # seconds of its UTC day (Skyfield's UTC calendar dates allow days
    # beyond the end of the month)
    days, second = np.divmod(
        utc.hour * 3600
        + utc.minute * 60
        + utc.second
        + utc.microsecond / 1e6
        + np.asarray(seconds, dtype=float),
        86400,
    )
    return constants.timescale.utc(
        utc.year, utc.month, utc.day + days.astype(int), 0, 0, second
    )


def _index_time(t: Time, index: npt.ArrayLike) -> Time:
    """
    Indexes a Skyfield `Time`, carrying over its sidereal time and
    precession-nutation matrix: Skyfield caches these costly per-instant
    quantities (used to convert between the inertial and Earth-fixed
    frames) on a `Time`, but indexing a `Time` does not carry them over.
    Computes them for all of `t` if not already cached.

    Args:
        t (skyfield.timelib.Time): The time(s).
        index (numpy.typing.ArrayLike): The index (an integer array or a
            boolean mask).

    Returns:
        skyfield.timelib.Time: the indexed time(s)
    """
    gast, precession_nutation = t.gast, t.M
    indexed = t[index]
    indexed.gast = gast[index]
    indexed.M = precession_nutation[:, :, index]
    if "_nutation_angles_radians" in vars(t):
        # and the nutation angles, if set or computed (see _interpolate_nutation)
        indexed._nutation_angles_radians = tuple(  # pylint: disable=protected-access
            np.asarray(angle)[index]
            for angle in t._nutation_angles_radians  # pylint: disable=protected-access
        )
    return indexed


def _index_orbit_track(orbit_track: Geocentric, index: npt.ArrayLike) -> Geocentric:
    """
    Indexes an orbit track, carrying over the per-instant quantities cached
    on its times (see `_index_time`).

    Args:
        orbit_track (skyfield.positionlib.Geocentric): The orbit track.
        index (numpy.typing.ArrayLike): The index (an integer array or a
            boolean mask).

    Returns:
        skyfield.positionlib.Geocentric: the indexed orbit track
    """
    return Geocentric(
        orbit_track.position.au[:, index],
        orbit_track.velocity.au_per_d[:, index],
        _index_time(orbit_track.t, index),
    )


TimeRequest = Generator[Time | list[Time], None, Any]
"""
A computation that yields each Skyfield `Time` (or list of `Time`s) it
creates before using it, so that the costly Earth orientation quantities of
the times of several computations can be computed together (see
`_run_together`), and returns its result.
"""


_NUTATION_TABLES: dict[tuple[int, float], tuple[npt.NDArray, npt.NDArray]] = {}
"""
The IAU 2000A nutation angles (radians) at each step of a TT day (with both
ends), by the number of steps per day and the day (a whole TT Julian date).
"""

_NUTATION_INTERPOLATION_VERIFIED: bool | None = None
"""
Whether interpolated nutation angles are verified to be used by Skyfield as
expected (see `_verify_nutation_interpolation`), once checked.
"""


def _interpolate_nutation(t: Time) -> None:
    """
    Sets the nutation angles of a Skyfield time, if not yet computed, by
    linear interpolation from a cached table of the IAU 2000A angles at
    steps of the `nutation_interpolation_minutes` runtime configuration
    (unless it is None). Skyfield computes the full IAU 2000A series for
    every time, at a cost of about 20 microseconds per time, which
    dominates the propagation of orbits at many times; the nutation angles
    change slowly (their shortest significant periods are days), so
    interpolating them at 15 minute steps changes them by about a
    microarcsecond, at a small fraction of the cost.

    This sets Skyfield's private `Time._nutation_angles_radians` attribute
    (as Skyfield's own `almanac` module does to use the IAU 2000B model), so
    it is verified once (see `_verify_nutation_interpolation`): if Skyfield
    no longer uses that attribute as expected, nutation angles are computed
    by Skyfield as usual, with a warning. To compute them for every time,
    set the runtime configuration to None.

    Args:
        t (skyfield.timelib.Time): The time(s).
    """
    minutes = config.get_rc().nutation_interpolation_minutes
    if (
        minutes is None
        or "_nutation_angles_radians" in vars(t)
        or ("gast" in vars(t) and "M" in vars(t))
        or np.size(t.tt) == 0
        or not _check_nutation_interpolation()
    ):
        return
    _set_interpolated_nutation(t, max(1, int(round(1440 / minutes))))


def _set_interpolated_nutation(t: Time, steps: int) -> None:
    """
    Sets the nutation angles of a Skyfield time by linear interpolation from
    a cached table at a number of steps per day (see `_interpolate_nutation`).

    Args:
        t (skyfield.timelib.Time): The time(s).
        steps (int): The number of steps per day.
    """
    whole = np.reshape(np.asarray(t.whole, dtype=float), -1)
    fraction = np.reshape(np.asarray(t.tt_fraction, dtype=float), -1)
    whole, fraction = np.broadcast_arrays(whole, fraction)
    day = np.floor(whole + fraction)
    # position within the day, in steps
    position = ((whole - day) + fraction) * steps
    index = np.clip(np.floor(position).astype(int), 0, steps - 1)
    weight = position - index
    days, row = np.unique(day, return_inverse=True)
    missing = [d for d in days if (steps, d) not in _NUTATION_TABLES]
    if len(missing) > 0:
        d_psi, d_eps = iau2000a_radians(
            constants.timescale.tt_jd(
                np.repeat(missing, steps + 1),
                np.tile(np.arange(steps + 1) / steps, len(missing)),
            )
        )
        for k, d in enumerate(missing):
            part = slice(k * (steps + 1), (k + 1) * (steps + 1))
            _NUTATION_TABLES[(steps, d)] = (d_psi[part], d_eps[part])
    table_psi = np.array([_NUTATION_TABLES[(steps, d)][0] for d in days])
    table_eps = np.array([_NUTATION_TABLES[(steps, d)][1] for d in days])
    shape = np.shape(t.tt)
    t._nutation_angles_radians = tuple(  # pylint: disable=protected-access
        np.reshape(
            table[row, index] * (1 - weight) + table[row, index + 1] * weight, shape
        )[()]
        for table in (table_psi, table_eps)
    )


def _check_nutation_interpolation() -> bool:
    """
    Checks, once, whether Skyfield uses interpolated nutation angles as
    expected (see `_verify_nutation_interpolation`), warning if not.

    Returns:
        bool: True, if it does
    """
    global _NUTATION_INTERPOLATION_VERIFIED  # pylint: disable=global-statement
    if _NUTATION_INTERPOLATION_VERIFIED is None:
        _NUTATION_INTERPOLATION_VERIFIED = _verify_nutation_interpolation()
        if not _NUTATION_INTERPOLATION_VERIFIED:
            warnings.warn(
                "Skyfield no longer uses the nutation angles set on a time "
                "(`Time._nutation_angles_radians`) as expected, so they are "
                "computed for every time rather than interpolated (see the "
                "`nutation_interpolation_minutes` runtime configuration).",
                stacklevel=3,
            )
    return _NUTATION_INTERPOLATION_VERIFIED


def _verify_nutation_interpolation() -> bool:
    """
    Verifies that Skyfield uses the nutation angles set on a time
    (`Time._nutation_angles_radians`) as expected: that setting them
    changes its sidereal time and precession-nutation matrix, and that
    interpolated angles (see `_set_interpolated_nutation`) reproduce those
    computed for every time, to well within a milliarcsecond.

    Returns:
        bool: True, if verified
    """
    try:
        jd = 2461041.5 + np.array([0.0, 0.2913, 0.5, 0.75, 0.9999])
        exact = constants.timescale.tt_jd(jd)
        # setting the angles must change the results
        perturbed = constants.timescale.tt_jd(jd)
        d_psi, d_eps = iau2000a_radians(perturbed)
        perturbed._nutation_angles_radians = (  # pylint: disable=protected-access
            d_psi + 1e-6,
            d_eps + 1e-6,
        )
        if (
            np.max(np.abs(perturbed.gast - exact.gast)) < 1e-9
            or np.max(np.abs(perturbed.M - exact.M)) < 1e-9
        ):
            return False
        # and interpolated angles must reproduce the results (to 1e-9 hours
        # of sidereal time and 1e-9 in the matrix, about 0.2 milliarcseconds)
        interpolated = constants.timescale.tt_jd(jd)
        _set_interpolated_nutation(interpolated, 96)
        return bool(
            np.max(np.abs(interpolated.gast - exact.gast)) < 1e-9
            and np.max(np.abs(interpolated.M - exact.M)) < 1e-9
        )
    except Exception:  # pylint: disable=broad-exception-caught
        return False


def _share_earth_orientation(times: list[Time]) -> None:
    """
    Computes the sidereal time and precession-nutation matrix of several
    Skyfield times together and caches them on each (see `_index_time`):
    Skyfield caches these costly per-instant quantities on each `Time`, but
    computes them separately for every `Time`, at a cost dominated by a large
    fixed overhead for each. Times that already have them are skipped.

    Args:
        times (list[skyfield.timelib.Time]): The times.
    """
    unique = {id(t): t for t in times if "gast" not in vars(t) or "M" not in vars(t)}
    times = [t for t in unique.values() if np.size(t.tt) > 0]
    if len(times) < 2:
        # computed as needed
        return
    combined = constants.timescale.tt_jd(
        np.concatenate([np.reshape(t.whole, -1) for t in times]),
        np.concatenate([np.reshape(t.tt_fraction, -1) for t in times]),
    )
    _interpolate_nutation(combined)
    gast, precession_nutation = combined.gast, combined.M
    offset = 0
    for t in times:
        size = np.size(t.tt)
        t.gast = np.reshape(gast[offset : offset + size], np.shape(t.tt))
        t.M = np.reshape(
            precession_nutation[:, :, offset : offset + size], (3, 3) + np.shape(t.tt)
        )
        offset += size


def _run_together(computations: list[TimeRequest]) -> list[Any]:
    """
    Runs several computations (see `TimeRequest`) together, in steps: at each
    step, the Earth orientation quantities of all times that the computations
    yield are computed together (see `_share_earth_orientation`), and each
    computation then continues to its next time.

    Args:
        computations (list[TimeRequest]): The computations.

    Returns:
        list[Any]: the result of each computation
    """
    results: list[Any] = [None] * len(computations)
    pending: dict[int, Time | list[Time]] = {}

    def advance(i: int, first: bool = False) -> None:
        try:
            pending[i] = next(computations[i]) if first else computations[i].send(None)
        except StopIteration as stop:
            results[i] = stop.value
            pending.pop(i, None)

    for i in range(len(computations)):
        advance(i, first=True)
    while len(pending) > 0:
        _share_earth_orientation(
            [
                t
                for request in pending.values()
                for t in (request if isinstance(request, list) else [request])
            ]
        )
        for i in list(pending):
            advance(i)
    return results


def _run(computation: TimeRequest) -> Any:
    """
    Runs a single computation (see `TimeRequest`).

    Args:
        computation (TimeRequest): The computation.

    Returns:
        Any: its result
    """
    return _run_together([computation])[0]


def _value(value: Any) -> TimeRequest:
    """
    Wraps a value as a computation (see `TimeRequest`) that yields no times.
    """
    return value
    yield  # pylint: disable=unreachable


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
    return _run(_bisect_steps(lambda x: _value(excess(x)), lower, upper))


def _bisect_steps(
    excess: Callable[[npt.NDArray], TimeRequest],
    lower: npt.NDArray,
    upper: npt.NDArray,
) -> TimeRequest:
    """
    Bisects brackets of a change of sign of a function to a millisecond, as
    a computation (see `_bisect` and `TimeRequest`) of a function that is
    itself a computation. Each bracket is bisected only until it is
    narrower than a millisecond.
    """
    lower = np.array(lower, dtype=float)
    upper = np.array(upper, dtype=float)
    if len(lower) == 0:
        return lower
    f_lower = np.array((yield from excess(lower)), dtype=float)
    for _ in range(64):
        active = np.flatnonzero(upper - lower >= 1e-3 / 86400)
        if len(active) == 0:
            break
        middle = (lower[active] + upper[active]) / 2
        f_middle = yield from excess(middle)
        same = np.sign(f_middle) == np.sign(f_lower[active])
        lower[active] = np.where(same, middle, lower[active])
        f_lower[active] = np.where(same, f_middle, f_lower[active])
        upper[active] = np.where(same, upper[active], middle)
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
    return _run(_find_crossing_steps(lambda x: _value(excess(x)), start, stop))


def _find_crossing_steps(
    excess: Callable[[npt.NDArray], TimeRequest], start: float, stop: float
) -> TimeRequest:
    """
    Finds the first crossing of a function, as a computation (see
    `_find_crossing` and `TimeRequest`) of a function that is itself a
    computation.
    """
    samples = np.linspace(start, stop, 26)[:-1]
    negative = np.flatnonzero((yield from excess(samples)) < 0)
    if len(negative) == 0:
        return None
    bracket = np.sort(samples[negative[0] - 1 : negative[0] + 1])
    return float((yield from _bisect_steps(excess, bracket[:1], bracket[1:]))[0])


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
    minimum, and usually within a minute before the reported time, so the
    narrower of these brackets that contains it is bisected to a
    millisecond.

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
    return _run(_find_events_steps(satellite, topos, t_0, t_1, min_elevation_angle))


def _find_events_steps(
    satellite: EarthSatellite,
    topos: GeographicPosition,
    t_0: Time,
    t_1: Time,
    min_elevation_angle: float,
) -> TimeRequest:
    """
    Finds the rise, culminate, and set events of a satellite with respect to
    a ground position, as a computation (see `_find_events` and
    `TimeRequest`). Skyfield's `find_events` itself is not shared.
    """
    times, events = satellite.find_events(topos, t_0, t_1, min_elevation_angle)
    jd = np.array(times.tt, dtype=float, ndmin=1)
    events = np.array(events, dtype=int, ndmin=1)
    relative = satellite - topos

    def excess(x: npt.NDArray) -> TimeRequest:
        t = constants.timescale.tt_jd(x)
        yield t
        return relative.at(t).altaz()[0].degrees - min_elevation_angle

    refine = np.flatnonzero(events != 1)
    if len(refine) > 0:
        lower = np.concatenate(([t_0.tt], jd))[refine]
        upper = jd[refine]
        # Skyfield reports the later end of its last bracket, so the change
        # is usually within a minute before the reported time: bracket it
        # there if possible, or otherwise from the preceding event
        near = np.maximum(lower, upper - 60 / 86400)
        f_near, f_upper = np.split(
            (yield from excess(np.concatenate([near, upper]))), 2
        )
        narrow = np.sign(f_near) != np.sign(f_upper)
        f_lower = np.array(f_near, dtype=float)
        wide = np.flatnonzero(~narrow)
        if len(wide) > 0:
            f_lower[wide] = yield from excess(lower[wide])
        lower = np.where(narrow, near, lower)
        bracketed = np.sign(f_lower) != np.sign(f_upper)
        jd[refine[bracketed]] = yield from _bisect_steps(
            excess, lower[bracketed], upper[bracketed]
        )
    crossing_jd, crossing_events = jd[events != 1], events[events != 1]
    f_0, f_1 = yield from excess(np.array([t_0.tt, t_1.tt]))
    added = []
    if f_0 >= 0 and (len(crossing_events) == 0 or crossing_events[0] == 0):
        if len(crossing_events) > 0 or f_1 < 0:
            stop = crossing_jd[0] if len(crossing_events) > 0 else t_1.tt
            added.append(((yield from _find_crossing_steps(excess, t_0.tt, stop)), 2))
    if f_1 >= 0 and (len(crossing_events) == 0 or crossing_events[-1] == 2):
        if len(crossing_events) > 0 or f_0 < 0:
            stop = crossing_jd[-1] if len(crossing_events) > 0 else t_0.tt
            added.append(((yield from _find_crossing_steps(excess, t_1.tt, stop)), 0))
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
        _interpolate_nutation(shifted)
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
        return _run(self.find_events_steps(topos, start, end, min_elevation_angle))

    def find_events_steps(
        self,
        topos: GeographicPosition,
        start: datetime,
        end: datetime,
        min_elevation_angle: float,
    ) -> TimeRequest:
        """
        Finds the observation events between `start` and `end`, as a
        computation (see `find_events` and `TimeRequest`).
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
        times, codes = yield from _find_events_steps(
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
