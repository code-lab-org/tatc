"""
Time conversion utility functions, and the conversion and indexing of
Skyfield times.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import overload

import numpy as np
import numpy.typing as npt
from skyfield.positionlib import Geocentric
from skyfield.timelib import Time

from .. import constants
from .earth_orientation import _interpolate_nutation


@overload
def to_datetime64_ns(value: datetime) -> np.datetime64: ...


@overload
def to_datetime64_ns(
    value: list[datetime] | npt.NDArray[np.datetime64],
) -> npt.NDArray[np.datetime64]: ...


@overload
def to_datetime64_ns(value: Time) -> np.datetime64 | npt.NDArray[np.datetime64]: ...


def to_datetime64_ns(
    value: datetime | list[datetime] | npt.NDArray[np.datetime64] | Time,
) -> np.datetime64 | npt.NDArray[np.datetime64]:
    """
    Converts a datetime (or a list/array of datetimes, or a Skyfield time)
    to numpy datetime64[ns]. Timezone-aware datetimes are normalized to
    naive UTC first because numpy deprecates implicit timezone-aware
    conversion.

    A Skyfield `Time` is converted from its UTC calendar fields in
    vectorized form, which is much faster for large arrays than converting
    each time to a Python `datetime` (e.g. via `Time.utc_datetime`). Times
    within a leap second carry over into the following minute.

    Args:
        value (datetime | list[datetime] | npt.NDArray[np.datetime64] | Time): the time(s) to convert

    Returns:
        np.datetime64 | npt.NDArray[np.datetime64]: the converted datetime64 value(s)
    """
    if isinstance(value, Time):
        year, month, day, hour, minute, second = value.utc
        months = (np.asarray(year) - 1970).astype("datetime64[Y]") + (
            np.asarray(month) - 1
        ).astype("timedelta64[M]")
        result = (
            months.astype("datetime64[D]").astype("datetime64[ns]")
            + ((np.asarray(day) - 1) * 86400 + np.asarray(hour) * 3600).astype(
                "timedelta64[s]"
            )
            + np.asarray(minute).astype("timedelta64[m]")
            + np.round(np.asarray(second) * 1e9).astype("timedelta64[ns]")
        )
        return result[()] if result.ndim == 0 else result
    if isinstance(value, datetime):
        return np.datetime64(value.astimezone(timezone.utc).replace(tzinfo=None), "ns")
    if isinstance(value, np.ndarray) and np.issubdtype(value.dtype, np.datetime64):
        return value.astype("datetime64[ns]")
    return np.array(
        [
            (
                v.astimezone(timezone.utc).replace(tzinfo=None)
                if isinstance(v, datetime)
                else v
            )
            for v in value
        ],
        dtype="datetime64[ns]",
    )


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
    Computes them for all of `t` if not already cached (with interpolated
    nutation angles, see `_interpolate_nutation`).

    Args:
        t (skyfield.timelib.Time): The time(s).
        index (numpy.typing.ArrayLike): The index (an integer array or a
            boolean mask).

    Returns:
        skyfield.timelib.Time: the indexed time(s)
    """
    _interpolate_nutation(t)
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
