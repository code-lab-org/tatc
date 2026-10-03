"""
Time conversion utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import overload

import numpy as np
import numpy.typing as npt
from skyfield.timelib import Time


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
