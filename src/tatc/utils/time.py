"""
Time conversion utility functions.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import overload

import numpy as np
import numpy.typing as npt


@overload
def to_datetime64_ns(value: datetime) -> np.datetime64: ...


@overload
def to_datetime64_ns(
    value: list[datetime] | npt.NDArray[np.datetime64],
) -> npt.NDArray[np.datetime64]: ...


def to_datetime64_ns(
    value: datetime | list[datetime] | npt.NDArray[np.datetime64],
) -> np.datetime64 | npt.NDArray[np.datetime64]:
    """
    Converts a datetime (or a list/array of datetimes) to numpy
    datetime64[ns]. Timezone-aware datetimes are normalized to naive UTC
    first because numpy deprecates implicit timezone-aware conversion.

    Args:
        value (datetime | list[datetime] | npt.NDArray[np.datetime64]): the datetime(s) to convert

    Returns:
        np.datetime64 | npt.NDArray[np.datetime64]: the converted datetime64 value(s)
    """
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
