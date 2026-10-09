"""
Earth orientation utilities: the IAU 2000A nutation angles of Skyfield
times, interpolated from a cached table (see the
`nutation_interpolation_minutes` runtime configuration), and the sidereal
time and precession-nutation matrix of several Skyfield times computed
together. This is the only module that sets Skyfield's private
`Time._nutation_angles_radians` attribute.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import warnings

import numpy as np
import numpy.typing as npt
from skyfield.nutationlib import iau2000a_radians
from skyfield.timelib import Time

from .. import config, constants

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
