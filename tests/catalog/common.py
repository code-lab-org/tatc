"""
Helpers shared by the tests of catalog fixtures.

The fixture tests compare values computed now with values recorded by the
validation scripts to a fraction of a millisecond, which is fine enough to see
the difference between UT1 tables. TAT-C's timescale uses the table bundled with
the installed skyfield, which changes between skyfield releases and is
extrapolated past its last date, so these tests compute with the table the
records were made with (skyfield 1.55's, saved as
`fixtures/earth_orientation.npz` in skyfield's own encoding).

Re-recording fixture metrics with a different table (a different skyfield, or a
validation day past 23 January 2027, where this table ends) needs a new snapshot,
written from the skyfield used to record, and every fixture record rescored.
"""

import pathlib
from unittest import mock

import numpy as np
from skyfield.timelib import build_delta_t

from tatc.constants import timescale

EARTH_ORIENTATION = pathlib.Path(__file__).parent / "fixtures" / "earth_orientation.npz"


def _load_recorded_delta_t() -> tuple:
    """
    Returns the daily TT dates (Julian days) and delta T values (seconds) of the
    recorded table, decoded as skyfield decodes its bundled `iers.npz`.
    """
    with np.load(EARTH_ORIENTATION) as arrays:
        daily_tt = arrays["tt_jd_minus_arange"] + np.arange(
            len(arrays["tt_jd_minus_arange"])
        )
        daily_delta_t = (arrays["delta_t_1e7"] / 1e7).round(7)
    return daily_tt, daily_delta_t


def recorded_earth_orientation() -> mock._patch:
    """
    Returns a patcher that makes TAT-C's shared timescale use the UT1 table the
    catalog records were made with. Start it in `setUpModule` and stop it in
    `tearDownModule`. It changes the timescale object in place, because modules
    import it by name.
    """
    table = _load_recorded_delta_t()
    return mock.patch.multiple(
        timescale, delta_t_table=table, delta_t_function=build_delta_t(table)
    )
