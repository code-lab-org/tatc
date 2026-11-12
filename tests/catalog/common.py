"""
Helpers shared by the tests of catalog fixtures.

The fixture tests compare values computed now with values recorded by the
validation scripts to a fraction of a millisecond, which is fine enough to see
the difference between UT1 tables. TAT-C's timescale uses the table bundled with
the installed skyfield, which changes between skyfield releases (predicted values
become measured ones) and is extrapolated past its last date. So each fixture
holds the delta T (TT - UT1) its records were made with, sampled around its days
(`delta_t_tt`, `delta_t`, written by `docs/validation/validation_io.py`), and its
tests compute with those samples.

Times outside a fixture's days use the installed table. Only the checks of
published footprint dimensions use such a time (a circular orbit on 1 October
2026), and they compare shares of dimensions within 3%, which UT1 moves by far
too little to change.
"""

from unittest import mock

import numpy as np

from tatc.constants import timescale

# samples further apart than this (days) belong to different days of a fixture
_MAX_SAMPLE_GAP = 0.25


def _build_delta_t_function(sample_tt, sample_delta_t, fallback):
    """
    Returns a delta T function that interpolates a fixture's samples linearly
    within the days they cover, as skyfield interpolates its daily table, and uses
    `fallback` elsewhere.
    """
    sample_tt = np.asarray(sample_tt, dtype=float)
    sample_delta_t = np.asarray(sample_delta_t, dtype=float)

    def delta_t_function(tt):
        tt_array = np.atleast_1d(np.asarray(tt, dtype=float))
        upper = np.searchsorted(sample_tt, tt_array).clip(1, len(sample_tt) - 1)
        inside = (
            (tt_array >= sample_tt[0])
            & (tt_array <= sample_tt[-1])
            & (sample_tt[upper] - sample_tt[upper - 1] <= _MAX_SAMPLE_GAP)
        )
        values = np.interp(tt_array, sample_tt, sample_delta_t)
        if not inside.all():
            values[~inside] = fallback(tt_array[~inside])
        return values if np.ndim(tt) else float(values[0])

    return delta_t_function


def use_fixture_earth_orientation(test_class, fixture) -> None:
    """
    Makes TAT-C's shared timescale use a fixture's delta T samples until the
    tests of a class end. Call it in `setUpClass` with the fixture's arrays. It
    changes the timescale in place, because modules import it by name.
    """
    patcher = mock.patch.object(
        timescale,
        "delta_t_function",
        _build_delta_t_function(
            fixture["delta_t_tt"], fixture["delta_t"], timescale.delta_t_function
        ),
    )
    patcher.start()
    test_class.addClassCleanup(patcher.stop)
