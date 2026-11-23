"""
Helpers shared by the tests of catalog fixtures.

The fixture tests compare values computed now with values recorded by the
validation scripts to a fraction of a millisecond, which is fine enough to see
the difference between UT1 tables. TAT-C's timescale uses the table bundled with
the installed skyfield, which changes between skyfield releases (predicted values
become measured ones) and is extrapolated past its last date. So each fixture
holds the delta T (TT - UT1) its records were made with, sampled around its days
(`delta_t_tt`, `delta_t`, written by `validation/catalog/validation_io/records.py`), and its
tests compute with those samples.

Times outside a fixture's days use the installed table. Only the checks of
published footprint dimensions use such a time (a circular orbit on 1 October
2026), and they compare shares of dimensions within 3%, which UT1 moves by far
too little to change.
"""

from datetime import timedelta
from unittest import mock

import numpy as np
import pandas as pd
import shapely.geometry as geo

from tatc.analysis import collect_observations
from tatc.constants import timescale
from tatc.schemas import Point
from tatc.utils.geometry import hash_geometry

# the window around each recorded pass in which its observations are found
PASS_WINDOW = timedelta(minutes=10)

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


def find_pass_differences(satellite, instrument_index, day, passes) -> np.ndarray:
    """
    Gets the difference (seconds) of the nearest epoch from each recorded pass,
    given as rows of latitude, longitude, and seconds from the start of `day`, or
    NaN if the pass is not observed, as `validation_io.recompute_passes` gets them
    for the window of 10 minutes around each pass. The observations of all passes
    are collected together (observations of a point do not depend on the other
    points), which is much faster than one window at a time.

    Args:
        satellite: The observing satellite.
        instrument_index: The index of the observing instrument, `None` for every
            instrument, or a sequence of one index for each pass.
        day: The start of the day of the passes.
        passes: The rows of latitude, longitude, and seconds of each pass.
    """
    passes = np.asarray(passes, dtype=float)
    times = pd.Timestamp(day) + pd.to_timedelta(passes[:, 2], unit="s")
    targets = np.array(
        [hash_geometry(geo.Point(lon, lat)) for lat, lon, _ in passes], dtype=object
    )
    if instrument_index is None or np.ndim(instrument_index) == 0:
        groups = {instrument_index: np.arange(len(passes))}
    else:
        indices = np.asarray(instrument_index, dtype=int)
        groups = {int(i): np.flatnonzero(indices == i) for i in np.unique(indices)}
    differences = np.full(len(passes), np.nan)
    for index, rows in groups.items():
        found = collect_observations(
            [Point(latitude=passes[r, 0], longitude=passes[r, 1]) for r in rows],
            satellite,
            times[rows].min() - PASS_WINDOW,
            times[rows].max() + PASS_WINDOW,
            index,
        )
        by_target = {key: group for key, group in found.groupby("target_hash")}
        for r in rows:
            observed = by_target.get(targets[r])
            if observed is None:
                continue
            # the observations overlapping the window of this pass
            observed = observed[
                (observed.start <= times[r] + PASS_WINDOW)
                & (observed.end >= times[r] - PASS_WINDOW)
            ]
            if not observed.empty:
                seconds = (observed.epoch - times[r]).dt.total_seconds().to_numpy()
                differences[r] = seconds[np.argmin(np.abs(seconds))]
    return differences
