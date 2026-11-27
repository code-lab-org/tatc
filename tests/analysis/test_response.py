"""
Unit tests for the tatc.analysis.response module.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from tatc.analysis import compute_event_responses, compute_response_statistics
from tatc.utils.geometry import hash_geometry

START = datetime(2026, 1, 1, tzinfo=timezone.utc)
END = START + timedelta(hours=10)
PLACE = Point(10, 20)
OTHER = Point(0, 0)


def _observations(looks):
    """
    Observations of places at hours from the start, each with the hours of
    its downlink (`None` for none): {place: [(hour, downlink hour), ...]}.
    """
    rows = [
        (place, hour, downlink)
        for place, items in looks.items()
        for hour, downlink in items
    ]
    return gpd.GeoDataFrame(
        {
            "target_hash": [hash_geometry(place) for place, _, _ in rows],
            "epoch": [START + timedelta(hours=hour) for _, hour, _ in rows],
            "downlinked": pd.to_datetime(
                [
                    None if downlink is None else START + timedelta(hours=downlink)
                    for _, _, downlink in rows
                ],
                utc=True,
            ),
        },
        geometry=[place for place, _, _ in rows],
        crs="EPSG:4326",
    )


def _brute_force(looks, delivered):
    """
    The waits (hours, NaN if not answered) of onsets every second of the
    period (at the middle of each), answered by the next look (or the
    earliest delivery of the looks at or after the onset).
    """
    times = np.array([hour for hour, _ in looks]) * 3600
    answers = (
        np.array([np.inf if d is None else d for _, d in looks]) * 3600
        if delivered
        else times
    )
    onsets = np.arange(0, (END - START).total_seconds()) + 0.5
    waits = np.full(len(onsets), np.nan)
    for i, onset in enumerate(onsets):
        later = answers[times >= onset]
        if len(later) and np.isfinite(later.min()):
            waits[i] = (later.min() - onset) / 3600
    return waits


class TestResponseStatistics(unittest.TestCase):
    """
    Tests the statistics of responses to events at uniform onsets.
    """

    LOOKS = [(1.0, 3.0), (1.5, 2.0), (4.0, 7.0), (9.0, None)]

    def _check(self, delivered):
        waits = _brute_force(self.LOOKS, delivered)
        answered = ~np.isnan(waits)
        statistics = compute_response_statistics(
            _observations({PLACE: self.LOOKS}),
            START,
            END,
            quantiles=(0.25, 0.5, 0.9),
            deadline=timedelta(hours=1),
            delivered=delivered,
        ).iloc[0]
        hours = timedelta(hours=1)
        self.assertAlmostEqual(statistics.responded, answered.mean(), places=6)
        self.assertAlmostEqual(
            statistics.response / hours, waits[answered].mean(), places=3
        )
        for q in (0.25, 0.5, 0.9):
            self.assertAlmostEqual(
                statistics[f"response p{round(100 * q)}"] / hours,
                np.quantile(waits[answered], q),
                delta=1e-3,
            )
        self.assertAlmostEqual(
            statistics.max_response / hours, np.nanmax(waits), delta=1e-3
        )
        self.assertAlmostEqual(
            statistics.within_deadline,
            np.mean(answered & (np.nan_to_num(waits, nan=np.inf) <= 1)),
            places=4,
        )

    def test_first_looks(self):
        """
        Tests the waits for the next look against brute force.
        """
        self._check(delivered=False)

    def test_deliveries(self):
        """
        Tests the waits for the earliest delivery against brute force.
        """
        self._check(delivered=True)

    def test_mean_of_gaps(self):
        """
        Tests that the mean wait for a look is the sum of the squares of the
        gaps over twice the period, for looks that end the period.
        """
        statistics = compute_response_statistics(
            _observations({PLACE: [(2.0, None), (5.0, None), (10.0, None)]}),
            START,
            END,
        ).iloc[0]
        self.assertAlmostEqual(statistics.responded, 1)
        self.assertAlmostEqual(
            statistics.response / timedelta(hours=1), (4 + 9 + 25) / (2 * 10)
        )

    def test_targets_without_observations(self):
        """
        Tests that targets without observations answer no onsets.
        """
        statistics = compute_response_statistics(
            _observations({PLACE: self.LOOKS}),
            START,
            END,
            deadline=timedelta(hours=1),
            targets=gpd.GeoSeries([PLACE, OTHER], crs="EPSG:4326"),
        )
        self.assertEqual(len(statistics), 2)
        other = statistics.iloc[1]
        self.assertEqual(other.responded, 0)
        self.assertEqual(other.within_deadline, 0)
        self.assertTrue(pd.isna(other.response))
        self.assertNotIn(
            "within_deadline",
            compute_response_statistics(
                _observations({PLACE: self.LOOKS}), START, END
            ).columns,
        )

    def test_invalid(self):
        """
        Tests that empty periods and deliveries without downlinks raise
        errors.
        """
        observations = _observations({PLACE: self.LOOKS})
        with self.assertRaises(ValueError):
            compute_response_statistics(observations, END, START)
        with self.assertRaises(ValueError):
            compute_response_statistics(
                observations.drop(columns="downlinked"), START, END, delivered=True
            )


class TestEventResponses(unittest.TestCase):
    """
    Tests the responses to given events.
    """

    def test_responses(self):
        """
        Tests the first looks and deliveries of events, before and after
        their expiry, and of places not observed.
        """
        observations = _observations(
            {PLACE: [(1.0, 3.0), (1.5, 2.0), (4.0, 7.0), (9.0, None)]}
        )
        events = gpd.GeoDataFrame(
            {
                "onset": [START + timedelta(hours=h) for h in (1.2, 8.5, 0, 1.2)],
                "expires": pd.to_datetime(
                    [None, None, None, START + timedelta(hours=1.4)], utc=True
                ),
                "name": ["a", "b", "c", "d"],
            },
            geometry=[PLACE, PLACE, OTHER, PLACE],
            crs="EPSG:4326",
        )
        responses = compute_event_responses(observations, events)
        self.assertEqual(list(responses["name"]), ["a", "b", "c", "d"])
        self.assertEqual(list(responses["looks"]), [3, 1, 0, 0])
        self.assertEqual(
            list(responses["response"].iloc[:2]),
            [timedelta(minutes=18), timedelta(minutes=30)],
        )
        # the earliest delivery is of a later look than the first
        self.assertEqual(responses["delivery"].iloc[0], timedelta(minutes=48))
        self.assertTrue(pd.isna(responses["delivered"].iloc[1]))
        self.assertTrue(responses["first_look"].iloc[2:].isna().all())

    def test_without_downlinks(self):
        """
        Tests that observations without downlinks give no deliveries, and
        that events require onsets.
        """
        observations = _observations({PLACE: [(1.0, None)]}).drop(columns="downlinked")
        events = gpd.GeoDataFrame({"onset": [START]}, geometry=[PLACE], crs="EPSG:4326")
        responses = compute_event_responses(observations, events)
        self.assertNotIn("delivery", responses.columns)
        self.assertEqual(responses["response"].iloc[0], timedelta(hours=1))
        with self.assertRaises(ValueError):
            compute_event_responses(observations, events.drop(columns="onset"))
