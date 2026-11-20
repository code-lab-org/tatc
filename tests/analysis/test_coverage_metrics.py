"""
Unit tests for counting observations over time.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import pandas as pd

from tatc.analysis import count_observations

START = datetime(2026, 1, 1, tzinfo=timezone.utc)


def _minutes(value):
    return pd.Timestamp(START + timedelta(minutes=value))


class TestCountObservations(unittest.TestCase):
    """
    Unit tests for `count_observations`.
    """

    def setUp(self):
        self.observations = pd.DataFrame(
            {
                "start": [_minutes(1), _minutes(2.5), _minutes(-5), _minutes(9)],
                "end": [_minutes(3), _minutes(2.5), _minutes(0.5), _minutes(20)],
                "kind": ["a", "b", "a", "b"],
                "plane": [0, 0, 1, 1],
            }
        )

    def test_counts(self):
        """
        Test that each step counts the observations that overlap it
        (inclusive of their ends), within the period, including those that
        start before it or end after it, and the last partial step.
        """
        counts = count_observations(
            self.observations,
            START,
            START + timedelta(minutes=9.5),
            timedelta(minutes=1),
        )
        self.assertListEqual(list(counts.columns), ["observations"])
        self.assertEqual(counts.index[0], pd.Timestamp(START))
        self.assertEqual(len(counts), 10)
        self.assertListEqual(list(counts.observations), [1, 1, 2, 1, 0, 0, 0, 0, 0, 1])

    def test_groups(self):
        """
        Test that the counts of each group (of one or more columns) are
        separate columns that sum to the counts.
        """
        end = START + timedelta(minutes=10)
        step = timedelta(minutes=1)
        total = count_observations(self.observations, START, end, step)
        by_kind = count_observations(self.observations, START, end, step, by="kind")
        self.assertListEqual(list(by_kind.columns), ["a", "b"])
        pd.testing.assert_series_equal(
            by_kind.sum(axis=1), total.observations, check_names=False
        )
        by_both = count_observations(
            self.observations, START, end, step, by=["kind", "plane"]
        )
        self.assertIn(("a", 1), by_both.columns)
        self.assertEqual(by_both[("a", 1)].iloc[0], 1)

    def test_invalid(self):
        """
        Test that a step that is not positive or a reversed period raises a
        `ValueError`, and that no observations give counts of zero.
        """
        with self.assertRaises(ValueError):
            count_observations(
                self.observations, START, START + timedelta(hours=1), timedelta(0)
            )
        with self.assertRaises(ValueError):
            count_observations(
                self.observations, START, START - timedelta(hours=1), timedelta(1)
            )
        empty = count_observations(
            self.observations.iloc[:0],
            START,
            START + timedelta(minutes=2),
            timedelta(minutes=1),
        )
        self.assertListEqual(list(empty.observations), [0, 0])


if __name__ == "__main__":
    unittest.main()
