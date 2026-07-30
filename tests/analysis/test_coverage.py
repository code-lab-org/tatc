"""
Unit tests for the coverage analysis functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
from datetime import datetime, timedelta, timezone

from tatc.analysis import (
    aggregate_observations,
    collect_multi_observations,
    collect_observations,
    reduce_observations,
)
from tatc.schemas import Point

from .common import IssConstellationTestCase


class TestCoverageAnalysis(IssConstellationTestCase):
    """
    Unit tests for the coverage analysis functions in tatc.analysis.
    """
    def setUp(self):
        super().setUp()
        self.point = Point(id=0, latitude=0, longitude=0)

    def test_collect_observations(self):
        """
        Test that observations can be collected for a single satellite and point.
        """
        collect_observations(
            self.point,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
            instrument_index=0,
            omit_solar=True,
        )

    def test_collect_observations_with_solar(self):
        """
        Test that observations can be collected for a single satellite and point with solar constraints.
        """
        collect_observations(
            self.point,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
            instrument_index=0,
            omit_solar=False,
        )

    def test_collect_observations_all_culminate(self):
        """
        Test that observations can be collected for a single satellite 
        and point when all observations culminate.
        """
        start = datetime(2022, 6, 1, 0, 43, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 45, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_all_culminate_no_instrument_index(self):
        """
        Test that observations can be collected for a single satellite and p
        oint when all observations culminate and no instrument index is provided.
        """
        start = datetime(2022, 6, 1, 0, 43, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 45, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_miss_first_rise(self):
        """
        Test that observations can be collected for a single satellite and point
        when the first observation is missed due to the start time.
        """
        start = datetime(2022, 6, 1, 0, 43, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 1, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)

    def test_collect_observations_miss_last_set(self):
        """
        Test that observations can be collected for a single satellite and point
        when the last observation is missed due to the end time.
        """
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 45, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_null(self):
        """
        Test that no observations are collected for a single satellite and point
        when the time window does not overlap with any observations.
        """
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertTrue(results.empty)

    def test_collect_observations_null_with_solar(self):
        """
        Test that no observations are collected for a single satellite and point
        when the time window does not overlap with any observations, even with solar constraints.
        """
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point, self.satellite, start, end, instrument_index=0, omit_solar=False
        )
        self.assertTrue(results.empty)

    def test_collect_multi_observations(self):
        """
        Test that observations can be collected for multiple satellites and a single point.
        """
        collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )

    def test_collect_multi_observations_null(self):
        """
        Test that no observations are collected for multiple satellites and a single point
        when the time window does not overlap with any observations.
        """
        start = datetime(2022, 6, 1, 0, 10, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc)
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            start,
            end,
        )
        self.assertTrue(results.empty)

    def test_aggregate_observations(self):
        """
        Test that observations can be aggregated for multiple satellites and a single point.
        """
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )
        results = aggregate_observations(results)
        for i in range(len(results.index)):
            self.assertEqual(
                results.iloc[i].end - results.iloc[i].start, results.iloc[i].access
            )
        for i in range(1, len(results.index)):
            self.assertEqual(
                results.iloc[i].start - results.iloc[i - 1].end, results.iloc[i].revisit
            )

    def test_aggregate_observations_null(self):
        """
        Test that no observations are aggregated for multiple satellites and a single point
        when the time window does not overlap with any observations.
        """
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, 0, 10, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
        )
        results = aggregate_observations(results)
        self.assertTrue(results.empty)

    def test_reduce_observations(self):
        """
        Test that observations can be reduced for multiple satellites and a single point.
        """
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 10, tzinfo=timezone.utc),
        )
        aggregated_results = aggregate_observations(results)
        reduced_results = reduce_observations(aggregated_results)
        self.assertEqual(len(reduced_results.index), 1)
        self.assertAlmostEqual(
            reduced_results.iloc[0].access,
            aggregated_results[
                aggregated_results.point_id == reduced_results.iloc[0].point_id
            ].access.mean(),
            delta=timedelta(seconds=0.01),
        )
        self.assertAlmostEqual(
            reduced_results.iloc[0].revisit,
            aggregated_results[
                aggregated_results.point_id == reduced_results.iloc[0].point_id
            ].revisit.mean(),
            delta=timedelta(seconds=0.01),
        )

    def test_reduce_observations_null(self):
        """
        Test that no observations are reduced for multiple satellites and a single point
        when the time window does not overlap with any observations.
        """
        results = collect_multi_observations(
            self.point,
            self.constellation.generate_members(),
            datetime(2022, 6, 1, 0, 10, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
        )
        aggregated_results = aggregate_observations(results)
        reduced_results = reduce_observations(aggregated_results)
        self.assertTrue(reduced_results.empty)
