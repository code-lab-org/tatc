"""
Unit tests for the coverage analysis functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import box

from tatc.analysis import (
    aggregate_observations,
    collect_multi_observations,
    collect_observations,
    grid_observations,
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

    def test_collect_observations_narrow_window_fully_visible(self):
        """
        Regression test: a narrow window entirely contained within a longer
        visible pass, with no rise, set, or culmination event inside it
        (elevation angle stays continuously above the threshold and never
        reaches a local maximum in-window), must still be reported as one
        continuous observation spanning the whole window -- not as no
        observation at all.
        """
        start = datetime(2022, 6, 1, 0, 43, 0, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 43, 30, tzinfo=timezone.utc)
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

    def test_collect_observations_narrow_window_fully_invisible(self):
        """
        Test that a narrow window entirely outside any visible pass (also
        producing no rise/set/culmination events) correctly yields no
        observations, distinguishing this from the fully-visible case in
        `test_collect_observations_narrow_window_fully_visible` (both
        produce zero events, but only one is a true miss).
        """
        start = datetime(2022, 6, 1, 0, 30, 0, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            start,
            end,
            instrument_index=0,
        )
        self.assertTrue(results.empty)

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

    def test_collect_multi_observations_no_satellites(self):
        """
        Regression test: an empty `satellites` list must return an empty
        DataFrame with the expected columns, not raise (there is nothing to
        concatenate when no satellite/instrument pair ever runs).
        """
        results = collect_multi_observations(
            self.point,
            [],
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )
        self.assertTrue(results.empty)
        self.assertIn("start", results.columns)

    @staticmethod
    def _make_observation(point_id, satellite, instrument, start, end):
        """
        Build a synthetic observation record matching the schema
        `collect_observations` produces, for direct, deterministic control
        over `aggregate_observations` inputs (bypassing orbital propagation).
        """
        return {
            "point_id": point_id,
            "geometry": ShapelyPoint(0, 0),
            "satellite": satellite,
            "instrument": instrument,
            "start": pd.Timestamp(start),
            "epoch": pd.Timestamp(start) + (pd.Timestamp(end) - pd.Timestamp(start)) / 2,
            "end": pd.Timestamp(end),
        }

    def test_aggregate_observations_merges_overlapping_and_nested_intervals(self):
        """
        Test the core interval-merging algorithm directly with synthetic,
        deliberately overlapping/nested/gapped observations: satellite B's
        window is fully nested inside A's, and C's window starts before A
        ends but after B ends (a case a naive "compare only to the previous
        row" check would wrongly split, since C.start > B.end even though
        C still overlaps A's still-ongoing window -- the running max via
        `.cummax()` is what keeps this correct). D is fully separate.
        """
        t0 = datetime(2022, 6, 1, tzinfo=timezone.utc)
        observations = gpd.GeoDataFrame(
            [
                self._make_observation(0, "A", "Test", t0, t0 + timedelta(minutes=10)),
                self._make_observation(
                    0, "B", "Test", t0 + timedelta(minutes=2), t0 + timedelta(minutes=4)
                ),
                self._make_observation(
                    0, "C", "Test", t0 + timedelta(minutes=9), t0 + timedelta(minutes=15)
                ),
                self._make_observation(
                    0, "D", "Test", t0 + timedelta(minutes=20), t0 + timedelta(minutes=25)
                ),
            ],
            crs="EPSG:4326",
        )
        results = aggregate_observations(observations)
        self.assertEqual(len(results.index), 2)
        self.assertEqual(results.iloc[0].satellite, "A, B, C")
        self.assertEqual(results.iloc[0].start, t0)
        self.assertEqual(results.iloc[0].end, t0 + timedelta(minutes=15))
        self.assertTrue(pd.isna(results.iloc[0].revisit))
        self.assertEqual(results.iloc[1].satellite, "D")
        self.assertEqual(results.iloc[1].start, t0 + timedelta(minutes=20))
        self.assertEqual(results.iloc[1].end, t0 + timedelta(minutes=25))
        self.assertEqual(results.iloc[1].revisit, timedelta(minutes=5))

    def test_aggregate_observations_isolates_point_ids(self):
        """
        Test that merging and revisit computation are scoped per point_id:
        a point_id=1 observation must not be merged with, or treated as a
        revisit predecessor for, a point_id=0 observation, even if their
        windows would otherwise overlap/abut.
        """
        t0 = datetime(2022, 6, 1, tzinfo=timezone.utc)
        observations = gpd.GeoDataFrame(
            [
                self._make_observation(0, "A", "Test", t0, t0 + timedelta(minutes=10)),
                self._make_observation(
                    1, "B", "Test", t0 + timedelta(minutes=5), t0 + timedelta(minutes=15)
                ),
            ],
            crs="EPSG:4326",
        )
        results = aggregate_observations(observations)
        self.assertEqual(len(results.index), 2)
        for i in range(len(results.index)):
            self.assertTrue(pd.isna(results.iloc[i].revisit))

    def test_aggregate_observations_epoch_is_merged_midpoint(self):
        """
        Test that the merged group's epoch is the midpoint of its (merged)
        start/end, not the mean of the constituent observations' own
        (pre-merge) epochs -- these differ here since B's epoch sits much
        earlier than the midpoint of the full A+B merged window.
        """
        t0 = datetime(2022, 6, 1, tzinfo=timezone.utc)
        observations = gpd.GeoDataFrame(
            [
                self._make_observation(0, "A", "Test", t0, t0 + timedelta(minutes=10)),
                self._make_observation(
                    0, "B", "Test", t0 + timedelta(minutes=1), t0 + timedelta(minutes=2)
                ),
            ],
            crs="EPSG:4326",
        )
        results = aggregate_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].epoch, t0 + timedelta(minutes=5))

    def test_aggregate_observations_drops_per_observation_columns(self):
        """
        Test that per-observation columns that lose their meaning once
        merged across satellites (sat_alt, sat_az, sat_sunlit, solar_alt,
        solar_az, solar_time) are dropped, even if present on the input.
        """
        t0 = datetime(2022, 6, 1, tzinfo=timezone.utc)
        record = self._make_observation(0, "A", "Test", t0, t0 + timedelta(minutes=10))
        record["sat_alt"] = 45.0
        record["sat_az"] = 180.0
        record["solar_alt"] = 10.0
        observations = gpd.GeoDataFrame([record], crs="EPSG:4326")
        results = aggregate_observations(observations)
        for column in ("sat_alt", "sat_az", "solar_alt"):
            self.assertNotIn(column, results.columns)

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

    @staticmethod
    def _make_aggregated_observation(point_id, access_minutes, revisit_minutes):
        """
        Build a synthetic aggregated-observation record (matching
        `aggregate_observations`'s output schema) for direct, deterministic
        control over `reduce_observations` inputs. `revisit_minutes=None`
        produces `pandas.NaT`, matching the first observation for a point.
        """
        return {
            "point_id": point_id,
            "geometry": ShapelyPoint(0, 0),
            "access": pd.Timedelta(minutes=access_minutes),
            "revisit": (
                pd.NaT if revisit_minutes is None else pd.Timedelta(minutes=revisit_minutes)
            ),
        }

    def test_reduce_observations_computes_mean_access_and_skips_first_revisit(self):
        """
        Test that access is averaged over every sample, but revisit is
        averaged only over the samples with a defined revisit -- the first
        sample's revisit is undefined (NaT, no prior observation), and must
        be skipped rather than treated as a zero-minute revisit (which
        would otherwise skew the mean down substantially).
        """
        observations = gpd.GeoDataFrame(
            [
                self._make_aggregated_observation(0, 5, None),
                self._make_aggregated_observation(0, 10, 55),
                self._make_aggregated_observation(0, 15, 50),
            ],
            crs="EPSG:4326",
        )
        results = reduce_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].samples, 3)
        self.assertEqual(results.iloc[0].access, timedelta(minutes=10))
        self.assertEqual(results.iloc[0].revisit, timedelta(minutes=52.5))

    def test_reduce_observations_isolates_point_ids(self):
        """
        Test that statistics are computed independently per point_id, not
        pooled across points.
        """
        observations = gpd.GeoDataFrame(
            [
                self._make_aggregated_observation(0, 5, None),
                self._make_aggregated_observation(1, 5, None),
                self._make_aggregated_observation(1, 7, 30),
            ],
            crs="EPSG:4326",
        )
        results = reduce_observations(observations)
        self.assertEqual(len(results.index), 2)
        point_0 = results[results.point_id == 0].iloc[0]
        point_1 = results[results.point_id == 1].iloc[0]
        self.assertEqual(point_0.samples, 1)
        self.assertEqual(point_0.access, timedelta(minutes=5))
        self.assertTrue(pd.isna(point_0.revisit))
        self.assertEqual(point_1.samples, 2)
        self.assertEqual(point_1.access, timedelta(minutes=6))
        self.assertEqual(point_1.revisit, timedelta(minutes=30))

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

    @staticmethod
    def _make_cell(cell_id, min_lon, min_lat, max_lon, max_lat):
        """
        Build a synthetic cell record matching the schema expected by
        `grid_observations` (a `cell_id` plus a polygon `geometry`).
        """
        return {"cell_id": cell_id, "geometry": box(min_lon, min_lat, max_lon, max_lat)}

    @staticmethod
    def _make_reduced_observation(point_id, lon, lat, access_seconds, revisit_seconds, samples):
        """
        Build a synthetic reduced-observation record (matching
        `reduce_observations`'s output schema) for direct, deterministic
        control over `grid_observations` inputs.
        """
        return {
            "point_id": point_id,
            "geometry": ShapelyPoint(lon, lat),
            "access": pd.Timedelta(seconds=access_seconds),
            "revisit": pd.Timedelta(seconds=revisit_seconds),
            "samples": samples,
        }

    def test_grid_observations_empty_reduced_observations(self):
        """
        Test that an empty `reduced_observations` yields every cell with
        zero samples and undefined access/revisit, rather than an empty
        result.
        """
        cells = gpd.GeoDataFrame(
            [self._make_cell(0, 0, 0, 1, 1), self._make_cell(1, 2, 2, 3, 3)],
            crs="EPSG:4326",
        )
        reduced = gpd.GeoDataFrame(
            columns=["point_id", "geometry", "access", "revisit", "samples"],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        self.assertEqual(len(result.index), 2)
        self.assertTrue((result.samples == 0).all())
        self.assertTrue(result.access.isna().all())
        self.assertTrue(result.revisit.isna().all())

    def test_grid_observations_single_point_passes_through_unchanged(self):
        """
        Test that a single point within a cell yields that point's own
        access/revisit/samples unchanged (the trivial case of a weighted
        mean over one value).
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        reduced = gpd.GeoDataFrame(
            [self._make_reduced_observation(0, 0.5, 0.5, 5, 10, 100)],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].samples, 100)
        self.assertEqual(result.iloc[0].access, timedelta(seconds=5))
        self.assertEqual(result.iloc[0].revisit, timedelta(seconds=10))

    def test_grid_observations_uses_weighted_arithmetic_mean_for_access(self):
        """
        Test that access is combined across points in the same cell using
        a sample-weighted arithmetic mean.
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        reduced = gpd.GeoDataFrame(
            [
                self._make_reduced_observation(0, 0.4, 0.4, 5, 10, 100),
                self._make_reduced_observation(1, 0.6, 0.6, 50, 100, 10),
            ],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].samples, 110)
        expected_access = (5 * 100 + 50 * 10) / 110
        self.assertAlmostEqual(
            result.iloc[0].access.total_seconds(), expected_access, places=6
        )

    def test_grid_observations_uses_reciprocal_of_summed_rate_for_revisit(self):
        """
        Test that revisit is combined across points in the same cell as the
        reciprocal of the summed per-point rates (1/revisit), not an
        arithmetic mean, nor the standard sample-weighted harmonic mean:
        revisit is a time-between-events (reciprocal-of-rate) quantity, and
        this specific combination is the one that keeps the cell's revisit
        exactly consistent with its summed sample count (100+10=110) under
        a shared mission duration -- e.g. duration/samples ~= revisit holds
        for point 0 (1000/100=10) and point 1 (1000/10=100) alike, and only
        this formula preserves that relationship at the cell level
        (1000/110 ~= 9.09s). The three candidate formulas all give clearly
        different results here (~9.09s reciprocal-of-summed-rate vs. ~10.89s
        standard weighted harmonic mean vs. 55s arithmetic mean), so this
        distinguishes all three concretely.
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        reduced = gpd.GeoDataFrame(
            [
                self._make_reduced_observation(0, 0.4, 0.4, 5, 10, 100),
                self._make_reduced_observation(1, 0.6, 0.6, 50, 100, 10),
            ],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        expected_revisit = 1 / (1 / 10 + 1 / 100)
        weighted_harmonic_mean_revisit = (100 + 10) / (100 / 10 + 10 / 100)
        naive_arithmetic_revisit = (10 * 100 + 100 * 10) / 110
        self.assertAlmostEqual(
            result.iloc[0].revisit.total_seconds(), expected_revisit, places=6
        )
        self.assertNotAlmostEqual(
            result.iloc[0].revisit.total_seconds(),
            weighted_harmonic_mean_revisit,
            places=1,
        )
        self.assertNotAlmostEqual(
            result.iloc[0].revisit.total_seconds(), naive_arithmetic_revisit, places=1
        )

    def test_grid_observations_cell_without_points_is_omitted(self):
        """
        Documents current behavior: unlike the fully-empty
        `reduced_observations` case (which zero-fills every cell), a cell
        with no points inside it is simply absent from a non-empty result,
        since the point-in-cell join is an inner join.
        """
        cells = gpd.GeoDataFrame(
            [self._make_cell(0, 0, 0, 1, 1), self._make_cell(1, 2, 2, 3, 3)],
            crs="EPSG:4326",
        )
        reduced = gpd.GeoDataFrame(
            [self._make_reduced_observation(0, 0.5, 0.5, 5, 10, 100)],
            crs="EPSG:4326",
        )
        result = grid_observations(reduced, cells)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].cell_id, 0)
