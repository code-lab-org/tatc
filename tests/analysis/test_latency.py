"""
Unit tests for latency analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import box

from tatc.analysis import (
    collect_downlinks,
    collect_observations,
    compute_latencies,
    grid_latencies,
    reduce_latencies,
)
from tatc.schemas import GroundStation, Point

from .common import IssConstellationTestCase


class TestLatencyAnalysis(IssConstellationTestCase):
    """
    Unit tests for latency analysis functions.
    """

    def setUp(self):
        super().setUp()
        self.point = Point(id=0, latitude=0, longitude=0, min_elevation_angle=10)
        self.station = GroundStation(
            name="Station 1", latitude=0, longitude=180, min_elevation_angle=10
        )
        self.stations = [
            self.station,
            GroundStation(
                name="Station 2", latitude=50, longitude=0, min_elevation_angle=10
            ),
            GroundStation(
                name="Station 3", latitude=50, longitude=90, min_elevation_angle=10
            ),
            GroundStation(
                name="Station 4", latitude=50, longitude=-90, min_elevation_angle=10
            ),
        ]

    def test_collect_downlinks(self):
        """
        Test that downlink collection works for a single station and satellite.
        """
        collect_downlinks(
            self.station,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 10, tzinfo=timezone.utc),
        )

    def test_collect_downlinks_empty(self):
        """
        Test that downlink collection returns an empty DataFrame when no downlinks are available.
        """
        results = collect_downlinks(
            self.station,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 1, tzinfo=timezone.utc),
        )
        self.assertTrue(results.empty)

    def test_collect_multi_downlinks(self):
        """
        Test that downlink collection works for multiple stations and a single satellite.
        """
        collect_downlinks(
            self.stations,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 10, tzinfo=timezone.utc),
        )

    def test_collect_multi_downlinks_empty(self):
        """
        Test that downlink collection returns an empty DataFrame when no downlinks are available.
        """
        results = collect_downlinks(
            self.stations,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 0, 10, tzinfo=timezone.utc),
        )
        self.assertTrue(results.empty)

    def test_compute_latency(self):
        """
        Test that latency computation works for a single point, satellite, and station.
        """
        compute_latencies(
            collect_observations(
                self.point,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
                instrument_index=0,
            ),
            collect_downlinks(
                self.station,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
        )

    def test_compute_latency_empty(self):
        """
        Test that latency computation returns an empty DataFrame when no latencies are available.
        """
        results = compute_latencies(
            collect_observations(
                self.point,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
                instrument_index=0,
            ),
            collect_downlinks(
                self.station,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
        )
        self.assertTrue(results.empty)

    def test_compute_latency_no_downlinks(self):
        """
        Test that latency computation returns an empty DataFrame when no downlinks are available.
        """
        results = compute_latencies(
            collect_observations(
                self.point,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
                instrument_index=0,
            ),
            collect_downlinks(
                self.station,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 1, 1, tzinfo=timezone.utc),
            ),
        )
        self.assertTrue(results.empty)

    def test_compute_latency_multi_station(self):
        """
        Test that latency computation works for a single point, satellite, and multiple stations.
        """
        results = compute_latencies(
            collect_observations(
                self.point,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
                instrument_index=0,
            ),
            collect_downlinks(
                self.stations,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
        )

    def test_reduce_latency(self):
        """
        Test that latency reduction works for a single point, satellite, and multiple stations.
        """
        reduce_latencies(
            compute_latencies(
                collect_observations(
                    self.point,
                    self.satellite,
                    datetime(2022, 6, 1, tzinfo=timezone.utc),
                    datetime(2022, 6, 10, tzinfo=timezone.utc),
                    instrument_index=0,
                ),
                collect_downlinks(
                    self.stations,
                    self.satellite,
                    datetime(2022, 6, 1, tzinfo=timezone.utc),
                    datetime(2022, 6, 10, tzinfo=timezone.utc),
                ),
            )
        )

    def test_reduce_latency_empty(self):
        """
        Test that latency reduction returns an empty DataFrame when no latencies are available.
        """
        results = reduce_latencies(
            compute_latencies(
                collect_observations(
                    self.point,
                    self.satellite,
                    datetime(2022, 6, 1, tzinfo=timezone.utc),
                    datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
                    instrument_index=0,
                ),
                collect_downlinks(
                    self.stations,
                    self.satellite,
                    datetime(2022, 6, 1, tzinfo=timezone.utc),
                    datetime(2022, 6, 10, tzinfo=timezone.utc),
                ),
            )
        )
        self.assertTrue(results.empty)

    def test_compute_latency_columns_present_and_correct(self):
        """
        Regression test: `sat_alt`/`sat_az`/`station` must survive
        `compute_latencies` with their real values -- these columns are
        never suffixed by the underlying `merge_asof` (each exists in only
        one of `observations`/`downlinks`), unlike `epoch`/`geometry`
        (which exist in both, and are suffixed and explicitly renamed).
        """
        observations = collect_observations(
            self.point,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 10, tzinfo=timezone.utc),
            instrument_index=0,
        )
        downlinks = collect_downlinks(
            self.station,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 10, tzinfo=timezone.utc),
        )
        results = compute_latencies(observations, downlinks)
        self.assertIn("sat_alt", results.columns)
        self.assertIn("sat_az", results.columns)
        self.assertIn("station", results.columns)
        for column in ("sat_alt", "sat_az"):
            self.assertTrue((results[column] == observations[column]).all())
        matched = results.station.notna()
        self.assertTrue(matched.any())
        self.assertTrue((results.station[matched] == self.station.name).all())

    def test_reduce_latencies_all_unmatched(self):
        """
        Test that observations with no matching downlink (NaT latency)
        still count toward `samples`, but are excluded from the mean
        latency (which stays undefined/NaT), rather than being treated as
        a zero-latency sample or crashing.
        """
        observations = gpd.GeoDataFrame(
            [
                {
                    "point_id": 0,
                    "geometry": ShapelyPoint(0, 0),
                    "satellite": "A",
                    "instrument": "I",
                    "start": pd.Timestamp("2022-06-01T10:00", tz="UTC"),
                    "epoch": pd.Timestamp("2022-06-01T10:01", tz="UTC"),
                    "end": pd.Timestamp("2022-06-01T10:02", tz="UTC"),
                    "sat_alt": 45.0,
                    "sat_az": 90.0,
                },
                {
                    "point_id": 0,
                    "geometry": ShapelyPoint(0, 0),
                    "satellite": "A",
                    "instrument": "I",
                    "start": pd.Timestamp("2022-06-01T20:00", tz="UTC"),
                    "epoch": pd.Timestamp("2022-06-01T20:01", tz="UTC"),
                    "end": pd.Timestamp("2022-06-01T20:02", tz="UTC"),
                    "sat_alt": 30.0,
                    "sat_az": 100.0,
                },
            ],
            crs="EPSG:4326",
        )
        downlinks = gpd.GeoDataFrame(
            [
                {
                    "station": "S1",
                    "geometry": ShapelyPoint(1, 1),
                    "satellite": "A",
                    "start": pd.Timestamp("2022-06-01T00:00", tz="UTC"),
                    "epoch": pd.Timestamp("2022-06-01T00:01", tz="UTC"),
                    "end": pd.Timestamp("2022-06-01T00:02", tz="UTC"),
                }
            ],
            crs="EPSG:4326",
        )
        latencies = compute_latencies(observations, downlinks)
        self.assertTrue(latencies.latency.isna().all())
        result = reduce_latencies(latencies)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].samples, 2)
        self.assertTrue(pd.isna(result.iloc[0].latency))

    @staticmethod
    def _make_cell(cell_id, min_lon, min_lat, max_lon, max_lat):
        """
        Build a synthetic cell record matching the schema expected by
        `grid_latencies` (a `cell_id` plus a polygon `geometry`).
        """
        return {"cell_id": cell_id, "geometry": box(min_lon, min_lat, max_lon, max_lat)}

    @staticmethod
    def _make_reduced_latency(point_id, lon, lat, latency_seconds, samples):
        """
        Build a synthetic reduced-latency record (matching
        `reduce_latencies`'s output schema) for direct, deterministic
        control over `grid_latencies` inputs.
        """
        return {
            "point_id": point_id,
            "geometry": ShapelyPoint(lon, lat),
            "latency": pd.Timedelta(seconds=latency_seconds),
            "samples": samples,
        }

    def test_grid_latencies_empty_reduced_latencies(self):
        """
        Test that an empty `reduced_latencies` yields every cell with zero
        samples and undefined latency, rather than an empty result.
        """
        cells = gpd.GeoDataFrame(
            [self._make_cell(0, 0, 0, 1, 1), self._make_cell(1, 2, 2, 3, 3)],
            crs="EPSG:4326",
        )
        reduced = gpd.GeoDataFrame(
            columns=["point_id", "geometry", "latency", "samples"], crs="EPSG:4326"
        )
        result = grid_latencies(reduced, cells)
        self.assertEqual(len(result.index), 2)
        self.assertTrue((result.samples == 0).all())
        self.assertTrue(result.latency.isna().all())

    def test_grid_latencies_single_point_passes_through_unchanged(self):
        """
        Test that a single point within a cell yields that point's own
        latency/samples unchanged (the trivial case of a weighted mean
        over one value).
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        reduced = gpd.GeoDataFrame(
            [self._make_reduced_latency(0, 0.5, 0.5, 5, 100)], crs="EPSG:4326"
        )
        result = grid_latencies(reduced, cells)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].samples, 100)
        self.assertEqual(result.iloc[0].latency, timedelta(seconds=5))

    def test_grid_latencies_uses_weighted_arithmetic_mean(self):
        """
        Test that latency is combined across points in the same cell using
        a sample-weighted arithmetic mean (latency is a per-observation
        duration, unlike revisit, so no harmonic-mean treatment applies).
        """
        cells = gpd.GeoDataFrame([self._make_cell(0, 0, 0, 1, 1)], crs="EPSG:4326")
        reduced = gpd.GeoDataFrame(
            [
                self._make_reduced_latency(0, 0.4, 0.4, 5, 100),
                self._make_reduced_latency(1, 0.6, 0.6, 50, 10),
            ],
            crs="EPSG:4326",
        )
        result = grid_latencies(reduced, cells)
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].samples, 110)
        expected_latency = (5 * 100 + 50 * 10) / 110
        self.assertAlmostEqual(
            result.iloc[0].latency.total_seconds(), expected_latency, places=6
        )
