"""
Unit tests for latency analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timezone

from tatc.analysis import (
    collect_downlinks,
    collect_observations,
    compute_latencies,
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
