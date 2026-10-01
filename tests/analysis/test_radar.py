"""
Unit tests for radar analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from shapely.geometry import box

from tatc.analysis import collect_radar_track, compute_radar_track
from tatc.schemas import RadarStation


class TestRadarAnalysis(unittest.TestCase):
    """
    Unit tests for radar analysis functions.
    """

    def setUp(self):
        self.station = RadarStation(name="Station 1", latitude=0, longitude=0)
        self.stations = [
            self.station,
            RadarStation(name="Station 2", latitude=0, longitude=1.5),
        ]

    def test_collect_radar_track_empty(self):
        """
        Test that an empty station list produces an empty data frame.
        """
        track = collect_radar_track([])
        self.assertTrue(track.empty)
        self.assertListEqual(
            list(track.columns),
            [
                "station",
                "elevation",
                "inner_ground_range",
                "outer_ground_range",
                "geometry",
            ],
        )

    def test_collect_radar_track_singular_matches_list(self):
        """
        Test that passing a single RadarStation produces the same result as
        passing a list containing that one station.
        """
        singular = collect_radar_track(self.station)
        as_list = collect_radar_track([self.station])
        self.assertEqual(len(singular), 1)
        self.assertTrue(singular.geometry.iloc[0].equals(as_list.geometry.iloc[0]))

    def test_collect_radar_track_multiple_stations(self):
        """
        Test that collecting multiple stations produces one row per station
        with ground ranges matching compute_ground_ranges directly.
        """
        track = collect_radar_track(self.stations, elevation=3048)
        self.assertEqual(len(track), 2)
        self.assertListEqual(list(track.station), ["Station 1", "Station 2"])
        for i, station in enumerate(self.stations):
            inner, outer = station.compute_ground_ranges(3048)
            self.assertAlmostEqual(track.inner_ground_range.iloc[i], inner, delta=1e-6)
            self.assertAlmostEqual(track.outer_ground_range.iloc[i], outer, delta=1e-6)

    def test_collect_radar_track_mask(self):
        """
        Test that a mask clips results to the intersecting station(s) only.
        """
        mask = box(-1, -1, 1, 1)
        track = collect_radar_track(self.stations, mask=mask)
        # only "Station 1" at (0, 0) falls within the mask's footprint overlap
        self.assertTrue((track.station == "Station 1").any())

    def test_compute_radar_track_dissolves_overlap(self):
        """
        Test that compute_radar_track merges overlapping footprints into a
        single geometry smaller than the sum of the individual footprints.
        """
        local_crs = "+proj=eqc +lat_ts=0 +datum=WGS84"
        individual = collect_radar_track(self.stations)
        total_area = individual.to_crs(local_crs).geometry.area.sum()
        merged = compute_radar_track(self.stations)
        self.assertEqual(len(merged), 1)
        self.assertListEqual(list(merged.columns), ["geometry"])
        self.assertLess(merged.to_crs(local_crs).geometry.iloc[0].area, total_area)

    def test_compute_radar_track_empty(self):
        """
        Test that an empty station list produces an empty data frame.
        """
        merged = compute_radar_track([])
        self.assertTrue(merged.empty)
