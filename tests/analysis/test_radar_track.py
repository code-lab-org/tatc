"""
Unit tests for radar analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import pandas as pd

from shapely.geometry import box

from tatc.analysis import collect_radar_track, compute_radar_track
from tatc.schemas import RadarBand, RadarStation


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
        track = collect_radar_track([], elevation=3048)
        self.assertTrue(track.empty)
        self.assertListEqual(
            list(track.columns),
            [
                "radar_station",
                "band",
                "elevation",
                "inner_ground_range",
                "outer_ground_range",
                "valid_obs",
                "geometry",
            ],
        )

    def test_collect_radar_track_singular_matches_list(self):
        """
        Test that passing a single RadarStation produces the same result as
        passing a list containing that one station.
        """
        singular = collect_radar_track(self.station, elevation=3048)
        as_list = collect_radar_track([self.station], elevation=3048)
        self.assertEqual(len(singular), 1)
        self.assertTrue(singular.geometry.iloc[0].equals(as_list.geometry.iloc[0]))

    def test_collect_radar_track_multiple_stations(self):
        """
        Test that collecting multiple stations produces one row per station
        with ground ranges matching compute_ground_ranges directly.
        """
        track = collect_radar_track(self.stations, elevation=3048)
        self.assertEqual(len(track), 2)
        self.assertListEqual(list(track.radar_station), ["Station 1", "Station 2"])
        for i, station in enumerate(self.stations):
            inner, outer = station.compute_ground_ranges(3048)
            self.assertAlmostEqual(track.inner_ground_range.iloc[i], inner, delta=1e-6)
            self.assertAlmostEqual(track.outer_ground_range.iloc[i], outer, delta=1e-6)

    def test_collect_radar_track_reports_band(self):
        """
        Test that the band column reports each station's tagged band (or
        missing if untagged).
        """
        tagged = RadarStation.from_band(
            RadarBand.X, name="X1", latitude=10, longitude=10
        )
        track = collect_radar_track([self.station, tagged], elevation=3048)
        self.assertTrue(pd.isna(track.band.iloc[0]))
        self.assertEqual(track.band.iloc[1], RadarBand.X)

    def test_collect_radar_track_mask(self):
        """
        Test that a mask clips results to the intersecting station(s) only.
        """
        mask = box(-1, -1, 1, 1)
        track = collect_radar_track(self.stations, elevation=3048, mask=mask)
        # only "Station 1" at (0, 0) falls within the mask's footprint overlap
        self.assertTrue((track.radar_station == "Station 1").any())

    def test_compute_radar_track_dissolves_overlap(self):
        """
        Test that compute_radar_track merges overlapping footprints into a
        single geometry smaller than the sum of the individual footprints.
        """
        local_crs = "+proj=eqc +lat_ts=0 +datum=WGS84"
        individual = collect_radar_track(self.stations, elevation=3048)
        total_area = individual.to_crs(local_crs).geometry.area.sum()
        merged = compute_radar_track(self.stations, elevation=3048)
        self.assertEqual(len(merged), 1)
        self.assertListEqual(list(merged.columns), ["geometry"])
        self.assertLess(merged.to_crs(local_crs).geometry.iloc[0].area, total_area)

    def test_compute_radar_track_empty(self):
        """
        Test that an empty station list produces an empty data frame.
        """
        merged = compute_radar_track([], elevation=3048)
        self.assertTrue(merged.empty)
        self.assertListEqual(list(merged.columns), ["geometry"])

    def test_collect_radar_track_valid_obs(self):
        """
        Test that a station whose beams do not reach the target elevation has
        an empty footprint, which is not a valid observation, and is not
        dissolved.
        """
        track = collect_radar_track(self.stations, elevation=0)
        self.assertFalse(track["valid_obs"].any())
        self.assertTrue(track.geometry.is_empty.all())
        self.assertTrue(compute_radar_track(self.stations, elevation=0).empty)
        self.assertTrue(collect_radar_track(self.stations, 3048)["valid_obs"].all())

    def test_collect_radar_track_mask_keeps_order(self):
        """
        Test that clipping to a mask keeps the order of the stations, and
        that a mask across the anti-meridian is split.
        """
        stations = [
            RadarStation(name="East", latitude=0, longitude=179.5),
            RadarStation(name="West", latitude=0, longitude=-179.5),
            RadarStation(name="Far", latitude=0, longitude=0),
        ]
        mask = box(178, -2, 182, 2)
        track = collect_radar_track(stations, elevation=3048, mask=mask)
        self.assertListEqual(list(track.radar_station), ["East", "West"])

    def test_collect_radar_track_type(self):
        """
        Test that each station must be a RadarStation.
        """
        with self.assertRaisesRegex(TypeError, r"stations\[1\] must be a RadarStation"):
            collect_radar_track([self.station, "KXXX"], elevation=3048)

    def test_collect_radar_track_requires_elevation(self):
        """
        Test that the target elevation must be specified (there is no
        meaningful default common to all stations).
        """
        with self.assertRaises(TypeError):
            collect_radar_track(self.stations)  # pylint: disable=no-value-for-parameter
        with self.assertRaises(TypeError):
            compute_radar_track(self.stations)  # pylint: disable=no-value-for-parameter

    def test_collect_radar_track_times(self):
        """
        Test that, with times, each station's static footprint is repeated
        at each time, sorted by time and then by station, as for surface
        platforms.
        """
        start = datetime(2024, 1, 1, tzinfo=timezone.utc)
        times = [start + timedelta(minutes=minutes) for minutes in (5, 0)]
        static = collect_radar_track(self.stations, elevation=3048)
        track = collect_radar_track(self.stations, elevation=3048, times=times)
        self.assertEqual(list(track.columns), ["time", *static.columns])
        self.assertListEqual(list(track.time), [times[1]] * 2 + [times[0]] * 2)
        self.assertListEqual(list(track.radar_station), ["Station 1", "Station 2"] * 2)
        for row in range(4):
            self.assertTrue(
                track.geometry.iloc[row].equals(static.geometry.iloc[row % 2])
            )
        empty = collect_radar_track(self.stations, elevation=3048, times=[])
        self.assertTrue(empty.empty)
        self.assertIn("time", empty.columns)
        with self.assertRaisesRegex(ValueError, "timezone-aware"):
            collect_radar_track(
                self.station, elevation=3048, times=[start.replace(tzinfo=None)]
            )
