"""
Unit tests for the coincident observation analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point as ShapelyPoint, box

from tatc.analysis import (
    aggregate_observations,
    collect_ground_track,
    collect_observations,
    find_coincident_ground_track,
    find_coincident_observations,
)
from tatc.schemas import Instrument, Satellite

from .common import IssConstellationTestCase

T0 = datetime(2022, 6, 1, tzinfo=timezone.utc)


def _minutes(value):
    return T0 + timedelta(minutes=value)


def _make_observation(satellite, start, end, instrument="Test", target_hash="0"):
    """
    Build a synthetic observation record matching the schema
    `collect_observations` produces, with start and end in minutes after T0.
    """
    return {
        "target_hash": target_hash,
        "geometry": ShapelyPoint(0, 0),
        "satellite": satellite,
        "instrument": instrument,
        "start": pd.Timestamp(_minutes(start)),
        "epoch": pd.Timestamp(_minutes((start + end) / 2)),
        "end": pd.Timestamp(_minutes(end)),
    }


def _make_footprint(satellite, time, geometry, instrument="Test"):
    """
    Build a synthetic footprint record matching the schema
    `collect_ground_track` produces, with time in minutes after T0.
    """
    return {
        "time": pd.Timestamp(_minutes(time)),
        "satellite": satellite,
        "instrument": instrument,
        "valid_obs": True,
        "geometry": geometry,
    }


class TestFindCoincidentObservations(IssConstellationTestCase):
    """
    Unit tests for `find_coincident_observations`.
    """

    def test_overlap(self):
        """
        Test that the coincident period of two satellites is the overlap of
        their observation periods.
        """
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 6, 15)],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite, "A, B")
        self.assertEqual(results.iloc[0].instrument, "Test, Test")
        self.assertEqual(results.iloc[0].start, _minutes(6))
        self.assertEqual(results.iloc[0].end, _minutes(10))
        self.assertEqual(results.iloc[0].epoch, _minutes(8))

    def test_no_overlap(self):
        """
        Test that separate observation periods are not coincident, giving an
        empty data frame with the result columns.
        """
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 12, 15)],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations)
        self.assertTrue(results.empty)
        self.assertListEqual(
            list(results.columns),
            [
                "target_hash",
                "geometry",
                "satellite",
                "instrument",
                "start",
                "epoch",
                "end",
            ],
        )

    def test_max_time_difference(self):
        """
        Test that observations coincide if the gap between their periods is
        at most the maximum time difference, widening each period by half.
        """
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 12, 15)],
            crs="EPSG:4326",
        )
        self.assertTrue(
            find_coincident_observations(
                observations, max_time_difference=timedelta(minutes=1)
            ).empty
        )
        results = find_coincident_observations(
            observations, max_time_difference=timedelta(minutes=4)
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].start, _minutes(10))
        self.assertEqual(results.iloc[0].end, _minutes(12))

    def test_instantaneous(self):
        """
        Test that instantaneous observations coincide if their epochs differ
        by at most the maximum time difference, regardless of their periods,
        and that each coincident observation spans the contributing epochs.
        """
        # periods overlap (minutes 8 to 10), epochs (minutes 5 and 19) do not
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 8, 30)],
            crs="EPSG:4326",
        )
        self.assertEqual(len(find_coincident_observations(observations).index), 1)
        self.assertTrue(
            find_coincident_observations(
                observations,
                max_time_difference=timedelta(minutes=10),
                instantaneous=True,
            ).empty
        )
        results = find_coincident_observations(
            observations, max_time_difference=timedelta(minutes=15), instantaneous=True
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite, "A, B")
        self.assertEqual(results.iloc[0].start, _minutes(5))
        self.assertEqual(results.iloc[0].end, _minutes(19))
        self.assertEqual(results.iloc[0].epoch, _minutes(12))

    def test_same_satellite(self):
        """
        Test that observations by different instruments on the same
        satellite are not coincident with each other.
        """
        observations = gpd.GeoDataFrame(
            [
                _make_observation("A", 0, 10, instrument="X"),
                _make_observation("A", 5, 15, instrument="Y"),
            ],
            crs="EPSG:4326",
        )
        self.assertTrue(find_coincident_observations(observations).empty)

    def test_maximal_period(self):
        """
        Test that a coincident period spans the time that at least two
        satellites observe, even as the satellites observing change, and
        records each satellite and instrument once.
        """
        observations = gpd.GeoDataFrame(
            [
                _make_observation("A", 0, 10),
                _make_observation("B", 2, 4, instrument="X"),
                _make_observation("B", 2, 4, instrument="X"),
                _make_observation("C", 3, 12),
                _make_observation("D", 20, 25),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite, "A, B, C")
        self.assertEqual(results.iloc[0].instrument, "Test, X, Test")
        self.assertEqual(results.iloc[0].start, _minutes(2))
        self.assertEqual(results.iloc[0].end, _minutes(10))

    def test_min_platforms(self):
        """
        Test that a minimum of three satellites requires all three to
        observe at once.
        """
        observations = gpd.GeoDataFrame(
            [
                _make_observation("A", 0, 10),
                _make_observation("B", 2, 8),
                _make_observation("C", 6, 12),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations, min_platforms=3)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite, "A, B, C")
        self.assertEqual(results.iloc[0].start, _minutes(6))
        self.assertEqual(results.iloc[0].end, _minutes(8))

    def test_isolates_targets(self):
        """
        Test that observations of different targets are not coincident.
        """
        observations = gpd.GeoDataFrame(
            [
                _make_observation("A", 0, 10, target_hash="0"),
                _make_observation("B", 5, 15, target_hash="1"),
                _make_observation("C", 5, 15, target_hash="1"),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].target_hash, "1")
        self.assertEqual(results.iloc[0].satellite, "B, C")

    def test_invalid_arguments(self):
        """
        Test that invalid arguments and region observations raise errors.
        """
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 6, 15)],
            crs="EPSG:4326",
        )
        with self.assertRaises(ValueError):
            find_coincident_observations(observations, min_platforms=1)
        with self.assertRaises(ValueError):
            find_coincident_observations(
                observations, max_time_difference=timedelta(minutes=-1)
            )
        observations.geometry = [box(0, 0, 1, 1), box(0, 0, 1, 1)]
        with self.assertRaises(ValueError):
            find_coincident_observations(observations)

    def test_collect_observations(self):
        """
        Test that the coincident observations of a constellation lie within
        observations by at least two of its satellites, and aggregate.
        """
        satellites = self.constellation.generate_members()
        observations = collect_observations(
            ShapelyPoint(0, 45), satellites, T0, T0 + timedelta(days=1)
        )
        # successive passes of the constellation are about 40 minutes apart
        max_time_difference = timedelta(minutes=45)
        results = find_coincident_observations(
            observations, max_time_difference=max_time_difference
        )
        self.assertFalse(results.empty)
        half = max_time_difference / 2
        for result in results.itertuples():
            contributors = observations[
                (observations.start - half <= result.end)
                & (observations.end + half >= result.start)
            ]
            self.assertGreaterEqual(contributors.satellite.nunique(), 2)
            self.assertSetEqual(
                set(result.satellite.split(", ")), set(contributors.satellite)
            )
        self.assertEqual(len(aggregate_observations(results).index), len(results.index))


class TestFindCoincidentGroundTrack(IssConstellationTestCase):
    """
    Unit tests for `find_coincident_ground_track`.
    """

    def test_overlap(self):
        """
        Test that the overlapping footprints of two satellites at the same
        time are coincident, with their overlap as the geometry.
        """
        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("A", 0, box(0, 0, 2, 2)),
                _make_footprint("B", 0, box(1, 1, 3, 3)),
                _make_footprint("B", 1, box(1, 1, 3, 3)),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_ground_track(ground_track)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite_1, "A")
        self.assertEqual(results.iloc[0].satellite_2, "B")
        self.assertEqual(results.iloc[0].time_1, _minutes(0))
        self.assertEqual(results.iloc[0].time_2, _minutes(0))
        self.assertEqual(results.iloc[0].time_difference, timedelta(0))
        self.assertTrue(results.iloc[0].valid_obs_1 and results.iloc[0].valid_obs_2)
        self.assertTrue(results.iloc[0].geometry.equals(box(1, 1, 2, 2)))
        self.assertEqual(results.crs, ground_track.crs)

    def test_max_time_difference(self):
        """
        Test that footprints coincide if their times differ by at most the
        maximum time difference, ordered by time.
        """
        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("B", 0, box(1, 1, 3, 3)),
                _make_footprint("A", 5, box(0, 0, 2, 2)),
            ],
            crs="EPSG:4326",
        )
        self.assertTrue(find_coincident_ground_track(ground_track).empty)
        results = find_coincident_ground_track(
            ground_track, max_time_difference=timedelta(minutes=5)
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite_1, "B")
        self.assertEqual(results.iloc[0].satellite_2, "A")
        self.assertEqual(results.iloc[0].time_difference, timedelta(minutes=5))

    def test_same_satellite_and_touching(self):
        """
        Test that footprints of the same satellite, or that only touch,
        are not coincident.
        """
        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("A", 0, box(0, 0, 2, 2), instrument="X"),
                _make_footprint("A", 0, box(1, 1, 3, 3), instrument="Y"),
                _make_footprint("B", 0, box(3, 3, 4, 4)),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_ground_track(ground_track)
        self.assertTrue(results.empty)
        self.assertIn("time_difference", results.columns)
        self.assertIn("satellite_2", results.columns)

    def test_empty(self):
        """
        Test that an empty ground track has no coincident footprints.
        """
        ground_track = collect_ground_track([], [])
        results = find_coincident_ground_track(ground_track)
        self.assertTrue(results.empty)
        with self.assertRaises(ValueError):
            find_coincident_ground_track(
                ground_track, max_time_difference=timedelta(minutes=-1)
            )

    def test_collect_ground_track(self):
        """
        Test that the coincident footprints of a constellation lie within
        the footprints of both satellites.
        """
        instrument = Instrument(name="Test", field_of_regard=60.0)
        satellites = [
            Satellite(name=member.name, orbit=member.orbit, instruments=[instrument])
            for member in self.constellation.generate_members()
        ]
        times = [T0 + timedelta(minutes=i) for i in range(0, 24 * 60, 2)]
        ground_track = collect_ground_track(satellites, times)
        results = find_coincident_ground_track(
            ground_track, max_time_difference=timedelta(minutes=10)
        )
        self.assertFalse(results.empty)
        self.assertTrue((results.satellite_1 != results.satellite_2).all())
        self.assertTrue((results.time_difference <= timedelta(minutes=10)).all())
        self.assertTrue((results.time_1 <= results.time_2).all())
        footprints = ground_track.set_index(["time", "satellite"]).geometry
        for result in results.itertuples():
            for time, satellite in (
                (result.time_1, result.satellite_1),
                (result.time_2, result.satellite_2),
            ):
                footprint = footprints.loc[(time, satellite)]
                self.assertTrue(footprint.buffer(1e-9).contains(result.geometry))
