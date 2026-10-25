"""
Unit tests for the feature observation functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import box
from skyfield.api import wgs84

from tatc.analysis import collect_feature_observations
from tatc.schemas import PointedInstrument, Satellite, SunSynchronousOrbit
from tatc.utils import compute_along_track_field_of_view
from tatc.utils.time import _to_time


class TestFeatureObservations(unittest.TestCase):
    """
    Unit tests for `collect_feature_observations`.
    """

    def setUp(self):
        self.start = datetime(2026, 3, 1, tzinfo=timezone.utc)
        orbit = SunSynchronousOrbit(
            altitude=700e3, equator_crossing_time="10:30", epoch=self.start
        )
        self.time_step = 30
        self.satellite = Satellite(
            name="sat",
            orbit=orbit,
            instruments=[
                PointedInstrument(
                    name="imager",
                    field_of_regard=40,
                    cross_track_field_of_view=40,
                    along_track_field_of_view=compute_along_track_field_of_view(
                        700e3, self.time_step, orbit.get_inclination()
                    ),
                    is_rectangular=True,
                )
            ],
        )
        self.times = [
            self.start + timedelta(seconds=self.time_step * (i + 0.5))
            for i in range(int(3 * 3600 / self.time_step))
        ]
        # the sub-satellite point at a sampled time about an hour after the start
        self.nadir_time = self.times[120]
        subpoint = wgs84.subpoint_of(
            orbit.to_gp_orbit().get_orbit_track_at_time(_to_time(self.nadir_time))
        )
        self.nadir = (subpoint.longitude.degrees, subpoint.latitude.degrees)

    def _features(self, geometries, start, end, index=None):
        return gpd.GeoDataFrame(
            {
                "start": pd.DatetimeIndex(
                    [start] * len(geometries), dtype="datetime64[ns, UTC]"
                ),
                "end": pd.DatetimeIndex(
                    [end] * len(geometries), dtype="datetime64[ns, UTC]"
                ),
            },
            geometry=geometries,
            crs="EPSG:4326",
            index=index,
        )

    def test_feature_at_nadir(self):
        """
        Test that a small feature beneath the satellite while it exists is
        viewed entirely, near nadir, at a time while it exists.
        """
        features = self._features(
            [ShapelyPoint(*self.nadir).buffer(0.2)],
            self.nadir_time - timedelta(minutes=5),
            self.nadir_time + timedelta(minutes=5),
        )
        observations = collect_feature_observations(
            features, self.satellite, self.times
        )
        self.assertEqual(len(observations), 1)
        observation = observations.iloc[0]
        self.assertEqual(observation.feature, 0)
        self.assertAlmostEqual(observation.viewed_fraction, 1)
        self.assertLess(
            abs(observation.time - pd.Timestamp(self.nadir_time)),
            pd.Timedelta(seconds=self.time_step),
        )
        self.assertGreater(observation.sat_alt, 85)
        # near nadir, the slant range is about the satellite's height
        self.assertLess(abs(observation.sat_range - observation.sat_height), 5e3)
        self.assertGreater(observation.viewed_area, 0)

    def test_feature_between_passes(self):
        """
        Test that a feature that exists only when no satellite passes over it
        is not observed.
        """
        features = self._features(
            [ShapelyPoint(*self.nadir).buffer(0.2)],
            self.nadir_time + timedelta(minutes=20),
            self.nadir_time + timedelta(minutes=30),
        )
        observations = collect_feature_observations(
            features, self.satellite, self.times
        )
        self.assertEqual(len(observations), 0)
        self.assertIn("viewed_fraction", observations.columns)

    def test_feature_wider_than_swath(self):
        """
        Test that a feature wider than the swath is viewed only in part, and
        that the viewed part is the feature's overlap with the swath.
        """
        longitude, latitude = self.nadir
        feature = box(longitude - 20, latitude - 1, longitude + 20, latitude + 1)
        features = self._features(
            [feature],
            self.nadir_time - timedelta(minutes=2),
            self.nadir_time + timedelta(minutes=2),
        )
        observations = collect_feature_observations(
            features, self.satellite, self.times
        )
        self.assertEqual(len(observations), 1)
        observation = observations.iloc[0]
        self.assertGreater(observation.viewed_fraction, 0)
        self.assertLess(observation.viewed_fraction, 0.5)
        self.assertTrue(feature.buffer(1e-9).contains(observation.geometry))

    def test_feature_index_labels(self):
        """
        Test that observations record each feature's index label.
        """
        features = self._features(
            [ShapelyPoint(*self.nadir).buffer(0.2), ShapelyPoint(0, -80).buffer(0.2)],
            self.nadir_time - timedelta(minutes=5),
            self.nadir_time + timedelta(minutes=5),
            index=[17, 42],
        )
        observations = collect_feature_observations(
            features, self.satellite, self.times
        )
        self.assertEqual(list(observations.feature), [17])

    def test_invalid_features(self):
        """
        Test that features without timezone-aware start and end times are
        rejected.
        """
        geometry = [ShapelyPoint(*self.nadir).buffer(0.2)]
        with self.assertRaises(ValueError):
            collect_feature_observations(
                gpd.GeoDataFrame(geometry=geometry, crs="EPSG:4326"),
                self.satellite,
                self.times,
            )
        naive = gpd.GeoDataFrame(
            {"start": [datetime(2026, 3, 1)], "end": [datetime(2026, 3, 2)]},
            geometry=geometry,
            crs="EPSG:4326",
        )
        with self.assertRaises(ValueError):
            collect_feature_observations(naive, self.satellite, self.times)
        with self.assertRaises(TypeError):
            collect_feature_observations(geometry, self.satellite, self.times)

    def test_empty_features(self):
        """
        Test that no features give an empty data frame of observations.
        """
        features = self._features([], self.start, self.start)
        observations = collect_feature_observations(
            features, self.satellite, self.times
        )
        self.assertEqual(len(observations), 0)
        self.assertTrue(np.isin(["sat_alt", "sat_range"], observations.columns).all())
