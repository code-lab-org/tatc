"""
Unit tests for the sample size functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point as ShapelyPoint

from tatc.analysis import (
    collect_feature_observations,
    collect_observations,
    compute_sample_sizes,
)
from tatc.schemas import Instrument, PointedInstrument, Satellite, SunSynchronousOrbit
from tatc.utils import VelocityFrame, compute_along_track_field_of_view


class TestSampleSizes(unittest.TestCase):
    """
    Unit tests for `compute_sample_sizes`.
    """

    def setUp(self):
        self.start = datetime(2026, 3, 1, tzinfo=timezone.utc)
        self.end = self.start + timedelta(days=1)
        # a cross-track scanner like ATMS: 2.2 deg beams at 1.1 deg spacing
        self.scanner = PointedInstrument(
            name="scanner",
            field_of_regard=105.6,
            cross_track_field_of_view=105.6,
            along_track_field_of_view=2.2,
            cross_track_pixels=96,
            cross_track_oversampling=0.5,
            is_rectangular=True,
            velocity_frame=VelocityFrame.INERTIAL,
        )
        self.simple = Instrument(name="simple", field_of_regard=100)
        self.satellite = Satellite(
            name="sat",
            orbit=SunSynchronousOrbit(
                altitude=830e3, equator_crossing_time="13:30", epoch=self.start
            ),
            instruments=[self.scanner, self.simple],
        )
        self.points = [ShapelyPoint(-100, 35), ShapelyPoint(-60, -10)]

    def test_point_observations(self):
        """
        Test that the pixels of a pointed instrument are larger across track
        than along track off nadir, that the spacing of overlapping pixels is
        half their size across track, and that the scan angle lies within
        the field of regard.
        """
        observations = collect_observations(
            self.points, self.satellite, self.start, self.end, instrument_index=0
        )
        sizes = compute_sample_sizes(observations, self.satellite)
        self.assertGreater(len(sizes), 0)
        self.assertEqual(list(sizes.index), list(observations.index))
        self.assertTrue((sizes.sample_along > 0).all())
        self.assertTrue((sizes.sample_cross >= sizes.sample_along).all())
        np.testing.assert_allclose(
            sizes.spacing_cross / sizes.sample_cross, 0.5, rtol=0.02
        )
        np.testing.assert_allclose(sizes.spacing_along, sizes.sample_along)
        self.assertTrue((sizes.scan_angle.abs() <= 105.6 / 2 + 0.1).all())
        # the slant range is at least the satellite's altitude
        self.assertTrue((sizes.sat_range > 800e3).all())

    def test_instrument_without_pixels(self):
        """
        Test that an instrument without pixels has undefined sample sizes, but
        a slant range and scan angle.
        """
        observations = collect_observations(
            self.points, self.satellite, self.start, self.end, instrument_index=1
        )
        sizes = compute_sample_sizes(observations, self.satellite)
        self.assertGreater(len(sizes), 0)
        self.assertTrue(sizes.sample_along.isna().all())
        self.assertTrue(sizes.spacing_cross.isna().all())
        self.assertTrue(sizes.sat_range.notna().all())
        self.assertTrue(sizes.scan_angle.notna().all())

    def test_feature_observations(self):
        """
        Test that the samples of feature observations are computed at the
        centroids of the viewed parts, at their times.
        """
        satellite = self.satellite.model_copy(
            update={
                "instruments": [
                    self.scanner.model_copy(
                        update={
                            "along_track_field_of_view": compute_along_track_field_of_view(
                                830e3, 30, 98.7
                            )
                        }
                    )
                ]
            }
        )
        times = [self.start + timedelta(seconds=30 * (i + 0.5)) for i in range(2880)]
        features = gpd.GeoDataFrame(
            {
                "start": pd.DatetimeIndex([self.start], dtype="datetime64[ns, UTC]"),
                "end": pd.DatetimeIndex([self.end], dtype="datetime64[ns, UTC]"),
            },
            geometry=[ShapelyPoint(-100, 35).buffer(0.5)],
            crs="EPSG:4326",
        )
        observations = collect_feature_observations(features, satellite, times)
        sizes = compute_sample_sizes(observations, satellite)
        self.assertGreater(len(sizes), 0)
        np.testing.assert_allclose(sizes.sat_range, observations.sat_range, rtol=1e-9)
        self.assertTrue((sizes.sample_cross > 0).all())

    def test_unknown_satellite(self):
        """
        Test that observations by a satellite not given raise an error.
        """
        observations = collect_observations(
            self.points, self.satellite, self.start, self.end
        )
        other = self.satellite.model_copy(update={"name": "other"})
        with self.assertRaises(ValueError):
            compute_sample_sizes(observations, other)

    def test_empty_observations(self):
        """
        Test that no observations give no sample sizes, with their columns.
        """
        observations = collect_observations(
            self.points[0],
            self.satellite,
            self.start,
            self.start + timedelta(minutes=1),
        )
        sizes = compute_sample_sizes(observations, self.satellite)
        self.assertEqual(len(sizes), 0)
        self.assertIn("sample_cross", sizes.columns)
