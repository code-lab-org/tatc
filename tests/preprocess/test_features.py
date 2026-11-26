"""
Unit tests for finding the features of a nature run's fields.
"""

import unittest
from datetime import timedelta

import numpy as np
import pandas as pd
import pytest
from pyproj import Geod

from tatc import constants
from tatc.preprocess import find_features

xr = pytest.importorskip("xarray")
pytest.importorskip("scipy")

GEOD = Geod(ellps="WGS84")
RADIUS = constants.EARTH_MEAN_RADIUS


def _distance(lon, lat, lon0, lat0):
    """The great-circle distance (meters) from a center."""
    lon, lat, lon0, lat0 = (np.radians(v) for v in (lon, lat, lon0, lat0))
    cosine = np.sin(lat) * np.sin(lat0) + np.cos(lat) * np.cos(lat0) * np.cos(
        lon - lon0
    )
    return RADIUS * np.arccos(np.clip(cosine, -1, 1))


def _discs(lon, lat, discs, times=2):
    """A mask (time, lat, lon) of discs ((lon, lat, radius, frames)) on a grid."""
    lon2, lat2 = np.meshgrid(lon, lat)
    data = np.zeros((times, len(lat), len(lon)), dtype=bool)
    for lon0, lat0, radius, frames in discs:
        inside = _distance(lon2, lat2, lon0, lat0) <= radius
        for frame in frames:
            data[frame] |= inside
    return xr.DataArray(
        data,
        dims=("time", "lat", "lon"),
        coords={
            "time": pd.date_range("2026-01-01", periods=times, freq="30min"),
            "lat": lat,
            "lon": lon,
        },
    )


class TestFindFeatures(unittest.TestCase):
    """
    Unit tests for `find_features`.
    """

    def setUp(self):
        self.lon = np.arange(0.0, 20.0, 0.1)
        self.lat = np.arange(-10.0, 10.0, 0.1)
        self.mask = _discs(
            self.lon,
            self.lat,
            [(5.0, 0.0, 100e3, [0]), (14.0, 3.0, 300e3, [0, 1])],
        )

    def test_discs(self):
        """Discs are found with their areas, centroids, frames, and polygons."""
        features = find_features(self.mask)
        self.assertEqual(list(features.frame), [0, 0, 1])
        self.assertEqual(list(features.index), [0, 1, 2])
        expected = np.pi * np.array([100e3, 300e3, 300e3]) ** 2
        np.testing.assert_allclose(features.area_m2, expected, rtol=0.04)
        np.testing.assert_allclose(
            features.diameter_m, np.sqrt(4 * features.area_m2 / np.pi)
        )
        np.testing.assert_allclose(features.lon, [5.0, 14.0, 14.0], atol=0.02)
        np.testing.assert_allclose(features.lat, [0.0, 3.0, 3.0], atol=0.02)
        geodesic = [abs(GEOD.geometry_area_perimeter(g)[0]) for g in features.geometry]
        np.testing.assert_allclose(geodesic, features.area_m2, rtol=0.01)
        self.assertTrue(features.geometry.is_valid.all())
        # frames of 30 minutes, starting at their times
        self.assertEqual(features.start[2], pd.Timestamp("2026-01-01 00:30", tz="UTC"))
        self.assertEqual(features.end[2] - features.start[2], timedelta(minutes=30))
        self.assertIsInstance(features.start.dtype, pd.DatetimeTZDtype)

    def test_reductions(self):
        """Variables are reduced over each feature's cells."""
        lon2 = np.broadcast_to(self.lon, (2, len(self.lat), len(self.lon)))
        variable = self.mask.copy(data=lon2.astype(float))
        constant = self.mask.copy(data=np.full(self.mask.shape, 2.0))
        features = find_features(self.mask, variables={"x": variable, "c": constant})
        np.testing.assert_allclose(features.mean_c, 2.0)
        np.testing.assert_allclose(features.integral_c, 2 * features.area_m2)
        self.assertLess(features.min_x[0], 5.0)
        self.assertGreater(features.max_x[0], 5.0)
        np.testing.assert_allclose(features.mean_x[0], 5.0, atol=0.02)
        features = find_features(
            self.mask, variables={"c": constant}, reductions=["max"]
        )
        self.assertEqual([c for c in features.columns if c.endswith("_c")], ["max_c"])

    def test_min_diameter_and_labels(self):
        """Small features are left out and unlabeled."""
        features, labels = find_features(
            self.mask, min_diameter=300e3, return_labels=True
        )
        self.assertEqual(len(features), 2)
        self.assertEqual(labels.dims, self.mask.dims)
        for row in features.itertuples():
            cells = int((labels.isel(time=row.frame) == row.label).sum())
            self.assertEqual(cells, row.n_cells)
        # the small disc's cells are unlabeled
        self.assertEqual(
            int(labels.isel(time=0).sel(lon=5.0, lat=0.0, method="nearest")), 0
        )

    def test_descending_latitudes(self):
        """Grids with descending latitudes give the same features and labels."""
        flipped = self.mask.isel(lat=slice(None, None, -1))
        features, labels = find_features(flipped, return_labels=True)
        expected = find_features(self.mask)
        np.testing.assert_allclose(features.area_m2, expected.area_m2)
        np.testing.assert_array_equal(labels.lat, flipped.lat)
        self.assertTrue(((labels > 0) == flipped).all())

    def test_curvilinear(self):
        """A curvilinear grid gives about the same features as a regular one."""
        lon2, lat2 = np.meshgrid(self.lon, self.lat)
        curvilinear = xr.DataArray(
            self.mask.values,
            dims=("time", "y", "x"),
            coords={
                "time": self.mask.time.values,
                "lat": (("y", "x"), lat2),
                "lon": (("y", "x"), lon2),
            },
        )
        features = find_features(curvilinear)
        expected = find_features(self.mask)
        self.assertEqual(len(features), len(expected))
        np.testing.assert_allclose(features.area_m2, expected.area_m2, rtol=0.01)
        np.testing.assert_allclose(features.lon, expected.lon, atol=0.01)
        geodesic = [abs(GEOD.geometry_area_perimeter(g)[0]) for g in features.geometry]
        np.testing.assert_allclose(geodesic, features.area_m2, rtol=0.02)

    def test_around_the_globe(self):
        """Features across the ends of a global grid are merged and split at 180 deg."""
        lon = np.arange(0.0, 360.0, 0.5)
        lat = np.arange(-30.0, 30.1, 0.5)
        mask = _discs(lon, lat, [(0.0, 0.0, 300e3, [0]), (180.0, 10.0, 300e3, [0])], 1)
        features = find_features(mask, frame_duration=timedelta(hours=1))
        self.assertEqual(len(features), 2)
        # cells of about 55 km quantize the discs' areas
        np.testing.assert_allclose(features.area_m2, np.pi * 300e3**2, rtol=0.08)
        np.testing.assert_allclose(
            np.abs((features.lon + 180) % 360 - 180), [0.0, 180.0], atol=0.05
        )
        for geometry in features.geometry:
            minx, _, maxx, _ = geometry.bounds
            self.assertGreaterEqual(minx, -180)
            self.assertLessEqual(maxx, 180)
        geodesic = [abs(GEOD.geometry_area_perimeter(g)[0]) for g in features.geometry]
        np.testing.assert_allclose(geodesic, features.area_m2, rtol=0.02)
        self.assertEqual(features.geometry[1].geom_type, "MultiPolygon")

    def test_alignment(self):
        """Frames may be centered on their times."""
        features = find_features(self.mask, alignment="center")
        self.assertEqual(features.start[0], pd.Timestamp("2025-12-31 23:45", tz="UTC"))

    def test_errors_and_empty(self):
        """Masks need a time; masks without features give none."""
        with self.assertRaisesRegex(ValueError, "time dimension"):
            find_features(self.mask.isel(time=0, drop=True))
        with self.assertRaisesRegex(ValueError, "frame_duration"):
            find_features(self.mask.isel(time=[0]))
        empty = find_features(self.mask & False, variables={"c": self.mask * 1.0})
        self.assertEqual(len(empty), 0)
        self.assertIn("mean_c", empty.columns)
        with self.assertRaisesRegex(ValueError, "reductions"):
            find_features(self.mask, reductions=["median"])

    def test_feature_observations(self):
        """Features are observed as collect_feature_observations takes them."""
        # pylint: disable=import-outside-toplevel
        from datetime import datetime, timezone

        from tatc.analysis import collect_feature_observations
        from tatc.schemas import Instrument, Satellite, SunSynchronousOrbit

        mask = self.mask.assign_coords(
            time=pd.date_range("2026-03-01", periods=2, freq="12h")
        )
        features = find_features(mask)
        satellite = Satellite(
            name="sat",
            orbit=SunSynchronousOrbit(
                altitude=700e3,
                equator_crossing_time="10:30",
                epoch=datetime(2026, 3, 1, tzinfo=timezone.utc),
            ),
            instruments=[Instrument(name="wide", field_of_regard=120)],
        )
        times = [
            datetime(2026, 3, 1, tzinfo=timezone.utc) + timedelta(minutes=m)
            for m in range(0, 24 * 60, 1)
        ]
        looks = collect_feature_observations(features, satellite, times)
        self.assertGreater(len(looks), 0)
