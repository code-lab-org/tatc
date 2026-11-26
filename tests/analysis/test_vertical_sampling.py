"""
Tests for sampling fields with a vertical dimension.
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import MultiPoint, Point, box

from tatc.analysis import collect_limb_observations, sample_field, sample_profiles
from tatc.schemas import Satellite, SunSynchronousOrbit

xr = pytest.importorskip("xarray")
pytest.importorskip("scipy")

START = datetime(2026, 1, 1, tzinfo=timezone.utc)
HEIGHTS = np.arange(0.0, 20001.0, 1000.0)


def _field(levels=HEIGHTS, name="level", function=None):
    """A global field, hourly, on levels: by default, its level's height (km)
    plus its longitude (degrees east, -180 to 180) / 100."""
    lat = np.arange(-90.0, 90.1, 1.0)
    lon = np.arange(-180.0, 180.0, 1.0)
    hour = np.arange(3.0)
    if function is None:
        values = (
            levels[None, :, None, None] / 1e3
            + lon[None, None, None, :] / 100
            + 0 * lat[None, None, :, None]
            + 0 * hour[:, None, None, None]
        )
    else:
        values = function(
            hour[:, None, None, None],
            levels[None, :, None, None],
            lat[None, None, :, None],
            lon[None, None, None, :],
        )
    return xr.DataArray(
        values,
        dims=("time", name, "lat", "lon"),
        coords={
            "time": pd.date_range("2026-01-01", periods=3, freq="h"),
            name: levels,
            "lat": lat,
            "lon": lon,
        },
        name="temperature",
    )


def _observations(geometries, hour=1.0):
    """Observations of geometries at an hour after the start."""
    return gpd.GeoDataFrame(
        {"epoch": [START + timedelta(hours=hour)] * len(geometries)},
        geometry=list(geometries),
        crs="EPSG:4326",
    )


class TestLevels(unittest.TestCase):
    """
    Unit tests for sampling fields at levels and by weights of levels.
    """

    def setUp(self):
        self.field = _field()
        self.observations = _observations([Point(10, 5), Point(-20, -5)])

    def test_level(self):
        """A level is selected, or interpolated between levels."""
        sampled = sample_field(self.observations, self.field, level=5000.0)
        np.testing.assert_allclose(sampled.temperature, [5.1, 4.8])
        sampled = sample_field(self.observations, self.field, level=5500.0)
        np.testing.assert_allclose(sampled.temperature, [5.6, 5.3])
        sampled = sample_field(self.observations, self.field, level=30000.0)
        self.assertTrue(sampled.temperature.isna().all())

    def test_weights(self):
        """Weights of the levels sum the levels, by any method."""
        weights = xr.DataArray(np.full(len(HEIGHTS), 1 / len(HEIGHTS)), dims="level")
        sampled = sample_field(self.observations, self.field, weights=weights)
        np.testing.assert_allclose(sampled.temperature, [10.1, 9.8])
        regions = _observations([box(9.5, 4.5, 10.5, 5.5)])
        area = sample_field(regions, self.field, weights=weights, method="area")
        self.assertAlmostEqual(area.temperature[0], 10.1, places=6)

    def test_weights_by_coordinate(self):
        """Weights with a vertical coordinate are aligned to the field's levels."""
        weights = xr.DataArray(
            [1.0, 1.0], dims="level", coords={"level": [2000.0, 1000.0]}
        ).reindex(level=HEIGHTS, fill_value=0.0)
        sampled = sample_field(self.observations, self.field, weights=weights)
        np.testing.assert_allclose(sampled.temperature, [3.2, 2.6])

    def test_weights_by_observation(self):
        """Weights that vary by observation weigh each observation's levels."""
        weights = np.zeros((2, len(HEIGHTS)))
        weights[0, 3] = 1.0
        weights[1, 10] = 0.5
        weights = xr.DataArray(weights, dims=("observation", "level"))
        sampled = sample_field(self.observations, self.field, weights=weights)
        np.testing.assert_allclose(sampled.temperature, [3.1, 0.5 * 9.8])
        with self.assertRaisesRegex(ValueError, "nearest"):
            sample_field(self.observations, self.field, weights=weights, method="area")

    def test_levels_without_values(self):
        """Levels without values (below the ground) are left out of the sum."""
        field = self.field.where(self.field.level > 0)
        weights = xr.DataArray(np.ones(len(HEIGHTS)), dims="level")
        sampled = sample_field(self.observations, field, weights=weights)
        expected = np.sum(HEIGHTS[1:] / 1e3) + 20 * np.array([0.1, -0.2])
        np.testing.assert_allclose(sampled.temperature, expected)

    def test_errors(self):
        """Levels and weights need a vertical dimension, and one of them."""
        with self.assertRaisesRegex(ValueError, "give a level or weights"):
            sample_field(self.observations, self.field)
        with self.assertRaisesRegex(ValueError, "not both"):
            sample_field(
                self.observations,
                self.field,
                level=1000.0,
                weights=xr.DataArray(np.ones(len(HEIGHTS)), dims="level"),
            )
        flat = self.field.isel(level=0, drop=True)
        with self.assertRaisesRegex(ValueError, "vertical dimension"):
            sample_field(self.observations, flat, level=1000.0)
        with self.assertRaisesRegex(ValueError, "weights must be over"):
            sample_field(
                self.observations,
                self.field,
                weights=xr.DataArray(np.ones(3), dims="other"),
            )
        with self.assertRaisesRegex(ValueError, "no vertical dimension named"):
            sample_field(
                self.observations,
                self.field,
                level=1000.0,
                dimensions={"level": "height"},
            )


class TestProfiles(unittest.TestCase):
    """
    Unit tests for `sample_profiles`.
    """

    def setUp(self):
        self.field = _field()

    def test_points_without_heights(self):
        """Points without heights are sampled at every level."""
        observations = _observations([Point(10, 5), Point(-20, -5)])
        profiles = sample_profiles(observations, self.field)
        self.assertEqual(len(profiles), 2 * len(HEIGHTS))
        self.assertEqual(list(profiles.observation[: len(HEIGHTS)]), [0] * len(HEIGHTS))
        np.testing.assert_allclose(profiles.level[: len(HEIGHTS)], HEIGHTS)
        np.testing.assert_allclose(
            profiles.temperature[: len(HEIGHTS)], HEIGHTS / 1e3 + 0.1
        )
        self.assertTrue(profiles.height.isna().all())

    def test_points_at_heights(self):
        """Points at heights are interpolated between levels, NaN beyond them."""
        observations = _observations(
            [Point(10, 5, 2500), MultiPoint([(10, 5, 7250), (11, 5, 25000)])]
        )
        observations["satellite"] = ["a", "b"]
        profiles = sample_profiles(observations, self.field)
        self.assertEqual(len(profiles), 3)
        np.testing.assert_allclose(profiles.height, [2500, 7250, 25000])
        np.testing.assert_allclose(profiles.temperature[:2], [2.6, 7.35])
        self.assertTrue(np.isnan(profiles.temperature[2]))
        self.assertEqual(list(profiles.observation), [0, 1, 1])
        self.assertEqual(list(profiles.satellite), ["a", "b", "b"])
        self.assertAlmostEqual(profiles.geometry[2].x, 11)

    def test_heights_of_pressure_levels(self):
        """Heights of levels (from the top down) place points at heights."""
        pressure = np.array([250.0, 500.0, 850.0, 1000.0])
        level_heights = np.array([10400.0, 5600.0, 1500.0, 100.0])
        field = _field(
            levels=pressure,
            name="pressure",
            function=lambda h, p, lat, lon: p + 0 * lat + 0 * lon + 0 * h,
        )
        heights = xr.DataArray(
            level_heights, dims="pressure", coords={"pressure": pressure}
        )
        observations = _observations([Point(0, 0, 3550), Point(0, 0)])
        profiles = sample_profiles(observations, field, heights=heights)
        self.assertAlmostEqual(profiles.temperature[0], (850 + 500) / 2)
        np.testing.assert_allclose(profiles.level[1:], pressure)
        np.testing.assert_allclose(profiles.height[1:], level_heights)

    def test_heights_varying_in_space(self):
        """Heights of levels may vary in time and space."""
        field = _field()
        heights = field.copy(
            data=np.broadcast_to(
                HEIGHTS[None, :, None, None]
                + 100 * field.lat.values[None, None, :, None],
                field.shape,
            ).copy()
        )
        observations = _observations([Point(0, 2, 3200)])
        profiles = sample_profiles(observations, field, heights=heights)
        # 3200 m is level 3000 m + 200 m there, so 3.0 km plus 0 longitude
        self.assertAlmostEqual(profiles.temperature[0], 3.0)

    def test_limb_observations(self):
        """The tangent points of limb scans are sampled at their heights."""
        satellite = Satellite(
            name="sat",
            orbit=SunSynchronousOrbit(
                altitude=600e3, equator_crossing_time="13:30", epoch=START
            ),
        )
        observations = collect_limb_observations(
            satellite,
            [START + timedelta(minutes=m) for m in range(0, 60, 10)],
            scan_azimuth=0,
            scan_elevations=[5e3, 10e3, 15e3],
            scan_duration=timedelta(seconds=10),
        )
        profiles = sample_profiles(observations, self.field)
        self.assertEqual(len(profiles), 3 * len(observations))
        expected = profiles.height / 1e3 + profiles.geometry.x / 100
        np.testing.assert_allclose(profiles.temperature, expected, atol=0.006)

    def test_without_vertical_dimension(self):
        """Profiles need a vertical dimension."""
        with self.assertRaisesRegex(ValueError, "vertical dimension"):
            sample_profiles(
                _observations([Point(0, 0)]), self.field.isel(level=0, drop=True)
            )

    def test_empty(self):
        """Empty observations give empty profiles."""
        profiles = sample_profiles(_observations([]), self.field)
        self.assertEqual(len(profiles), 0)
        self.assertIn("temperature", profiles.columns)
