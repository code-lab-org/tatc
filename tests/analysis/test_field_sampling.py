"""
Tests for sampling fields on grids.
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point

from tatc.analysis import compare_sampled_means, sample_field
from tatc.utils import hash_geometry

xr = pytest.importorskip("xarray")
pytest.importorskip("scipy")

START = datetime(2026, 1, 1, tzinfo=timezone.utc)


def _observations(places, hours):
    """Observations of places (longitude, latitude), each at hours after the start."""
    rows = [
        {
            "geometry": Point(lon, lat),
            "target_hash": hash_geometry(Point(lon, lat)),
            "epoch": START + timedelta(hours=hour),
        }
        for (lon, lat), hour in zip(places, hours)
    ]
    return gpd.GeoDataFrame(rows, crs="EPSG:4326")


def _regular():
    """A regional field, hourly over a day, linear in time, latitude, and longitude."""
    times = pd.date_range("2026-01-01", periods=25, freq="h")
    lat = np.arange(30.0, 41.0, 0.5)
    lon = np.arange(-110.0, -94.5, 0.5)
    values = (
        np.arange(25.0)[:, None, None]
        + 10 * lat[None, :, None]
        + 0.1 * lon[None, None, :]
    )
    return xr.DataArray(
        values,
        dims=("time", "lat", "lon"),
        coords={"time": times, "lat": lat, "lon": lon},
        name="temperature",
    )


def _curvilinear(field):
    """The field on a curvilinear grid: two-dimensional coordinates over (y, x)."""
    lon, lat = np.meshgrid(field.lon.values, field.lat.values)
    return xr.DataArray(
        field.values,
        dims=("time", "y", "x"),
        coords={
            "time": field.time.values,
            "lat": (("y", "x"), lat),
            "lon": (("y", "x"), lon),
        },
        name=field.name,
    )


def _unstructured(field):
    """The field on an unstructured grid: coordinates of a dimension of cells."""
    curvilinear = _curvilinear(field)
    return xr.DataArray(
        curvilinear.values.reshape(len(field.time), -1),
        dims=("time", "cell"),
        coords={
            "time": field.time.values,
            "lat": ("cell", curvilinear.lat.values.ravel()),
            "lon": ("cell", curvilinear.lon.values.ravel()),
        },
        name=field.name,
    )


class TestGrids(unittest.TestCase):
    """
    Unit tests for sampling fields on regular, curvilinear, and unstructured grids.
    """

    def setUp(self):
        self.field = _regular()
        rng = np.random.default_rng(1)
        # places near the centers of cells, at times between steps
        lon = rng.choice(self.field.lon.values[1:-1], 40) + rng.uniform(-0.2, 0.2, 40)
        lat = rng.choice(self.field.lat.values[1:-1], 40) + rng.uniform(-0.2, 0.2, 40)
        self.observations = _observations(list(zip(lon, lat)), rng.uniform(0, 24, 40))

    def test_same_nearest_values(self):
        """Every kind of grid gives the same nearest values."""
        regular = sample_field(self.observations, self.field)
        self.assertFalse(regular.temperature.isna().any())
        for other in (_curvilinear(self.field), _unstructured(self.field)):
            sampled = sample_field(self.observations, other)
            np.testing.assert_allclose(sampled.temperature, regular.temperature)

    def test_outside_extent(self):
        """Places beyond the edge of a grid by more than half a cell have no value."""
        observations = _observations(
            [(-102.0, 40.7), (-102.0, 41.0), (-94.85, 35.0), (-94.4, 35.0)],
            [3, 3, 3, 3],
        )
        for field in (_curvilinear(self.field), _unstructured(self.field)):
            sampled = sample_field(observations, field)
            self.assertFalse(np.isnan(sampled.temperature[0]))
            self.assertTrue(np.isnan(sampled.temperature[1]))
            self.assertFalse(np.isnan(sampled.temperature[2]))
            self.assertTrue(np.isnan(sampled.temperature[3]))

    def test_outside_period(self):
        """Times up to half a step beyond the field's period have the nearest value."""
        observations = _observations([(-102.0, 35.0)] * 4, [-0.4, -0.6, 24.4, 24.6])
        for field in (self.field, _unstructured(self.field)):
            sampled = sample_field(observations, field)
            self.assertAlmostEqual(sampled.temperature[0], 350 - 10.2)
            self.assertTrue(np.isnan(sampled.temperature[1]))
            self.assertAlmostEqual(sampled.temperature[2], 24 + 350 - 10.2)
            self.assertTrue(np.isnan(sampled.temperature[3]))

    def test_corners_inside(self):
        """Places at the corners of cells inside a grid of oblong cells have values."""
        lon = self.field.lon.values[1:-1] + 0.25
        lat = np.full(len(lon), 40.25)
        observations = _observations(list(zip(lon, lat)), [3] * len(lon))
        sampled = sample_field(observations, _unstructured(self.field))
        self.assertFalse(sampled.temperature.isna().any())

    def test_linear_requires_regular(self):
        """Linear interpolation needs a regular grid."""
        with self.assertRaisesRegex(ValueError, "regular grid"):
            sample_field(self.observations, _unstructured(self.field), method="linear")

    def test_antimeridian(self):
        """Cells across the antimeridian are nearest to places on either side."""
        lon = np.array([178.0, 179.0, -180.0, -179.0, -178.0])
        field = xr.DataArray(
            np.arange(5.0)[None, :],
            dims=("time", "cell"),
            coords={
                "time": [np.datetime64("2026-01-01")],
                "lat": ("cell", np.zeros(5)),
                "lon": ("cell", lon),
            },
            name="index",
        )
        observations = _observations([(179.9, 0), (-179.9, 0), (180.6, 0)], [0, 0, 0])
        sampled = sample_field(observations, field)
        np.testing.assert_allclose(sampled["index"], [2, 2, 3])

    def test_radians(self):
        """Coordinates in radians (as of an MPAS mesh) are converted."""
        field = _unstructured(self.field)
        field = field.assign_coords(
            lat=("cell", np.radians(field.lat.values), {"units": "radian"}),
            lon=("cell", np.radians(field.lon.values), {"units": "radians"}),
        )
        sampled = sample_field(self.observations, field)
        regular = sample_field(self.observations, self.field)
        np.testing.assert_allclose(sampled.temperature, regular.temperature)

    def test_coordinates_in_time(self):
        """Coordinates that vary in time (as of a WRF run) are taken at the first time."""
        curvilinear = _curvilinear(self.field)
        times = len(curvilinear.time)
        field = curvilinear.assign_coords(
            lat=(("time", "y", "x"), np.repeat(curvilinear.lat.values[None], times, 0)),
            lon=(("time", "y", "x"), np.repeat(curvilinear.lon.values[None], times, 0)),
        )
        sampled = sample_field(self.observations, field)
        regular = sample_field(self.observations, self.field)
        np.testing.assert_allclose(sampled.temperature, regular.temperature)

    def test_named_coordinates(self):
        """Coordinates of other names are given by the dimensions."""
        field = _unstructured(self.field).rename(lat="latCell", lon="lonCell")
        field = field.rename(latCell="la", lonCell="lo")
        sampled = sample_field(
            self.observations, field, dimensions={"lat": "la", "lon": "lo"}
        )
        regular = sample_field(self.observations, self.field)
        np.testing.assert_allclose(sampled.temperature, regular.temperature)
        with self.assertRaisesRegex(ValueError, "no latitude and longitude"):
            sample_field(self.observations, field)

    def test_compare_means(self):
        """The true means at targets are the same on every kind of grid."""
        observations = _observations([(-102.0, 35.0), (-100.0, 33.0)], [6, 18])
        targets = gpd.GeoSeries(
            [Point(-102.0, 35.0), Point(-100.0, 33.0), Point(-90.0, 35.0)]
        )
        results = [
            compare_sampled_means(
                sample_field(observations, field),
                field,
                START,
                START + timedelta(days=1),
                targets=targets,
            )
            for field in (
                self.field,
                _curvilinear(self.field),
                _unstructured(self.field),
            )
        ]
        for result in results:
            np.testing.assert_allclose(result.true_mean, results[0].true_mean)
            self.assertTrue(np.isnan(result.true_mean[2]))
