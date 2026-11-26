"""
Tests for sampling fields over the footprints and regions of observations.
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Point, box

from tatc import constants
from tatc.analysis import (
    collect_ground_pixels,
    compute_sample_sizes,
    sample_field,
)
from tatc.schemas import PointedInstrument, Satellite, SunSynchronousOrbit

xr = pytest.importorskip("xarray")
pytest.importorskip("scipy")

START = datetime(2026, 1, 1, tzinfo=timezone.utc)
RADIUS = constants.EARTH_MEAN_RADIUS
SPACING = 0.02  # degrees, about 2.2 km


def _field(function, lon=None, lat=None, hours=3):
    """A field on a regular grid near the equator, from a function of
    (hours, longitude, latitude), hourly."""
    lon = np.arange(0.0, 4.0 + SPACING / 2, SPACING) if lon is None else lon
    lat = np.arange(-2.0, 2.0 + SPACING / 2, SPACING) if lat is None else lat
    hour = np.arange(float(hours))
    values = function(
        hour[:, None, None], lon[None, None, :], lat[None, :, None]
    ) * np.ones((len(hour), len(lat), len(lon)))
    return xr.DataArray(
        values,
        dims=("time", "lat", "lon"),
        coords={
            "time": pd.date_range("2026-01-01", periods=hours, freq="h"),
            "lat": lat,
            "lon": lon,
        },
        name="value",
    )


def _samples(places, along, cross, azimuth, hour=1.0):
    """Observations at places (longitude, latitude), with sample sizes (meters)
    and azimuths (degrees)."""
    count = len(places)
    return gpd.GeoDataFrame(
        {
            "epoch": [START + timedelta(hours=hour)] * count,
            "sample_along": np.broadcast_to(along, count).astype(float),
            "sample_cross": np.broadcast_to(cross, count).astype(float),
            "sample_azimuth": np.broadcast_to(azimuth, count).astype(float),
        },
        geometry=[Point(lon, lat) for lon, lat in places],
        crs="EPSG:4326",
    )


def _unstructured(field):
    """The field on an unstructured grid of the same cells."""
    lon, lat = np.meshgrid(field.lon.values, field.lat.values)
    return xr.DataArray(
        field.values.reshape(len(field.time), -1),
        dims=("time", "cell"),
        coords={
            "time": field.time.values,
            "lat": ("cell", lat.ravel()),
            "lon": ("cell", lon.ravel()),
        },
        name=field.name,
    )


def _east(lon):
    """The distance (meters) east of longitude 2 at the equator."""
    return RADIUS * np.radians(lon - 2.0)


class TestFootprintSampling(unittest.TestCase):
    """
    Unit tests for sampling fields over footprints.
    """

    def test_constant(self):
        """A constant field is sampled exactly, with nearly all its weight valid."""
        field = _field(lambda h, lon, lat: 7.0)
        samples = _samples([(2, 0), (1.5, 0.7)], 20e3, 10e3, [0, 37])
        for footprint in ("gaussian", "uniform"):
            sampled = sample_field(
                samples, field, method="footprint", footprint=footprint
            )
            np.testing.assert_allclose(sampled.value, 7.0)
            self.assertTrue(sampled.value_valid_fraction.between(0.97, 1.0).all())

    def test_linear(self):
        """A symmetric footprint returns a linear field's value at its center."""
        field = _field(lambda h, lon, lat: _east(lon) / 1e3 + 3 * lat)
        samples = _samples([(2.013, 0.004), (1.7, -0.5)], 30e3, 12e3, [20, 110])
        sampled = sample_field(samples, field, method="footprint")
        expected = _east(np.array([2.013, 1.7])) / 1e3 + 3 * np.array([0.004, -0.5])
        np.testing.assert_allclose(sampled.value, expected, atol=0.02)

    def test_sinusoid_attenuation(self):
        """A Gaussian footprint attenuates a sinusoid as its width along the wave."""
        wavelength = 60e3
        field = _field(lambda h, lon, lat: np.cos(2 * np.pi * _east(lon) / wavelength))
        along, cross = 30e3, 10e3
        for azimuth, width in ((90, along), (0, cross), (270, along)):
            sampled = sample_field(
                _samples([(2, 0)], along, cross, azimuth), field, method="footprint"
            )
            sigma = width / (2 * np.sqrt(2 * np.log(2)))
            expected = np.exp(-0.5 * (2 * np.pi * sigma / wavelength) ** 2)
            self.assertAlmostEqual(sampled.value[0], expected, delta=0.02)

    def test_edge(self):
        """A footprint centered on the edge of a regional field (half a cell beyond
        its outer cells) has half its weight."""
        field = _field(lambda h, lon, lat: 1.0)
        samples = _samples([(-SPACING / 2, 0.0)], 20e3, 20e3, 0)
        sampled = sample_field(samples, field, method="footprint")
        self.assertAlmostEqual(sampled.value_valid_fraction[0], 0.5, delta=0.05)
        sampled = sample_field(
            samples, field, method="footprint", min_valid_fraction=0.6
        )
        self.assertTrue(np.isnan(sampled.value[0]))

    def test_missing_cells(self):
        """Cells without values (such as land) reduce the valid share."""
        field = _field(lambda h, lon, lat: np.where(lon < 2.0, np.nan, 5.0))
        sampled = sample_field(
            _samples([(2.0, 0.0), (2.3, 0.0)], 20e3, 20e3, 0),
            field,
            method="footprint",
            min_valid_fraction=0.3,
        )
        np.testing.assert_allclose(sampled.value, 5.0)
        self.assertLess(sampled.value_valid_fraction[0], 0.7)
        self.assertGreater(sampled.value_valid_fraction[1], 0.97)

    def test_unstructured(self):
        """An unstructured grid gives the same footprint means as a regular one."""
        field = _field(lambda h, lon, lat: np.sin(lon * 3) + np.cos(lat * 5))
        samples = _samples([(2, 0), (1.2, 1.1), (3.1, -0.4)], 25e3, 8e3, [0, 45, 300])
        regular = sample_field(samples, field, method="footprint")
        other = sample_field(samples, _unstructured(field), method="footprint")
        np.testing.assert_allclose(other.value, regular.value, rtol=1e-3)
        np.testing.assert_allclose(
            other.value_valid_fraction, regular.value_valid_fraction, rtol=2e-2
        )

    def test_antimeridian(self):
        """Footprints across the antimeridian are weighted on either side."""
        lon = np.arange(0.0, 360.0, 0.25)
        lat = np.arange(-3.0, 3.01, 0.25)
        # a field linear in the signed distance east of the antimeridian
        field = _field(
            lambda h, lon, lat: lon % 360 - 180 + 0 * lat,
            lon=lon,
            lat=lat,
        )
        field = field.where(np.abs(field) < 10, 0.0)
        samples = _samples([(180.0, 0.0), (-179.9, 0.0)], 100e3, 100e3, 0)
        sampled = sample_field(samples, field, method="footprint")
        np.testing.assert_allclose(sampled.value, [0.0, 0.1], atol=0.01)
        self.assertTrue(sampled.value_valid_fraction.gt(0.97).all())
        unstructured = sample_field(samples, _unstructured(field), method="footprint")
        np.testing.assert_allclose(unstructured.value, sampled.value, atol=0.01)

    def test_small_footprint(self):
        """A footprint smaller than a cell takes the value of its nearest cell."""
        field = _field(lambda h, lon, lat: lon)
        sampled = sample_field(
            _samples([(2.011, 0.0)], 500.0, 500.0, 0), field, method="footprint"
        )
        self.assertAlmostEqual(sampled.value[0], 2.02)
        self.assertEqual(sampled.value_valid_fraction[0], 1.0)

    def test_without_sizes(self):
        """Footprints need sample sizes; samples without them take the nearest cell."""
        field = _field(lambda h, lon, lat: lon)
        samples = _samples([(2.011, 0.0)], np.nan, np.nan, np.nan)
        sampled = sample_field(samples, field, method="footprint")
        self.assertAlmostEqual(sampled.value[0], 2.02)
        self.assertTrue(np.isnan(sampled.value_valid_fraction[0]))
        with self.assertRaisesRegex(ValueError, "compute_sample_sizes"):
            sample_field(
                samples.drop(columns="sample_azimuth"), field, method="footprint"
            )

    def test_time(self):
        """Footprints take the nearest time step, or interpolate between steps."""
        field = _field(lambda h, lon, lat: 10 * h + 0 * lon)
        samples = _samples([(2, 0)], 20e3, 20e3, 0, hour=1.3)
        nearest = sample_field(samples, field, method="footprint")
        linear = sample_field(samples, field, method="footprint", time_method="linear")
        self.assertAlmostEqual(nearest.value[0], 10.0)
        self.assertAlmostEqual(linear.value[0], 13.0)
        late = _samples([(2, 0)], 20e3, 20e3, 0, hour=2.4)
        self.assertAlmostEqual(
            sample_field(late, field, method="footprint").value[0], 20
        )
        self.assertTrue(
            np.isnan(
                sample_field(
                    late, field, method="footprint", time_method="linear"
                ).value[0]
            )
        )

    def test_extra_dimensions(self):
        """Fields with a vertical dimension need a level or weights."""
        field = _field(lambda h, lon, lat: lon).expand_dims(level=[1000.0, 500.0])
        with self.assertRaisesRegex(ValueError, "give a level or weights"):
            sample_field(_samples([(2, 0)], 20e3, 20e3, 0), field, method="footprint")

    def test_observed_pixels(self):
        """The pixels of a scanner are sampled over their footprints."""
        scanner = PointedInstrument(
            name="scanner",
            field_of_regard=60,
            cross_track_field_of_view=60,
            along_track_field_of_view=1.0,
            cross_track_pixels=40,
            is_rectangular=True,
        )
        satellite = Satellite(
            name="sat",
            orbit=SunSynchronousOrbit(
                altitude=700e3, equator_crossing_time="10:30", epoch=START
            ),
            instruments=[scanner],
        )
        times = [START + timedelta(seconds=10 * k) for k in range(8640)]
        pixels = collect_ground_pixels(satellite, times, mask=box(0.2, -1.8, 3.8, 1.8))
        pixels = compute_sample_sizes(
            pixels[pixels.valid_obs].reset_index(drop=True), satellite
        )
        self.assertGreater(len(pixels), 0)
        field = _field(lambda h, lon, lat: 2.0 + 0 * lon, hours=24)
        sampled = sample_field(pixels, field, method="footprint")
        self.assertGreater(sampled.value.notna().sum(), len(pixels) // 2)
        np.testing.assert_allclose(sampled.value.dropna(), 2.0)


class TestAreaSampling(unittest.TestCase):
    """
    Unit tests for sampling fields inside regions.
    """

    def test_region_mean(self):
        """A region takes the area-weighted mean of the cells inside it."""
        field = _field(lambda h, lon, lat: lon + 0 * lat)
        regions = gpd.GeoDataFrame(
            {"epoch": [START + timedelta(hours=1)] * 2},
            geometry=[box(1.0, -0.5, 2.0, 0.5), box(-0.5, -0.5, 0.5, 0.5)],
            crs="EPSG:4326",
        )
        sampled = sample_field(regions, field, method="area", min_valid_fraction=0.3)
        self.assertAlmostEqual(sampled.value[0], 1.5, delta=0.01)
        self.assertGreater(sampled.value_valid_fraction[0], 0.97)
        # half outside the field
        self.assertAlmostEqual(sampled.value_valid_fraction[1], 0.5, delta=0.05)
        self.assertAlmostEqual(sampled.value[1], 0.25, delta=0.02)

    def test_points(self):
        """Points (without area) take the value of their nearest cell."""
        field = _field(lambda h, lon, lat: lon + 0 * lat)
        points = gpd.GeoDataFrame(
            {"epoch": [START + timedelta(hours=1)]},
            geometry=[Point(2.011, 0.0)],
            crs="EPSG:4326",
        )
        sampled = sample_field(points, field, method="area")
        self.assertAlmostEqual(sampled.value[0], 2.02)
        self.assertEqual(sampled.value_valid_fraction[0], 1.0)

    def test_empty(self):
        """Empty observations get empty columns."""
        field = _field(lambda h, lon, lat: lon)
        empty = gpd.GeoDataFrame(
            {"epoch": pd.to_datetime([], utc=True)}, geometry=[], crs="EPSG:4326"
        )
        sampled = sample_field(empty, field, method="area")
        self.assertIn("value_valid_fraction", sampled.columns)
        with self.assertRaisesRegex(ValueError, "method must be"):
            sample_field(empty, field, method="mean")
