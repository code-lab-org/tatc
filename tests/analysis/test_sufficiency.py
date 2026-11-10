"""
Unit tests for methods to assess the statistical sufficiency of observations.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from tatc.analysis import (
    compare_sampled_means,
    compute_local_time_coverage,
    compute_revisit_statistics,
    estimate_sampling_error,
    sample_field,
    sampling_error_variance,
)
from tatc.utils import hash_geometry

HAS_XARRAY = importlib.util.find_spec("xarray") is not None
START = datetime(2026, 1, 1, tzinfo=timezone.utc)
DAY = timedelta(days=1)


def _observations(samples):
    """Observations of points (longitude, latitude) at hours after the start."""
    rows = [
        {
            "geometry": Point(lon, lat),
            "target_hash": hash_geometry(Point(lon, lat)),
            "epoch": START + timedelta(hours=hour),
        }
        for (lon, lat), hours in samples.items()
        for hour in hours
    ]
    return gpd.GeoDataFrame(rows, crs="EPSG:4326")


class TestSamplingErrorVariance(unittest.TestCase):
    """
    Unit tests for the variance of sampling errors.
    """

    def test_limits(self):
        """Independent samples give 1/n; dense samples give no error."""
        length = 1000.0
        # samples far apart relative to the decorrelation time
        sparse = sampling_error_variance(np.linspace(50, 950, 10), length, 1.0)
        self.assertAlmostEqual(sparse, 1 / 10, delta=0.01)
        dense = sampling_error_variance(np.arange(0.5, 1000, 1.0), length, 10.0)
        self.assertLess(dense, 0.001)
        # a field that hardly varies over the period is estimated by one sample
        self.assertLess(sampling_error_variance([500.0], length, 1e6), 0.01)
        self.assertTrue(np.isnan(sampling_error_variance([], length, 1.0)))
        # samples closer than the decorrelation time repeat each other
        clustered = sampling_error_variance(np.linspace(0, 1, 10), length, 10.0)
        self.assertGreater(clustered, 0.5)
        with self.assertRaises(ValueError):
            sampling_error_variance([1.0], 0, 1.0)

    def test_monte_carlo(self):
        """The variance matches simulations of a first-order Markov process."""
        rng = np.random.default_rng(0)
        dt, length, scale = 0.1, 48.0, 3.0
        steps = int(length / dt)
        rho = np.exp(-dt / scale)
        n = 4000
        x = np.empty((n, steps))
        x[:, 0] = rng.normal(size=n)
        noise = rng.normal(size=(n, steps)) * np.sqrt(1 - rho**2)
        for k in range(1, steps):
            x[:, k] = rho * x[:, k - 1] + noise[:, k]
        grid = (np.arange(steps) + 0.5) * dt
        for times in ([10.0, 22.0, 34.0, 46.0], [12.0, 13.0, 14.0], [24.0]):
            index = [int(t / dt) for t in times]
            errors = x[:, index].mean(axis=1) - x.mean(axis=1)
            expected = sampling_error_variance(grid[index], length, scale)
            self.assertAlmostEqual(np.var(errors) / expected, 1, delta=0.08)


class TestEstimateSamplingError(unittest.TestCase):
    """
    Unit tests for estimates of sampling errors.
    """

    def test_targets(self):
        """Each target's error follows its samples, and unobserved targets have none."""
        observations = _observations(
            {(0, 0): [2, 8, 14, 20], (10, 0): [2, 3], (20, 0): [12]}
        )
        targets = gpd.GeoSeries([Point(0, 0), Point(10, 0), Point(20, 0), Point(30, 0)])
        error = estimate_sampling_error(
            observations,
            START,
            START + DAY,
            timedelta(hours=1),
            standard_deviation=2.0,
            targets=targets,
        )
        self.assertEqual(list(error.samples), [4, 2, 1, 0])
        self.assertTrue(
            error.relative_error[0] < error.relative_error[1] < error.relative_error[2]
        )
        np.testing.assert_allclose(
            error.sampling_error[:3], 2 * error.relative_error[:3]
        )
        self.assertTrue(np.isnan(error.relative_error[3]))
        self.assertEqual(error.effective_samples[3], 0)
        # four samples six hours apart are nearly independent; two an hour apart are not
        self.assertAlmostEqual(error.effective_samples[0], 4, delta=0.5)
        self.assertLess(error.effective_samples[1], 2)

    def test_standard_deviation_by_target(self):
        """The standard deviation may differ between targets."""
        observations = _observations({(0, 0): [6, 18], (10, 0): [6, 18]})
        sigma = {hash_geometry(Point(0, 0)): 1.0, hash_geometry(Point(10, 0)): 3.0}
        error = estimate_sampling_error(
            observations, START, START + DAY, timedelta(hours=2), sigma
        )
        self.assertAlmostEqual(error.sampling_error[1], 3 * error.sampling_error[0])


class TestRevisitStatistics(unittest.TestCase):
    """
    Unit tests for revisit statistics.
    """

    def test_gaps(self):
        """Revisits are between samples; gaps include the edges of the period."""
        observations = _observations({(0, 0): [1, 3, 11], (10, 0): [12]})
        targets = gpd.GeoSeries([Point(0, 0), Point(10, 0), Point(20, 0)])
        stats = compute_revisit_statistics(
            observations, START, START + DAY, quantiles=(0.5,), targets=targets
        )
        self.assertEqual(list(stats.samples), [3, 1, 0])
        self.assertEqual(stats.revisit[0], timedelta(hours=5))
        self.assertEqual(stats["revisit p50"][0], timedelta(hours=5))
        self.assertEqual(stats.max_gap[0], timedelta(hours=13))
        self.assertEqual(stats.max_gap[1], timedelta(hours=12))
        self.assertTrue(pd.isna(stats.revisit[1]))
        self.assertEqual(stats.max_gap[2], DAY)


class TestLocalTimeCoverage(unittest.TestCase):
    """
    Unit tests for local time coverage.
    """

    def test_bins(self):
        """Samples are binned by local mean solar time."""
        # at 90 degrees east, local time is six hours after universal time
        observations = _observations({(90, 0): [0, 1, 12], (-90, 0): [10, 22]})
        coverage = compute_local_time_coverage(observations, bins=4)
        east = coverage.iloc[0]
        self.assertEqual(east["local 06:00"], 2)
        self.assertEqual(east["local 18:00"], 1)
        self.assertEqual(east.local_bins_sampled, 2)
        west = coverage.iloc[1]
        # 04:00 and 16:00 local time
        self.assertEqual((west["local 00:00"], west["local 12:00"]), (1, 1))
        self.assertEqual(west.local_time_coverage, 0.5)
        with self.assertRaises(ValueError):
            compute_local_time_coverage(observations, bins=0)


@unittest.skipUnless(HAS_XARRAY, "requires xarray")
class TestSampleField(unittest.TestCase):
    """
    Unit tests for sampling fields.
    """

    def setUp(self):
        import xarray as xr  # pylint: disable=import-outside-toplevel

        times = pd.date_range("2026-01-01", periods=25, freq="h")
        lat = np.arange(-10.0, 11.0, 1.0)
        lon = np.arange(0.0, 360.0, 1.0)
        hours = np.arange(25, dtype=float)
        # a field linear in time, latitude, and longitude (0-360)
        values = (
            hours[:, None, None] + 10 * lat[None, :, None] + 0.1 * lon[None, None, :]
        )
        self.field = xr.DataArray(
            values,
            dims=("time", "lat", "lon"),
            coords={"time": times, "lat": lat, "lon": lon},
            name="temperature",
        )

    def test_linear(self):
        """Linear interpolation recovers a linear field, across longitude conventions."""
        observations = _observations({(10.5, 2.25): [3.5], (-20, -1): [6]})
        sampled = sample_field(observations, self.field, method="linear")
        expected = [3.5 + 22.5 + 1.05, 6 - 10 + 0.1 * 340]
        np.testing.assert_allclose(sampled.temperature, expected)

    def test_nearest_and_outside(self):
        """Nearest values snap to the grid; times outside the field have no value."""
        observations = _observations({(10.4, 2.2): [3.4, 30]})
        sampled = sample_field(observations, self.field, name="t")
        self.assertAlmostEqual(sampled.t[0], 3 + 20 + 1.0)
        self.assertTrue(np.isnan(sampled.t[1]))

    def test_compare_means(self):
        """Sampled means are compared with the means of the field over the period."""
        observations = _observations({(10, 2): [0, 24], (20, 0): [12]})
        sampled = sample_field(observations, self.field)
        comparison = compare_sampled_means(
            sampled,
            self.field,
            START,
            START + DAY,
            targets=gpd.GeoSeries([Point(10, 2), Point(20, 0), Point(30, 0)]),
        )
        # the field's mean over hours 0 to 24 is its value at hour 12
        np.testing.assert_allclose(comparison.true_mean, [12 + 20 + 1, 12 + 2, 12 + 3])
        np.testing.assert_allclose(comparison.error[:2], [0, 0], atol=1e-9)
        self.assertEqual(list(comparison.samples), [2, 1, 0])
        self.assertTrue(np.isnan(comparison.sampled_mean[2]))
