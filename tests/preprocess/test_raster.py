"""
Unit tests for the tatc.preprocess.raster module, using small synthetic
GeoTIFF rasters written to a temporary directory.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

from tatc.generation import generate_points_random
from tatc.preprocess import read_raster_weights


def _write_raster(path, data, west, north, resolution, crs="EPSG:4326", nodata=-1):
    """
    Writes a single-band float32 GeoTIFF with its north-west corner at
    (`west`, `north`) and square cells of `resolution` degrees.
    """
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype="float32",
        crs=crs,
        transform=from_origin(west, north, resolution, resolution),
        nodata=nodata,
    ) as dst:
        dst.write(data.astype(np.float32), 1)


class TestReadRasterWeights(unittest.TestCase):
    """
    Unit tests for the tatc.preprocess.read_raster_weights function.
    """

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = str(Path(self.directory.name) / "weights.tif")
        # 1-degree cells over 0 to 10 longitude, 0 to 6 latitude, with
        # values numbering the cells and one nodata cell
        self.data = np.arange(60, dtype=np.float64).reshape(6, 10)
        self.data[0, 0] = -1
        _write_raster(self.path, self.data, 0, 6, 1)

    def tearDown(self):
        self.directory.cleanup()

    def test_read_all(self):
        """
        Test that the full raster is read with nodata as zero weight.
        """
        weights, bounds = read_raster_weights(self.path)
        expected = self.data.copy()
        expected[0, 0] = 0
        np.testing.assert_array_equal(weights, expected)
        self.assertEqual(bounds, (0, 0, 10, 6))

    def test_read_window(self):
        """
        Test that only cells overlapping the mask bounds are read.
        """
        weights, bounds = read_raster_weights(self.path, mask=box(2.5, 1.5, 4.5, 3))
        np.testing.assert_array_equal(weights, self.data[3:5, 2:5])
        self.assertEqual(bounds, (2, 1, 5, 3))

    def test_read_window_beyond_180(self):
        """
        Test that a raster extending slightly beyond 180 degrees longitude
        (as some global rasters do) is still read only within the mask
        bounds, while a mask beyond 180 degrees reads all longitudes.
        """
        path = str(Path(self.directory.name) / "beyond.tif")
        _write_raster(path, np.ones((18, 37)), -180.5, 90, 10)
        weights, bounds = read_raster_weights(path, mask=box(0, 0, 20, 10))
        self.assertEqual(weights.shape, (1, 3))
        self.assertEqual(bounds, (-0.5, 0, 29.5, 10))
        weights, bounds = read_raster_weights(path, mask=box(170, 0, 190, 10))
        self.assertEqual(weights.shape, (1, 37))
        self.assertEqual(bounds, (-180.5, 0, 189.5, 10))

    def test_aggregate_sum(self):
        """
        Test that aggregated counts are summed, with blocks shifted to stay
        within the raster.
        """
        weights, bounds = read_raster_weights(
            self.path, mask=box(7.5, 0, 10, 2), aggregate=2
        )
        # the window 7 to 10 expands to 4 columns, shifted to 6 to 10
        self.assertEqual(bounds, (6, 0, 10, 2))
        np.testing.assert_array_equal(
            weights, [[self.data[4:6, 6:8].sum(), self.data[4:6, 8:10].sum()]]
        )

    def test_aggregate_density_padded(self):
        """
        Test that aggregated densities are averaged over cells within the
        raster, padding a raster smaller than a whole number of blocks.
        """
        weights, bounds = read_raster_weights(self.path, aggregate=4, density=True)
        self.assertEqual(bounds, (0, -2, 12, 6))
        data = self.data.copy()
        data[0, 0] = 0
        np.testing.assert_allclose(weights[0, 0], data[0:4, 0:4].mean())
        np.testing.assert_allclose(weights[1, 2], data[4:6, 8:10].mean())

    def test_aggregate_beyond_pole(self):
        """
        Test that an aggregate extending a global raster beyond the poles
        raises a ValueError.
        """
        path = str(Path(self.directory.name) / "global.tif")
        _write_raster(path, np.ones((18, 36)), -180, 90, 10)
        with self.assertRaises(ValueError):
            read_raster_weights(path, aggregate=4)
        weights, bounds = read_raster_weights(path, aggregate=3)
        self.assertEqual(weights.shape, (6, 12))
        self.assertEqual(bounds, (-180, -90, 180, 90))

    def test_invalid(self):
        """
        Test that a projected raster, a raster outside the mask, and an
        aggregate below 1 raise a ValueError.
        """
        path = str(Path(self.directory.name) / "projected.tif")
        _write_raster(path, self.data, 500000, 4000000, 1000, crs="EPSG:32617")
        with self.assertRaises(ValueError):
            read_raster_weights(path)
        with self.assertRaises(ValueError):
            read_raster_weights(self.path, mask=box(20, 20, 30, 30))
        with self.assertRaises(ValueError):
            read_raster_weights(self.path, aggregate=0)

    def test_generate_points(self):
        """
        Test generating points weighted by a raster read within a mask.
        """
        mask = box(2.5, 1.5, 4.5, 3)
        weights, bounds = read_raster_weights(self.path, mask=mask)
        points = generate_points_random(
            100, mask=mask, weights=weights, weights_bounds=bounds, seed=0
        )
        self.assertEqual(len(points), 100)
        self.assertTrue(points.intersects(mask).all())


if __name__ == "__main__":
    unittest.main()
