"""
Unit tests for the tatc.preprocess.terrain module.

These tests use a small synthetic GeoTIFF written to a temporary
directory and never access the network, so they remain fast and
deterministic regardless of connectivity.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from tatc.preprocess import (
    compute_terrain_mask,
    compute_terrain_mask_for_station,
    get_copernicus_dem_tile_urls,
    sample_dem_elevation,
)
from tatc.schemas import RadarStation


def _write_ridge_dem(path: str) -> None:
    """
    Writes a small synthetic DEM, centered near (0, 0): flat terrain at 0
    meters everywhere, except a 2000-meter "ridge" directly east of the
    origin at roughly 20-36 km range, used to test that terrain masking
    identifies a single, specific blocked azimuth (90 degrees) and leaves
    all others clear.
    """
    resolution = 0.001  # ~111 m/pixel at the equator
    width = height = 2000
    transform = from_origin(-1.0, 1.0, resolution, resolution)
    data = np.zeros((height, width), dtype=np.float32)
    col_min = int((0.18 - (-1.0)) / resolution)
    col_max = int((0.32 - (-1.0)) / resolution)
    row_center = int((1.0 - 0.0) / resolution)
    data[row_center - 5 : row_center + 5, col_min:col_max] = 2000.0
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=-9999,
    ) as dst:
        dst.write(data, 1)


class TestGetCopernicusDemTileUrls(unittest.TestCase):
    """
    Unit tests for get_copernicus_dem_tile_urls.
    """

    def test_single_tile_northeast(self):
        """
        Test tile naming for a point in the northern, eastern hemisphere.
        """
        urls = get_copernicus_dem_tile_urls(10.5, 20.5)
        self.assertEqual(len(urls), 1)
        self.assertIn("N20_00_E010_00", urls[0])

    def test_single_tile_southwest(self):
        """
        Test tile naming for a point in the southern, western hemisphere.
        Tiles are named by their lower-left (southernmost, westernmost)
        corner, so (-70.5, -30.5) falls in tile S31/W071, which spans
        [-31, -30] latitude and [-71, -70] longitude (verified against the
        real dataset's published tile bounds).
        """
        urls = get_copernicus_dem_tile_urls(-70.5, -30.5)
        self.assertEqual(len(urls), 1)
        self.assertIn("S31_00_W071_00", urls[0])

    def test_zero_radius_is_single_tile(self):
        """
        Test that a zero radius (the default) returns exactly one tile.
        """
        urls = get_copernicus_dem_tile_urls(-111.670, 33.289)
        self.assertEqual(len(urls), 1)
        self.assertIn("N33_00_W112_00", urls[0])

    def test_positive_radius_spans_multiple_tiles(self):
        """
        Test that a large enough radius spans more than one tile.
        """
        urls = get_copernicus_dem_tile_urls(-111.670, 33.289, radius=100000)
        self.assertGreater(len(urls), 1)

    def test_urls_are_unique(self):
        """
        Test that no tile URL is duplicated in the result.
        """
        urls = get_copernicus_dem_tile_urls(-111.670, 33.289, radius=150000)
        self.assertEqual(len(urls), len(set(urls)))


class TestSampleDemElevation(unittest.TestCase):
    """
    Unit tests for sample_dem_elevation.
    """

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.dem_path = str(Path(self.tmpdir.name) / "ridge.tif")
        _write_ridge_dem(self.dem_path)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_samples_flat_terrain(self):
        """
        Test sampling a point on the flat (0 m) part of the synthetic DEM.
        """
        self.assertAlmostEqual(
            sample_dem_elevation(self.dem_path, 0.0, 0.0), 0.0, delta=1e-6
        )

    def test_samples_ridge(self):
        """
        Test sampling a point on the synthetic ridge (2000 m).
        """
        self.assertAlmostEqual(
            sample_dem_elevation(self.dem_path, 0.25, 0.0), 2000.0, delta=1e-6
        )

    def test_accepts_list_of_paths(self):
        """
        Test that a list containing the DEM path works the same as a bare
        string path.
        """
        self.assertAlmostEqual(
            sample_dem_elevation([self.dem_path], 0.25, 0.0), 2000.0, delta=1e-6
        )

    def test_raises_outside_coverage(self):
        """
        Test that a point outside every provided DEM's extent raises a
        ValueError.
        """
        with self.assertRaises(ValueError):
            sample_dem_elevation(self.dem_path, 50.0, 50.0)


class TestComputeTerrainMask(unittest.TestCase):
    """
    Unit tests for compute_terrain_mask and compute_terrain_mask_for_station.
    """

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.dem_path = str(Path(self.tmpdir.name) / "ridge.tif")
        _write_ridge_dem(self.dem_path)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_identifies_blocked_azimuth(self):
        """
        Test that the synthetic ridge (due east of the station) is
        identified as a significant obstruction near azimuth 90 degrees,
        with a plausible (few-degree) blocking angle.
        """
        mask = compute_terrain_mask(
            self.dem_path,
            0.0,
            0.0,
            station_elevation=0,
            search_radius=50000,
            number_azimuths=36,
            number_range_samples=50,
        )
        angles = dict(zip(mask.azimuth, mask.min_elevation_angle))
        blocked_azimuth = min(angles, key=lambda az: abs(az - 90))
        self.assertAlmostEqual(blocked_azimuth, 90.0, delta=1e-6)
        self.assertGreater(angles[blocked_azimuth], 2.0)
        self.assertLess(angles[blocked_azimuth], 10.0)

    def test_clear_azimuths_report_low_angle(self):
        """
        Test that azimuths far from the ridge report a low (near-zero or
        negative) blocking angle, well below the ridge's azimuth.
        """
        mask = compute_terrain_mask(
            self.dem_path,
            0.0,
            0.0,
            station_elevation=0,
            search_radius=50000,
            number_azimuths=36,
            number_range_samples=50,
        )
        angles = dict(zip(mask.azimuth, mask.min_elevation_angle))
        clear_azimuth = min(angles, key=lambda az: abs(az - 270))
        self.assertLess(angles[clear_azimuth], 1.0)

    def test_result_is_valid_terrain_mask(self):
        """
        Test that the result is a usable TerrainMask with matching-length
        azimuth/min_elevation_angle lists.
        """
        mask = compute_terrain_mask(
            self.dem_path, 0.0, 0.0, station_elevation=0, number_azimuths=8
        )
        self.assertEqual(len(mask.azimuth), 8)
        self.assertEqual(len(mask.min_elevation_angle), 8)

    def test_raises_when_no_dem_overlaps(self):
        """
        Test that a station far outside the DEM's coverage raises a
        ValueError rather than silently returning an empty/meaningless mask.
        """
        with self.assertRaises(ValueError):
            compute_terrain_mask(
                self.dem_path, 50.0, 50.0, station_elevation=0, search_radius=10000
            )

    def test_compute_terrain_mask_for_station_matches_direct_call(self):
        """
        Test that the RadarStation convenience wrapper produces the same
        result as calling compute_terrain_mask directly with the station's
        own longitude/latitude/elevation.
        """
        station = RadarStation(name="test", latitude=0.0, longitude=0.0, elevation=0)
        direct = compute_terrain_mask(
            self.dem_path,
            station.longitude,
            station.latitude,
            station.elevation,
            search_radius=50000,
            number_azimuths=36,
            number_range_samples=50,
        )
        wrapped = compute_terrain_mask_for_station(
            self.dem_path,
            station,
            search_radius=50000,
            number_azimuths=36,
            number_range_samples=50,
        )
        self.assertEqual(direct.azimuth, wrapped.azimuth)
        self.assertEqual(direct.min_elevation_angle, wrapped.min_elevation_angle)

    def test_terrain_mask_usable_by_radar_station(self):
        """
        Test that the resulting TerrainMask can be attached to a
        RadarStation and used to compute a (smaller, terrain-shaped)
        footprint relative to the unmasked symmetric case.
        """
        mask = compute_terrain_mask(
            self.dem_path,
            0.0,
            0.0,
            station_elevation=0,
            search_radius=50000,
            number_azimuths=36,
            number_range_samples=50,
        )
        masked_station = RadarStation(
            name="test", latitude=0.0, longitude=0.0, elevation=0, terrain_mask=mask
        )
        unmasked_station = RadarStation(
            name="test", latitude=0.0, longitude=0.0, elevation=0
        )
        masked_footprint = masked_station.compute_footprint(
            elevation=1000, number_points=72
        )
        unmasked_footprint = unmasked_station.compute_footprint(elevation=1000)
        self.assertLess(masked_footprint.area, unmasked_footprint.area)
