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

import geopandas as gpd
import numpy as np
import rasterio
import rasterio.windows
from rasterio.transform import from_origin
from shapely.geometry import Point, box

from tatc import constants
from tatc.preprocess import (
    add_local_incidence_angles,
    compute_terrain_mask,
    compute_terrain_mask_for_platform,
    compute_terrain_mask_for_station,
    compute_terrain_slope_aspect,
    get_copernicus_dem_tile_urls,
    sample_dem_elevation,
)
from tatc.schemas import RadarStation, SurfacePlatform


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

    def test_linear_and_geometric_range_spacing_find_ridge(self):
        """
        Test that both range spacings identify the ridge, with the same
        blocking angle to within the sampling of its crest.
        """
        masks = {
            range_spacing: compute_terrain_mask(
                self.dem_path,
                0.0,
                0.0,
                station_elevation=0,
                search_radius=50000,
                number_azimuths=36,
                range_spacing=range_spacing,
            )
            for range_spacing in ("linear", "geometric")
        }
        angles = {
            key: dict(zip(m.azimuth, m.min_elevation_angle)) for key, m in masks.items()
        }
        self.assertGreater(angles["geometric"][90.0], 2.0)
        self.assertAlmostEqual(
            angles["linear"][90.0], angles["geometric"][90.0], delta=0.5
        )

    def test_geometric_range_spacing_samples_densely_near_station(self):
        """
        Test that geometric range spacing captures a narrow obstruction
        close to the station that sparse, evenly spaced samples skip.
        """
        resolution = 0.0001  # ~11 m/pixel
        data = np.zeros((400, 400), dtype=np.float32)
        # a 50 m tall wall due east, from about 1.40 to 1.48 km
        data[195:205, 326:333] = 50.0
        path = str(Path(self.tmpdir.name) / "wall.tif")
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            height=400,
            width=400,
            count=1,
            dtype="float32",
            crs="EPSG:4326",
            transform=from_origin(-0.02, 0.02, resolution, resolution),
        ) as dst:
            dst.write(data, 1)
        # evenly spaced samples (every 192 m) fall at 1.38 and 1.58 km,
        # on either side of the wall; geometric samples are ~43 m apart there
        kwargs = {
            "search_radius": 20000,
            "number_azimuths": 4,
            "number_range_samples": 100,
        }
        linear = compute_terrain_mask(
            path, 0.0, 0.0, 0, range_spacing="linear", **kwargs
        )
        geometric = compute_terrain_mask(path, 0.0, 0.0, 0, **kwargs)
        self.assertLess(linear.min_elevation_angle[1], 1.0)
        self.assertGreater(geometric.min_elevation_angle[1], 1.0)

    def test_bad_range_spacing(self):
        """
        Test that an unknown range spacing, or geometric spacing starting
        at zero range, raises a ValueError.
        """
        with self.assertRaises(ValueError):
            compute_terrain_mask(self.dem_path, 0.0, 0.0, 0, range_spacing="log")
        with self.assertRaises(ValueError):
            compute_terrain_mask(self.dem_path, 0.0, 0.0, 0, min_range=0)

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

    def test_compute_terrain_mask_for_platform_matches_direct_call(self):
        """
        Test that the SurfacePlatform convenience wrapper places the
        instruments a platform height above the DEM ground, without
        refraction by default.
        """
        platform = SurfacePlatform(
            name="test", latitude=0.0, longitude=0.0, elevation=-30, height=30
        )
        kwargs = {"search_radius": 50000, "number_azimuths": 36}
        kwargs["number_range_samples"] = 50
        direct = compute_terrain_mask(
            self.dem_path, 0.0, 0.0, 30, effective_earth_radius_factor=1, **kwargs
        )
        wrapped = compute_terrain_mask_for_platform(self.dem_path, platform, **kwargs)
        self.assertEqual(direct.azimuth, wrapped.azimuth)
        self.assertEqual(direct.min_elevation_angle, wrapped.min_elevation_angle)
        # refraction (4/3) raises the apparent terrain
        refracted = compute_terrain_mask_for_platform(
            self.dem_path, platform, effective_earth_radius_factor=4 / 3, **kwargs
        )
        self.assertTrue(
            all(
                a > b
                for a, b in zip(
                    refracted.min_elevation_angle, wrapped.min_elevation_angle
                )
            )
        )

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


def _write_plane_dem(path: str, slope: float, aspect: float) -> None:
    """
    Writes a small synthetic DEM near (0, 0) of a tilted plane with a given
    slope (degrees) facing a given aspect (degrees), with meters per degree
    of the WGS 84 ellipsoid at the equator, and a band of nodata cells
    along its northern edge.
    """
    resolution = 0.0003  # about 33 m/pixel
    size = 400
    transform = from_origin(-0.06, 0.06, resolution, resolution)
    centers = -0.06 + resolution * (np.arange(size) + 0.5)
    east = np.radians(centers) * constants.EARTH_EQUATORIAL_RADIUS
    north = (
        np.radians(-centers)
        * constants.EARTH_EQUATORIAL_RADIUS
        * (1 - constants.EARTH_ECCENTRICITY**2)
    )
    east_grid, north_grid = np.meshgrid(east, north)
    # elevation rises opposite the downslope (aspect) direction
    rise = -np.tan(np.radians(slope)) * (
        east_grid * np.sin(np.radians(aspect)) + north_grid * np.cos(np.radians(aspect))
    )
    data = (1000 + rise).astype(np.float32)
    data[:20, :] = -9999
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=size,
        width=size,
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=transform,
        nodata=-9999,
    ) as dst:
        dst.write(data, 1)


class TestTerrainSlopeAspect(unittest.TestCase):
    """
    Unit tests for compute_terrain_slope_aspect and add_local_incidence_angles.
    """

    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.west_path = str(Path(self.tmp_dir.name) / "west.tif")
        _write_plane_dem(self.west_path, 20, 270)
        self.northeast_path = str(Path(self.tmp_dir.name) / "northeast.tif")
        _write_plane_dem(self.northeast_path, 10, 45)
        self.flat_path = str(Path(self.tmp_dir.name) / "flat.tif")
        _write_plane_dem(self.flat_path, 0, 0)

    def tearDown(self):
        self.tmp_dir.cleanup()

    def test_slope_and_aspect_of_planes(self):
        """
        Test the slope and aspect of tilted planes at several points,
        including one between cell centers.
        """
        longitudes = np.array([0.0, 0.01, -0.023456])
        latitudes = np.array([0.0, -0.02, 0.011111])
        for path, slope, aspect in [
            (self.west_path, 20, 270),
            (self.northeast_path, 10, 45),
        ]:
            with self.subTest(aspect=aspect):
                slopes, aspects = compute_terrain_slope_aspect(
                    path, longitudes, latitudes
                )
                np.testing.assert_allclose(slopes, slope, atol=0.05)
                np.testing.assert_allclose(aspects, aspect, atol=0.1)

    def test_slope_across_tile_seam(self):
        """
        Test that points whose samples fall between the cells of two
        adjacent tiles are interpolated across the seam.
        """
        paths = []
        with rasterio.open(self.northeast_path) as source:
            for name, rows in [("north", (0, 200)), ("south", (200, 400))]:
                window = rasterio.windows.Window(0, rows[0], source.width, 200)
                profile = source.profile.copy()
                profile.update(height=200, transform=source.window_transform(window))
                paths.append(str(Path(self.tmp_dir.name) / f"{name}.tif"))
                with rasterio.open(paths[-1], "w", **profile) as dst:
                    dst.write(source.read(1, window=window), 1)
        # the seam is at latitude 0; samples lie 30 m north and south
        longitudes = np.array([0.0, 0.0001, -0.0002])
        latitudes = np.array([0.0, 0.0001, -0.00005])
        slopes, aspects = compute_terrain_slope_aspect(paths, longitudes, latitudes)
        np.testing.assert_allclose(slopes, 10, atol=0.05)
        np.testing.assert_allclose(aspects, 45, atol=0.1)

    def test_flat_terrain(self):
        """
        Test that flat terrain has zero slope and zero aspect, for a scalar
        point and a larger baseline.
        """
        slope, aspect = compute_terrain_slope_aspect(
            self.flat_path, 0.0, 0.0, baseline=500
        )
        self.assertEqual(slope.shape, ())
        self.assertAlmostEqual(float(slope), 0)
        self.assertAlmostEqual(float(aspect), 0)

    def test_nodata_and_outside_are_nan(self):
        """
        Test that points with nodata or outside the DEM have NaN slopes.
        """
        slopes, aspects = compute_terrain_slope_aspect(
            self.west_path, np.array([0.0, 0.0, 0.0599]), np.array([0.0, 0.056, 0.0])
        )
        self.assertFalse(np.isnan(slopes[0]))
        self.assertTrue(np.isnan(slopes[1]))
        self.assertTrue(np.isnan(slopes[2]))
        self.assertTrue(np.isnan(aspects[1]))

    def test_invalid_inputs(self):
        """
        Test that a non-positive baseline or DEMs that do not overlap the
        points raise a ValueError.
        """
        with self.assertRaises(ValueError):
            compute_terrain_slope_aspect(self.west_path, 0.0, 0.0, baseline=0)
        with self.assertRaises(ValueError):
            compute_terrain_slope_aspect(self.west_path, 10.0, 10.0)

    def test_add_local_incidence_angles(self):
        """
        Test the local incidence angles of a satellite to the west of a
        west-facing slope, and on flat terrain.
        """
        observations = gpd.GeoDataFrame(
            {"sat_alt": [40.0, 40.0], "sat_az": [270.0, 90.0]},
            geometry=[Point(0, 0), Point(0.01, 0.01)],
            crs="EPSG:4326",
        )
        result = add_local_incidence_angles(observations, self.west_path)
        self.assertNotIn("local_incidence_angle", observations.columns)
        np.testing.assert_allclose(result.terrain_slope, 20, atol=0.05)
        np.testing.assert_allclose(result.local_incidence_angle, [30, 70], atol=0.1)
        flat = add_local_incidence_angles(observations, self.flat_path)
        np.testing.assert_allclose(flat.local_incidence_angle, 50, atol=1e-6)

    def test_add_local_incidence_angles_invalid(self):
        """
        Test that observations without satellite angles or with geometries
        other than points raise a ValueError, and that empty observations
        gain the columns.
        """
        with self.assertRaises(ValueError):
            add_local_incidence_angles(
                gpd.GeoDataFrame(geometry=[Point(0, 0)], crs="EPSG:4326"),
                self.west_path,
            )
        with self.assertRaises(ValueError):
            add_local_incidence_angles(
                gpd.GeoDataFrame(
                    {"sat_alt": [40.0], "sat_az": [0.0]},
                    geometry=[box(0, 0, 0.01, 0.01)],
                    crs="EPSG:4326",
                ),
                self.west_path,
            )
        empty = add_local_incidence_angles(
            gpd.GeoDataFrame(
                {"sat_alt": [], "sat_az": []}, geometry=[], crs="EPSG:4326"
            ),
            self.west_path,
        )
        self.assertIn("local_incidence_angle", empty.columns)
        self.assertEqual(len(empty.index), 0)
