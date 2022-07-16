import unittest

import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import geopandas as gpd

from tatc.utils import (
    mean_anomaly_to_true_anomaly,
    true_anomaly_to_mean_anomaly,
    compute_number_samples,
    swath_width_to_field_of_regard,
    field_of_regard_to_swath_width,
    compute_field_of_regard,
    compute_min_elevation_angle,
    compute_orbit_period,
    compute_max_access_time,
    split_polygon,
    normalize_geometry,
)

from tatc import constants


class TestUtils(unittest.TestCase):
    def test_mean_anomaly_to_true_anomaly(self):
        self.assertAlmostEqual(
            mean_anomaly_to_true_anomaly(78.940629, 0.0001492), 78.95065818, delta=0.01
        )

    def test_true_anomaly_to_mean_anomaly(self):
        self.assertAlmostEqual(
            true_anomaly_to_mean_anomaly(78.95065818, 0.0001492), 78.940629, delta=0.01
        )

    def test_compute_number_samples(self):
        # rough approximation based on flat sample areas
        sample_distance = 10000
        num_samples = int(
            constants.earth_surface_area / (np.pi * (sample_distance / 2) ** 2)
        )
        self.assertEqual(compute_number_samples(sample_distance), num_samples)

    def test_swath_width_to_field_of_regard(self):
        self.assertAlmostEqual(
            swath_width_to_field_of_regard(705000, 185815), 15.0, delta=0.001
        )

    def test_field_of_regard_to_swath_width(self):
        self.assertAlmostEqual(
            field_of_regard_to_swath_width(705000, 15.0), 185815, delta=1.0
        )

    def test_compute_field_of_regard(self):
        self.assertAlmostEqual(
            compute_field_of_regard(705000, 81.66446), 15.0, delta=0.001
        )

    def test_compute_min_elevation_angle(self):
        self.assertAlmostEqual(
            compute_min_elevation_angle(705000, 15.0), 81.66446, delta=0.001
        )

    def test_compute_min_elevation_angle_saturated(self):
        self.assertEqual(compute_min_elevation_angle(30000000, 180.0), 0.0)

    def test_compute_max_access_time(self):
        self.assertAlmostEqual(
            compute_max_access_time(705000, 81.66446), 27.33097, delta=0.001
        )

    def test_split_polygon_nominal(self):
        polygon = Polygon([(170, 10), (-170, 10), (-170, -10), (170, -10), (170, 10)])
        result = MultiPolygon(
            [
                Polygon([(170, -10), (170, 10), (180, 10), (180, -10), (170, -10)]),
                Polygon(
                    [(-180, -10), (-180, 10), (-170, 10), (-170, -10), (-180, -10)]
                ),
            ]
        )
        self.assertEqual(split_polygon(polygon), result)

    def test_normalize_geometry_polygon(self):
        polygon = Polygon([(-150, 10), (150, 10), (150, -10), (-150, -10), (-150, 10)])
        result = normalize_geometry(polygon)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(result.iloc[0].geometry, split_polygon(polygon))

    def test_normalize_geometry_polygon_invalid(self):
        polygon = Polygon([(-150, 0), (-50, 0), (0, 0), (50, 0), (150, 0)])
        with self.assertRaises(ValueError):
            normalize_geometry(polygon)

    def test_normalize_geometry_multipolygon(self):
        geometries = [
            [(-150, 10), (150, 10), (150, -10), (-150, -10), (-150, 10)],
            [(-150, 30), (150, 30), (150, 20), (-150, 20), (-150, 30)],
        ]
        multipolygon = MultiPolygon([[geometry, []] for geometry in geometries])
        result = normalize_geometry(multipolygon)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.iloc[0].geometry,
            split_polygon(multipolygon),
        )

    def test_normalize_geometry_geoseries(self):
        polygon = Polygon([(-150, 10), (150, 10), (150, -10), (-150, -10), (-150, 10)])
        gs = gpd.GeoSeries(polygon, crs="EPSG:4326")
        result = normalize_geometry(gs)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.iloc[0].geometry,
            split_polygon(polygon),
        )

    def test_normalize_geometry_geodataframe(self):
        polygon = Polygon([(-150, 10), (150, 10), (150, -10), (-150, -10), (-150, 10)])
        df = gpd.GeoDataFrame(geometry=[polygon], index=[0], crs="EPSG:4326")
        result = normalize_geometry(df)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.iloc[0].geometry,
            split_polygon(polygon),
        )
