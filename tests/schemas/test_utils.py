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
    wrap_coordinates_antimeridian,
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

    def test_wrap_coordinates_antimeridian_nominal(self):
        test_data = [(-150, 0), (-50, 0), (0, 0), (50, 0), (150, 0)]
        test_result = [(210, 0), (310, 0), (0, 0), (50, 0), (150, 0)]
        self.assertEqual(wrap_coordinates_antimeridian(test_data), test_result)

    def test_wrap_coordinates_antimeridian_no_op_all_non_negative(self):
        test_data = [(0, 0), (50, 0), (150, 0)]
        self.assertEqual(wrap_coordinates_antimeridian(test_data), test_data)

    def test_wrap_coordinates_antimeridian_no_op_all_non_positive(self):
        test_data = [(-150, 0), (-50, 0), (0, 0)]
        self.assertEqual(wrap_coordinates_antimeridian(test_data), test_data)

    def test_wrap_coordinates_antimeridian_no_op_small_range(self):
        test_data = [(-50, 0), (0, 0), (50, 0)]
        self.assertEqual(wrap_coordinates_antimeridian(test_data), test_data)

    def test_normalize_geometry_polygon(self):
        geometry = [(-150, 10), (150, 10), (150, -10), (-150, -10), (-150, 10)]
        result = normalize_geometry(Polygon(geometry))
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.geometry[0], Polygon(wrap_coordinates_antimeridian(geometry))
        )

    def test_normalize_geometry_polygon_invalid(self):
        geometry = [(-150, 0), (-50, 0), (0, 0), (50, 0), (150, 0)]
        with self.assertRaises(ValueError):
            normalize_geometry(Polygon(geometry))

    def test_normalize_geometry_multipolygon(self):
        geometries = [
            [(-150, 10), (150, 10), (150, -10), (-150, -10), (-150, 10)],
            [(-150, 30), (150, 30), (150, 20), (-150, 20), (-150, 30)],
        ]
        result = normalize_geometry(
            MultiPolygon([[geometry, []] for geometry in geometries])
        )
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.geometry[0],
            MultiPolygon(
                [
                    [wrap_coordinates_antimeridian(geometry), []]
                    for geometry in geometries
                ]
            ),
        )

    def test_normalize_geometry_geoseries(self):
        geometry = [(-150, 10), (150, 10), (150, -10), (-150, -10), (-150, 10)]
        gs = gpd.GeoSeries(Polygon(geometry), crs="EPSG:4326")
        result = normalize_geometry(gs)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.geometry[0], Polygon(wrap_coordinates_antimeridian(geometry))
        )

    def test_normalize_geometry_geodataframe(self):
        geometry = [(-150, 10), (150, 10), (150, -10), (-150, -10), (-150, 10)]
        df = gpd.GeoDataFrame(geometry=[Polygon(geometry)], index=[0], crs="EPSG:4326")
        result = normalize_geometry(df)
        self.assertIsInstance(result, gpd.GeoDataFrame)
        self.assertEqual(result.crs, "EPSG:4326")
        self.assertEqual(len(result.index), 1)
        self.assertEqual(
            result.geometry[0], Polygon(wrap_coordinates_antimeridian(geometry))
        )
