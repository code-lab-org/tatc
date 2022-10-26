import unittest

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import transform
import pyproj

from tatc.generation import (
    generate_fibonacci_lattice_points,
    generate_cubed_sphere_points,
    _generate_cubed_sphere_points,
)
from tatc import constants, utils


class TestPointGenerators(unittest.TestCase):
    def test_generate_fibonacci_lattice_points_no_mask_count(self):
        distance = 2000000
        points = generate_fibonacci_lattice_points(distance)
        num_samples = utils.compute_number_samples(distance)
        self.assertEqual(len(points), num_samples)

    def test_generate_cubed_sphere_points_no_mask_count(self):
        points = _generate_cubed_sphere_points(10, 10)
        num_samples = (360 / 10) * (180 / 10)
        self.assertEqual(len(points), num_samples)

    def test_generate_fibonacci_lattice_points_no_mask_distance(self):
        distance = 2000000
        points = generate_fibonacci_lattice_points(distance)
        num_samples = utils.compute_number_samples(distance)
        test_points_m = points.to_crs("EPSG:32663")
        test_points_m["closest_neighbor"] = test_points_m.apply(
            lambda r: test_points_m[test_points_m.point_id != r.point_id]
            .distance(r.geometry)
            .min(),
            axis=1,
        )
        # mean should be within 10% of target sample distance
        self.assertAlmostEqual(
            test_points_m.closest_neighbor.mean(),
            distance,
            delta=0.10 * distance,
        )

    def test_generate_fibonacci_lattice_points_mask_count(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-50, -25], [-100, -25], [-100, 25]])
        points = generate_fibonacci_lattice_points(distance, mask=mask)
        num_samples = utils.compute_number_samples(distance)
        transformer = pyproj.Transformer.from_crs(
            pyproj.CRS("EPSG:4326"), pyproj.CRS("EPSG:32663"), always_xy=True
        ).transform
        mask_m = transform(transformer, mask)
        self.assertAlmostEqual(
            len(points),
            int(num_samples * mask_m.area / constants.EARTH_SURFACE_AREA),
            delta=1,
        )

    def test_generate_fibonacci_lattice_points_mask_contains(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-50, -25], [-100, -25], [-100, 25]])
        points = generate_fibonacci_lattice_points(distance, mask=mask)
        points.apply(lambda r: self.assertTrue(mask.contains(r.geometry)), axis=1)

    def test_generate_fibonacci_lattice_points_mask_distance(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-50, -25], [-100, -25], [-100, 25]])
        points = generate_fibonacci_lattice_points(distance, mask=mask)
        transformer = pyproj.Transformer.from_crs(
            pyproj.CRS("EPSG:4326"), pyproj.CRS("EPSG:32663"), always_xy=True
        ).transform
        mask_m = transform(transformer, mask)
        test_points_m = points.to_crs("EPSG:32663")
        test_points_m["closest_neighbor"] = test_points_m.apply(
            lambda r: test_points_m[test_points_m.point_id != r.point_id]
            .distance(r.geometry)
            .min(),
            axis=1,
        )
        # mean should be within 20% of target sample distance
        self.assertAlmostEqual(
            test_points_m.closest_neighbor.mean(),
            distance,
            delta=0.20 * distance,
        )

    def test_generate_cubed_sphere_points_no_mask_distance(self):
        distance = 2000000
        points = generate_cubed_sphere_points(distance)
        num_samples = utils.compute_number_samples(distance)
        test_points_m = points.to_crs("EPSG:32663")
        test_points_m["closest_neighbor"] = test_points_m.apply(
            lambda r: test_points_m[test_points_m.point_id != r.point_id]
            .distance(r.geometry)
            .min(),
            axis=1,
        )
        # mean should be within 10% of target sample distance
        self.assertAlmostEqual(
            test_points_m.closest_neighbor.mean(),
            distance,
            delta=0.10 * distance,
        )

    def test_generate_cubed_sphere_points_mask_contains(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-50, -25], [-100, -25], [-100, 25]])
        points = generate_cubed_sphere_points(distance, mask=mask)
        points.apply(lambda r: self.assertTrue(mask.contains(r.geometry)), axis=1)

    def test_generate_cubed_sphere_points_mask_distance(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-50, -25], [-100, -25], [-100, 25]])
        points = generate_cubed_sphere_points(distance, mask=mask)
        transformer = pyproj.Transformer.from_crs(
            pyproj.CRS("EPSG:4326"), pyproj.CRS("EPSG:32663"), always_xy=True
        ).transform
        mask_m = transform(transformer, mask)
        test_points_m = points.to_crs("EPSG:32663")
        test_points_m["closest_neighbor"] = test_points_m.apply(
            lambda r: test_points_m[test_points_m.point_id != r.point_id]
            .distance(r.geometry)
            .min(),
            axis=1,
        )
        # mean should be within 20% of target sample distance
        self.assertAlmostEqual(
            test_points_m.closest_neighbor.mean(),
            distance,
            delta=0.20 * distance,
        )

    def test_generate_fibonacci_lattice_points_mask_small(self):
        distance = 2000000
        mask = Polygon([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
        points = generate_fibonacci_lattice_points(distance, mask=mask)
        self.assertEqual(len(points), 0)

    def test_generate_cubed_sphere_points_mask_small(self):
        distance = 2000000
        mask = Polygon([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
        points = generate_cubed_sphere_points(distance, mask=mask)
        self.assertEqual(len(points), 0)

    def test_generate_fibonacci_lattice_points_mask_invalid(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-100, -25], [-50, -25], [-100, 25]])
        with self.assertRaises(ValueError):
            generate_fibonacci_lattice_points(distance, mask=mask)

    def test_generate_cubed_sphere_points_mask_invalid(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-100, -25], [-50, -25], [-100, 25]])
        with self.assertRaises(ValueError):
            generate_cubed_sphere_points(distance, mask=mask)
