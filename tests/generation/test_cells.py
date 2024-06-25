import unittest

from tatc.generation import generate_equally_spaced_cells, _generate_equally_spaced_cells
from shapely.geometry import Polygon


class TestCellGenerators(unittest.TestCase):
    def test_generate_equally_spaced_cells_no_mask_count(self):
        points = _generate_equally_spaced_cells(10, 10)
        num_samples = (360 / 10) * (180 / 10)
        self.assertEqual(len(points), num_samples)

    def test_generate_equally_spaced_cells_mask_contains(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-50, -25], [-100, -25], [-100, 25]])
        points = generate_equally_spaced_cells(distance, mask=mask)
        points.apply(lambda r: self.assertTrue(mask.contains(r.geometry)), axis=1)

    def test_generate_equally_spaced_cells_mask_small(self):
        distance = 2000000
        mask = Polygon([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
        points = generate_equally_spaced_cells(distance, mask=mask)
        self.assertEqual(len(points), 0)

    def test_generate_equally_spaced_cells_mask_invalid(self):
        distance = 2000000
        mask = Polygon([[-100, 25], [-50, 25], [-100, -25], [-50, -25], [-100, 25]])
        with self.assertRaises(ValueError):
            generate_equally_spaced_cells(distance, mask=mask)

    def test_generate_equally_spaced_cells_no_mask_lat(self):
        points = _generate_equally_spaced_cells(10, 10, strips="lat")
        num_samples = 180 / 10
        self.assertEqual(len(points), num_samples)

    def test_generate_equally_spaced_cells_no_mask_lon(self):
        points = _generate_equally_spaced_cells(10, 10, strips="lon")
        num_samples = 360 / 10
        self.assertEqual(len(points), num_samples)
