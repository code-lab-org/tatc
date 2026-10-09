"""
Unit tests for the tatc.generation._grid module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from shapely.geometry import Polygon

from tatc.generation._grid import generate_indices_uniform_spacing


class TestGenerateEquallySpacedIndices(unittest.TestCase):
    """
    Unit tests for the tatc.generation._grid.generate_indices_uniform_spacing
    function.
    """

    def test_global_grid_index_count(self):
        """
        Test that a global (no mask) grid produces the expected number of
        two-dimensional (longitude, latitude) index pairs.
        """
        indices = generate_indices_uniform_spacing(10, 10)
        self.assertEqual(len(indices), (360 / 10) * (180 / 10))

    def test_latitude_strips_produce_one_dimensional_indices(self):
        """
        Test that strips="lat" produces one index per latitude row, each
        with longitude index fixed at 0.
        """
        indices = generate_indices_uniform_spacing(10, 10, strips="lat")
        self.assertEqual(len(indices), 180 / 10)
        self.assertTrue(all(i == 0 for i, j in indices))

    def test_longitude_strips_produce_one_dimensional_indices(self):
        """
        Test that strips="lon" produces one index per longitude column,
        each with latitude index fixed at 0.
        """
        indices = generate_indices_uniform_spacing(10, 10, strips="lon")
        self.assertEqual(len(indices), 360 / 10)
        self.assertTrue(all(j == 0 for i, j in indices))

    def test_mask_restricts_indices_to_its_bounds(self):
        """
        Test that supplying a mask restricts generated indices to the
        mask's bounding box rather than the full globe. Uses bounds that
        land on exact grid-index boundaries (rather than a half-step
        offset) to avoid np.round's round-half-to-even tie-breaking.
        """
        mask = Polygon([[-100, 20], [-50, 20], [-50, -20], [-100, -20], [-100, 20]])
        indices = generate_indices_uniform_spacing(10, 10, mask=mask)
        expected_count = ((-50) - (-100)) / 10 * (20 - (-20)) / 10
        self.assertEqual(len(indices), expected_count)

    def test_indices_are_unique(self):
        """
        Test that a global grid produces no duplicate (i, j) index pairs.
        """
        indices = generate_indices_uniform_spacing(10, 10)
        self.assertEqual(len(indices), len(set(map(tuple, indices))))


if __name__ == "__main__":
    unittest.main()
