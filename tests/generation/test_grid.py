"""
Unit tests for the tatc.generation._grid module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

from shapely.geometry import Polygon

from tatc.generation._grid import (
    compute_point_id_uniform_spacing,
    generate_indices_uniform_spacing,
)


class TestComputeEquallySpacedPointId(unittest.TestCase):
    """
    Unit tests for the tatc.generation._grid.compute_point_id_uniform_spacing
    function.
    """

    def test_first_point_has_id_zero(self):
        """
        Test that the first grid point (i=0, j=0) is assigned id 0.
        """
        self.assertEqual(compute_point_id_uniform_spacing(0, 0, 10), 0)

    def test_ids_increment_west_to_east_within_a_latitude_row(self):
        """
        Test that the longitude index i increments the id by 1 within a
        fixed latitude row j, matching the documented west-to-east
        ordering.
        """
        ids = [compute_point_id_uniform_spacing(i, 0, 10) for i in range(5)]
        self.assertEqual(ids, [0, 1, 2, 3, 4])

    def test_ids_increment_by_full_row_width_south_to_north(self):
        """
        Test that the latitude index j increments the id by a full row's
        worth of longitude bins (360/theta_i), matching the documented
        south-to-north ordering.
        """
        row_width = int(360 / 10)
        self.assertEqual(compute_point_id_uniform_spacing(0, 1, 10), row_width)

    def test_longitude_index_wraps_around_row_width(self):
        """
        Test that a longitude index at the row width wraps back to the
        start of the row (np.mod behavior) rather than overflowing into
        the next latitude row's id range.
        """
        row_width = int(360 / 10)
        self.assertEqual(compute_point_id_uniform_spacing(row_width, 0, 10), 0)

    def test_row_width_scales_with_theta_i(self):
        """
        Test the id formula with a longitude step that doesn't evenly
        divide into the default 10-degree examples above, confirming the
        row width (360/theta_i) scales correctly.
        """
        # 360/20 = 18 longitude bins per row
        self.assertEqual(compute_point_id_uniform_spacing(0, 1, 20), 18)
        self.assertEqual(compute_point_id_uniform_spacing(5, 1, 20), 23)

    def test_ids_are_unique_for_unequal_longitude_and_latitude_steps(self):
        """
        Regression test: ensure no two distinct grid points collide on the
        same id when the grid's longitude and latitude angular steps
        differ (compute_point_id_uniform_spacing only takes the longitude
        step; a prior bug used the latitude step for the row multiplier
        instead, causing collisions in exactly this scenario).
        """
        theta_longitude, theta_latitude = 10, 20
        indices = generate_indices_uniform_spacing(theta_longitude, theta_latitude)
        ids = [
            compute_point_id_uniform_spacing(i, j, theta_longitude)
            for (i, j) in indices
        ]
        self.assertEqual(len(ids), len(set(ids)))


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
        self.assertEqual(len(indices), len(set(indices)))


if __name__ == "__main__":
    unittest.main()
