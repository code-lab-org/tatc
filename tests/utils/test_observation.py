"""
Unit tests for the tatc.utils.observation module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

from tatc.utils import (
    compute_field_of_regard,
    compute_max_access_time,
    compute_min_elevation_angle,
    field_of_regard_to_swath_width,
    swath_width_to_field_of_regard,
)


class TestObservation(unittest.TestCase):
    """
    Unit tests for the tatc.utils.observation module.
    """
    def test_swath_width_to_field_of_regard(self):
        """
        Test that the swath width can be converted to field of regard.
        """
        self.assertAlmostEqual(
            swath_width_to_field_of_regard(705000, 185815), 15.0, delta=0.001
        )

    def test_field_of_regard_to_swath_width(self):
        """
        Test that the field of regard can be converted to swath width.
        """
        self.assertAlmostEqual(
            field_of_regard_to_swath_width(705000, 15.0), 185815, delta=1.0
        )

    def test_compute_field_of_regard(self):
        """
        Test that the field of regard can be computed for a given altitude and minimum elevation angle.
        """
        self.assertAlmostEqual(
            compute_field_of_regard(705000, 81.66446), 15.0, delta=0.001
        )

    def test_compute_min_elevation_angle(self):
        """
        Test that the minimum elevation angle can be computed for a given altitude and field of regard.
        """
        self.assertAlmostEqual(
            compute_min_elevation_angle(705000, 15.0), 81.66446, delta=0.001
        )

    def test_compute_min_elevation_angle_saturated(self):
        """
        Test that the minimum elevation angle is saturated at 0 degrees for a given altitude and field of regard.
        """
        self.assertEqual(compute_min_elevation_angle(30000000, 180.0), 0.0)

    def test_compute_max_access_time(self):
        """
        Test that the maximum access time can be computed for a given altitude and minimum elevation angle.
        """
        self.assertAlmostEqual(
            compute_max_access_time(705000, 81.66446), 28, delta=1
        )
