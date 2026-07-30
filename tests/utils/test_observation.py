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
    swath_width_to_field_of_view,
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

    def test_swath_width_to_field_of_regard_zero_swath(self):
        """
        Test that a zero swath width (nadir-only observation) requires no
        field of regard.
        """
        self.assertEqual(swath_width_to_field_of_regard(705000, 0), 0.0)

    def test_swath_width_to_field_of_regard_elevation(self):
        """
        Test that a positive elevation (observing a point above the mean
        Earth radius) changes the field of regard required for the same
        nominal swath width and altitude.
        """
        self.assertAlmostEqual(
            swath_width_to_field_of_regard(705000, 185815, 5000),
            15.105822,
            delta=1e-5,
        )

    def test_swath_width_to_field_of_regard_inverts_field_of_regard_to_swath_width(
        self,
    ):
        """
        Test that swath_width_to_field_of_regard is the inverse of
        field_of_regard_to_swath_width across a range of altitudes and
        field of regard values (as fractions of the horizon-limited maximum
        field of regard at each altitude, so every combination stays within
        the physically achievable range).
        """
        for altitude in (500000, 705000, 800000, 20000000):
            max_field_of_regard = compute_field_of_regard(altitude, 0)
            for fraction in (0.05, 0.2, 0.5, 0.9):
                field_of_regard = max_field_of_regard * fraction
                swath_width = field_of_regard_to_swath_width(altitude, field_of_regard)
                self.assertAlmostEqual(
                    swath_width_to_field_of_regard(altitude, swath_width),
                    field_of_regard,
                    delta=1e-6,
                )

    def test_swath_width_to_field_of_view_matches_field_of_regard_at_nadir(self):
        """
        Test that, at a look angle of zero (swath centered on nadir), the
        field of view matches the symmetric field of regard exactly.
        """
        self.assertAlmostEqual(
            swath_width_to_field_of_view(705000, 185815, 0),
            swath_width_to_field_of_regard(705000, 185815),
            delta=1e-9,
        )

    def test_swath_width_to_field_of_view_zero_swath(self):
        """
        Test that a zero swath width requires no field of view, regardless
        of look angle.
        """
        self.assertEqual(swath_width_to_field_of_view(705000, 0, 20), 0.0)

    def test_swath_width_to_field_of_view_decreases_with_look_angle(self):
        """
        Test that, for a fixed altitude and swath width, the field of view
        decreases as the off-nadir look angle increases (a fixed ground
        swath subtends a smaller angle as viewing becomes more oblique).
        """
        field_of_views = [
            swath_width_to_field_of_view(705000, 185815, look_angle)
            for look_angle in (0, 5, 10, 20, 30)
        ]
        self.assertEqual(field_of_views, sorted(field_of_views, reverse=True))

    def test_swath_width_to_field_of_view_off_nadir(self):
        """
        Test that the field of view can be computed for a specified
        off-nadir look angle.
        """
        self.assertAlmostEqual(
            swath_width_to_field_of_view(705000, 185815, 20),
            12.992488,
            delta=1e-5,
        )

    def test_swath_width_to_field_of_view_saturates_beyond_horizon(self):
        """
        Test that a look angle at or beyond the horizon-limited maximum
        saturates to the same field of view rather than raising an error.
        """
        max_look_angle = compute_field_of_regard(705000, 0) / 2
        self.assertAlmostEqual(
            swath_width_to_field_of_view(705000, 185815, max_look_angle),
            swath_width_to_field_of_view(705000, 185815, max_look_angle + 10),
            delta=1e-12,
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
        Test that the field of regard can be computed for a given altitude
        and minimum elevation angle.
        """
        self.assertAlmostEqual(
            compute_field_of_regard(705000, 81.66446), 15.0, delta=0.001
        )

    def test_compute_min_elevation_angle(self):
        """
        Test that the minimum elevation angle can be computed for a given
        altitude and field of regard.
        """
        self.assertAlmostEqual(
            compute_min_elevation_angle(705000, 15.0), 81.66446, delta=0.001
        )

    def test_compute_min_elevation_angle_saturated(self):
        """
        Test that the minimum elevation angle is saturated at 0 degrees for
        a given altitude and field of regard.
        """
        self.assertEqual(compute_min_elevation_angle(30000000, 180.0), 0.0)

    def test_compute_max_access_time(self):
        """
        Test that the maximum access time can be computed for a given
        altitude and minimum elevation angle.
        """
        self.assertAlmostEqual(
            compute_max_access_time(705000, 81.66446), 28, delta=1
        )
