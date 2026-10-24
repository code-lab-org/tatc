"""
Unit tests for the TerrainMask schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from pydantic import ValidationError

from tatc.schemas import TerrainMask


class TestTerrainMask(unittest.TestCase):
    """
    Unit tests for the TerrainMask schema.
    """

    def test_good_data(self):
        """
        Test that the TerrainMask schema correctly initializes with valid data.
        """
        mask = TerrainMask(
            azimuth=[0, 90, 180, 270], min_elevation_angle=[0.5, 5, 0.5, 2]
        )
        self.assertEqual(mask.azimuth, [0, 90, 180, 270])
        self.assertEqual(mask.min_elevation_angle, [0.5, 5, 0.5, 2])

    def test_bad_mismatched_lengths(self):
        """
        Test that mismatched azimuth/min_elevation_angle lengths are rejected.
        """
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[0, 90, 180], min_elevation_angle=[0.5, 5])

    def test_bad_azimuth_out_of_range(self):
        """
        Test that azimuth values outside [0, 360) are rejected.
        """
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[-1, 90], min_elevation_angle=[0.5, 5])
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[0, 360], min_elevation_angle=[0.5, 5])

    def test_bad_azimuth_not_increasing(self):
        """
        Test that non-strictly-increasing azimuth values are rejected.
        """
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[90, 0], min_elevation_angle=[0.5, 5])
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[0, 0], min_elevation_angle=[0.5, 5])

    def test_bad_too_few_samples(self):
        """
        Test that a single-sample mask is rejected (at least two required).
        """
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[0], min_elevation_angle=[0.5])

    def test_get_min_elevation_angle_at_samples(self):
        """
        Test that the interpolated value exactly at a sample azimuth
        matches that sample.
        """
        mask = TerrainMask(
            azimuth=[0, 90, 180, 270], min_elevation_angle=[0.5, 5, 0.5, 2]
        )
        for azimuth, expected in zip([0, 90, 180, 270], [0.5, 5, 0.5, 2]):
            self.assertAlmostEqual(
                mask.get_min_elevation_angle(azimuth), expected, delta=1e-9
            )

    def test_get_min_elevation_angle_interpolates_between_samples(self):
        """
        Test that the interpolated value midway between two samples is the
        linear average.
        """
        mask = TerrainMask(azimuth=[0, 90], min_elevation_angle=[0, 10])
        self.assertAlmostEqual(mask.get_min_elevation_angle(45), 5, delta=1e-9)

    def test_get_min_elevation_angle_wraps_around_zero(self):
        """
        Test that interpolation wraps around the 0/360 degree boundary
        (e.g. between the last sample and the first).
        """
        mask = TerrainMask(azimuth=[0, 270], min_elevation_angle=[0, 10])
        # halfway from 270 to 360 (0), wrapping, should be ~5
        self.assertAlmostEqual(mask.get_min_elevation_angle(315), 5, delta=1e-9)

    def test_get_min_elevation_angle_accepts_azimuth_outside_0_360(self):
        """
        Test that an azimuth outside [0, 360) is wrapped before lookup.
        """
        mask = TerrainMask(
            azimuth=[0, 90, 180, 270], min_elevation_angle=[0.5, 5, 0.5, 2]
        )
        self.assertAlmostEqual(
            mask.get_min_elevation_angle(90),
            mask.get_min_elevation_angle(450),
            delta=1e-9,
        )
        self.assertAlmostEqual(
            mask.get_min_elevation_angle(0),
            mask.get_min_elevation_angle(-360),
            delta=1e-9,
        )


if __name__ == "__main__":
    unittest.main()
