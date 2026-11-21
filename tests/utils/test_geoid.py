"""
Unit tests for the geoid utilities.

The reference values are PROJ's transformation of WGS 84 heights to EGM2008
heights (EPSG:4979 to EPSG:4326+3855) with the 2.5-minute EGM2008 grid
(`us_nga_egm08_25.tif`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np

from tatc.utils import compute_geoid_height

# longitude, latitude (degrees), and the geoid height (meters)
REFERENCES = [
    (-112.074, 33.448, -30.268),  # Phoenix
    (0.0, 0.0, 17.225),
    (78.0, 5.0, -104.557),  # the low of the Indian Ocean
    (147.0, -6.0, 73.022),  # the high of New Guinea
    (179.99, -17.0, 51.313),
    (-179.99, -17.0, 51.251),
    (0.0, 90.0, 14.899),
    (30.0, -90.0, -30.15),
]


class TestGeoid(unittest.TestCase):
    """
    Unit tests for the geoid utilities.
    """

    def test_reference_values(self):
        """
        Tests the geoid heights at points against PROJ's.
        """
        for longitude, latitude, height in REFERENCES:
            with self.subTest(longitude=longitude, latitude=latitude):
                value = compute_geoid_height(longitude, latitude)
                self.assertIsInstance(value, float)
                self.assertAlmostEqual(value, height, delta=0.05)

    def test_arrays(self):
        """
        Tests that coordinates are broadcast and that longitudes wrap.
        """
        longitude, latitude, height = np.array(REFERENCES).T
        np.testing.assert_allclose(
            compute_geoid_height(longitude, latitude), height, atol=0.05
        )
        np.testing.assert_allclose(
            compute_geoid_height(longitude + 360, latitude),
            compute_geoid_height(longitude, latitude),
            atol=1e-9,
        )
        self.assertEqual(compute_geoid_height([0.0, 10.0], 90.0).shape, (2,))
        self.assertAlmostEqual(
            compute_geoid_height(180.0, -17.0), compute_geoid_height(-180.0, -17.0)
        )

    def test_latitude_range(self):
        """
        Tests that latitudes beyond the poles are not accepted.
        """
        with self.assertRaises(ValueError):
            compute_geoid_height(0.0, 91.0)


if __name__ == "__main__":
    unittest.main()
