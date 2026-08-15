"""
Unit tests for the tatc.utils.atmosphere module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np

from tatc.utils import altitude_to_pressure, pressure_to_altitude


class TestAtmosphere(unittest.TestCase):
    """
    Unit tests for the tatc.utils.atmosphere module.
    """

    def test_altitude_to_pressure_at_sea_level(self):
        """
        Test that sea level (0 m) maps to standard atmospheric pressure.
        """
        self.assertAlmostEqual(altitude_to_pressure(0.0), 101325.0, delta=0.01)

    def test_altitude_to_pressure_decreases_with_altitude(self):
        """
        Test that pressure decreases monotonically with increasing
        altitude across the full valid range of the model.
        """
        altitudes = np.linspace(0, 86000, 200)
        pressures = [altitude_to_pressure(h) for h in altitudes]
        self.assertTrue(all(p1 > p2 for p1, p2 in zip(pressures, pressures[1:])))

    def test_altitude_to_pressure_matches_known_reference_values(self):
        """
        Test that pressure at several published U.S. Standard Atmosphere
        layer boundaries matches commonly cited reference values.
        """
        # (altitude in m, expected pressure in Pa, tolerance in Pa)
        references = [
            (11000.0, 22632.1, 1.0),  # tropopause
            (20000.0, 5474.9, 1.0),
            (32000.0, 868.0, 0.5),
        ]
        for altitude, expected, tolerance in references:
            with self.subTest(altitude=altitude):
                self.assertAlmostEqual(
                    altitude_to_pressure(altitude), expected, delta=tolerance
                )

    def test_altitude_to_pressure_saturates_below_sea_level(self):
        """
        Test that altitudes below sea level saturate to the sea-level
        pressure rather than extrapolating.
        """
        self.assertEqual(altitude_to_pressure(-5000.0), altitude_to_pressure(0.0))

    def test_altitude_to_pressure_saturates_above_model_top(self):
        """
        Test that altitudes above the model's 86 km upper bound saturate
        to the pressure at that bound rather than extrapolating.
        """
        self.assertEqual(
            altitude_to_pressure(200000.0), altitude_to_pressure(86000.0)
        )

    def test_pressure_to_altitude_saturates_above_sea_level_pressure(self):
        """
        Test that pressures above the sea-level value saturate to sea
        level (0 m) rather than extrapolating to a negative altitude.
        """
        self.assertEqual(pressure_to_altitude(150000.0), pressure_to_altitude(101325.0))
        self.assertAlmostEqual(pressure_to_altitude(101325.0), 0.0, delta=0.01)

    def test_pressure_to_altitude_saturates_below_model_bottom_pressure(self):
        """
        Test that pressures below the model's minimum (the pressure at
        86 km) saturate to 86 km rather than extrapolating.
        """
        min_pressure = altitude_to_pressure(86000.0)
        self.assertEqual(
            pressure_to_altitude(min_pressure / 100), pressure_to_altitude(min_pressure)
        )
        self.assertAlmostEqual(pressure_to_altitude(min_pressure), 86000.0, delta=0.01)

    def test_pressure_to_altitude_inverts_altitude_to_pressure(self):
        """
        Test that pressure_to_altitude exactly inverts altitude_to_pressure
        across the full valid range, including every layer boundary.
        """
        altitudes = list(np.linspace(0, 86000, 50)) + [
            0.0,
            11000.0,
            20000.0,
            32000.0,
            47000.0,
            51000.0,
            71000.0,
            86000.0,
        ]
        for altitude in altitudes:
            with self.subTest(altitude=altitude):
                pressure = altitude_to_pressure(altitude)
                self.assertAlmostEqual(
                    pressure_to_altitude(pressure), altitude, delta=1e-6
                )

    def test_pressure_to_altitude_matches_common_science_reference_levels(self):
        """
        Test that standard atmospheric pressure levels (in hPa, converted
        to Pa) used throughout atmospheric science map to their commonly
        cited approximate altitudes.
        """
        # (pressure level in hPa, expected altitude in km, tolerance in km)
        references = [
            (100, 16.2, 0.5),  # near the tropopause
            (10, 31.1, 0.5),  # stratosphere
            (1, 47.8, 0.5),  # stratopause
        ]
        for hpa, expected_km, tolerance_km in references:
            with self.subTest(hpa=hpa):
                altitude_km = pressure_to_altitude(hpa * 100) / 1e3
                self.assertAlmostEqual(altitude_km, expected_km, delta=tolerance_km)
