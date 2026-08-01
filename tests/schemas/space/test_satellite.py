"""
Unit tests for the Satellite schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from pydantic import ValidationError

from tatc.schemas import CircularOrbit, Instrument, Satellite


class TestSatellite(unittest.TestCase):
    """
    Unit tests for the Satellite schema.
    """
    def setUp(self):
        self.test_data = {
            "name": "Test Satellite",
            "orbit": {
                "type": "circular",
                "mean_altitude": 400000,
                "inclination": 51.6,
                "epoch": "2000-01-01T00:00:00Z",
            },
            "instruments": [{"name": "Test Instrument", "field_of_regard": 25.0}],
        }
        self.test_sat = Satellite(**self.test_data)

    def test_good_data(self):
        """
        Test that the Satellite schema correctly initializes with valid data.
        """
        self.assertEqual(self.test_sat.name, self.test_data.get("name"))
        self.assertEqual(
            self.test_sat.orbit, CircularOrbit(**self.test_data.get("orbit"))
        )
        self.assertEqual(len(self.test_sat.instruments), 1)
        self.assertEqual(
            self.test_sat.instruments[0],
            Instrument(**self.test_data.get("instruments")[0]),
        )

    def test_type_defaults_to_satellite(self):
        """
        Test that omitting type defaults to the "satellite" discriminator.
        """
        self.assertEqual(self.test_sat.type, "satellite")

    def test_type_rejects_other_values(self):
        """
        Test that type only accepts the "satellite" literal, rejecting
        other space system type discriminators (e.g. from a constellation).
        """
        with self.assertRaises(ValidationError):
            Satellite(**{**self.test_data, "type": "walker"})

    def test_orbit_required(self):
        """
        Test that orbit is required.
        """
        with self.assertRaises(ValidationError):
            Satellite(name="Test Satellite")

    def test_instruments_default_inherited_from_space_system(self):
        """
        Test that omitting instruments falls back to SpaceSystem's default
        of a single generic Instrument.
        """
        sat = Satellite(name="Test Satellite", orbit=self.test_data["orbit"])
        self.assertEqual(len(sat.instruments), 1)
        self.assertEqual(sat.instruments[0], Instrument())
