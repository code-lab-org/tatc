"""
Unit tests for the Satellite schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

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

    def test_generate_members(self):
        """
        Test that the Satellite schema correctly generates members.
        """
        members = self.test_sat.generate_members()
        self.assertEqual(len(members), 1)
        self.assertEqual(members[0], self.test_sat)
