"""
Unit tests for the MolniyaOrbit schema.

@author: Paul T. Grogan <paul.t.grogan@asu.edu>
"""

import unittest

from tatc.schemas import MolniyaOrbit


class TestMolniyaOrbit(unittest.TestCase):
    """
    Unit tests for the MolniyaOrbit schema.
    """

    def setUp(self):
        self.test_data = {
            "perigee_altitude": 1199100,
            "right_ascension_ascending_node": 11.0394,
        }
        self.test_orbit = MolniyaOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the MolniyaOrbit schema correctly initializes with valid data.
        """
        self.assertEqual(
            self.test_orbit.perigee_altitude, self.test_data.get("perigee_altitude")
        )
        self.assertEqual(
            self.test_orbit.right_ascension_ascending_node,
            self.test_data.get("right_ascension_ascending_node"),
        )
        self.assertAlmostEqual(self.test_orbit.get_inclination(), 63.4, delta=0.1)
        self.assertAlmostEqual(
            self.test_orbit.get_orbit_period().total_seconds(),
            718*60, delta=60
        )
        self.assertEqual(
            self.test_orbit.get_perigee_argument(),
            270
        )
        self.assertAlmostEqual(
            self.test_orbit.get_eccentricity(),
            0.74, delta=0.1
        )

    def test_get_derived_orbit(self):
        """
        Test that the MolniyaOrbit schema correctly derives a new orbit with the specified parameters.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.right_ascension_ascending_node,
            self.test_orbit.right_ascension_ascending_node + 10,
            delta=0.001,
        )

    def test_to_gp_orbit(self):
        """
        Test that the MolniyaOrbit schema correctly converts to a general perturbations representation.
        """
        gp_orbit = self.test_orbit.to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(),
            self.test_data.get("right_ascension_ascending_node"),
            delta=0.1
        )
