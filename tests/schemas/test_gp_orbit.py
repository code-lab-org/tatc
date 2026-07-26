"""
Unit tests for the GeneralPerturbationsOrbit schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest
from datetime import datetime, timezone

from pydantic import ValidationError

from tatc.schemas import GeneralPerturbationsOrbit


class TestGPOrbit(unittest.TestCase):
    """
    Unit tests for the GeneralPerturbationsOrbit schema.
    """
    def setUp(self):
        self.test_tle = [
            "1 25544U 98067A   21156.30527927  .00003432  00000-0  70541-4 0  9993",
            "2 25544  51.6455  41.4969 0003508  68.0432  78.3395 15.48957534286754",
        ]

    def test_get_catalog_number(self):
        """
        Test that the catalog number can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.elements[0].norad_cat_id, 25544)

    def test_get_epoch(self):
        """
        Test that the epoch can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(
            gp_orbit.elements[0].epoch,
            datetime(2021, 6, 5, 7, 19, 36, 128928, tzinfo=timezone.utc),
        )

    def test_get_first_derivative_mean_motion(self):
        """
        Test that the first derivative of the mean motion can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.elements[0].mean_motion_dot * (24*60*60)**2/360, 0.00003432)

    def test_get_second_derivative_mean_motion(self):
        """
        Test that the second derivative of the mean motion can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.elements[0].mean_motion_ddot * (24*60*60)**3/360, 0.0)

    def test_get_b_star(self):
        """
        Test that the B* drag term can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.elements[0].bstar, 0.000070541)

    def test_get_inclination(self):
        """
        Test that the inclination can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.get_inclination(), 51.6455)

    def test_get_right_ascension_ascending_node(self):
        """
        Test that the right ascension of the ascending node can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.get_right_ascension_ascending_node(), 41.4969)

    def test_get_eccentricity(self):
        """
        Test that the eccentricity can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.get_eccentricity(), 0.0003508)

    def test_get_perigee_argument(self):
        """
        Test that the argument of perigee can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.get_perigee_argument(), 68.0432)

    def test_get_mean_anomaly(self):
        """
        Test that the mean anomaly can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertEqual(gp_orbit.get_mean_anomaly(), 78.3395)

    def test_get_mean_motion(self):
        """
        Test that the mean motion can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertAlmostEqual(gp_orbit.get_mean_motion() * (24*60*60) / 360, 15.48957534)

    def test_get_semimajor_axis(self):
        """
        Test that the semi-major axis can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertAlmostEqual(gp_orbit.get_semimajor_axis(), 6797911, delta=1.0)

    def test_get_mean_altitude(self):
        """
        Test that the mean altitude can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertAlmostEqual(gp_orbit.get_mean_altitude(), 426902, delta=1.0)

    def test_get_true_anomaly(self):
        """
        Test that the true anomaly can be retrieved from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        self.assertAlmostEqual(gp_orbit.elements[0].get_true_anomaly(), 78.3788725993742)

    def test_get_derived_orbit(self):
        """
        Test that a derived orbit can be created from the GeneralPerturbationsOrbit object.
        """
        gp_orbit = GeneralPerturbationsOrbit.from_tle(self.test_tle)
        derived_orbit = gp_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.get_mean_anomaly(),
            gp_orbit.get_mean_anomaly() + 20,
            delta=0.001,
        )
        self.assertAlmostEqual(
            derived_orbit.get_right_ascension_ascending_node(),
            gp_orbit.get_right_ascension_ascending_node() + 10,
            delta=0.001,
        )
