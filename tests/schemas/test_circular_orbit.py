"""
Unit tests for the CircularOrbit schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest
from datetime import datetime, timedelta, timezone

import numpy as np

from tatc.constants import EARTH_MEAN_RADIUS, EARTH_MU
from tatc.schemas import CircularOrbit
from tatc.utils.orbital import semimajor_axis_to_mean_motion


class TestCircularOrbit(unittest.TestCase):
    """
    Unit tests for the CircularOrbit schema.
    """
    def setUp(self):
        self.test_data = {
            "mean_altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, tzinfo=timezone.utc),
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
        }
        self.test_orbit = CircularOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the CircularOrbit schema correctly initializes with valid data.
        """
        self.assertEqual(self.test_orbit.mean_altitude, self.test_data.get("mean_altitude"))
        self.assertEqual(
            self.test_orbit.true_anomaly, self.test_data.get("true_anomaly")
        )
        self.assertEqual(self.test_orbit.epoch, self.test_data.get("epoch"))
        self.assertEqual(self.test_orbit.inclination, self.test_data.get("inclination"))
        self.assertEqual(
            self.test_orbit.right_ascension_ascending_node,
            self.test_data.get("right_ascension_ascending_node"),
        )

    def test_good_data_iso8601_datetime(self):
        """
        Test that the CircularOrbit schema correctly initializes with valid data
        when the epoch is provided as an ISO 8601 string.
        """
        good_data = {
            "mean_altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": "2022-01-01T12:00:00Z",
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
        }
        o = CircularOrbit(**good_data)
        self.assertEqual(o.mean_altitude, good_data.get("mean_altitude"))
        self.assertEqual(o.true_anomaly, good_data.get("true_anomaly"))
        self.assertEqual(o.epoch, datetime(2022, 1, 1, 12, tzinfo=timezone.utc))
        self.assertEqual(o.inclination, good_data.get("inclination"))
        self.assertEqual(
            o.right_ascension_ascending_node,
            good_data.get("right_ascension_ascending_node"),
        )

    def test_get_semimajor_axis(self):
        """
        Test that the CircularOrbit class correctly calculates the semimajor axis.
        """
        self.assertEqual(
            self.test_orbit.get_semimajor_axis(),
            self.test_data.get("mean_altitude") + EARTH_MEAN_RADIUS,
        )

    def test_get_mean_anomaly(self):
        """
        Test that the CircularOrbit class correctly calculates the mean anomaly.
        """
        self.assertEqual(
            self.test_orbit.get_mean_anomaly(), self.test_data.get("true_anomaly")
        )

    def test_get_mean_motion(self):
        """
        Test that the CircularOrbit class correctly calculates the mean motion.
        """
        self.assertAlmostEqual(
            self.test_orbit.get_mean_motion(), 
            semimajor_axis_to_mean_motion(EARTH_MEAN_RADIUS + self.test_orbit.mean_altitude),
            delta=0.001
        )

    def test_get_orbit_period(self):
        """
        Test that the CircularOrbit class correctly calculates the orbit period.
        """
        orbit_period = (
            2
            * np.pi
            * np.sqrt(
                np.power(EARTH_MEAN_RADIUS + self.test_orbit.mean_altitude, 3) / EARTH_MU
            )
        )
        self.assertAlmostEqual(
            self.test_orbit.get_orbit_period(),
            timedelta(seconds=orbit_period),
            delta=1.0,
        )

    def test_get_derived_orbit(self):
        """
        Test that the CircularOrbit schema correctly derives a new orbit.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.get_mean_anomaly(),
            self.test_orbit.get_mean_anomaly() + 20,
            delta=0.001,
        )
        self.assertAlmostEqual(
            derived_orbit.right_ascension_ascending_node,
            self.test_orbit.right_ascension_ascending_node + 10,
            delta=0.001,
        )

    def test_to_gp_orbit(self):
        """
        Test that the CircularOrbit schema correctly converts to a general perturbations orbit.
        """
        gp_orbit = self.test_orbit.to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_mean_altitude(), self.test_data.get("mean_altitude"), delta=1.0
        )
        self.assertAlmostEqual(
            gp_orbit.get_true_anomaly(), self.test_data.get("true_anomaly"), delta=0.001
        )
        self.assertAlmostEqual(
            gp_orbit.get_epoch().timestamp(),
            self.test_data.get("epoch").timestamp(),
            delta=1,
        )
        self.assertEqual(
            gp_orbit.get_inclination(),
            self.test_data.get("inclination"),
        )
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(),
            self.test_data.get("right_ascension_ascending_node"),
        )
