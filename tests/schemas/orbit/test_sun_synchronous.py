"""
Unit tests for the SunSynchronousOrbit schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, time, timezone

from tatc.schemas import SunSynchronousOrbit


class TestSunSynchronousOrbit(unittest.TestCase):
    """
    Unit tests for the SunSynchronousOrbit schema.
    """

    def setUp(self):
        self.test_data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(10, 30),
            "equator_crossing_ascending": True,
        }
        self.test_orbit = SunSynchronousOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the SunSynchronousOrbit schema correctly initializes with valid data.
        """
        good_data = {
            "mean_altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(10, 30),
            "equator_crossing_ascending": True,
        }
        o = SunSynchronousOrbit(**good_data)
        self.assertEqual(o.mean_altitude, good_data.get("mean_altitude"))
        self.assertEqual(o.true_anomaly, good_data.get("true_anomaly"))
        self.assertEqual(
            o.equator_crossing_time, good_data.get("equator_crossing_time")
        )
        self.assertEqual(
            o.equator_crossing_ascending, good_data.get("equator_crossing_ascending")
        )

    def test_get_derived_orbit(self):
        """
        Test that the SunSynchronousOrbit schema correctly derives a new 
        orbit with specified mean anomaly and RAAN offsets.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.get_mean_anomaly(),
            self.test_orbit.get_mean_anomaly() + 20,
            delta=0.001,
        )
        self.assertAlmostEqual(
            derived_orbit.get_right_ascension_ascending_node(),
            self.test_orbit.get_right_ascension_ascending_node() + 10,
            delta=0.001,
        )

    def test_to_gp_orbit(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts 
        to a general perturbations orbit.
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
        self.assertAlmostEqual(gp_orbit.get_inclination(), 97.7, delta=0.1)

    def test_to_gp_orbit_raan_ascending_equinox(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts 
        to a general perturbations orbit with RAAN at ascending equinox.
        """
        data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 3, 20, 3, 49, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": True,
        }
        gp_orbit = SunSynchronousOrbit(**data).to_gp_orbit()
        self.assertAlmostEqual(
            min(
                gp_orbit.get_right_ascension_ascending_node(),
                360.0 - gp_orbit.get_right_ascension_ascending_node(),
            ),
            0.0,
            delta=0.25,
        )

    def test_to_gp_orbit_raan_descending_equinox(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts 
        to a general perturbations orbit with RAAN at descending equinox.
        """
        data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 3, 20, 3, 49, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": False,
        }
        gp_orbit = SunSynchronousOrbit(**data).to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(), 180.0, delta=0.25
        )

    def test_to_gp_orbit_raan_ascending_solstice(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts 
        to a general perturbations orbit with RAAN at ascending solstice.
        """
        data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 6, 21, 9, 14, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": True,
        }
        gp_orbit = SunSynchronousOrbit(**data).to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(), 90.0, delta=0.25
        )

    def test_to_gp_orbit_raan_descending_solstice(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts 
        to a general perturbations orbit with RAAN at descending solstice.
        """
        data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 6, 21, 9, 14, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": False,
        }
        gp_orbit = SunSynchronousOrbit(**data).to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(), 270.0, delta=0.25
        )
