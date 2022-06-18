import unittest

from datetime import datetime, timedelta, timezone
import numpy as np

from tatc.schemas import CircularOrbit
from tatc.constants import earth_mean_radius, earth_mu


class TestOrbitBase(unittest.TestCase):
    def setUp(self):
        self.test_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
        }
        self.test_orbit = CircularOrbit(**self.test_data)

    def test_get_semimajor_axis(self):
        self.assertEqual(
            self.test_orbit.get_semimajor_axis(),
            self.test_data.get("altitude") + earth_mean_radius,
        )

    def test_get_mean_anomaly(self):
        self.assertEqual(
            self.test_orbit.get_mean_anomaly(), self.test_data.get("true_anomaly")
        )

    def test_get_mean_motion(self):
        orbit_period = (
            2
            * np.pi
            * np.sqrt(
                np.power(earth_mean_radius + self.test_orbit.altitude, 3) / earth_mu
            )
        )
        self.assertAlmostEqual(
            self.test_orbit.get_mean_motion(), 1 / (orbit_period / 86400), delta=0.001
        )

    def test_get_orbit_period(self):
        orbit_period = (
            2
            * np.pi
            * np.sqrt(
                np.power(earth_mean_radius + self.test_orbit.altitude, 3) / earth_mu
            )
        )
        self.assertAlmostEqual(
            self.test_orbit.get_orbit_period(),
            timedelta(seconds=orbit_period),
            delta=1.0,
        )
