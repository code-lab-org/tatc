import unittest

from datetime import datetime, timedelta, timezone
import numpy as np

from tatc.schemas import CircularOrbit
from tatc.constants import earth_mean_radius, earth_mu


class TestOrbitBase(unittest.TestCase):
    def test_get_semimajor_axis(self):
        good_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
        }
        o = CircularOrbit(**good_data)
        self.assertEqual(
            o.get_semimajor_axis(), good_data.get("altitude") + earth_mean_radius
        )

    def test_get_mean_motion(self):
        good_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
        }
        o = CircularOrbit(**good_data)
        orbit_period = (
            2 * np.pi * np.sqrt(np.power(earth_mean_radius + o.altitude, 3) / earth_mu)
        )
        self.assertAlmostEqual(
            o.get_mean_motion(), 1 / (orbit_period / 86400), delta=0.001
        )

    def test_get_orbit_period(self):
        good_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
        }
        o = CircularOrbit(**good_data)
        orbit_period = (
            2 * np.pi * np.sqrt(np.power(earth_mean_radius + o.altitude, 3) / earth_mu)
        )
        self.assertAlmostEqual(
            o.get_orbit_period(), timedelta(seconds=orbit_period), delta=1.0
        )
