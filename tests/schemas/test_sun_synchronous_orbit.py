import unittest

from tatc.schemas import SunSynchronousOrbit
from datetime import datetime, time, timezone


class TestSunSynchronousOrbit(unittest.TestCase):
    def test_good_data(self):
        good_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(10, 30),
            "equator_crossing_ascending": True,
        }
        o = SunSynchronousOrbit(**good_data)
        self.assertEqual(o.altitude, good_data.get("altitude"))
        self.assertEqual(o.true_anomaly, good_data.get("true_anomaly"))
        self.assertEqual(
            o.equator_crossing_time, good_data.get("equator_crossing_time")
        )
        self.assertEqual(
            o.equator_crossing_ascending, good_data.get("equator_crossing_ascending")
        )

    def test_to_tle(self):
        good_data = {
            "altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(10, 30),
            "equator_crossing_ascending": True,
        }
        o = SunSynchronousOrbit(**good_data)
        self.assertAlmostEqual(
            o.to_tle().get_altitude(), good_data.get("altitude"), delta=1.0
        )
        self.assertAlmostEqual(
            o.to_tle().get_true_anomaly(), good_data.get("true_anomaly"), delta=0.001
        )
        self.assertAlmostEqual(
            o.to_tle().get_epoch().timestamp(),
            good_data.get("epoch").timestamp(),
            delta=1,
        )
        self.assertAlmostEqual(o.to_tle().get_inclination(), 97.7, delta=0.1)

    def test_to_tle_raan_ascending_equinox(self):
        good_data = {
            "altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 3, 20, 3, 49, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": True,
        }
        o = SunSynchronousOrbit(**good_data)
        self.assertAlmostEqual(
            abs(o.to_tle().get_right_ascension_ascending_node() - 360.0),
            0.0,
            delta=0.25,
        )

    def test_to_tle_raan_descending_equinox(self):
        good_data = {
            "altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 3, 20, 3, 49, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": False,
        }
        o = SunSynchronousOrbit(**good_data)
        self.assertAlmostEqual(
            o.to_tle().get_right_ascension_ascending_node(), 180.0, delta=0.25
        )

    def test_to_tle_raan_ascending_solstice(self):
        good_data = {
            "altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 6, 21, 9, 14, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": True,
        }
        o = SunSynchronousOrbit(**good_data)
        self.assertAlmostEqual(
            o.to_tle().get_right_ascension_ascending_node(), 90.0, delta=0.25
        )

    def test_to_tle_raan_descending_solstice(self):
        good_data = {
            "altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 6, 21, 9, 14, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": False,
        }
        o = SunSynchronousOrbit(**good_data)
        self.assertAlmostEqual(
            o.to_tle().get_right_ascension_ascending_node(), 270.0, delta=0.25
        )
