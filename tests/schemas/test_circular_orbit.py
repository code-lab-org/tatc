import unittest

from tatc.schemas import CircularOrbit
from datetime import datetime, timezone


class TestCircularOrbit(unittest.TestCase):
    def test_good_data(self):
        good_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
        }
        o = CircularOrbit(**good_data)
        self.assertEqual(o.altitude, good_data.get("altitude"))
        self.assertEqual(o.true_anomaly, good_data.get("true_anomaly"))
        self.assertEqual(o.epoch, good_data.get("epoch"))
        self.assertEqual(o.inclination, good_data.get("inclination"))
        self.assertEqual(
            o.right_ascension_ascending_node,
            good_data.get("right_ascension_ascending_node"),
        )

    def test_good_data_iso8601_datetime(self):
        good_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": "2022-01-01T12:00:00Z",
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
        }
        o = CircularOrbit(**good_data)
        self.assertEqual(o.altitude, good_data.get("altitude"))
        self.assertEqual(o.true_anomaly, good_data.get("true_anomaly"))
        self.assertEqual(o.epoch, datetime(2022, 1, 1, 12, tzinfo=timezone.utc))
        self.assertEqual(o.inclination, good_data.get("inclination"))
        self.assertEqual(
            o.right_ascension_ascending_node,
            good_data.get("right_ascension_ascending_node"),
        )

    def test_to_tle(self):
        good_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, tzinfo=timezone.utc),
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
        }
        o = CircularOrbit(**good_data)
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
        self.assertEqual(
            o.to_tle().get_inclination(),
            good_data.get("inclination"),
        )
        self.assertAlmostEqual(
            o.to_tle().get_right_ascension_ascending_node(),
            good_data.get("right_ascension_ascending_node"),
        )
