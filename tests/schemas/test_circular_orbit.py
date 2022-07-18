import unittest

from datetime import datetime, timezone

from tatc.schemas import CircularOrbit


class TestCircularOrbit(unittest.TestCase):
    def setUp(self):
        self.test_data = {
            "altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, tzinfo=timezone.utc),
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
        }
        self.test_orbit = CircularOrbit(**self.test_data)

    def test_good_data(self):
        self.assertEqual(self.test_orbit.altitude, self.test_data.get("altitude"))
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

    def test_get_derived_orbit(self):
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

    def test_to_tle(self):
        tle = self.test_orbit.to_tle()
        self.assertAlmostEqual(
            tle.get_altitude(), self.test_data.get("altitude"), delta=1.0
        )
        self.assertAlmostEqual(
            tle.get_true_anomaly(), self.test_data.get("true_anomaly"), delta=0.001
        )
        self.assertAlmostEqual(
            tle.get_epoch().timestamp(),
            self.test_data.get("epoch").timestamp(),
            delta=1,
        )
        self.assertEqual(
            tle.get_inclination(),
            self.test_data.get("inclination"),
        )
        self.assertAlmostEqual(
            tle.get_right_ascension_ascending_node(),
            self.test_data.get("right_ascension_ascending_node"),
        )
