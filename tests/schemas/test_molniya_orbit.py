import unittest

from datetime import datetime, timezone

from tatc.schemas import MolniyaOrbit


class TestMolniyaOrbit(unittest.TestCase):
    def setUp(self):
        self.test_data = {
            "altitude": 1199100,
            "true_anomaly": 53.42,
            "epoch": datetime(2024, 11, 20, 7, 11, 20, 742432, tzinfo=timezone.utc),
            "inclination": 64.2583,
            "right_ascension_ascending_node": 11.0394,
            "perigee_argument": 289.3942,
            "orbit_period": 86400,
        }
        self.test_orbit = MolniyaOrbit(**self.test_data)

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
        self.assertEqual(
            self.test_orbit.perigee_argument, self.test_data.get("perigee_argument")
        )

    def test_get_derived_orbit(self):
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.right_ascension_ascending_node,
            self.test_orbit.right_ascension_ascending_node + 10,
            delta=0.001,
        )

    def test_to_tle(self):
        tle = self.test_orbit.to_tle()
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
        self.assertAlmostEqual(
            tle.get_perigee_argument(),
            self.test_data.get("perigee_argument"),
        )
