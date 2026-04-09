import unittest

from datetime import datetime, timezone

from tatc import constants
from tatc.schemas import GeosynchronousOrbit


class TestGeosynchronousOrbit(unittest.TestCase):
    def setUp(self):
        self.test_data = {
            "altitude": 35786000,
            "epoch": datetime(2022, 1, 1, 12, tzinfo=timezone.utc),
            "inclination": 0.0,
            "right_ascension_ascending_node": 0.0,
            "longitude": 10.0,
        }
        self.test_orbit = GeosynchronousOrbit(**self.test_data)

    def test_good_data(self):
        self.assertEqual(self.test_orbit.altitude, self.test_data.get("altitude"))
        self.assertEqual(self.test_orbit.epoch, self.test_data.get("epoch"))
        self.assertEqual(self.test_orbit.inclination, self.test_data.get("inclination"))
        self.assertEqual(
            self.test_orbit.right_ascension_ascending_node,
            self.test_data.get("right_ascension_ascending_node"),
        )
        self.assertEqual(self.test_orbit.longitude, self.test_data.get("longitude"))

    def test_default_altitude(self):
        o = GeosynchronousOrbit(longitude=0.0)
        self.assertEqual(o.altitude, 35786000)

    def test_good_data_iso8601_datetime(self):
        good_data = {
            "altitude": 35786000,
            "epoch": "2022-01-01T12:00:00Z",
            "inclination": 0.0,
            "right_ascension_ascending_node": 0.0,
            "longitude": 10.0,
        }
        o = GeosynchronousOrbit(**good_data)
        self.assertEqual(o.altitude, good_data.get("altitude"))
        self.assertEqual(o.epoch, datetime(2022, 1, 1, 12, tzinfo=timezone.utc))
        self.assertEqual(o.inclination, good_data.get("inclination"))
        self.assertEqual(
            o.right_ascension_ascending_node,
            good_data.get("right_ascension_ascending_node"),
        )
        self.assertEqual(o.longitude, good_data.get("longitude"))

    def test_get_true_anomaly_matches_sidereal_time(self):
        t = constants.timescale.from_datetime(self.test_data.get("epoch"))
        expected = (self.test_data.get("longitude") + t.gmst * 15) % 360
        self.assertAlmostEqual(
            self.test_orbit.get_true_anomaly(), expected, delta=1e-6
        )

    def test_get_derived_orbit(self):
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.longitude,
            (self.test_orbit.longitude + 20) % 360,
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
            delta=0.001,
        )

    def test_to_tle_subpoint_matches_longitude(self):
        tle = self.test_orbit.to_tle()
        sat = tle.as_skyfield()
        t = constants.timescale.from_datetime(self.test_data.get("epoch"))
        subpoint = sat.at(t).subpoint()
        self.assertAlmostEqual(
            subpoint.longitude.degrees, self.test_data.get("longitude"), delta=0.1
        )
        self.assertAlmostEqual(subpoint.latitude.degrees, 0.0, delta=0.1)

    def test_to_tle_subpoint_multiple_longitudes(self):
        for longitude in (0.0, 45.0, 137.5, 270.0):
            with self.subTest(longitude=longitude):
                o = GeosynchronousOrbit(
                    epoch=self.test_data.get("epoch"),
                    longitude=longitude,
                    inclination=0.0,
                )
                tle = o.to_tle()
                sat = tle.as_skyfield()
                t = constants.timescale.from_datetime(self.test_data.get("epoch"))
                subpoint = sat.at(t).subpoint()
                sub_lon = subpoint.longitude.degrees % 360
                self.assertAlmostEqual(sub_lon, longitude, delta=0.1)
                self.assertAlmostEqual(subpoint.latitude.degrees, 0.0, delta=0.1)

    def test_bad_longitude_negative(self):
        with self.assertRaises(Exception):
            GeosynchronousOrbit(longitude=-1.0)

    def test_bad_longitude_too_large(self):
        with self.assertRaises(Exception):
            GeosynchronousOrbit(longitude=360.0)
