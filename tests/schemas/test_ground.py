import unittest

from tatc.schemas import GroundStation
from datetime import timedelta


class TestGroundStation(unittest.TestCase):
    def test_good_data(self):
        good_data = {
            "name": "test",
            "latitude": 40.74259,
            "longitude": -74.02686,
            "min_elevation_angle": 20.0,
            "min_access_time": timedelta(20),
        }
        o = GroundStation(**good_data)
        self.assertEqual(o.name, good_data.get("name"))
        self.assertEqual(o.latitude, good_data.get("latitude"))
        self.assertEqual(o.longitude, good_data.get("longitude"))
        self.assertEqual(o.min_elevation_angle, good_data.get("min_elevation_angle"))
        self.assertEqual(o.min_access_time, good_data.get("min_access_time"))

    def test_good_data_timedelta_seconds(self):
        good_data = {
            "name": "test",
            "latitude": 40.74259,
            "longitude": -74.02686,
            "min_elevation_angle": 20.0,
            "min_access_time": 20,
        }
        o = GroundStation(**good_data)
        self.assertEqual(o.name, good_data.get("name"))
        self.assertEqual(o.latitude, good_data.get("latitude"))
        self.assertEqual(o.longitude, good_data.get("longitude"))
        self.assertEqual(o.min_elevation_angle, good_data.get("min_elevation_angle"))
        self.assertEqual(
            o.min_access_time, timedelta(seconds=good_data.get("min_access_time"))
        )
