import unittest

from pydantic import ValidationError

from tatc.schemas import Point


class TestPoint(unittest.TestCase):
    def test_good_data(self):
        good_data = {"id": 42, "latitude": 40.74259, "longitude": -74.02686}
        o = Point(**good_data)
        self.assertEqual(o.id, good_data.get("id"))
        self.assertEqual(o.latitude, good_data.get("latitude"))
        self.assertEqual(o.longitude, good_data.get("longitude"))

    def test_bad_latitude_too_big(self):
        bad_data = {"id": 0, "latitude": 100.0, "longitude": -74.02686}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_latitude_too_small(self):
        bad_data = {"id": 0, "latitude": -90.1, "longitude": -74.02686}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_latitude_missing(self):
        bad_data = {"id": 0, "longitude": -74.02686}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_too_big(self):
        bad_data = {"id": 0, "latitude": 40.74259, "longitude": 180.1}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_too_small(self):
        bad_data = {"id": 0, "latitude": 40.74259, "longitude": -180.1}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_missing(self):
        bad_data = {"id": 0, "latitude": 40.74259}
        with self.assertRaises(ValidationError):
            Point(**bad_data)
