import unittest

from tatc.schemas.point import Point
from pydantic import ValidationError


class TestPoint(unittest.TestCase):
    def test_good_data(self):
        good_data = {"id": 0, "latitude": 0.0, "longitude": 0.0}
        point = Point(**good_data)
        self.assertEqual(point.id, good_data.get("id"))
        self.assertEqual(point.latitude, good_data.get("latitude"))
        self.assertEqual(point.longitude, good_data.get("longitude"))

    def test_bad_latitude_too_big(self):
        bad_data = {"id": 0, "latitude": 100.0, "longitude": 0.0}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_latitude_too_small(self):
        bad_data = {"id": 0, "latitude": -90.1, "longitude": 0.0}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_latitude_missing(self):
        bad_data = {"id": 0, "longitude": 0.0}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_too_big(self):
        bad_data = {"id": 0, "latitude": 0.0, "longitude": 180.1}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_too_small(self):
        bad_data = {"id": 0, "latitude": 0.0, "longitude": -180.1}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_missing(self):
        bad_data = {"id": 0, "latitude": 0.0}
        with self.assertRaises(ValidationError):
            Point(**bad_data)
