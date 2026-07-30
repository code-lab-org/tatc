"""
Unit tests for the Point schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from pydantic import ValidationError

from tatc.schemas import Point


class TestPoint(unittest.TestCase):
    """
    Unit tests for the Point schema.
    """
    def test_good_data(self):
        """
        Test that the Point schema correctly initializes with valid data.
        """
        good_data = {"id": 42, "latitude": 40.74259, "longitude": -74.02686}
        o = Point(**good_data)
        self.assertEqual(o.id, good_data.get("id"))
        self.assertEqual(o.latitude, good_data.get("latitude"))
        self.assertEqual(o.longitude, good_data.get("longitude"))

    def test_bad_latitude_too_big(self):
        """
        Test that the Point schema raises a ValidationError when the latitude is too large.
        """
        bad_data = {"id": 0, "latitude": 100.0, "longitude": -74.02686}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_latitude_too_small(self):
        """
        Test that the Point schema raises a ValidationError when the latitude is too small.
        """
        bad_data = {"id": 0, "latitude": -90.1, "longitude": -74.02686}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_latitude_missing(self):
        """
        Test that the Point schema raises a ValidationError when the latitude is missing.
        """
        bad_data = {"id": 0, "longitude": -74.02686}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_too_big(self):
        """
        Test that the Point schema raises a ValidationError when the longitude is too large.
        """
        bad_data = {"id": 0, "latitude": 40.74259, "longitude": 180.1}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_too_small(self):
        """
        Test that the Point schema raises a ValidationError when the longitude is too small.
        """
        bad_data = {"id": 0, "latitude": 40.74259, "longitude": -180.1}
        with self.assertRaises(ValidationError):
            Point(**bad_data)

    def test_bad_longitude_missing(self):
        """
        Test that the Point schema raises a ValidationError when the longitude is missing.
        """
        bad_data = {"id": 0, "latitude": 40.74259}
        with self.assertRaises(ValidationError):
            Point(**bad_data)
