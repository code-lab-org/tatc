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

    def test_latitude_boundary_values(self):
        """
        Test that latitude values exactly at the poles (-90, 90 degrees)
        are accepted rather than rejected by an off-by-one bound.
        """
        self.assertEqual(Point(latitude=90, longitude=0).latitude, 90)
        self.assertEqual(Point(latitude=-90, longitude=0).latitude, -90)

    def test_longitude_boundary_values(self):
        """
        Test that longitude values exactly at the antimeridian (-180, 180
        degrees) are accepted rather than rejected by an off-by-one bound.
        """
        self.assertEqual(Point(latitude=0, longitude=180).longitude, 180)
        self.assertEqual(Point(latitude=0, longitude=-180).longitude, -180)

    def test_id_defaults_to_zero(self):
        """
        Test that omitting id defaults to 0.
        """
        self.assertEqual(Point(latitude=0, longitude=0).id, 0)

    def test_id_rejects_negative(self):
        """
        Test that a negative id is rejected, since id is a
        NonNegativeInt.
        """
        with self.assertRaises(ValidationError):
            Point(id=-1, latitude=0, longitude=0)

    def test_id_rejects_non_integer(self):
        """
        Test that a fractional id (not cleanly convertible to int) is
        rejected rather than silently truncated.
        """
        with self.assertRaises(ValidationError):
            Point(id=5.5, latitude=0, longitude=0)

    def test_elevation_defaults_to_zero(self):
        """
        Test that omitting elevation defaults to 0.
        """
        self.assertEqual(Point(latitude=0, longitude=0).elevation, 0)

    def test_elevation_accepts_negative_values(self):
        """
        Test that elevation accepts negative values (e.g. below-sea-level
        locations like Death Valley or the Dead Sea), since it has no
        lower-bound constraint.
        """
        self.assertEqual(
            Point(latitude=0, longitude=0, elevation=-430.5).elevation, -430.5
        )
