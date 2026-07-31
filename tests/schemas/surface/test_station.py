"""
Unit tests for the GroundStation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest
from datetime import timedelta

from pydantic import ValidationError

from tatc.schemas import GroundStation


class TestGroundStation(unittest.TestCase):
    """
    Unit tests for the GroundStation schema.
    """
    def test_good_data(self):
        """
        Test that the GroundStation schema correctly initializes with valid data.
        """
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
        """
        Test that the GroundStation schema correctly initializes with valid data
        when min_access_time is provided as seconds.
        """
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

    def test_defaults(self):
        """
        Test that min_elevation_angle defaults to 0 and min_access_time
        defaults to a zero timedelta when omitted.
        """
        o = GroundStation(name="test", latitude=40.74259, longitude=-74.02686)
        self.assertEqual(o.min_elevation_angle, 0)
        self.assertEqual(o.min_access_time, timedelta(0))

    def test_bad_name_missing(self):
        """
        Test that the GroundStation schema raises a ValidationError when
        the required name field is missing.
        """
        bad_data = {"latitude": 40.74259, "longitude": -74.02686}
        with self.assertRaises(ValidationError):
            GroundStation(**bad_data)

    def test_bad_min_elevation_angle_negative(self):
        """
        Test that a negative min_elevation_angle is rejected.
        """
        with self.assertRaises(ValidationError):
            GroundStation(
                name="test", latitude=0, longitude=0, min_elevation_angle=-0.1
            )

    def test_bad_min_elevation_angle_too_large(self):
        """
        Test that a min_elevation_angle above 90 degrees is rejected.
        """
        with self.assertRaises(ValidationError):
            GroundStation(
                name="test", latitude=0, longitude=0, min_elevation_angle=90.1
            )

    def test_min_elevation_angle_boundary_values(self):
        """
        Test that min_elevation_angle values exactly at its bounds (0 and
        90 degrees) are accepted.
        """
        self.assertEqual(
            GroundStation(
                name="test", latitude=0, longitude=0, min_elevation_angle=0
            ).min_elevation_angle,
            0,
        )
        self.assertEqual(
            GroundStation(
                name="test", latitude=0, longitude=0, min_elevation_angle=90
            ).min_elevation_angle,
            90,
        )

    def test_inherits_point_latitude_validation(self):
        """
        Test that GroundStation inherits Point's latitude validation
        (out-of-range latitude is still rejected), confirming the
        inheritance relationship is functionally in effect and not just
        structural.
        """
        with self.assertRaises(ValidationError):
            GroundStation(name="test", latitude=100, longitude=0)
