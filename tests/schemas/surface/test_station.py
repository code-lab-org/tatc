"""
Unit tests for the GroundStation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest
from datetime import timedelta

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
