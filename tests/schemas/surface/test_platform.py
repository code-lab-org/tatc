"""
Unit tests for the SurfacePlatform schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import numpy as np
from pydantic import ValidationError
from pyproj import Geod
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc.schemas import (
    ConicalInstrument,
    Instrument,
    PointedInstrument,
    SurfacePlatform,
    TerrainMask,
)
from tatc.utils import VelocityFrame


class TestSurfacePlatform(unittest.TestCase):
    """
    Unit tests for the SurfacePlatform schema.
    """

    def test_good_data(self):
        """
        Test that the SurfacePlatform schema correctly initializes with valid data.
        """
        good_data = {
            "name": "Tower",
            "latitude": 33.45,
            "longitude": -112.07,
            "elevation": 300,
            "height": 30,
            "heading": 90,
            "instruments": [
                {"name": "Camera", "field_of_regard": 60},
                {
                    "name": "Pointed",
                    "cross_track_field_of_view": 40,
                    "along_track_field_of_view": 30,
                    "pitch_angle": 60,
                },
            ],
        }
        o = SurfacePlatform(**good_data)
        self.assertEqual(o.name, "Tower")
        self.assertEqual(o.height, 30)
        self.assertEqual(o.heading, 90)
        self.assertIsInstance(o.instruments[0], Instrument)
        self.assertIsInstance(o.instruments[1], PointedInstrument)
        self.assertEqual(o.get_altitude(), 330)

    def test_default_instrument(self):
        """
        Test that a platform hosts a default instrument.
        """
        o = SurfacePlatform(name="Tower", latitude=0, longitude=0, height=30)
        self.assertEqual(len(o.instruments), 1)
        self.assertEqual(o.heading, 0)
        self.assertEqual(o.elevation, 0)

    def test_bad_data(self):
        """
        Test that invalid heights and headings raise validation errors.
        """
        for bad in [
            {"height": 0},
            {"height": -10},
            {"height": 30, "heading": 360},
            {"height": 30, "heading": -1},
            {"height": 30, "instruments": []},
        ]:
            with self.assertRaises(ValidationError):
                SurfacePlatform(name="Tower", latitude=0, longitude=0, **bad)
        with self.assertRaises(ValidationError):
            SurfacePlatform(name="Tower", latitude=0, longitude=0)

    def test_unsupported_instruments(self):
        """
        Test that instruments whose views depend on orbital motion raise
        validation errors.
        """
        for instrument in [
            ConicalInstrument(cone_angle=45, along_track_field_of_view=1),
            PointedInstrument(
                cross_track_field_of_view=10,
                along_track_field_of_view=10,
                roll_angle_profile=[(0, 10)],
            ),
            PointedInstrument(
                cross_track_field_of_view=10,
                along_track_field_of_view=10,
                pitch_angle_profile=[(0, 10)],
            ),
            PointedInstrument(
                cross_track_field_of_view=10,
                along_track_field_of_view=10,
                velocity_frame=VelocityFrame.INERTIAL,
            ),
        ]:
            with self.assertRaises(ValidationError):
                SurfacePlatform(
                    name="Tower",
                    latitude=0,
                    longitude=0,
                    height=30,
                    instruments=[instrument],
                )

    def test_get_track(self):
        """
        Test that the track is at the instruments' Earth-fixed position, at
        one or more times.
        """
        o = SurfacePlatform(
            name="Balloon", latitude=33.45, longitude=-112.07, elevation=300, height=500
        )
        expected = wgs84.latlon(33.45, -112.07, 800).itrs_xyz.m
        position, _ = o.get_track().frame_xyz_and_velocity(itrs)
        np.testing.assert_allclose(position.m, expected, atol=1e-3)
        times = [datetime(2020, 6, 1, hour, tzinfo=timezone.utc) for hour in (0, 12)]
        track = o.get_track(times)
        position, _ = track.frame_xyz_and_velocity(itrs)
        np.testing.assert_allclose(
            position.m, np.stack([expected] * 2, axis=1), atol=1e-3
        )
        track = o.get_track(times[0])
        self.assertEqual(track.t.utc_datetime(), times[0])

    def test_valid_observation(self):
        """
        Test that an instrument's sunlit requirement varies with time.
        """
        o = SurfacePlatform(
            name="Tower",
            latitude=0,
            longitude=0,
            height=30,
            instruments=[Instrument(req_target_sunlit=True)],
        )
        times = [datetime(2020, 3, 20, hour, tzinfo=timezone.utc) for hour in (0, 12)]
        valid = o.instruments[0].is_valid_observation(
            o.get_track(times), wgs84.latlon(0, 0)
        )
        self.assertListEqual(list(valid), [False, True])

    def test_compute_visible_region(self):
        """
        Test that the region visible from a platform at a target elevation
        above it reaches the horizon of the ground without a terrain mask,
        and only as far as the lines of sight clearing a terrain mask.
        """
        o = SurfacePlatform(
            name="Tower", latitude=33.45, longitude=-112.07, elevation=300, height=30
        )
        geod = Geod(ellps="WGS84")

        def distances(region):
            lon, lat = np.array(region.exterior.coords).T
            return geod.inv(
                np.full(len(lon), o.longitude), np.full(len(lat), o.latitude), lon, lat
            )

        # lines of sight tangent to the ground, extended to 1 km above it
        _, _, d = distances(o.compute_visible_region(1300))
        horizon = np.sqrt(2 * 6371e3 * 30) + np.sqrt(2 * 6371e3 * 1000)
        np.testing.assert_allclose(d, horizon, rtol=0.005)
        # a 2 degree ridge to the west
        o.terrain_mask = TerrainMask(
            azimuth=[0, 240, 250, 290, 300], min_elevation_angle=[-1, -1, 2, 2, -1]
        )
        azimuth, _, d = distances(o.compute_visible_region(1300))
        west = (azimuth % 360 >= 250) & (azimuth % 360 <= 290)
        # lines of sight 2 degrees above horizontal rise 970 m (with the
        # curvature drop) at the distance d: d tan(2) + d^2 / 2R = 970
        slope, curvature = np.tan(np.radians(2)), 1 / (2 * 6371e3)
        reach = (np.sqrt(slope**2 + 4 * curvature * 970) - slope) / (2 * curvature)
        np.testing.assert_allclose(d[west], reach, rtol=0.005)
        east = (azimuth > 0) & (azimuth < 180)
        np.testing.assert_allclose(d[east], horizon, rtol=0.005)
        with self.assertRaises(ValueError):
            o.compute_visible_region(330)
