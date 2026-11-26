"""
Unit tests for the azimuths of samples on the ground.
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
from pyproj import Geod
from shapely.geometry import Point
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc.analysis import collect_observations, compute_sample_sizes
from tatc.schemas import (
    ConicalInstrument,
    GeosynchronousOrbit,
    PointedInstrument,
    Satellite,
    SunSynchronousOrbit,
)
from tatc.utils import (
    VelocityFrame,
    compute_ellipsoid_intersection,
    compute_ground_sample_azimuth,
    rectangular_to_geodetic,
)
from tatc.utils.projection import _compute_ground_axes, _compute_view_frame

GEOD = Geod(ellps="WGS84")
START = datetime(2026, 3, 1, 12, tzinfo=timezone.utc)


def _subpoint(satellite, time):
    """The point on the ellipsoid below a satellite at a time."""
    return wgs84.subpoint_of(satellite.get_track(time))


class TestGroundSampleAzimuth(unittest.TestCase):
    """
    Unit tests for `compute_ground_sample_azimuth`.
    """

    def setUp(self):
        self.satellite = Satellite(
            name="sat",
            orbit=SunSynchronousOrbit(
                altitude=700e3, equator_crossing_time="10:30", epoch=START
            ),
        )
        self.times = [START + timedelta(minutes=m) for m in range(0, 100, 7)]

    def test_nadir_along_ground_track(self):
        """At nadir, the along-track axis lies along the ground track."""
        for time in self.times:
            track = self.satellite.get_track(time)
            target = _subpoint(self.satellite, time)
            before = _subpoint(self.satellite, time - timedelta(seconds=0.5))
            after = _subpoint(self.satellite, time + timedelta(seconds=0.5))
            bearing, _, _ = GEOD.inv(
                before.longitude.degrees,
                before.latitude.degrees,
                after.longitude.degrees,
                after.latitude.degrees,
            )
            azimuth = compute_ground_sample_azimuth(track, target)
            difference = (azimuth - bearing + 180) % 360 - 180
            self.assertLess(abs(difference), 0.1)

    def test_yaw(self):
        """A yawed view turns its along-track axis to the left (counterclockwise)."""
        time = self.times[3]
        track = self.satellite.get_track(time)
        target = _subpoint(self.satellite, time)
        base = compute_ground_sample_azimuth(track, target)
        yawed = compute_ground_sample_azimuth(track, target, yaw_angle=30)
        self.assertAlmostEqual((base - yawed) % 360, 30, delta=0.01)

    def test_projection_matches_intersections(self):
        """The azimuth is that between the points where nearby rays meet the ground."""
        times = self.times[:5]
        track = self.satellite.get_track(times)
        for roll, pitch in ((30, 0), (-20, 15), (0, -25)):
            # targets at the center of the rolled and pitched view
            p_m, boresight, _, _ = _compute_view_frame(track, roll, pitch)
            center, _ = compute_ellipsoid_intersection(p_m, boresight, 0)
            lon, lat, _ = rectangular_to_geodetic(center)
            target = wgs84.latlon(lat, lon)
            azimuth = compute_ground_sample_azimuth(
                track, target, roll_angle=roll, pitch_angle=pitch
            )
            _, los, along, _ = _compute_ground_axes(
                track, target, roll_angle=roll, pitch_angle=pitch
            )
            half = np.radians(0.01)
            ends = [
                rectangular_to_geodetic(
                    compute_ellipsoid_intersection(
                        p_m, np.cos(half) * los + sign * np.sin(half) * along, 0
                    )[0]
                )
                for sign in (-1, 1)
            ]
            expected, _, _ = GEOD.inv(ends[0][0], ends[0][1], ends[1][0], ends[1][1])
            difference = (azimuth - expected + 180) % 360 - 180
            np.testing.assert_array_less(np.abs(difference), 0.05)

    def test_geostationary_inertial(self):
        """A geostationary satellite's along-track axis is east in the inertial frame."""
        satellite = Satellite(
            name="geo", orbit=GeosynchronousOrbit(longitude=-75, epoch=START)
        )
        track = satellite.get_track(START)
        target = wgs84.latlon(10, -70)
        azimuth = compute_ground_sample_azimuth(
            track, target, velocity_frame=VelocityFrame.INERTIAL
        )
        self.assertAlmostEqual(float(azimuth), 90, delta=2)


class TestSampleAzimuths(unittest.TestCase):
    """
    Unit tests for the azimuths of instruments' samples and `compute_sample_sizes`.
    """

    def setUp(self):
        self.end = START + timedelta(days=1)
        self.scanner = PointedInstrument(
            name="scanner",
            field_of_regard=105.6,
            cross_track_field_of_view=105.6,
            along_track_field_of_view=2.2,
            cross_track_pixels=96,
            is_rectangular=True,
        )
        self.conical = ConicalInstrument(
            name="conical",
            cone_angle=47.5,
            along_track_field_of_view=0.82,
            scan_half_width=75,
            beam_width=1.8,
        )
        self.satellite = Satellite(
            name="sat",
            orbit=SunSynchronousOrbit(
                altitude=830e3, equator_crossing_time="13:30", epoch=START
            ),
            instruments=[self.scanner, self.conical],
        )

    def test_sample_sizes_azimuth(self):
        """Observations have the azimuth of their samples, near that of the ground track."""
        points = [Point(-100, 35), Point(-60, -10), Point(20, 60)]
        observations = collect_observations(points, self.satellite, START, self.end)
        sizes = compute_sample_sizes(observations, self.satellite)
        self.assertGreater(len(sizes), 0)
        self.assertTrue(sizes.sample_azimuth.between(0, 360).all())
        # against the bearing of the ground track at each observation
        for row in sizes.itertuples():
            before = _subpoint(self.satellite, row.epoch - timedelta(seconds=0.5))
            after = _subpoint(self.satellite, row.epoch + timedelta(seconds=0.5))
            bearing, _, _ = GEOD.inv(
                before.longitude.degrees,
                before.latitude.degrees,
                after.longitude.degrees,
                after.latitude.degrees,
            )
            difference = (row.sample_azimuth - bearing + 180) % 360 - 180
            # the along-track axis turns from the ground track off nadir
            self.assertLess(abs(difference), 15)

    def test_conical_instrument(self):
        """A conical instrument's beams have azimuths too."""
        track = self.satellite.get_track(START)
        target = _subpoint(self.satellite, START)
        instrument = self.satellite.get_instrument(1)
        azimuth = instrument.compute_sample_azimuth(track, target)
        self.assertTrue(0 <= float(azimuth) < 360)
        self.assertAlmostEqual(
            float(azimuth),
            float(compute_ground_sample_azimuth(track, target)),
            places=6,
        )
