"""
Unit tests for the attitude schemas in tatc.schemas.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import ValidationError
from shapely.geometry import Point as ShapelyPoint

from tatc.analysis import collect_ground_track, collect_observations
from tatc.schemas import (
    Attitude,
    AttitudeSchedule,
    AttitudeSegment,
    CircularOrbit,
    ConicalInstrument,
    PointedInstrument,
    Satellite,
    SunSynchronousOrbit,
    WalkerConstellation,
)
from tatc.utils import VelocityFrame
from tatc.utils.time import _to_time


class TestAttitude(unittest.TestCase):
    """
    Unit tests for the Attitude and AttitudeSchedule schemas.
    """

    def setUp(self):
        self.start = datetime(2026, 1, 1, tzinfo=timezone.utc)
        self.orbit = SunSynchronousOrbit(
            altitude=700e3, equator_crossing_time="10:30", epoch=self.start
        )
        self.times = [self.start + timedelta(minutes=13 * i) for i in range(7)]
        self.track = self.orbit.to_gp_orbit().get_orbit_track_at_time(
            _to_time(self.times)
        )

    def test_angle_bounds_and_profiles(self):
        """
        Test that angles must lie in [-180, 180] and that profiles are
        checked and sorted by argument of latitude.
        """
        with self.assertRaises(ValidationError):
            Attitude(yaw_angle=190)
        with self.assertRaises(ValidationError):
            Attitude(roll_angle_profile=[(360, 0)])
        with self.assertRaises(ValidationError):
            Attitude(pitch_angle_profile=[(0, 200)])
        attitude = Attitude(yaw_angle_profile=[(180, 5), (0, -5)])
        self.assertEqual(attitude.yaw_angle_profile, [(0, -5), (180, 5)])

    def test_schedule_segments(self):
        """
        Test that a schedule sorts its segments, requires timezone-aware
        start times, and applies the first segment before its start.
        """
        flip = self.start + timedelta(hours=1)
        schedule = AttitudeSchedule.from_attitudes(
            [flip, self.start], [Attitude(yaw_angle=0), Attitude(yaw_angle=180)]
        )
        self.assertEqual([s.start for s in schedule.segments], [self.start, flip])
        self.assertEqual(
            schedule.get_attitude(self.start - timedelta(days=1)).yaw_angle, 180
        )
        self.assertEqual(schedule.get_attitude(flip).yaw_angle, 0)
        with self.assertRaises(ValidationError):
            AttitudeSegment(start=datetime(2026, 1, 1), attitude=Attitude())
        with self.assertRaises(ValueError):
            AttitudeSchedule.from_attitudes([self.start], [])

    def test_body_angles_of_schedule(self):
        """
        Test that a schedule gives each time the angles and velocity frame of
        its segment, for a vector or a single time.
        """
        flip = self.start + timedelta(minutes=40)
        schedule = AttitudeSchedule.from_attitudes(
            [self.start, flip],
            [
                Attitude(yaw_angle=180, velocity_frame=VelocityFrame.INERTIAL),
                Attitude(roll_angle=10),
            ],
        )
        yaw, roll, pitch, inertial = schedule.get_body_angles(self.track)
        after = np.array([t >= flip for t in self.times])
        np.testing.assert_array_equal(yaw, np.where(after, 0, 180))
        np.testing.assert_array_equal(roll, np.where(after, 10, 0))
        np.testing.assert_array_equal(pitch, 0)
        np.testing.assert_array_equal(inertial, ~after)
        single = self.orbit.to_gp_orbit().get_orbit_track_at_time(_to_time(self.start))
        yaw, _, _, inertial = schedule.get_body_angles(single)
        self.assertEqual(np.shape(yaw), ())
        self.assertTrue(inertial)

    def _center(self, satellite, index=0):
        center = satellite.get_instrument(index).compute_footprint_center(self.track)
        return np.array([center.latitude.degrees, center.longitude.degrees])

    def _assert_same_centers(self, a, b):
        np.testing.assert_allclose(self._center(a), self._center(b), atol=1e-9)

    def test_composition_with_instrument_pointing(self):
        """
        Test that instruments point relative to the body: an identity
        attitude changes nothing, an inertial attitude equals an instrument's
        inertial velocity frame, a yaw of 180 deg mirrors roll and pitch, and
        a body roll adds to an instrument's roll.
        """

        def instrument(**kwargs):
            return PointedInstrument(
                name="p",
                cross_track_field_of_view=2,
                along_track_field_of_view=2,
                **kwargs,
            )

        def satellite(attitude=None, **kwargs):
            return Satellite(
                name="s",
                orbit=self.orbit,
                instruments=[instrument(**kwargs)],
                attitude=attitude,
            )

        self._assert_same_centers(
            satellite(Attitude(), roll_angle=20, pitch_angle=5),
            satellite(roll_angle=20, pitch_angle=5),
        )
        self._assert_same_centers(
            satellite(Attitude(velocity_frame=VelocityFrame.INERTIAL), roll_angle=20),
            satellite(roll_angle=20, velocity_frame=VelocityFrame.INERTIAL),
        )
        self._assert_same_centers(
            satellite(Attitude(yaw_angle=180), roll_angle=20, pitch_angle=5),
            satellite(roll_angle=-20, pitch_angle=-5),
        )
        self._assert_same_centers(
            satellite(Attitude(roll_angle=15), roll_angle=10),
            satellite(roll_angle=25),
        )
        profile = [(0, 10), (90, 30), (180, 10), (270, -10)]
        self._assert_same_centers(
            satellite(Attitude(roll_angle_profile=profile)),
            satellite(roll_angle_profile=profile),
        )

    def test_conical_yaw(self):
        """
        Test that a yaw of 180 deg turns a forward conical scan aft.
        """

        def satellite(center, attitude=None):
            return Satellite(
                name="s",
                orbit=self.orbit,
                instruments=[
                    ConicalInstrument(
                        name="c",
                        cone_angle=45,
                        along_track_field_of_view=1,
                        scan_half_width=60,
                        scan_center_azimuth=center,
                    )
                ],
                attitude=attitude,
            )

        self._assert_same_centers(satellite(0, Attitude(yaw_angle=180)), satellite(180))

    def test_velocity_frame_conflict(self):
        """
        Test that an instrument may not set a velocity frame if the attitude
        sets the body's.
        """
        with self.assertRaises(ValidationError):
            Satellite(
                name="s",
                orbit=self.orbit,
                instruments=[
                    PointedInstrument(
                        name="p",
                        cross_track_field_of_view=2,
                        along_track_field_of_view=2,
                        velocity_frame=VelocityFrame.INERTIAL,
                    )
                ],
                attitude=Attitude(),
            )

    def test_serialization_and_constellations(self):
        """
        Test that a satellite with a schedule round-trips through JSON, and
        that the members of a constellation carry its attitude.
        """
        schedule = AttitudeSchedule.from_attitudes(
            [self.start, self.start + timedelta(days=40)],
            [Attitude(yaw_angle=180), Attitude(yaw_angle=0)],
        )
        satellite = Satellite(name="s", orbit=self.orbit, attitude=schedule)
        copy = Satellite.model_validate_json(satellite.model_dump_json())
        self.assertEqual(copy.attitude, schedule)
        self.assertIsNone(
            Satellite(name="s", orbit=self.orbit).get_instrument().platform_attitude
        )
        self.assertEqual(satellite.get_instrument().platform_attitude, schedule)
        constellation = WalkerConstellation(
            name="w",
            orbit=CircularOrbit(altitude=700e3, inclination=60, epoch=self.start),
            number_satellites=4,
            number_planes=2,
            attitude=Attitude(yaw_angle=90),
        )
        for member in constellation.generate_members():
            self.assertEqual(member.attitude.yaw_angle, 90)

    def test_observations_across_a_schedule(self):
        """
        Test that the observations of an analysis across a scheduled flip
        equal those of each attitude before and after it.
        """
        conical = ConicalInstrument(
            name="c", cone_angle=45, along_track_field_of_view=1, scan_half_width=60
        )
        flip = self.start + timedelta(hours=12)
        end = self.start + timedelta(days=1)
        points = [ShapelyPoint(-100, 35), ShapelyPoint(20, -10), ShapelyPoint(140, 50)]
        schedule = AttitudeSchedule.from_attitudes(
            [self.start, flip], [Attitude(yaw_angle=180), Attitude()]
        )

        def observe(attitude, start, stop):
            satellite = Satellite(
                name="s", orbit=self.orbit, instruments=[conical], attitude=attitude
            )
            return collect_observations(points, satellite, start, stop)

        scheduled = observe(schedule, self.start, end)
        before = observe(Attitude(yaw_angle=180), self.start, flip)
        after = observe(Attitude(), flip, end)
        self.assertGreater(len(before), 0)
        self.assertGreater(len(after), 0)
        expected = sorted(list(before.epoch) + list(after.epoch))
        self.assertEqual(len(scheduled), len(expected))
        for observed, wanted in zip(sorted(scheduled.epoch), expected):
            self.assertLess(abs((observed - wanted).total_seconds()), 1e-3)

    def test_ground_track_follows_yaw(self):
        """
        Test that the footprint of an instrument rolled to the left of the
        body moves to the right of the ground track while the body flies
        backward.
        """
        instrument = PointedInstrument(
            name="p",
            cross_track_field_of_view=5,
            along_track_field_of_view=1,
            roll_angle=30,
        )
        times = [self.start + timedelta(minutes=10)]
        tracks = {
            yaw: collect_ground_track(
                Satellite(
                    name="s",
                    orbit=self.orbit,
                    instruments=[instrument],
                    attitude=Attitude(yaw_angle=yaw),
                ),
                times,
            )
            .geometry.iloc[0]
            .centroid
            for yaw in (0, 180)
        }
        mirrored = (
            collect_ground_track(
                Satellite(
                    name="s",
                    orbit=self.orbit,
                    instruments=[instrument.model_copy(update={"roll_angle": -30})],
                ),
                times,
            )
            .geometry.iloc[0]
            .centroid
        )
        self.assertGreater(tracks[0].distance(tracks[180]), 5)
        self.assertLess(tracks[180].distance(mirrored), 1e-6)
