"""
Unit tests for the attitude schemas in tatc.schemas.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import ValidationError
from shapely.geometry import Point as ShapelyPoint
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc.analysis import collect_ground_track, collect_observations
from tatc.schemas import (
    Attitude,
    AttitudeSchedule,
    AttitudeSegment,
    CircularOrbit,
    ConicalInstrument,
    GlintAttitude,
    Point,
    PointedInstrument,
    Satellite,
    SunSynchronousOrbit,
    TargetAttitude,
    WalkerConstellation,
)
from tatc.utils import VelocityFrame, compute_glint_angle, geodesic_distance
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


class TestTargetAttitude(unittest.TestCase):
    """
    Unit tests for the TargetAttitude schema.
    """

    def setUp(self):
        self.start = datetime(2026, 1, 1, tzinfo=timezone.utc)
        self.orbit = SunSynchronousOrbit(
            altitude=700e3, equator_crossing_time="10:30", epoch=self.start
        )
        self.times = [self.start + timedelta(minutes=30 + 0.5 * i) for i in range(12)]
        self.track = self.orbit.to_gp_orbit().get_orbit_track_at_time(
            _to_time(self.times)
        )
        subpoint = wgs84.subpoint_of(self.track)
        self.subpoint = subpoint
        # a target about 250 km east of the ground track near the middle time
        self.target = Point(
            latitude=float(subpoint.latitude.degrees[6]),
            longitude=float(subpoint.longitude.degrees[6]) + 2.5,
        )
        self.instrument = PointedInstrument(
            name="p",
            field_of_regard=80,
            cross_track_field_of_view=1,
            along_track_field_of_view=1,
        )

    def _satellite(self, attitude):
        return Satellite(
            name="s", orbit=self.orbit, instruments=[self.instrument], attitude=attitude
        )

    def test_validation(self):
        """
        Test that a target attitude needs targets and a maximum off-nadir
        angle in (0, 90).
        """
        with self.assertRaises(ValidationError):
            TargetAttitude(targets=[])
        with self.assertRaises(ValidationError):
            TargetAttitude(targets=[self.target], max_off_nadir_angle=90)

    def test_points_at_target_within_reach(self):
        """
        Test that the body points at the target while it is within the
        maximum off-nadir angle (the view center on the target), and along
        nadir otherwise, choosing the target closest to nadir.
        """
        far = Point(
            latitude=-self.target.latitude, longitude=self.target.longitude + 120
        )
        attitude = TargetAttitude(targets=[far, self.target], max_off_nadir_angle=40)
        index = attitude.get_target_index(self.track)
        self.assertTrue(np.all(np.isin(index, [-1, 1])))
        self.assertGreater(np.sum(index == 1), 0)
        self.assertGreater(np.sum(index == -1), 0)
        center = (
            self._satellite(attitude)
            .get_instrument()
            .compute_footprint_center(self.track)
        )
        n = len(self.times)
        to_target = geodesic_distance(
            np.full(n, self.target.longitude),
            np.full(n, self.target.latitude),
            center.longitude.degrees,
            center.latitude.degrees,
        )
        to_nadir = geodesic_distance(
            self.subpoint.longitude.degrees,
            self.subpoint.latitude.degrees,
            center.longitude.degrees,
            center.latitude.degrees,
        )
        np.testing.assert_array_less(to_target[index == 1], 1)
        np.testing.assert_array_less(to_nadir[index == -1], 1)

    def test_yaw_keeps_pointing(self):
        """
        Test that a yaw of the body does not move a boresight view off the
        target, but mirrors the views of instruments rolled from it.
        """
        attitude = TargetAttitude(
            targets=[self.target], max_off_nadir_angle=40, yaw_angle=180
        )
        index = attitude.get_target_index(self.track)
        center = (
            self._satellite(attitude)
            .get_instrument()
            .compute_footprint_center(self.track)
        )
        pointing = index == 0
        n = int(np.sum(pointing))
        distance = geodesic_distance(
            np.full(n, self.target.longitude),
            np.full(n, self.target.latitude),
            center.longitude.degrees[pointing],
            center.latitude.degrees[pointing],
        )
        np.testing.assert_array_less(distance, 1)

    def test_observes_off_track_target(self):
        """
        Test that a narrow instrument observes a target off its ground track
        only while the body points at it.
        """
        end = self.start + timedelta(hours=2)
        point = ShapelyPoint(self.target.longitude, self.target.latitude)
        self.assertEqual(
            len(collect_observations(point, self._satellite(None), self.start, end)), 0
        )
        observations = collect_observations(
            point,
            self._satellite(
                TargetAttitude(targets=[self.target], max_off_nadir_angle=40)
            ),
            self.start,
            end,
        )
        self.assertEqual(len(observations), 1)

    def test_scheduled_and_serialized(self):
        """
        Test that a schedule may switch to target pointing, and that a
        satellite with it round-trips through JSON.
        """
        target = TargetAttitude(targets=[self.target], max_off_nadir_angle=40)
        schedule = AttitudeSchedule.from_attitudes(
            [self.start, self.times[5]], [Attitude(), target]
        )
        _, roll, _, _ = schedule.get_body_angles(self.track)
        index = target.get_target_index(self.track)
        self.assertTrue(np.all(roll[:5] == 0))
        self.assertTrue(np.any(roll[5:][index[5:] == 0] != 0))
        satellite = self._satellite(schedule)
        copy = Satellite.model_validate_json(satellite.model_dump_json())
        self.assertEqual(copy.attitude, schedule)
        self.assertIsInstance(copy.attitude.segments[1].attitude, TargetAttitude)


class TestGlintAttitude(unittest.TestCase):
    """
    Unit tests for the GlintAttitude schema.
    """

    def setUp(self):
        self.start = datetime(2026, 3, 20, 12, tzinfo=timezone.utc)
        self.orbit = CircularOrbit(altitude=705e3, inclination=98.2, epoch=self.start)
        self.times = [self.start + timedelta(minutes=m) for m in range(0, 100, 4)]
        self.track = self.orbit.to_gp_orbit().get_orbit_track_at_time(
            _to_time(self.times)
        )
        self.instrument = PointedInstrument(
            name="p",
            field_of_regard=120,
            cross_track_field_of_view=1,
            along_track_field_of_view=1,
            max_target_glint_angle=1,
        )

    def _satellite(self, attitude):
        return Satellite(
            name="s", orbit=self.orbit, instruments=[self.instrument], attitude=attitude
        )

    def test_validation(self):
        """
        Test the bounds of the maximum off-nadir and minimum solar elevation
        angles.
        """
        with self.assertRaises(ValidationError):
            GlintAttitude(max_off_nadir_angle=90)
        with self.assertRaises(ValidationError):
            GlintAttitude(min_solar_elevation=-1)

    def test_points_at_glint(self):
        """
        Test that the view center lies on the specular point (a zero glint
        angle) while the body points at it, and at nadir otherwise.
        """
        attitude = GlintAttitude(max_off_nadir_angle=50)
        point = attitude.get_aim_point(self.track)
        pointing = np.isfinite(point.latitude.degrees)
        self.assertGreater(np.sum(pointing), 0)
        self.assertGreater(np.sum(~pointing), 0)
        center = (
            self._satellite(attitude).get_instrument().compute_footprint_center(self.track)
        )
        np.testing.assert_array_less(
            compute_glint_angle(self.track, center)[pointing], 1e-3
        )
        subpoint = wgs84.subpoint_of(self.track)
        to_nadir = geodesic_distance(
            subpoint.longitude.degrees,
            subpoint.latitude.degrees,
            center.longitude.degrees,
            center.latitude.degrees,
        )
        np.testing.assert_array_less(to_nadir[~pointing], 1)
        # a smaller maximum off-nadir angle or a larger minimum solar elevation points less
        fewer = GlintAttitude(max_off_nadir_angle=20).get_aim_point(self.track)
        self.assertLess(np.sum(np.isfinite(fewer.latitude.degrees)), np.sum(pointing))
        higher = GlintAttitude(max_off_nadir_angle=50, min_solar_elevation=60)
        higher_point = higher.get_aim_point(self.track)
        self.assertLess(np.sum(np.isfinite(higher_point.latitude.degrees)), np.sum(pointing))

    def test_offset_toward_nadir(self):
        """
        Test that an offset angle points the body short of the specular point
        toward nadir, on the line from nadir to it, by the offset (constant
        or from a profile of the specular point's off-nadir angle).
        """
        glint = GlintAttitude(max_off_nadir_angle=80)
        specular = glint.get_aim_point(self.track)
        pointing = np.flatnonzero(np.isfinite(specular.latitude.degrees))
        track = self.track[pointing]
        specular = wgs84.latlon(
            specular.latitude.degrees[pointing], specular.longitude.degrees[pointing]
        )
        subpoint = wgs84.subpoint_of(track)

        def off_nadir(point):
            # the angle from nadir: 90 deg less the angle of the line of sight below the horizon
            return np.degrees(
                np.arccos(
                    np.sum(
                        np.array((point.at(track.t) - track).frame_xyz(itrs).m)
                        * -np.array(subpoint.itrs_xyz.m)
                        / np.linalg.norm(np.array(subpoint.itrs_xyz.m), axis=0),
                        axis=0,
                    )
                    / np.linalg.norm(
                        np.array((point.at(track.t) - track).frame_xyz(itrs).m), axis=0
                    )
                )
            )

        glint_angle = off_nadir(specular)
        for attitude, expected in [
            (GlintAttitude(max_off_nadir_angle=80, offset_angle=5), glint_angle - 5),
            (
                GlintAttitude(
                    max_off_nadir_angle=80, offset_angle_profile=[(60, 10), (20, 2)]
                ),
                glint_angle - np.interp(glint_angle, [20, 60], [2, 10]),
            ),
        ]:
            aim = attitude.get_aim_point(track)
            # (at nadir where the offset exceeds the specular point's angle from nadir)
            np.testing.assert_allclose(off_nadir(aim), np.maximum(expected, 0), atol=0.3)
            # on the line from nadir to the specular point
            detour = (
                geodesic_distance(
                    subpoint.longitude.degrees, subpoint.latitude.degrees,
                    aim.longitude.degrees, aim.latitude.degrees,
                )
                + geodesic_distance(
                    aim.longitude.degrees, aim.latitude.degrees,
                    specular.longitude.degrees, specular.latitude.degrees,
                )
                - geodesic_distance(
                    subpoint.longitude.degrees, subpoint.latitude.degrees,
                    specular.longitude.degrees, specular.latitude.degrees,
                )
            )
            np.testing.assert_array_less(detour, 100)
        self.assertEqual(
            GlintAttitude(offset_angle_profile=[(60, 10), (20, 2)]).offset_angle_profile,
            [(20, 2), (60, 10)],
        )
        with self.assertRaises(ValidationError):
            GlintAttitude(offset_angle_profile=[(95, 1)])

    def test_glint_observations(self):
        """
        Test that an instrument requiring a small glint angle observes the
        specular point only with a glint attitude, and that its observations
        report the glint angle.
        """
        attitude = GlintAttitude(max_off_nadir_angle=50)
        k = int(np.flatnonzero(np.isfinite(attitude.get_aim_point(self.track).latitude.degrees))[2])
        point = attitude.get_aim_point(self.track[k])
        target = ShapelyPoint(float(point.longitude.degrees), float(point.latitude.degrees))
        start, end = self.times[k] - timedelta(minutes=10), self.times[k] + timedelta(minutes=10)
        self.assertEqual(len(collect_observations(target, self._satellite(None), start, end)), 0)
        observations = collect_observations(
            target, self._satellite(attitude), start, end, omit_solar=False
        )
        self.assertGreater(len(observations), 0)
        np.testing.assert_array_less(observations.glint_angle, 1)

    def test_scheduled_and_serialized(self):
        """
        Test that a schedule may switch to glint pointing, and that a
        satellite with it round-trips through JSON.
        """
        schedule = AttitudeSchedule.from_attitudes(
            [self.start, self.times[10]], [Attitude(), GlintAttitude()]
        )
        satellite = self._satellite(schedule)
        copy = Satellite.model_validate_json(satellite.model_dump_json())
        self.assertEqual(copy.attitude, schedule)
        self.assertIsInstance(copy.attitude.segments[1].attitude, GlintAttitude)
        self.assertIsInstance(
            Satellite.model_validate_json(self._satellite(GlintAttitude()).model_dump_json()).attitude,
            GlintAttitude,
        )
