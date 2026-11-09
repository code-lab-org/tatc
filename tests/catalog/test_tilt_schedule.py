"""
Unit tests for tilt schedules of catalog configurations: a pitch of the body that
changes by day, as OCI's tilt does.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import date, datetime, timedelta, timezone

import numpy as np
from pydantic import ValidationError

from tatc.catalog import TiltSchedule, get_entry, get_subsolar_latitude
from tatc.schemas import Attitude, AttitudeSchedule, GeneralPerturbationsOrbit
from tatc.utils import VelocityFrame
from tatc.utils.orbital import compute_argument_of_latitude

PACE_TLE = [
    "1 58928U 24025A   26313.27483190  .00000501  00000+0  10091-3 0  9993",
    "2 58928  98.1040 242.4857 0002666 129.8006 230.3433 14.65644467147180",
]
MEASURED = {
    "kind": "measured",
    "source": "a test",
    "platform": "pace",
    "period": "2026-10-04",
}


def _schedule(**changes) -> TiltSchedule:
    values = {
        "fore_pitch_angle": 20.0,
        "aft_pitch_angle": -20.0,
        "odd_day_offset": 5.0,
        "even_day_offset": -1.0,
        "transition": 2.0,
        "return_argument_of_latitude": 180,
        "blind_before": 15.0,
        "blind_after": 20.0,
    }
    values.update(changes)
    return TiltSchedule(**values, provenance={key: MEASURED for key in values})


class TestSubsolarLatitude(unittest.TestCase):
    """
    Unit tests for the subsolar latitude.
    """

    def test_solstices_and_equinox(self):
        """
        Test that the subsolar latitude is the obliquity at the solstices and near
        zero at an equinox.
        """
        june = get_subsolar_latitude(datetime(2026, 6, 21, 12, tzinfo=timezone.utc))
        december = get_subsolar_latitude(
            datetime(2026, 12, 21, 12, tzinfo=timezone.utc)
        )
        march = get_subsolar_latitude(datetime(2026, 3, 20, 15, tzinfo=timezone.utc))
        self.assertAlmostEqual(june, 23.44, delta=0.01)
        self.assertAlmostEqual(december, -23.44, delta=0.01)
        self.assertAlmostEqual(march, 0, delta=0.01)


class TestTiltSchedule(unittest.TestCase):
    """
    Unit tests for TiltSchedule.
    """

    def setUp(self):
        self.orbit = GeneralPerturbationsOrbit.from_tle(PACE_TLE)
        self.schedule = _schedule()

    def test_needs_provenance(self):
        """
        Test that each value needs its provenance.
        """
        with self.assertRaises(ValidationError):
            TiltSchedule(
                fore_pitch_angle=20,
                aft_pitch_angle=-20,
                odd_day_offset=5,
                even_day_offset=-1,
                transition=2,
                return_argument_of_latitude=180,
                blind_before=15,
                blind_after=20,
                provenance={"fore_pitch_angle": MEASURED},
            )

    def test_offset_alternates_by_day(self):
        """
        Test that the latitude of the change follows the subsolar latitude at noon,
        with the offset of the parity of the day of the year.
        """
        for day, offset in ((date(2026, 11, 7), 5.0), (date(2026, 11, 8), -1.0)):
            noon = datetime(day.year, day.month, day.day, 12, tzinfo=timezone.utc)
            self.assertAlmostEqual(
                self.schedule.get_switch_latitude(day),
                get_subsolar_latitude(noon) + offset,
            )

    def test_profile_changes_after_ascending_node(self):
        """
        Test that the profile changes from aft to fore over the transition, centered
        where the orbit passes over the latitude of the change, and back at the
        return argument of latitude.
        """
        day = date(2026, 8, 15)
        profile = self.schedule.get_pitch_angle_profile(day, self.orbit)
        self.assertEqual(len(profile), 4)
        arguments = [argument for argument, _ in profile]
        self.assertEqual(arguments, sorted(arguments))
        angles = dict(profile)
        fore = [a for a, angle in angles.items() if angle == 20.0]
        aft = [a for a, angle in angles.items() if angle == -20.0]
        self.assertAlmostEqual(min(fore) - max(a for a in aft if a < 90), 2.0)
        center = (min(fore) + max(a for a in aft if a < 90)) / 2
        latitude = np.degrees(
            np.arcsin(
                np.sin(np.radians(center))
                * np.sin(np.radians(self.orbit.get_inclination()))
            )
        )
        # the geocentric latitude of the change is below its geodetic latitude
        self.assertLess(latitude, self.schedule.get_switch_latitude(day))
        self.assertAlmostEqual(
            latitude, self.schedule.get_switch_latitude(day), delta=0.15
        )
        self.assertAlmostEqual(max(fore), 179.0)
        self.assertIn(181.0, angles)

    def test_blind_around_changes(self):
        """
        Test that the instruments do not observe from before to after the middle of
        each change, at the orbit's mean rate, and that no blind time gives no
        intervals.
        """
        day = date(2026, 8, 15)
        profile = self.schedule.get_pitch_angle_profile(day, self.orbit)
        # the change to fore ends a half transition (1 deg) after its middle
        switch = min(a for a, angle in profile if angle == 20.0) - 1.0
        rate = 360 / self.orbit.get_orbit_period().total_seconds()
        intervals = self.schedule.get_blind_intervals(day, self.orbit)
        np.testing.assert_allclose(
            intervals,
            [
                (switch - 15 * rate, switch + 20 * rate),
                (180 - 15 * rate, 180 + 20 * rate),
            ],
        )
        none = _schedule(blind_before=0.0, blind_after=0.0)
        self.assertIsNone(none.get_blind_intervals(day, self.orbit))

    def test_schedule_has_a_segment_a_day(self):
        """
        Test that a schedule has an attitude for each day of the analysis, from
        midnight (UTC), with the frame and other fields of the configuration.
        """
        start = datetime(2026, 11, 7, 18, tzinfo=timezone.utc)
        end = datetime(2026, 11, 9, 6, tzinfo=timezone.utc)
        schedule = self.schedule.build_schedule(
            self.orbit, start, end, VelocityFrame.EARTH_FIXED, {"yaw_angle": 0.1}
        )
        self.assertIsInstance(schedule, AttitudeSchedule)
        self.assertEqual(
            [s.start for s in schedule.segments],
            [datetime(2026, 11, d, tzinfo=timezone.utc) for d in (7, 8, 9)],
        )
        for segment in schedule.segments:
            day = segment.start.date()
            self.assertEqual(segment.attitude.velocity_frame, VelocityFrame.EARTH_FIXED)
            self.assertEqual(segment.attitude.yaw_angle, 0.1)
            self.assertEqual(
                segment.attitude.pitch_angle_profile,
                self.schedule.get_pitch_angle_profile(day, self.orbit),
            )
            self.assertEqual(
                segment.attitude.blind_intervals,
                self.schedule.get_blind_intervals(day, self.orbit),
            )
        with self.assertRaisesRegex(ValueError, "precedes"):
            self.schedule.build_schedule(
                self.orbit, end, start, VelocityFrame.EARTH_FIXED, {}
            )

    def test_body_pitch_follows_schedule(self):
        """
        Test that a satellite's body is pitched fore after the change and aft before
        it, on each day.
        """
        day = datetime(2026, 11, 8, tzinfo=timezone.utc)
        schedule = self.schedule.build_schedule(
            self.orbit, day, day, VelocityFrame.EARTH_FIXED, {}
        )
        times = [day + timedelta(seconds=s) for s in range(0, 86400, 30)]
        track = self.orbit.get_orbit_track(times)
        argument = compute_argument_of_latitude(
            track.position.m, track.velocity.m_per_s
        )
        pitch = schedule.get_body_angles(track)[2]
        profile = self.schedule.get_pitch_angle_profile(day.date(), self.orbit)
        switch = profile[-1][0]
        fore = (argument > switch) | (argument < 179)
        aft = (argument > 181) & (argument < switch - 2)
        np.testing.assert_allclose(pitch[fore], 20.0)
        np.testing.assert_allclose(pitch[aft], -20.0)

    def test_unreachable_latitude(self):
        """
        Test that a latitude the orbit does not reach is an error.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 39574U 14009C   26277.54227627  .00004528  00000+0  10004-3 0  9997",
                "2 39574  10.0000 301.6724 0010793 268.6572  91.3332 15.45217995714704",
            ]
        )
        with self.assertRaisesRegex(ValueError, "does not reach"):
            _schedule(odd_day_offset=20.0).get_pitch_angle_profile(
                date(2026, 6, 21), orbit
            )


class TestOciSchedule(unittest.TestCase):
    """
    Unit tests for the OCI entry's scheduled tilt.
    """

    def setUp(self):
        self.entry = get_entry("oci")
        self.orbit = GeneralPerturbationsOrbit.from_tle(PACE_TLE)

    def test_needs_period(self):
        """
        Test that building the scheduled configuration needs the period of the
        analysis, unless an attitude replaces the schedule.
        """
        with self.assertRaisesRegex(ValueError, "start and end"):
            self.entry.build(orbit=self.orbit, time_step=10)
        payload = self.entry.build(
            orbit=self.orbit,
            time_step=10,
            attitude=Attitude(velocity_frame=VelocityFrame.EARTH_FIXED),
        )
        self.assertIsInstance(payload.attitude, Attitude)

    def test_builds_schedule(self):
        """
        Test that the scheduled configuration builds a schedule of attitudes, and
        an instrument without a pitch of its own.
        """
        start = datetime(2026, 11, 7, tzinfo=timezone.utc)
        satellite = self.entry.to_satellite(
            self.orbit, time_step=10, start=start, end=start + timedelta(days=2)
        )
        self.assertIsInstance(satellite.attitude, AttitudeSchedule)
        self.assertEqual(len(satellite.attitude.segments), 3)
        # each change of tilt lies within an interval in which OCI does not observe
        for segment in satellite.attitude.segments:
            profile = segment.attitude.pitch_angle_profile
            for (first, before), (second, after) in zip(profile, profile[1:]):
                if before != after:
                    self.assertTrue(
                        any(
                            start < first and second < end
                            for start, end in segment.attitude.blind_intervals
                        )
                    )
        self.assertIsNone(satellite.instruments[0].pitch_angle_profile)
        self.assertEqual(satellite.instruments[0].pitch_angle, 0)

    def test_schedule_is_fitted_and_digested(self):
        """
        Test that the schedule's measured values make the configuration fitted, and
        that changing a value of the schedule changes the digest.
        """
        configuration = self.entry.configurations["swath"]
        self.assertTrue(configuration.is_fitted)
        changed = configuration.model_copy(
            update={
                "attitude_schedule": configuration.attitude_schedule.model_copy(
                    update={"odd_day_offset": 5.0}
                )
            }
        )
        self.assertNotEqual(changed.digest, configuration.digest)

    def test_schedule_sets_pitch(self):
        """
        Test that a configuration with a schedule cannot also set a pitch.
        """
        data = self.entry.configurations["swath"].model_dump(mode="json")
        for attitude in ({"pitch_angle": 1.0}, {"blind_intervals": [(0, 1)]}):
            data["attitude"] = attitude
            with self.assertRaisesRegex(ValidationError, "sets the pitch"):
                type(self.entry.configurations["swath"])(**data)


if __name__ == "__main__":
    unittest.main()
