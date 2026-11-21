"""
Unit tests for the RepeatGroundTrackOrbit schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import TypeAdapter, ValidationError

from tatc.schemas import (
    GeneralPerturbationsOrbit,
    RepeatGroundTrackOrbit,
    Satellite,
)
from tatc.schemas.orbit import AllOrbits, AscendingNode
from tatc.schemas.surface import Point

# NOAA-20 element set (NOAA's operational TLE-AUX, 2026-10-02)
NOAA_20_TLE = [
    "1 43013U 17073A   26275.00000000  .00000083  00000-0  39562-4 0  2455",
    "2 43013 098.7849 213.5412 0000651 045.7112 319.4843 14.19528658459636",
]
# ICESat-2 element set (2026-10-31), on its 91-day repeat ground track
ICESAT_2_TLE = [
    "1 43613U 18070A   26304.15488612  .00001662  00000+0  60207-4 0  9996",
    "2 43613  91.9951 327.6490 0002117  21.1753 338.9583 15.28276307453220",
]
# an ascending node of the NOAA-20 element set, as propagated
NOAA_20_NODE = AscendingNode(
    time=datetime(2026, 10, 2, 1, 40, 4, 139002, tzinfo=timezone.utc),
    longitude=177.795746,
)


def _nearest_node(orbit, time):
    """Gets the ascending node of an orbit nearest a time, within an hour."""
    nodes = orbit.get_ascending_nodes(
        time - timedelta(hours=1), time + timedelta(hours=1)
    )
    return nodes.iloc[int(np.argmin(np.abs((nodes.time - time).dt.total_seconds())))]


class TestRepeatGroundTrackOrbit(unittest.TestCase):
    """
    Unit tests for the RepeatGroundTrackOrbit schema.
    """

    def setUp(self):
        self.orbit = RepeatGroundTrackOrbit(
            revolutions=227, days=16, ascending_node=NOAA_20_NODE
        )

    def test_good_data(self):
        """
        Test that the orbit initializes and serializes with its discriminator.
        """
        self.assertEqual(self.orbit.type, "repeat")
        self.assertTrue(self.orbit.is_sun_synchronous())
        data = Satellite(name="NOAA-20", orbit=self.orbit).model_dump(mode="json")
        orbit = TypeAdapter(AllOrbits).validate_python(data["orbit"])
        self.assertEqual(orbit, self.orbit)

    def test_bad_repeats(self):
        """
        Test that repeats beyond geosynchronous orbit or below the Earth's
        surface are rejected, as are unknown fields.
        """
        with self.assertRaises(ValidationError):
            RepeatGroundTrackOrbit(revolutions=1, days=2, ascending_node=NOAA_20_NODE)
        with self.assertRaises(ValidationError):
            RepeatGroundTrackOrbit(revolutions=20, days=1, ascending_node=NOAA_20_NODE)
        with self.assertRaises(ValidationError):
            RepeatGroundTrackOrbit(
                revolutions=227,
                days=16,
                ascending_node=NOAA_20_NODE,
                epoch=NOAA_20_NODE.time,
            )

    def test_derived_values(self):
        """
        Test that the altitude, inclination, and local time of the ascending
        node are derived as those of NOAA-20.
        """
        tle = GeneralPerturbationsOrbit.from_tle(NOAA_20_TLE)
        self.assertAlmostEqual(
            self.orbit.get_mean_altitude(), tle.get_mean_altitude(), delta=1000
        )
        # NOAA-20's inclination is biased to steer its local time
        self.assertAlmostEqual(
            self.orbit.get_inclination(), tle.get_inclination(), delta=0.1
        )
        self.assertAlmostEqual(self.orbit.get_local_time_ascending_node(), 13.52, 2)
        self.assertEqual(self.orbit.get_eccentricity(), 0)
        self.assertEqual(self.orbit.get_epoch(), NOAA_20_NODE.time)

    def test_crosses_its_node(self):
        """
        Test that the satellite crosses the ascending node, as propagated.
        """
        node = _nearest_node(self.orbit.to_gp_orbit(), NOAA_20_NODE.time)
        self.assertAlmostEqual(
            (node.time - NOAA_20_NODE.time).total_seconds(), 0, delta=1e-3
        )
        self.assertAlmostEqual(node.longitude, NOAA_20_NODE.longitude, delta=1e-6)

    def test_agrees_with_element_set_near_node(self):
        """
        Test that the nodes agree with those of the element set over two days
        (with drag and within NOAA-20's maintenance).
        """
        start = NOAA_20_NODE.time
        tle = GeneralPerturbationsOrbit.from_tle(NOAA_20_TLE).get_ascending_nodes(
            start, start + timedelta(days=2)
        )
        model = self.orbit.to_gp_orbit().get_ascending_nodes(
            start, start + timedelta(days=2)
        )
        self.assertEqual(len(tle), len(model))
        delta_time = (tle.time - model.time).dt.total_seconds()
        delta_longitude = (tle.longitude - model.longitude + 180) % 360 - 180
        self.assertLess(np.max(np.abs(delta_time)), 5)
        self.assertLess(np.max(np.abs(delta_longitude)), 0.01)

    def test_repeats_years_from_node(self):
        """
        Test that the satellite crosses the same ascending node whole repeat
        cycles later (within the model, whose cycles are of Terrestrial Time),
        up to five years, by its track and by its observation events.
        """
        gp_orbit = self.orbit.to_gp_orbit()
        point = Point(id=0, latitude=33.4, longitude=-111.9)
        start = NOAA_20_NODE.time
        reference_times, reference_codes = gp_orbit.get_observation_events(
            point, start, start + timedelta(days=1), 10
        )
        for days in (16, 16 * 23, 16 * 114):
            time = start + timedelta(days=days)
            node = _nearest_node(gp_orbit, time)
            self.assertAlmostEqual((node.time - time).total_seconds(), 0, delta=0.01)
            self.assertAlmostEqual(node.longitude, NOAA_20_NODE.longitude, delta=1e-4)
            times, codes = gp_orbit.get_observation_events(
                point, time, time + timedelta(days=1), 10
            )
            np.testing.assert_array_equal(codes, reference_codes)
            np.testing.assert_allclose(
                (times.tt - reference_times.tt) * 86400, days * 86400, atol=0.01
            )

    def test_inclined_repeat(self):
        """
        Test that an inclined (not sun-synchronous) repeat with a long cycle,
        as of ICESat-2 (1387 revolutions in 91 days), completes its
        revolutions in its nodal days.
        """
        orbit = RepeatGroundTrackOrbit(
            revolutions=1387, days=91, inclination=91.9951, ascending_node=NOAA_20_NODE
        )
        self.assertFalse(orbit.is_sun_synchronous())
        self.assertAlmostEqual(orbit.get_inclination(), 91.9951, delta=1e-9)
        element = orbit.to_gp_orbit().elements[0]
        nodal_period, nodal_day = element.get_nodal_period_and_day()
        self.assertAlmostEqual(1387 * nodal_period, 91 * nodal_day, delta=1e-3)
        tle = GeneralPerturbationsOrbit.from_tle(ICESAT_2_TLE)
        self.assertAlmostEqual(
            orbit.get_mean_altitude(), tle.get_mean_altitude(), delta=100
        )

    def test_derived_orbit(self):
        """
        Test that a derived orbit is ahead in mean anomaly in a plane rotated
        by the right ascension of ascending node, at the same time.
        """
        for delta_mean_anomaly, delta_raan in ((90, 0), (0, 30), (-45, -20)):
            derived = self.orbit.get_derived_orbit(delta_mean_anomaly, delta_raan)
            time = NOAA_20_NODE.time + timedelta(hours=3)
            original = self.orbit.to_gp_orbit().get_orbit_track(time)
            other = derived.to_gp_orbit().get_orbit_track(time)
            normal = np.cross(original.position.m, original.velocity.m_per_s)
            other_normal = np.cross(other.position.m, other.velocity.m_per_s)
            # right ascension of each orbit normal (90 degrees from its node)
            raan = np.degrees(np.arctan2(normal[0], -normal[1]))
            other_raan = np.degrees(np.arctan2(other_normal[0], -other_normal[1]))
            self.assertAlmostEqual(
                (other_raan - raan - delta_raan + 180) % 360 - 180, 0, delta=0.01
            )
            # angle between the positions, within the (nearly) common plane
            angle = np.degrees(
                np.arctan2(
                    np.dot(np.cross(original.position.m, other.position.m), normal)
                    / np.linalg.norm(normal),
                    np.dot(original.position.m, other.position.m),
                )
            )
            if delta_raan == 0:
                # within the short-period (J2) oscillation of the osculating
                # argument of latitude, as each orbit is pinned at its node
                self.assertAlmostEqual(angle, delta_mean_anomaly, delta=0.15)


if __name__ == "__main__":
    unittest.main()
