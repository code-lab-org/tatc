"""
Unit tests for the KeplerianOrbit schema.
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

from pydantic import ValidationError

from tatc.schemas import KeplerianOrbit


class TestKeplerianOrbit(unittest.TestCase):
    """
    Unit tests for the KeplerianOrbit schema.
    """
    def setUp(self):
        self.test_data = {
            "semimajor_axis": 400000 + 6371000, # 400 km mean altitude
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
            "eccentricity": 0.01,
            "perigee_argument": 100.0,
        }
        self.test_orbit = KeplerianOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the KeplerianOrbit schema correctly initializes with valid data.
        """
        self.assertEqual(self.test_orbit.semimajor_axis, self.test_data.get("semimajor_axis"))
        self.assertEqual(
            self.test_orbit.true_anomaly, self.test_data.get("true_anomaly")
        )
        self.assertEqual(self.test_orbit.epoch, self.test_data.get("epoch"))
        self.assertEqual(self.test_orbit.inclination, self.test_data.get("inclination"))
        self.assertEqual(
            self.test_orbit.right_ascension_ascending_node,
            self.test_data.get("right_ascension_ascending_node"),
        )
        self.assertEqual(
            self.test_orbit.eccentricity, self.test_data.get("eccentricity")
        )
        self.assertEqual(
            self.test_orbit.perigee_argument, self.test_data.get("perigee_argument")
        )

    def test_defaults(self):
        """
        Test that inclination, right_ascension_ascending_node,
        eccentricity, and perigee_argument each default to 0 when
        omitted, leaving semimajor_axis as the only required field.
        """
        o = KeplerianOrbit(semimajor_axis=7000000)
        self.assertEqual(o.inclination, 0)
        self.assertEqual(o.right_ascension_ascending_node, 0)
        self.assertEqual(o.eccentricity, 0)
        self.assertEqual(o.perigee_argument, 0)

    def test_bad_semimajor_axis_missing(self):
        """
        Test that the KeplerianOrbit schema raises a ValidationError when
        the required semimajor_axis field is missing.
        """
        with self.assertRaises(ValidationError):
            KeplerianOrbit()

    def test_bad_semimajor_axis_negative(self):
        """
        Test that a negative semimajor_axis is rejected.
        """
        with self.assertRaises(ValidationError):
            KeplerianOrbit(semimajor_axis=-1)

    def test_bad_semimajor_axis_zero(self):
        """
        Test that a zero semimajor_axis is rejected, since it is
        physically meaningless (and previously produced nan downstream in
        get_mean_motion/get_orbit_period).
        """
        with self.assertRaises(ValidationError):
            KeplerianOrbit(semimajor_axis=0)

    def test_bad_eccentricity_parabolic(self):
        """
        Test that an eccentricity of exactly 1 (parabolic trajectory) is
        rejected, since KeplerianOrbit represents elliptical motion only
        (per its own docstring).
        """
        with self.assertRaises(ValidationError):
            KeplerianOrbit(semimajor_axis=7000000, eccentricity=1.0)

    def test_bad_eccentricity_hyperbolic(self):
        """
        Test that an eccentricity greater than 1 (hyperbolic escape
        trajectory) is rejected.
        """
        with self.assertRaises(ValidationError):
            KeplerianOrbit(semimajor_axis=7000000, eccentricity=1.5)

    def test_eccentricity_boundary_zero(self):
        """
        Test that an eccentricity of exactly 0 (the ge=0 boundary,
        a circular orbit) is accepted.
        """
        self.assertEqual(
            KeplerianOrbit(semimajor_axis=7000000, eccentricity=0).eccentricity, 0
        )

    def test_get_orbit_period_and_mean_motion_match_published_gps_orbit(self):
        """
        Test get_orbit_period and get_mean_motion against the well-known
        GPS constellation orbit: a semimajor axis of about 26,560 km
        gives an orbital period of half a sidereal day (~11h58m, by
        design, so that ground tracks repeat daily) and a mean motion of
        about 2 revolutions/day.
        """
        gps_orbit = KeplerianOrbit(
            semimajor_axis=26560000, inclination=55, eccentricity=0.01
        )
        self.assertAlmostEqual(
            gps_orbit.get_orbit_period(),
            timedelta(hours=11, minutes=58),
            delta=timedelta(minutes=1),
        )
        self.assertAlmostEqual(
            gps_orbit.get_mean_motion() * 86400 / 360, 2.0, delta=0.01
        )

    def test_to_gp_orbit_period_matches_published_gps_orbit(self):
        """
        Test that converting the same GPS-like orbit to a general
        perturbations representation preserves the published ~11h58m
        orbital period end-to-end.
        """
        gps_orbit = KeplerianOrbit(
            semimajor_axis=26560000, inclination=55, eccentricity=0.01
        )
        self.assertAlmostEqual(
            gps_orbit.to_gp_orbit().get_orbit_period(),
            timedelta(hours=11, minutes=58),
            delta=timedelta(minutes=1),
        )

    def test_getters_return_the_underlying_fields(self):
        """
        Test that the OrbitBase-interface getter methods each return the
        corresponding field value directly.
        """
        self.assertEqual(
            self.test_orbit.get_semimajor_axis(), self.test_data.get("semimajor_axis")
        )
        self.assertEqual(
            self.test_orbit.get_inclination(), self.test_data.get("inclination")
        )
        self.assertEqual(
            self.test_orbit.get_right_ascension_ascending_node(),
            self.test_data.get("right_ascension_ascending_node"),
        )
        self.assertEqual(
            self.test_orbit.get_eccentricity(), self.test_data.get("eccentricity")
        )
        self.assertEqual(
            self.test_orbit.get_perigee_argument(),
            self.test_data.get("perigee_argument"),
        )

    def test_get_derived_orbit(self):
        """
        Test that a derived orbit can be computed from the base orbit.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.get_mean_anomaly(),
            self.test_orbit.get_mean_anomaly() + 20,
            delta=0.001,
        )
        self.assertAlmostEqual(
            derived_orbit.right_ascension_ascending_node,
            self.test_orbit.right_ascension_ascending_node + 10,
            delta=0.001,
        )

    def test_to_gp_orbit(self):
        """
        Test that the KeplerianOrbit can be converted to a GP orbit.
        """
        gp_orbit = self.test_orbit.to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_mean_altitude(),
            self.test_data.get("semimajor_axis") - 6371000,
            delta=10
        )
        self.assertAlmostEqual(
            gp_orbit.get_true_anomaly(),
            self.test_data.get("true_anomaly"),
            delta=0.001
        )
        self.assertAlmostEqual(
            gp_orbit.get_epoch().timestamp(),
            self.test_data.get("epoch").timestamp(),
            delta=1,
        )
        self.assertEqual(
            gp_orbit.get_inclination(),
            self.test_data.get("inclination"),
        )
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(),
            self.test_data.get("right_ascension_ascending_node"),
        )
        self.assertAlmostEqual(
            gp_orbit.get_eccentricity(),
            self.test_data.get("eccentricity"),
        )
        self.assertAlmostEqual(
            gp_orbit.get_perigee_argument(),
            self.test_data.get("perigee_argument"),
        )
