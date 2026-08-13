"""
Unit tests for the CircularOrbit schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import ValidationError

from tatc.constants import EARTH_MEAN_RADIUS, EARTH_MU
from tatc.schemas import CircularOrbit
from tatc.utils.orbital import semimajor_axis_to_mean_motion


class TestCircularOrbit(unittest.TestCase):
    """
    Unit tests for the CircularOrbit schema.
    """

    def setUp(self):
        self.test_data = {
            "mean_altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, tzinfo=timezone.utc),
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
        }
        self.test_orbit = CircularOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the CircularOrbit schema correctly initializes with valid data.
        """
        self.assertEqual(
            self.test_orbit.mean_altitude, self.test_data.get("mean_altitude")
        )
        self.assertEqual(
            self.test_orbit.true_anomaly, self.test_data.get("true_anomaly")
        )
        self.assertEqual(self.test_orbit.epoch, self.test_data.get("epoch"))
        self.assertEqual(self.test_orbit.inclination, self.test_data.get("inclination"))
        self.assertEqual(
            self.test_orbit.right_ascension_ascending_node,
            self.test_data.get("right_ascension_ascending_node"),
        )

    def test_good_data_iso8601_datetime(self):
        """
        Test that the CircularOrbit schema correctly initializes with valid data
        when the epoch is provided as an ISO 8601 string.
        """
        good_data = {
            "mean_altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": "2022-01-01T12:00:00Z",
            "inclination": 45.0,
            "right_ascension_ascending_node": 50.0,
        }
        o = CircularOrbit(**good_data)
        self.assertEqual(o.mean_altitude, good_data.get("mean_altitude"))
        self.assertEqual(o.true_anomaly, good_data.get("true_anomaly"))
        self.assertEqual(o.epoch, datetime(2022, 1, 1, 12, tzinfo=timezone.utc))
        self.assertEqual(o.inclination, good_data.get("inclination"))
        self.assertEqual(
            o.right_ascension_ascending_node,
            good_data.get("right_ascension_ascending_node"),
        )

    def test_defaults(self):
        """
        Test that inclination, right_ascension_ascending_node, and the
        type discriminator each default correctly when omitted.
        """
        o = CircularOrbit(mean_altitude=500000)
        self.assertEqual(o.inclination, 0)
        self.assertEqual(o.right_ascension_ascending_node, 0)
        self.assertEqual(o.type, "circular")

    def test_bad_mean_altitude_negative(self):
        """
        Test that a negative mean_altitude is rejected, since a circular
        orbit cannot have a semimajor axis below Earth's mean radius.
        """
        with self.assertRaises(ValidationError):
            CircularOrbit(mean_altitude=-1)

    def test_mean_altitude_boundary_zero(self):
        """
        Test that a mean_altitude of exactly 0 (the ge=0 boundary) is
        accepted.
        """
        self.assertEqual(CircularOrbit(mean_altitude=0).mean_altitude, 0)

    def test_bad_mean_altitude_missing(self):
        """
        Test that the CircularOrbit schema raises a ValidationError when
        the required mean_altitude field is missing.
        """
        with self.assertRaises(ValidationError):
            CircularOrbit()

    def test_bad_inclination_negative(self):
        """
        Test that a negative inclination is rejected.
        """
        with self.assertRaises(ValidationError):
            CircularOrbit(mean_altitude=500000, inclination=-0.1)

    def test_bad_inclination_too_large(self):
        """
        Test that an inclination of 180 degrees or more is rejected
        (inclination is conventionally bounded to [0, 180)).
        """
        with self.assertRaises(ValidationError):
            CircularOrbit(mean_altitude=500000, inclination=180)

    def test_inclination_boundary_zero(self):
        """
        Test that an inclination of exactly 0 degrees (the ge=0 boundary)
        is accepted.
        """
        self.assertEqual(
            CircularOrbit(mean_altitude=500000, inclination=0).inclination, 0
        )

    def test_bad_right_ascension_ascending_node_negative(self):
        """
        Test that a negative right_ascension_ascending_node is rejected.
        """
        with self.assertRaises(ValidationError):
            CircularOrbit(mean_altitude=500000, right_ascension_ascending_node=-0.1)

    def test_bad_right_ascension_ascending_node_too_large(self):
        """
        Test that a right_ascension_ascending_node of 360 degrees or more
        is rejected.
        """
        with self.assertRaises(ValidationError):
            CircularOrbit(mean_altitude=500000, right_ascension_ascending_node=360)

    def test_get_semimajor_axis(self):
        """
        Test that the CircularOrbit class correctly calculates the semimajor axis.
        """
        self.assertEqual(
            self.test_orbit.get_semimajor_axis(),
            self.test_data.get("mean_altitude") + EARTH_MEAN_RADIUS,
        )

    def test_get_mean_altitude(self):
        """
        Test that get_mean_altitude returns the mean_altitude field
        directly.
        """
        self.assertEqual(
            self.test_orbit.get_mean_altitude(), self.test_data.get("mean_altitude")
        )

    def test_get_inclination(self):
        """
        Test that get_inclination returns the inclination field directly.
        """
        self.assertEqual(
            self.test_orbit.get_inclination(), self.test_data.get("inclination")
        )

    def test_get_right_ascension_ascending_node(self):
        """
        Test that get_right_ascension_ascending_node returns the
        right_ascension_ascending_node field directly.
        """
        self.assertEqual(
            self.test_orbit.get_right_ascension_ascending_node(),
            self.test_data.get("right_ascension_ascending_node"),
        )

    def test_get_eccentricity(self):
        """
        Test that get_eccentricity returns 0, since a CircularOrbit has no
        eccentricity.
        """
        self.assertEqual(self.test_orbit.get_eccentricity(), 0)

    def test_get_perigee_argument(self):
        """
        Test that get_perigee_argument returns 0, since a CircularOrbit
        has no perigee to reference.
        """
        self.assertEqual(self.test_orbit.get_perigee_argument(), 0)

    def test_get_mean_anomaly(self):
        """
        Test that the CircularOrbit class correctly calculates the mean anomaly.
        """
        self.assertEqual(
            self.test_orbit.get_mean_anomaly(), self.test_data.get("true_anomaly")
        )

    def test_get_mean_motion(self):
        """
        Test that the CircularOrbit class correctly calculates the mean motion.
        """
        self.assertAlmostEqual(
            self.test_orbit.get_mean_motion(),
            semimajor_axis_to_mean_motion(
                EARTH_MEAN_RADIUS + self.test_orbit.mean_altitude
            ),
            delta=0.001,
        )

    def test_get_orbit_period(self):
        """
        Test that the CircularOrbit class correctly calculates the orbit period.
        """
        orbit_period = (
            2
            * np.pi
            * np.sqrt(
                np.power(EARTH_MEAN_RADIUS + self.test_orbit.mean_altitude, 3)
                / EARTH_MU
            )
        )
        self.assertAlmostEqual(
            self.test_orbit.get_orbit_period(),
            timedelta(seconds=orbit_period),
            delta=1.0,
        )

    def test_get_derived_orbit(self):
        """
        Test that the CircularOrbit schema correctly derives a new orbit.
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

    def test_get_derived_orbit_preserves_other_fields(self):
        """
        Test that get_derived_orbit preserves mean_altitude, inclination,
        and epoch unchanged, only perturbing mean anomaly and RAAN.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertEqual(derived_orbit.mean_altitude, self.test_orbit.mean_altitude)
        self.assertEqual(derived_orbit.inclination, self.test_orbit.inclination)
        self.assertEqual(derived_orbit.epoch, self.test_orbit.epoch)

    def test_get_derived_orbit_wraps_past_360_degrees(self):
        """
        Test that get_derived_orbit wraps both mean anomaly and RAAN back
        into [0, 360) when the perturbation pushes them past 360 degrees.
        """
        orbit = CircularOrbit(
            mean_altitude=400000,
            true_anomaly=350.0,
            right_ascension_ascending_node=350.0,
        )
        derived_orbit = orbit.get_derived_orbit(20, 20)
        self.assertAlmostEqual(derived_orbit.get_mean_anomaly(), 10.0, delta=0.001)
        self.assertAlmostEqual(
            derived_orbit.right_ascension_ascending_node, 10.0, delta=0.001
        )

    def test_to_gp_orbit(self):
        """
        Test that the CircularOrbit schema correctly converts to a general perturbations orbit.
        """
        gp_orbit = self.test_orbit.to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_mean_altitude(), self.test_data.get("mean_altitude"), delta=1.0
        )
        self.assertAlmostEqual(
            gp_orbit.get_true_anomaly(), self.test_data.get("true_anomaly"), delta=0.001
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
        self.assertAlmostEqual(gp_orbit.get_eccentricity(), 0, delta=0.001)
        self.assertAlmostEqual(gp_orbit.get_perigee_argument(), 0, delta=0.001)
