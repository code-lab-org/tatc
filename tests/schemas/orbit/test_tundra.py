"""
Unit tests for the TundraOrbit schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from pydantic import ValidationError

from tatc.constants import EARTH_MEAN_RADIUS, EARTH_SIDEREAL_DAY_S
from tatc.schemas import TundraOrbit


class TestTundraOrbit(unittest.TestCase):
    """
    Unit tests for the TundraOrbit schema.
    """

    def setUp(self):
        self.test_data = {
            "perigee_altitude": 24480000,
            "right_ascension_ascending_node": 11.0394,
        }
        self.test_orbit = TundraOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the TundraOrbit schema correctly initializes with valid data.
        """
        self.assertEqual(
            self.test_orbit.perigee_altitude, self.test_data.get("perigee_altitude")
        )
        self.assertEqual(
            self.test_orbit.right_ascension_ascending_node,
            self.test_data.get("right_ascension_ascending_node"),
        )
        self.assertAlmostEqual(self.test_orbit.get_inclination(), 63.4, delta=0.1)
        self.assertAlmostEqual(
            self.test_orbit.get_orbit_period().total_seconds(), 1436 * 60, delta=60
        )
        self.assertEqual(self.test_orbit.get_perigee_argument(), 270)
        self.assertAlmostEqual(self.test_orbit.get_eccentricity(), 0.24, delta=0.03)

    def test_defaults(self):
        """
        Test that northern_coverage and right_ascension_ascending_node
        default correctly when omitted.
        """
        o = TundraOrbit(perigee_altitude=24480000)
        self.assertTrue(o.northern_coverage)
        self.assertEqual(o.right_ascension_ascending_node, 0)

    def test_bad_perigee_altitude_missing(self):
        """
        Test that the TundraOrbit schema raises a ValidationError when
        the required perigee_altitude field is missing.
        """
        with self.assertRaises(ValidationError):
            TundraOrbit()

    def test_bad_perigee_altitude_negative(self):
        """
        Test that a negative perigee_altitude is rejected.
        """
        with self.assertRaises(ValidationError):
            TundraOrbit(perigee_altitude=-1)

    def test_bad_perigee_altitude_exceeds_implied_apogee(self):
        """
        Test that a perigee_altitude above the Kepler-derived ceiling
        implied by Tundra's fixed (one sidereal day) orbit period is
        rejected, since it would otherwise silently produce a negative
        (physically invalid) eccentricity.
        """
        with self.assertRaises(ValidationError):
            TundraOrbit(perigee_altitude=36000000)

    def test_get_semimajor_axis_matches_geostationary_altitude(self):
        """
        Test get_semimajor_axis against a published reference: a Tundra
        orbit's period is defined as one sidereal day (same as a
        geostationary orbit), so its semimajor axis should match the
        well-known geostationary value of ~42,164 km.
        """
        self.assertAlmostEqual(
            self.test_orbit.get_semimajor_axis(),
            42164000,
            delta=1000,
        )
        self.assertAlmostEqual(
            self.test_orbit.get_semimajor_axis() - EARTH_MEAN_RADIUS,
            35786000,
            delta=10000,
        )

    def test_get_orbit_period_is_j2_corrected(self):
        """
        Regression test: get_orbit_period must not simply return the
        naive, uncorrected one sidereal day -- it should be shifted by a
        small (order 1-10 second), nonzero correction accounting for
        Earth's J2 oblateness perturbation to the true rate of mean
        anomaly advance. This exercises the fix for the historical
        "TODO this needs to be corrected to account for J2 effects".
        """
        naive_period_s = EARTH_SIDEREAL_DAY_S
        corrected_period_s = self.test_orbit.get_orbit_period().total_seconds()
        self.assertNotAlmostEqual(corrected_period_s, naive_period_s, delta=1e-6)
        self.assertAlmostEqual(corrected_period_s, naive_period_s, delta=10)

    def test_get_apogee_altitude_matches_published_sirius_xm_tundra_orbit(self):
        """
        Test the implied apogee altitude against the published Sirius/XM
        satellite radio Tundra orbit: a perigee altitude of ~24,480 km
        (this class's own test data) should yield an apogee altitude
        close to the published ~46,983 km, since both describe the same
        one-sidereal-day Tundra orbit family.

        See: https://en.wikipedia.org/wiki/Tundra_orbit
        """
        apogee_altitude = (
            2 * self.test_orbit.get_semimajor_axis()
            - (EARTH_MEAN_RADIUS + self.test_orbit.perigee_altitude)
            - EARTH_MEAN_RADIUS
        )
        self.assertAlmostEqual(apogee_altitude, 46983000, delta=500000)

    def test_get_mean_anomaly_at_perigee_is_zero(self):
        """
        Test that a true anomaly of 0 (perigee) always maps to a mean
        anomaly of 0, regardless of eccentricity.
        """
        o = TundraOrbit(perigee_altitude=24480000, true_anomaly=0)
        self.assertAlmostEqual(o.get_mean_anomaly(), 0, delta=1e-6)

    def test_get_derived_orbit_preserves_other_fields(self):
        """
        Test that get_derived_orbit preserves perigee_altitude,
        northern_coverage, and epoch unchanged.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertEqual(
            derived_orbit.perigee_altitude, self.test_orbit.perigee_altitude
        )
        self.assertEqual(
            derived_orbit.northern_coverage, self.test_orbit.northern_coverage
        )
        self.assertEqual(derived_orbit.epoch, self.test_orbit.epoch)

    def test_get_derived_orbit(self):
        """
        Test that the TundraOrbit schema correctly derives a new orbit with the specified parameters.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.right_ascension_ascending_node,
            self.test_orbit.right_ascension_ascending_node + 10,
            delta=0.001,
        )

    def test_to_gp_orbit(self):
        """
        Test that the TundraOrbit schema correctly converts to a general perturbations representation.
        """
        gp_orbit = self.test_orbit.to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(),
            self.test_data.get("right_ascension_ascending_node"),
            delta=0.1,
        )
        self.assertAlmostEqual(
            gp_orbit.get_inclination(), self.test_orbit.get_inclination(), delta=0.01
        )
        self.assertAlmostEqual(
            gp_orbit.get_eccentricity(), self.test_orbit.get_eccentricity(), delta=0.001
        )
        self.assertAlmostEqual(
            gp_orbit.get_perigee_argument(),
            self.test_orbit.get_perigee_argument(),
            delta=0.01,
        )
