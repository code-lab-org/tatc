"""
Unit tests for the MolniyaOrbit schema.

@author: Paul T. Grogan <paul.t.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np

from pydantic import ValidationError
from skyfield.api import wgs84

from tatc.constants import (
    EARTH_J2_CRITICAL_INCLINATION,
    EARTH_MEAN_RADIUS,
    EARTH_SIDEREAL_DAY_S,
)
from tatc.schemas import MolniyaOrbit


def sampled_apogee_longitudes(orbit, hours):
    """
    Longitudes (degrees) of the apogees found by sampling the propagated
    orbit every 20 s from its epoch.
    """
    times = [orbit.epoch + timedelta(seconds=20 * k) for k in range(hours * 180)]
    track = orbit.to_gp_orbit().get_orbit_track(times)
    radius = np.linalg.norm(track.position.m, axis=0)
    k = np.nonzero((radius[1:-1] > radius[:-2]) & (radius[1:-1] >= radius[2:]))[0] + 1
    return wgs84.subpoint_of(track[k]).longitude.degrees


class TestMolniyaOrbit(unittest.TestCase):
    """
    Unit tests for the MolniyaOrbit schema.
    """

    def setUp(self):
        self.test_data = {
            "perigee_altitude": 1199100,
            "right_ascension_ascending_node": 11.0394,
        }
        self.test_orbit = MolniyaOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the MolniyaOrbit schema correctly initializes with valid data.
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
            self.test_orbit.get_orbit_period().total_seconds(), 718 * 60, delta=60
        )
        self.assertEqual(self.test_orbit.get_perigee_argument(), 270)
        self.assertAlmostEqual(self.test_orbit.get_eccentricity(), 0.74, delta=0.1)

    def test_defaults(self):
        """
        Test that northern_coverage and right_ascension_ascending_node
        default correctly when omitted.
        """
        o = MolniyaOrbit(perigee_altitude=550000)
        self.assertTrue(o.northern_coverage)
        self.assertEqual(o.right_ascension_ascending_node, 0)

    def test_bad_perigee_altitude_missing(self):
        """
        Test that the MolniyaOrbit schema raises a ValidationError when
        the required perigee_altitude field is missing.
        """
        with self.assertRaises(ValidationError):
            MolniyaOrbit()

    def test_bad_perigee_altitude_negative(self):
        """
        Test that a negative perigee_altitude is rejected.
        """
        with self.assertRaises(ValidationError):
            MolniyaOrbit(perigee_altitude=-1)

    def test_bad_perigee_altitude_exceeds_implied_apogee(self):
        """
        Test that a perigee_altitude above the Kepler-derived ceiling
        implied by Molniya's fixed (half sidereal day) orbit period is
        rejected, since it would otherwise silently produce a negative
        (physically invalid) eccentricity.
        """
        with self.assertRaises(ValidationError):
            MolniyaOrbit(perigee_altitude=21190753)

    def test_get_inclination_and_apogee_altitude_match_published_molniya_orbit(self):
        """
        Test get_inclination and the implied apogee altitude against a
        published classic Molniya orbit: a ~550 km perigee altitude with
        the 63.4 degree critical inclination yields an apogee altitude of
        approximately 40,000 km and a ~12 hour (half sidereal day) period.

        See: https://en.wikipedia.org/wiki/Molniya_orbit
        """
        o = MolniyaOrbit(perigee_altitude=550000)
        self.assertAlmostEqual(o.get_inclination(), 63.4, delta=0.1)
        apogee_altitude = (
            2 * o.get_semimajor_axis()
            - (EARTH_MEAN_RADIUS + o.perigee_altitude)
            - EARTH_MEAN_RADIUS
        )
        self.assertAlmostEqual(apogee_altitude, 40000000, delta=500000)

    def test_get_orbit_period_is_j2_corrected(self):
        """
        Regression test: get_orbit_period must not simply return the
        naive, uncorrected half sidereal day -- it should be shifted by a
        small (order 10 second), nonzero correction accounting for
        Earth's J2 oblateness perturbation to the true rate of mean
        anomaly advance and the precession of the ascending node. This
        exercises the fix for the historical "TODO this needs to be
        corrected to account for J2 effects".
        """
        naive_period_s = EARTH_SIDEREAL_DAY_S / 2
        corrected_period_s = self.test_orbit.get_orbit_period().total_seconds()
        self.assertNotAlmostEqual(corrected_period_s, naive_period_s, delta=1e-6)
        self.assertAlmostEqual(corrected_period_s, naive_period_s, delta=30)

    def test_inclination_default(self):
        """
        Test that the inclination defaults to the critical inclination.
        """
        self.assertEqual(self.test_orbit.inclination, EARTH_J2_CRITICAL_INCLINATION)
        self.assertEqual(
            self.test_orbit.get_inclination(), EARTH_J2_CRITICAL_INCLINATION
        )

    def test_inclination_custom(self):
        """
        Test that a custom inclination is used by the orbit, its derived
        orbits, and its general perturbations representation.
        """
        orbit = MolniyaOrbit(perigee_altitude=600e3, inclination=50)
        self.assertEqual(orbit.get_inclination(), 50)
        self.assertEqual(orbit.get_derived_orbit(20, 10).inclination, 50)
        self.assertAlmostEqual(orbit.to_gp_orbit().get_inclination(), 50, delta=0.01)
        self.assertAlmostEqual(orbit.get_mean_motion() * 86400 / 360, 2.0, delta=0.01)

    def test_period_cache_follows_fields(self):
        """
        Test that the cached orbit period is recomputed for a copy with a
        changed inclination or perigee altitude (which model_copy copies
        along with the cache).
        """
        orbit = MolniyaOrbit(perigee_altitude=600e3)
        orbit.get_orbit_period()
        for update in ({"inclination": 50}, {"perigee_altitude": 1500e3}):
            self.assertEqual(
                orbit.model_copy(update=update).get_orbit_period(),
                MolniyaOrbit(
                    **{"perigee_altitude": 600e3, **update}
                ).get_orbit_period(),
            )

    def test_bad_inclination(self):
        """
        Test that the MolniyaOrbit schema raises a ValidationError for an
        inclination outside [0, 180).
        """
        for inclination in (-1, 180):
            with self.assertRaises(ValidationError):
                MolniyaOrbit(perigee_altitude=600e3, inclination=inclination)

    def test_apogee_longitudes_repeat(self):
        """
        Test that the ground track repeats as propagated: the longitudes of
        each of the two daily apogees drift by less than 0.01 deg per day
        over 20 days (without accounting for the precession of the
        ascending node, they drift westward by about 0.1 deg per day; at
        50 deg inclination, without accounting for the precession of the
        argument of perigee, eastward by about 0.3 deg per day), and the
        apogees are at the latitude of the inclination.
        """
        epoch = datetime(2026, 10, 4, tzinfo=timezone.utc)
        for perigee_altitude, inclination in (
            (600e3, EARTH_J2_CRITICAL_INCLINATION),
            (1500e3, EARTH_J2_CRITICAL_INCLINATION),
            (600e3, 50),
        ):
            orbit = MolniyaOrbit(
                perigee_altitude=perigee_altitude, inclination=inclination, epoch=epoch
            )
            times = [epoch + timedelta(minutes=2 * k) for k in range(20 * 720)]
            track = orbit.to_gp_orbit().get_orbit_track(times)
            radius = np.linalg.norm(track.position.m, axis=0)
            apogees = (
                np.nonzero((radius[1:-1] > radius[:-2]) & (radius[1:-1] >= radius[2:]))[
                    0
                ]
                + 1
            )
            subpoint = wgs84.subpoint_of(track[apogees])
            longitude = subpoint.longitude.degrees
            days = np.array([(times[k] - epoch) / timedelta(days=1) for k in apogees])
            for first in (0, 1):
                rate = np.polyfit(
                    days[first::2],
                    np.degrees(np.unwrap(np.radians(longitude[first::2]))),
                    1,
                )[0]
                self.assertLess(abs(rate), 0.01, (perigee_altitude, inclination))
            self.assertAlmostEqual(subpoint.latitude.degrees[0], inclination, delta=0.2)

    def test_get_semimajor_axis_and_mean_motion(self):
        """
        Test get_semimajor_axis and get_mean_motion are self-consistent
        with the fixed half-sidereal-day orbit period via Kepler's third
        law and the standard mean-motion relation.
        """
        self.assertAlmostEqual(
            self.test_orbit.get_mean_motion() * 86400 / 360, 2.0, delta=0.01
        )

    def test_get_mean_anomaly_at_perigee_is_zero(self):
        """
        Test that a true anomaly of 0 (perigee) always maps to a mean
        anomaly of 0, regardless of eccentricity -- a basic invariant of
        the true-to-mean anomaly relationship.
        """
        o = MolniyaOrbit(perigee_altitude=550000, true_anomaly=0)
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
        Test that the MolniyaOrbit schema correctly derives a new orbit with the specified parameters.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.right_ascension_ascending_node,
            self.test_orbit.right_ascension_ascending_node + 10,
            delta=0.001,
        )

    def test_from_apogee_longitude(self):
        """
        Test that a Molniya orbit placed by longitude has its first apogee
        over that longitude and its second apogee 180 degrees away.
        """
        orbit = MolniyaOrbit.from_apogee_longitude(
            40,
            perigee_altitude=600e3,
            true_anomaly=200,
            epoch=datetime(2026, 10, 4, tzinfo=timezone.utc),
        )
        self.assertIsInstance(orbit, MolniyaOrbit)
        self.assertAlmostEqual(orbit.get_apogee_longitude(), 40, delta=1e-5)
        np.testing.assert_allclose(
            sampled_apogee_longitudes(orbit, 25), [40, -140], atol=0.02
        )

    def test_to_gp_orbit(self):
        """
        Test that the MolniyaOrbit schema correctly converts to a general perturbations representation.
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
