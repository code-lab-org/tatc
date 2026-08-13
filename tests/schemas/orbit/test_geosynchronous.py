"""
Unit tests for the GeosynchronousOrbit schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

from pydantic import ValidationError

from tatc import constants
from tatc.constants import EARTH_SIDEREAL_DAY_S
from tatc.schemas import GeosynchronousOrbit


class TestGeosynchronousOrbit(unittest.TestCase):
    """
    Unit tests for the GeosynchronousOrbit schema.
    """

    def setUp(self):
        self.test_data = {
            "longitude": -75.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
        }
        self.test_orbit = GeosynchronousOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the GeosynchronousOrbit schema correctly initializes
        with valid data.
        """
        self.assertEqual(self.test_orbit.longitude, self.test_data.get("longitude"))
        self.assertEqual(self.test_orbit.epoch, self.test_data.get("epoch"))
        self.assertEqual(self.test_orbit.type, "geosynchronous")

    def test_defaults(self):
        """
        Test that mean_altitude defaults to the geosynchronous altitude
        (the value yielding an orbit period of exactly one sidereal day)
        and inclination defaults to 0, when omitted.
        """
        o = GeosynchronousOrbit(longitude=0)
        self.assertAlmostEqual(o.mean_altitude, 35786000, delta=10000)
        self.assertEqual(o.inclination, 0)

    def test_bad_longitude_missing(self):
        """
        Test that the GeosynchronousOrbit schema raises a ValidationError
        when the required longitude field is missing.
        """
        with self.assertRaises(ValidationError):
            GeosynchronousOrbit()

    def test_bad_longitude_too_large(self):
        """
        Test that a longitude above 180 degrees is rejected.
        """
        with self.assertRaises(ValidationError):
            GeosynchronousOrbit(longitude=180.1)

    def test_bad_longitude_too_small(self):
        """
        Test that a longitude below -180 degrees is rejected.
        """
        with self.assertRaises(ValidationError):
            GeosynchronousOrbit(longitude=-180.1)

    def test_longitude_boundary_values(self):
        """
        Test that longitude values exactly at the antimeridian (-180, 180
        degrees) are accepted.
        """
        self.assertEqual(GeosynchronousOrbit(longitude=180).longitude, 180)
        self.assertEqual(GeosynchronousOrbit(longitude=-180).longitude, -180)

    def test_bad_inclination_negative(self):
        """
        Test that a negative inclination is rejected.
        """
        with self.assertRaises(ValidationError):
            GeosynchronousOrbit(longitude=0, inclination=-0.1)

    def test_bad_inclination_too_large(self):
        """
        Test that an inclination of 180 degrees or more is rejected.
        """
        with self.assertRaises(ValidationError):
            GeosynchronousOrbit(longitude=0, inclination=180)

    def test_get_eccentricity_and_perigee_argument_inherited_from_circular_base(self):
        """
        Test that get_eccentricity and get_perigee_argument are inherited
        from CircularOrbitBase (both 0, since a geosynchronous orbit is
        circular).
        """
        self.assertEqual(self.test_orbit.get_eccentricity(), 0)
        self.assertEqual(self.test_orbit.get_perigee_argument(), 0)

    def test_get_semimajor_axis_matches_published_geosynchronous_altitude(self):
        """
        Test get_semimajor_axis against the well-known published
        geosynchronous altitude of ~35,786 km (semimajor axis ~42,164 km).
        """
        self.assertAlmostEqual(
            self.test_orbit.get_semimajor_axis(), 42164000, delta=10000
        )

    def test_get_orbit_period_is_one_sidereal_day(self):
        """
        Test that get_orbit_period matches Earth's sidereal day, the
        defining property of a geosynchronous orbit.
        """
        self.assertAlmostEqual(
            self.test_orbit.get_orbit_period().total_seconds(),
            EARTH_SIDEREAL_DAY_S,
            delta=1.0,
        )

    def test_get_right_ascension_ascending_node_at_greenwich_matches_gast(self):
        """
        Test that a satellite at longitude 0 (the Greenwich meridian) has
        a right ascension of ascending node exactly equal to Earth's
        rotation angle (Greenwich Apparent Sidereal Time) at epoch -- the
        defining relationship between Earth-fixed longitude and inertial
        right ascension.
        """
        epoch = datetime(2022, 6, 15, 6, 0, 0, tzinfo=timezone.utc)
        o = GeosynchronousOrbit(longitude=0, epoch=epoch)
        expected_raan = constants.timescale.from_datetime(epoch).gast * 15 % 360
        self.assertAlmostEqual(
            o.get_right_ascension_ascending_node(), expected_raan, delta=1e-6
        )

    def test_get_right_ascension_ascending_node_wraps_to_0_360(self):
        """
        Test that the right ascension of ascending node always falls
        within [0, 360), even when the raw longitude+rotation-angle sum
        would otherwise exceed that range.
        """
        raan = self.test_orbit.get_right_ascension_ascending_node()
        self.assertGreaterEqual(raan, 0)
        self.assertLess(raan, 360)

    def test_get_right_ascension_ascending_node_tracks_earth_rotation(self):
        """
        Test that a satellite parked at a fixed longitude has nearly the
        same right ascension of ascending node exactly one sidereal day
        later, since it remains above the same Earth-fixed point (the
        defining behavior of a geosynchronous orbit).
        """
        later_orbit = GeosynchronousOrbit(
            longitude=self.test_orbit.longitude,
            epoch=self.test_orbit.epoch + timedelta(seconds=EARTH_SIDEREAL_DAY_S),
        )
        self.assertAlmostEqual(
            later_orbit.get_right_ascension_ascending_node(),
            self.test_orbit.get_right_ascension_ascending_node(),
            delta=0.01,
        )

    def test_get_derived_orbit_shifts_longitude(self):
        """
        Test that get_derived_orbit shifts longitude by delta_raan
        degrees directly, since Earth's rotation angle at a fixed epoch
        does not change.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(0, 10)
        self.assertAlmostEqual(derived_orbit.longitude, -65.0, delta=1e-9)

    def test_get_derived_orbit_wraps_longitude_across_antimeridian(self):
        """
        Test that get_derived_orbit wraps longitude correctly when the
        shift crosses the +/-180 degree antimeridian.
        """
        o = GeosynchronousOrbit(longitude=170)
        derived_orbit = o.get_derived_orbit(0, 20)
        self.assertAlmostEqual(derived_orbit.longitude, -170.0, delta=1e-9)

    def test_get_derived_orbit_preserves_other_fields(self):
        """
        Test that get_derived_orbit preserves mean_altitude, inclination,
        and epoch unchanged.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertEqual(derived_orbit.mean_altitude, self.test_orbit.mean_altitude)
        self.assertEqual(derived_orbit.inclination, self.test_orbit.inclination)
        self.assertEqual(derived_orbit.epoch, self.test_orbit.epoch)

    def test_to_gp_orbit(self):
        """
        Test that the GeosynchronousOrbit schema correctly converts to a
        general perturbations orbit.
        """
        gp_orbit = self.test_orbit.to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_mean_altitude(),
            self.test_orbit.mean_altitude,
            delta=1.0,
        )
        self.assertAlmostEqual(
            gp_orbit.get_inclination(), self.test_orbit.get_inclination(), delta=0.001
        )
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(),
            self.test_orbit.get_right_ascension_ascending_node(),
            delta=0.001,
        )
        self.assertAlmostEqual(gp_orbit.get_eccentricity(), 0, delta=0.001)
        self.assertAlmostEqual(
            gp_orbit.get_orbit_period().total_seconds(),
            EARTH_SIDEREAL_DAY_S,
            delta=1.0,
        )


if __name__ == "__main__":
    unittest.main()
