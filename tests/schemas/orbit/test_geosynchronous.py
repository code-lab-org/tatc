"""
Unit tests for the GeosynchronousOrbit schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
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

    def test_get_right_ascension_ascending_node_near_gmst(self):
        """
        Test that a satellite at longitude 0 (the Greenwich meridian) has a right
        ascension of ascending node near Earth's rotation angle that orients the
        TEME frame of SGP4 (Greenwich Mean Sidereal Time) at epoch, within the
        correction for SGP4's lunar and solar terms (up to about 0.04 deg).
        """
        epoch = datetime(2022, 6, 15, 6, 0, 0, tzinfo=timezone.utc)
        o = GeosynchronousOrbit(longitude=0, epoch=epoch)
        gmst = constants.timescale.from_datetime(epoch).gmst * 15 % 360
        self.assertAlmostEqual(o.get_right_ascension_ascending_node(), gmst, delta=0.05)

    def test_longitude_at_epoch(self):
        """
        Test that the sub-satellite point at epoch is at the orbit's longitude on
        every month of a year, whose SGP4 lunar and solar terms differ (they
        displaced it by up to about 0.04 deg before the node was corrected).
        """
        for month in range(1, 13):
            for longitude in (-75.2, 0.0, 140.7):
                epoch = datetime(2026, month, 1, tzinfo=timezone.utc)
                orbit = GeosynchronousOrbit(longitude=longitude, epoch=epoch)
                position = orbit.to_gp_orbit().get_geographic_position(epoch)
                self.assertAlmostEqual(
                    position.longitude.degrees, longitude, delta=1e-6, msg=epoch
                )
                # the lunar and solar terms also move its latitude, by up to about
                # 0.004 deg, which the node does not correct
                self.assertAlmostEqual(position.latitude.degrees, 0, delta=0.01)

    def test_longitude_at_epoch_with_argument_of_latitude(self):
        """
        Test that a satellite offset in true anomaly (as in a constellation) and
        inclined lies at the longitude its argument of latitude gives, from the
        orbit's longitude at its node.
        """
        epoch = datetime(2026, 11, 8, 12, tzinfo=timezone.utc)
        orbit = GeosynchronousOrbit(
            longitude=10.0, inclination=5.0, true_anomaly=30.0, epoch=epoch
        )
        expected = 10.0 + np.degrees(
            np.arctan2(
                np.cos(np.radians(5)) * np.sin(np.radians(30)), np.cos(np.radians(30))
            )
        )
        position = orbit.to_gp_orbit().get_geographic_position(epoch)
        self.assertAlmostEqual(position.longitude.degrees, expected, delta=1e-6)

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

    def test_longitude_maintained(self):
        """
        Test that the orbit is propagated as maintained at its longitude, as
        by station-keeping: without it, the Earth's oblateness and the Moon
        and Sun move the longitude by tenths of a degree within a month.
        """
        gp_orbit = self.test_orbit.to_gp_orbit()
        self.assertTrue(gp_orbit.remove_drag)
        self.assertEqual(gp_orbit.repeat_cycle, timedelta(days=1))
        epoch = self.test_data["epoch"]
        times = [epoch + timedelta(days=d) for d in (0, 10.25, 30.5, 365, 3650)]
        longitudes = [
            gp_orbit.get_geographic_position(t).longitude.degrees for t in times
        ]
        for longitude in longitudes:
            self.assertAlmostEqual(longitude, self.test_data["longitude"], delta=0.01)


if __name__ == "__main__":
    unittest.main()
