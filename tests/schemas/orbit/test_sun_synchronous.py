"""
Unit tests for the SunSynchronousOrbit schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, time, timezone

from pydantic import ValidationError

from tatc.constants import EARTH_MEAN_RADIUS
from tatc.schemas import SunSynchronousOrbit
from tatc.schemas.orbit.gp import GeneralPerturbationsOrbit


class TestSunSynchronousOrbit(unittest.TestCase):
    """
    Unit tests for the SunSynchronousOrbit schema.
    """

    def setUp(self):
        self.test_data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(10, 30),
            "equator_crossing_ascending": True,
        }
        self.test_orbit = SunSynchronousOrbit(**self.test_data)

    def test_good_data(self):
        """
        Test that the SunSynchronousOrbit schema correctly initializes with valid data.
        """
        good_data = {
            "mean_altitude": 400000,
            "true_anomaly": 10.0,
            "epoch": datetime(2022, 1, 1, 12, 0, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(10, 30),
            "equator_crossing_ascending": True,
        }
        o = SunSynchronousOrbit(**good_data)
        self.assertEqual(o.mean_altitude, good_data.get("mean_altitude"))
        self.assertEqual(o.true_anomaly, good_data.get("true_anomaly"))
        self.assertEqual(
            o.equator_crossing_time, good_data.get("equator_crossing_time")
        )
        self.assertEqual(
            o.equator_crossing_ascending, good_data.get("equator_crossing_ascending")
        )

    def test_equator_crossing_ascending_defaults_to_true(self):
        """
        Test that equator_crossing_ascending defaults to True when
        omitted.
        """
        o = SunSynchronousOrbit(mean_altitude=500000, equator_crossing_time=time(12))
        self.assertTrue(o.equator_crossing_ascending)

    def test_bad_mean_altitude_missing(self):
        """
        Test that the SunSynchronousOrbit schema raises a ValidationError
        when the required mean_altitude field is missing.
        """
        with self.assertRaises(ValidationError):
            SunSynchronousOrbit(equator_crossing_time=time(12))

    def test_bad_equator_crossing_time_missing(self):
        """
        Test that the SunSynchronousOrbit schema raises a ValidationError
        when the required equator_crossing_time field is missing.
        """
        with self.assertRaises(ValidationError):
            SunSynchronousOrbit(mean_altitude=500000)

    def test_bad_mean_altitude_too_large(self):
        """
        Test that a mean_altitude at or above 12,352,000 - EARTH_MEAN_RADIUS
        meters is rejected, since the sun-synchronous inclination formula's
        arccos argument goes out of its valid [-1, 1] domain at or beyond
        that reference semimajor axis.
        """
        with self.assertRaises(ValidationError):
            SunSynchronousOrbit(
                mean_altitude=12352000 - EARTH_MEAN_RADIUS,
                equator_crossing_time=time(12),
            )

    def test_get_derived_orbit_shifts_equator_crossing_time(self):
        """
        Test that get_derived_orbit shifts equator_crossing_time by 1 hour
        per 15 degrees of delta_raan (the documented conversion rate),
        in both directions.
        """
        later = self.test_orbit.get_derived_orbit(0, 15)
        self.assertEqual(later.equator_crossing_time, time(11, 30))
        earlier = self.test_orbit.get_derived_orbit(0, -15)
        self.assertEqual(earlier.equator_crossing_time, time(9, 30))

    def test_get_derived_orbit_wraps_equator_crossing_time_across_midnight(self):
        """
        Test that get_derived_orbit wraps equator_crossing_time correctly
        when the RAAN shift pushes it past midnight.
        """
        o = SunSynchronousOrbit(mean_altitude=500000, equator_crossing_time=time(23, 0))
        derived = o.get_derived_orbit(0, 30)
        self.assertEqual(derived.equator_crossing_time, time(1, 0))

    def test_get_eccentricity_and_perigee_argument_inherited_from_circular_base(self):
        """
        Test that get_eccentricity and get_perigee_argument are inherited
        from CircularOrbitBase (both 0, since a sun-synchronous orbit is
        circular), confirming the inheritance is functionally in effect
        for this class specifically.
        """
        self.assertEqual(self.test_orbit.get_eccentricity(), 0)
        self.assertEqual(self.test_orbit.get_perigee_argument(), 0)

    def test_get_inclination_matches_landsat8_published_value(self):
        """
        Test get_inclination against Landsat 8's published orbit: a
        705 km mean altitude sun-synchronous orbit is documented (USGS)
        as having a 98.2 degree inclination.

        See: https://www.usgs.gov/landsat-missions/landsat-8
        """
        o = SunSynchronousOrbit(
            mean_altitude=705000,
            equator_crossing_time=time(10, 0),
            equator_crossing_ascending=False,
        )
        self.assertAlmostEqual(o.get_inclination(), 98.2, delta=0.05)

    def test_get_inclination_matches_noaa20_published_value(self):
        """
        Test get_inclination against NOAA-20/JPSS-1's published orbit: an
        824 km mean altitude sun-synchronous orbit is documented (NOAA
        eoPortal) as having a 98.7 degree inclination.

        See: https://www.eoportal.org/satellite-missions/noaa-20
        """
        o = SunSynchronousOrbit(
            mean_altitude=824000,
            equator_crossing_time=time(13, 30),
            equator_crossing_ascending=True,
        )
        self.assertAlmostEqual(o.get_inclination(), 98.7, delta=0.05)

    def test_get_right_ascension_ascending_node_matches_real_noaa20_tle(self):
        """
        Test get_right_ascension_ascending_node against a real NOAA-20
        two-line element set (the same TLE used in
        docs/examples/ComputeCoverage.ipynb), rather than only a
        synthetic/derived example. NOAA-20's documented local time of
        ascending node is 13:30 +/- 10 minutes (see eoPortal), which
        corresponds to a RAAN tolerance of about +/- 2.5 degrees (15
        degrees per hour of local time); the actual TLE RAAN falls
        comfortably within that documented station-keeping deadband of
        the nominal 13:30 crossing time.
        """
        tle = [
            "1 43013U 17073A   22195.78278435  .00000038  00000+0  38919-4 0  9996",
            "2 43013  98.7169 133.9110 0001202  63.8768 296.2532 14.19561306241107",
        ]
        real_gp_orbit = GeneralPerturbationsOrbit.from_tle(tle)
        o = SunSynchronousOrbit(
            mean_altitude=824000,
            equator_crossing_time=time(13, 30),
            equator_crossing_ascending=True,
            epoch=real_gp_orbit.get_epoch(),
        )
        self.assertAlmostEqual(
            o.get_right_ascension_ascending_node(),
            real_gp_orbit.get_right_ascension_ascending_node(),
            delta=3.0,
        )

    def test_get_derived_orbit(self):
        """
        Test that the SunSynchronousOrbit schema correctly derives a new
        orbit with specified mean anomaly and RAAN offsets.
        """
        derived_orbit = self.test_orbit.get_derived_orbit(20, 10)
        self.assertAlmostEqual(
            derived_orbit.get_mean_anomaly(),
            self.test_orbit.get_mean_anomaly() + 20,
            delta=0.001,
        )
        self.assertAlmostEqual(
            derived_orbit.get_right_ascension_ascending_node(),
            self.test_orbit.get_right_ascension_ascending_node() + 10,
            delta=0.001,
        )

    def test_to_gp_orbit(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts
        to a general perturbations orbit.
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
        self.assertAlmostEqual(gp_orbit.get_inclination(), 97.7, delta=0.1)

    def test_to_gp_orbit_raan_ascending_equinox(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts
        to a general perturbations orbit with RAAN at ascending equinox.
        """
        data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 3, 20, 3, 49, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": True,
        }
        gp_orbit = SunSynchronousOrbit(**data).to_gp_orbit()
        self.assertAlmostEqual(
            min(
                gp_orbit.get_right_ascension_ascending_node(),
                360.0 - gp_orbit.get_right_ascension_ascending_node(),
            ),
            0.0,
            delta=0.25,
        )

    def test_to_gp_orbit_raan_descending_equinox(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts
        to a general perturbations orbit with RAAN at descending equinox.
        """
        data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 3, 20, 3, 49, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": False,
        }
        gp_orbit = SunSynchronousOrbit(**data).to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(), 180.0, delta=0.25
        )

    def test_to_gp_orbit_raan_ascending_solstice(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts
        to a general perturbations orbit with RAAN at ascending solstice.
        """
        data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 6, 21, 9, 14, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": True,
        }
        gp_orbit = SunSynchronousOrbit(**data).to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(), 90.0, delta=0.25
        )

    def test_to_gp_orbit_raan_descending_solstice(self):
        """
        Test that the SunSynchronousOrbit schema correctly converts
        to a general perturbations orbit with RAAN at descending solstice.
        """
        data = {
            "mean_altitude": 567000,
            "true_anomaly": 0.0,
            "epoch": datetime(2020, 6, 21, 9, 14, 0, tzinfo=timezone.utc),
            "equator_crossing_time": time(12),
            "equator_crossing_ascending": False,
        }
        gp_orbit = SunSynchronousOrbit(**data).to_gp_orbit()
        self.assertAlmostEqual(
            gp_orbit.get_right_ascension_ascending_node(), 270.0, delta=0.25
        )
