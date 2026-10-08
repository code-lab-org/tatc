"""
Unit tests for the Earth orientation utilities: interpolated nutation
angles.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

# pylint: disable=protected-access

import unittest
import warnings

import numpy as np
from skyfield.nutationlib import iau2000a_radians

from tatc import config
from tatc.constants import timescale
from tatc.utils import earth_orientation


class TestNutationInterpolation(unittest.TestCase):
    """
    Unit tests for interpolated nutation angles.
    """

    def setUp(self):
        self.setting = config.get_rc().nutation_interpolation_minutes
        self.verified = earth_orientation._NUTATION_INTERPOLATION_VERIFIED

    def tearDown(self):
        config.get_rc().nutation_interpolation_minutes = self.setting
        earth_orientation._NUTATION_INTERPOLATION_VERIFIED = self.verified

    def test_interpolated_angles_match_iau2000a(self):
        """
        Test that nutation angles interpolated at 15 minute steps are within
        two microarcseconds of the IAU 2000A angles over a year.
        """
        config.get_rc().nutation_interpolation_minutes = 15
        jd = 2460676.5 + np.random.default_rng(0).uniform(0, 365, 20000)
        t = timescale.tt_jd(jd)
        earth_orientation._interpolate_nutation(t)
        self.assertIn("_nutation_angles_radians", vars(t))
        microarcseconds = np.degrees(1) * 3600e6
        for interpolated, exact in zip(
            t._nutation_angles_radians, iau2000a_radians(timescale.tt_jd(jd))
        ):
            self.assertLess(np.max(np.abs(interpolated - exact)) * microarcseconds, 2)

    def test_scalar_time(self):
        """
        Test that the nutation angles of a scalar time are interpolated as
        scalars, giving its sidereal time to well within a milliarcsecond.
        """
        config.get_rc().nutation_interpolation_minutes = 15
        t = timescale.tt_jd(2461041.123456)
        earth_orientation._interpolate_nutation(t)
        self.assertEqual(np.shape(t._nutation_angles_radians[0]), ())
        self.assertAlmostEqual(
            t.gast, timescale.tt_jd(2461041.123456).gast, delta=1e-10
        )

    def test_disabled(self):
        """
        Test that nutation angles are not interpolated when the runtime
        configuration is None.
        """
        config.get_rc().nutation_interpolation_minutes = None
        t = timescale.tt_jd(2461041.5 + np.arange(10) / 10)
        earth_orientation._interpolate_nutation(t)
        self.assertNotIn("_nutation_angles_radians", vars(t))

    def test_verified(self):
        """
        Test that Skyfield is verified to use the nutation angles set on a time.
        """
        self.assertTrue(earth_orientation._verify_nutation_interpolation())

    def test_not_used_if_not_verified(self):
        """
        Test that nutation angles are computed as usual, with a warning, if
        Skyfield is not verified to use those set on a time.
        """
        config.get_rc().nutation_interpolation_minutes = 15
        earth_orientation._NUTATION_INTERPOLATION_VERIFIED = None
        verify = earth_orientation._verify_nutation_interpolation
        earth_orientation._verify_nutation_interpolation = lambda: False
        try:
            t = timescale.tt_jd(2461041.5 + np.arange(10) / 10)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                earth_orientation._interpolate_nutation(t)
            self.assertNotIn("_nutation_angles_radians", vars(t))
            self.assertTrue(any("nutation" in str(w.message) for w in caught), caught)
        finally:
            earth_orientation._verify_nutation_interpolation = verify


if __name__ == "__main__":
    unittest.main()
