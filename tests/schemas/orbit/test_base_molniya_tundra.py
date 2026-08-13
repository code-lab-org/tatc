"""
Unit tests for the tatc.schemas.orbit.base_molniya_tundra module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from tatc.constants import EARTH_J2_CRITICAL_INCLINATION
from tatc.schemas.orbit.base_molniya_tundra import MolniyaTundraOrbitBase


class TestMolniyaTundraOrbitBase(unittest.TestCase):
    """
    Unit tests for the
    tatc.schemas.orbit.base_molniya_tundra.MolniyaTundraOrbitBase schema.
    """

    def test_get_inclination_is_the_j2_critical_inclination(self):
        """
        Test that get_inclination returns Earth's J2 critical inclination
        (~63.4 degrees), the inclination at which J2 perturbation of the
        argument of perigee vanishes -- the defining property of a frozen
        Molniya/Tundra orbit.
        """
        o = MolniyaTundraOrbitBase(perigee_altitude=1000000)
        self.assertEqual(o.get_inclination(), EARTH_J2_CRITICAL_INCLINATION)

    def test_get_perigee_argument_northern_coverage_true(self):
        """
        Test that get_perigee_argument returns 270 degrees (apogee over
        the northern hemisphere) when northern_coverage is True (the
        default).
        """
        o = MolniyaTundraOrbitBase(perigee_altitude=1000000, northern_coverage=True)
        self.assertEqual(o.get_perigee_argument(), 270)

    def test_get_perigee_argument_northern_coverage_false(self):
        """
        Test that get_perigee_argument returns 90 degrees (apogee over the
        southern hemisphere) when northern_coverage is False.
        """
        o = MolniyaTundraOrbitBase(perigee_altitude=1000000, northern_coverage=False)
        self.assertEqual(o.get_perigee_argument(), 90)

    def test_get_right_ascension_ascending_node(self):
        """
        Test that get_right_ascension_ascending_node returns the
        right_ascension_ascending_node field directly.
        """
        o = MolniyaTundraOrbitBase(
            perigee_altitude=1000000, right_ascension_ascending_node=123.4
        )
        self.assertEqual(o.get_right_ascension_ascending_node(), 123.4)

    def test_get_orbit_period_not_implemented_on_base(self):
        """
        Test that get_orbit_period raises NotImplementedError on the base
        class, since concrete Molniya/Tundra subclasses must each define
        their own fixed orbit period (half vs. full sidereal day).
        """
        o = MolniyaTundraOrbitBase(perigee_altitude=1000000)
        with self.assertRaises(NotImplementedError):
            o.get_orbit_period()

    def test_get_semimajor_axis_not_implemented_on_base(self):
        """
        Test that get_semimajor_axis also surfaces NotImplementedError on
        the base class, since it depends on get_orbit_period().
        """
        o = MolniyaTundraOrbitBase(perigee_altitude=1000000)
        with self.assertRaises(NotImplementedError):
            o.get_semimajor_axis()

    def test_eccentricity_validator_is_a_no_op_on_bare_base(self):
        """
        Test that the shared eccentricity model_validator does not raise
        on the bare base class, even though get_eccentricity() itself
        would raise NotImplementedError there (since get_orbit_period is
        abstract) -- the validator must catch and skip this case rather
        than blocking construction of the bare base class, which the
        other tests in this file rely on.
        """
        # construction succeeding at all is the assertion; no concrete
        # period exists yet for the validator to check against
        MolniyaTundraOrbitBase(perigee_altitude=1000000)


if __name__ == "__main__":
    unittest.main()
