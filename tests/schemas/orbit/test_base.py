"""
Unit tests for the tatc.schemas.orbit.base module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

from tatc import config
from tatc.constants import EARTH_MEAN_RADIUS
from tatc.schemas import CircularOrbit
from tatc.schemas.orbit.base import OrbitBase
from tatc.utils.orbital import (
    semimajor_axis_to_mean_motion,
    semimajor_axis_to_orbit_period,
)


class TestOrbitBase(unittest.TestCase):
    """
    Unit tests for the tatc.schemas.orbit.base.OrbitBase schema.
    """

    def test_epoch_defaults_to_fixed_reference_timestamp(self):
        """
        Test that omitting epoch defaults to the fixed reference timestamp
        2020-01-01T00:00:00Z, not the current time. A prior version used
        datetime.now(), which pydantic evaluates once at class-definition
        time (module import), freezing every instance that omits epoch to
        the same, increasingly stale timestamp.
        """
        self.assertEqual(OrbitBase().epoch, datetime(2020, 1, 1, tzinfo=timezone.utc))

    def test_true_anomaly_defaults_to_zero(self):
        """
        Test that omitting true_anomaly defaults to 0 degrees.
        """
        self.assertEqual(OrbitBase().true_anomaly, 0)

    def test_get_mean_anomaly_treats_orbit_as_circular(self):
        """
        Test that get_mean_anomaly assumes zero eccentricity by default,
        so mean anomaly equals true anomaly exactly.
        """
        self.assertEqual(OrbitBase(true_anomaly=123.4).get_mean_anomaly(), 123.4)

    def test_get_semimajor_axis_not_implemented_on_base(self):
        """
        Test that get_semimajor_axis raises NotImplementedError on the
        bare base class, since OrbitBase has no universal way to derive
        it.
        """
        with self.assertRaises(NotImplementedError):
            OrbitBase().get_semimajor_axis()

    def test_get_inclination_not_implemented_on_base(self):
        """
        Test that get_inclination raises NotImplementedError on the bare
        base class.
        """
        with self.assertRaises(NotImplementedError):
            OrbitBase().get_inclination()

    def test_get_right_ascension_ascending_node_not_implemented_on_base(self):
        """
        Test that get_right_ascension_ascending_node raises
        NotImplementedError on the bare base class.
        """
        with self.assertRaises(NotImplementedError):
            OrbitBase().get_right_ascension_ascending_node()

    def test_get_eccentricity_not_implemented_on_base(self):
        """
        Test that get_eccentricity raises NotImplementedError on the bare
        base class.
        """
        with self.assertRaises(NotImplementedError):
            OrbitBase().get_eccentricity()

    def test_get_perigee_argument_not_implemented_on_base(self):
        """
        Test that get_perigee_argument raises NotImplementedError on the
        bare base class.
        """
        with self.assertRaises(NotImplementedError):
            OrbitBase().get_perigee_argument()

    def test_compute_gp_orbit_not_implemented_on_base(self):
        """
        Test that _compute_gp_orbit raises NotImplementedError on the bare
        base class, since concrete orbit subclasses must supply their own
        conversion logic.
        """
        with self.assertRaises(NotImplementedError):
            OrbitBase()._compute_gp_orbit()

    def test_to_gp_orbit_surfaces_not_implemented_on_base(self):
        """
        Test that to_gp_orbit() on the bare base class surfaces the same
        NotImplementedError, since it only adds caching around whatever
        _compute_gp_orbit produces.
        """
        with self.assertRaises(NotImplementedError):
            OrbitBase().to_gp_orbit()


class _OrbitWithFixedSemimajorAxis(OrbitBase):
    """
    Minimal OrbitBase subclass implementing only get_semimajor_axis, used
    to test OrbitBase's generic get_mean_altitude/get_mean_motion/
    get_orbit_period defaults in isolation from any production subclass's
    own overrides.
    """

    semimajor_axis: float

    def get_semimajor_axis(self) -> float:
        return self.semimajor_axis


class TestOrbitBaseDerivedDefaults(unittest.TestCase):
    """
    Unit tests for OrbitBase's generic get_mean_altitude/get_mean_motion/
    get_orbit_period defaults, each derived solely from get_semimajor_axis.
    """

    def setUp(self):
        self.orbit = _OrbitWithFixedSemimajorAxis(semimajor_axis=7000000)

    def test_get_mean_altitude_derived_from_semimajor_axis(self):
        """
        Test that get_mean_altitude defaults to semimajor axis minus
        Earth's mean radius.
        """
        self.assertAlmostEqual(
            self.orbit.get_mean_altitude(), 7000000 - EARTH_MEAN_RADIUS, delta=0.01
        )

    def test_get_mean_motion_derived_from_semimajor_axis(self):
        """
        Test that get_mean_motion defaults to the standard
        semimajor-axis-to-mean-motion conversion.
        """
        self.assertAlmostEqual(
            self.orbit.get_mean_motion(),
            semimajor_axis_to_mean_motion(7000000),
            delta=1e-9,
        )

    def test_get_orbit_period_derived_from_semimajor_axis(self):
        """
        Test that get_orbit_period defaults to the standard
        semimajor-axis-to-orbit-period conversion.
        """
        self.assertAlmostEqual(
            self.orbit.get_orbit_period().total_seconds(),
            semimajor_axis_to_orbit_period(7000000),
            delta=1e-6,
        )


class TestOrbitBaseToGpOrbitCaching(unittest.TestCase):
    """
    Unit tests for the lazy-load caching behavior of
    OrbitBase.to_gp_orbit, shared by every concrete orbit subclass
    (CircularOrbit, KeplerianOrbit, SunSynchronousOrbit, MolniyaOrbit,
    TundraOrbit). Uses CircularOrbit as a concrete vehicle to exercise the
    shared implementation, since OrbitBase itself has no _compute_gp_orbit
    to cache.
    """

    def setUp(self):
        self.orbit = CircularOrbit(mean_altitude=500000)
        self._original_lazy_load = config.rc.gp_orbit_lazy_load

    def tearDown(self):
        config.rc.gp_orbit_lazy_load = self._original_lazy_load

    def test_repeat_call_reuses_cached_result(self):
        """
        Test that calling to_gp_orbit() twice with lazy_load=True (the
        default) returns the identical cached object rather than
        recomputing.
        """
        first = self.orbit.to_gp_orbit()
        second = self.orbit.to_gp_orbit()
        self.assertIs(first, second)

    def test_lazy_load_false_forces_recomputation(self):
        """
        Test that lazy_load=False always recomputes a fresh gp orbit, even
        if a cached result already exists.
        """
        first = self.orbit.to_gp_orbit()
        second = self.orbit.to_gp_orbit(lazy_load=False)
        self.assertIsNot(first, second)

    def test_lazy_load_false_still_updates_the_cache(self):
        """
        Test that a lazy_load=False call still stores its result in the
        cache, so a subsequent lazy_load=True call reuses that new result
        rather than the original.
        """
        first = self.orbit.to_gp_orbit()
        second = self.orbit.to_gp_orbit(lazy_load=False)
        third = self.orbit.to_gp_orbit()
        self.assertIs(third, second)
        self.assertIsNot(third, first)

    def test_lazy_load_none_follows_config_rc_true(self):
        """
        Test that lazy_load=None (the default) follows config.rc's
        gp_orbit_lazy_load setting when it is True: repeat calls reuse the
        cached result.
        """
        config.rc.gp_orbit_lazy_load = True
        first = self.orbit.to_gp_orbit()
        second = self.orbit.to_gp_orbit()
        self.assertIs(first, second)

    def test_lazy_load_none_follows_config_rc_false(self):
        """
        Test that lazy_load=None (the default) follows config.rc's
        gp_orbit_lazy_load setting when it is False: repeat calls always
        recompute.
        """
        config.rc.gp_orbit_lazy_load = False
        first = self.orbit.to_gp_orbit()
        second = self.orbit.to_gp_orbit()
        self.assertIsNot(first, second)


if __name__ == "__main__":
    unittest.main()
