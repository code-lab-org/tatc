"""
Unit tests for the tatc.utils.orbital module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

import numpy as np

from tatc import constants
from tatc.utils import mean_anomaly_to_true_anomaly, true_anomaly_to_mean_anomaly
from tatc.utils.orbital import (
    compute_apoapsis_radius,
    compute_ground_inertial_velocity,
    compute_ground_surface_velocity,
    compute_j2_aop_rate,
    compute_j2_mean_motion_rate,
    compute_j2_raan_rate,
    compute_orbit_inertial_velocity,
    mean_motion_to_orbit_period,
    mean_motion_to_semimajor_axis,
    semimajor_axis_to_mean_motion,
    semimajor_axis_to_orbit_period,
)


class TestOrbital(unittest.TestCase):  # pylint: disable=too-many-public-methods
    """
    Unit tests for the tatc.utils.orbital module.
    """
    def test_compute_orbit_inertial_velocity_iss(self):
        """
        Test against the commonly cited ISS orbital speed of
        approximately 7.66 km/s at its typical ~408 km altitude
        (per NASA ISS fact sheets).
        """
        self.assertAlmostEqual(
            compute_orbit_inertial_velocity(408000), 7660, delta=10
        )

    def test_compute_orbit_inertial_velocity_geo(self):
        """
        Test against the well-known geostationary orbital speed of
        approximately 3.07 km/s at ~35,786 km altitude.
        """
        self.assertAlmostEqual(
            compute_orbit_inertial_velocity(35786000), 3075, delta=5
        )

    def test_compute_orbit_inertial_velocity_surface(self):
        """
        Test against the classic "first cosmic velocity" of
        approximately 7.9 km/s for a circular orbit at the Earth's
        surface (zero altitude).
        """
        self.assertAlmostEqual(
            compute_orbit_inertial_velocity(0), 7910, delta=5
        )

    def test_compute_orbit_inertial_velocity_decreases_with_altitude(self):
        """
        Test that the inertial orbit velocity decreases monotonically as
        altitude increases (higher, slower orbits).
        """
        velocities = [
            compute_orbit_inertial_velocity(altitude)
            for altitude in (0, 400000, 705000, 20000000, 35786000)
        ]
        self.assertEqual(velocities, sorted(velocities, reverse=True))

    def test_compute_ground_inertial_velocity_matches_orbit_at_same_radius(self):
        """
        Test that the ground point velocity equals the orbital velocity
        exactly when projected to the satellite's own altitude (the
        projected point coincides with the satellite).
        """
        self.assertAlmostEqual(
            compute_ground_inertial_velocity(408000, 408000),
            compute_orbit_inertial_velocity(408000),
            delta=1e-9,
        )

    def test_compute_ground_inertial_velocity_slower_than_orbit_at_sea_level(self):
        """
        Test that the ground point velocity projected to sea level (zero
        elevation) is slower than the satellite's own orbital velocity,
        since it shares the same angular velocity at a smaller radius.
        """
        self.assertLess(
            compute_ground_inertial_velocity(408000, 0),
            compute_orbit_inertial_velocity(408000),
        )

    def test_compute_ground_inertial_velocity_increases_with_elevation(self):
        """
        Test that, for a fixed altitude, the ground point velocity
        increases monotonically with elevation.
        """
        velocities = [
            compute_ground_inertial_velocity(408000, elevation)
            for elevation in (-1000, 0, 100000, 300000, 408000)
        ]
        self.assertEqual(velocities, sorted(velocities))

    def test_compute_ground_inertial_velocity_decreases_with_altitude(self):
        """
        Test that, for a fixed (sea-level) elevation, the ground point
        velocity decreases monotonically as altitude increases.
        """
        velocities = [
            compute_ground_inertial_velocity(altitude, 0)
            for altitude in (200000, 408000, 705000, 20000000)
        ]
        self.assertEqual(velocities, sorted(velocities, reverse=True))

    def test_compute_ground_surface_velocity_equatorial_orbit_at_equator(self):
        """
        Test that an equatorial orbit's ground track, which moves purely
        eastward, has its ground-relative speed reduced by exactly Earth's
        eastward rotation speed at the equator.
        """
        v_inertial = compute_ground_inertial_velocity(408000, 0)
        v_earth_equator = (
            2 * np.pi / constants.EARTH_SIDEREAL_DAY_S * constants.EARTH_MEAN_RADIUS
        )
        self.assertAlmostEqual(
            compute_ground_surface_velocity(408000, 0, 0, 0),
            v_inertial - v_earth_equator,
            delta=1.0,
        )

    def test_compute_ground_surface_velocity_polar_orbit_at_equator(self):
        """
        Test that a polar orbit's ground track, which moves purely
        northward as it crosses the equator, combines with Earth's
        (perpendicular) eastward rotation speed in quadrature rather than
        by simple addition or subtraction.
        """
        v_inertial = compute_ground_inertial_velocity(408000, 0)
        v_earth_equator = (
            2 * np.pi / constants.EARTH_SIDEREAL_DAY_S * constants.EARTH_MEAN_RADIUS
        )
        self.assertAlmostEqual(
            compute_ground_surface_velocity(408000, 0, 90, 0),
            (v_inertial**2 + v_earth_equator**2) ** 0.5,
            delta=1.0,
        )

    def test_compute_ground_surface_velocity_landsat(self):
        """
        Test against the commonly cited Landsat ground track velocity of
        approximately 6.75 km/s (per USGS Landsat documentation), using
        Landsat's published altitude (705 km) and inclination (98.2
        degrees, sun-synchronous).
        """
        self.assertAlmostEqual(
            compute_ground_surface_velocity(705000, 0, 98.2, 0), 6750, delta=100
        )

    def test_compute_ground_surface_velocity_increases_with_inclination(self):
        """
        Test that, at the equator, the ground surface velocity increases
        monotonically with inclination from the equatorial (slowest) to
        the polar (fastest) case.
        """
        velocities = [
            compute_ground_surface_velocity(408000, 0, inclination, 0)
            for inclination in (0, 15, 30, 45, 60, 75, 90)
        ]
        self.assertEqual(velocities, sorted(velocities))

    def test_compute_ground_surface_velocity_pole_does_not_raise(self):
        """
        Test that evaluating at latitude = 90 degrees (a pole, where the
        azimuth formula's denominator approaches zero) does not raise an
        error, and instead saturates via clipping.
        """
        compute_ground_surface_velocity(408000, 0, 45, 90)

    def test_semimajor_axis_to_mean_motion_geo(self):
        """
        Test against the geostationary mean motion: 360 degrees per
        sidereal day, at the well-known geostationary geocentric
        semimajor axis of approximately 42,164 km.
        """
        self.assertAlmostEqual(
            semimajor_axis_to_mean_motion(42164000),
            360 / constants.EARTH_SIDEREAL_DAY_S,
            delta=1e-6,
        )

    def test_semimajor_axis_to_mean_motion_inverts_mean_motion_to_semimajor_axis(
        self,
    ):
        """
        Test that semimajor_axis_to_mean_motion is the exact inverse of
        mean_motion_to_semimajor_axis across a range of orbit altitudes.
        """
        for semimajor_axis in (6779008, 7076008, 8071008, 42164000):
            mean_motion = semimajor_axis_to_mean_motion(semimajor_axis)
            self.assertAlmostEqual(
                mean_motion_to_semimajor_axis(mean_motion),
                semimajor_axis,
                delta=1e-3,
            )

    def test_mean_motion_to_orbit_period_is_its_own_inverse(self):
        """
        Test that mean_motion_to_orbit_period is an involution: applying
        it twice recovers the original value, since a mean motion
        (degrees/second) and its orbital period (seconds) are each 360
        divided by the other.
        """
        mean_motion = 0.0648
        self.assertAlmostEqual(
            mean_motion_to_orbit_period(mean_motion_to_orbit_period(mean_motion)),
            mean_motion,
            delta=1e-12,
        )

    def test_semimajor_axis_to_orbit_period_geo(self):
        """
        Test against the geostationary orbital period, which by
        definition equals Earth's sidereal rotation period, at the
        well-known geostationary geocentric semimajor axis of
        approximately 42,164 km.
        """
        self.assertAlmostEqual(
            semimajor_axis_to_orbit_period(42164000),
            constants.EARTH_SIDEREAL_DAY_S,
            delta=1.0,
        )

    def test_semimajor_axis_to_orbit_period_iss(self):
        """
        Test against the commonly cited ISS orbital period of
        approximately 92.68 minutes, at its typical ~408 km altitude.
        """
        self.assertAlmostEqual(
            semimajor_axis_to_orbit_period(constants.EARTH_MEAN_RADIUS + 408000),
            92.68 * 60,
            delta=30,
        )

    def test_semimajor_axis_to_orbit_period_increases_with_semimajor_axis(self):
        """
        Test that the orbital period increases monotonically with the
        semimajor axis.
        """
        periods = [
            semimajor_axis_to_orbit_period(semimajor_axis)
            for semimajor_axis in (6779008, 7076008, 8071008, 42164000)
        ]
        self.assertEqual(periods, sorted(periods))

    def test_compute_apoapsis_radius_circular_orbit(self):
        """
        Test that a circular orbit's (zero eccentricity) apoapsis radius
        equals its semimajor axis exactly.
        """
        semimajor_axis = constants.EARTH_MEAN_RADIUS + 705000
        self.assertEqual(compute_apoapsis_radius(semimajor_axis, 0), semimajor_axis)

    def test_compute_apoapsis_radius_matches_perigee_apogee_definition(self):
        """
        Test that compute_apoapsis_radius recovers the apogee radius that
        defines the semimajor axis and eccentricity in the first place:
        semimajor axis is the mean of the perigee/apogee radii, and
        eccentricity is their normalized difference.
        """
        perigee_radius = constants.EARTH_MEAN_RADIUS + 550000
        apogee_radius = constants.EARTH_MEAN_RADIUS + 39900000
        semimajor_axis = (perigee_radius + apogee_radius) / 2
        eccentricity = (apogee_radius - perigee_radius) / (apogee_radius + perigee_radius)
        self.assertAlmostEqual(
            compute_apoapsis_radius(semimajor_axis, eccentricity),
            apogee_radius,
            delta=1e-6,
        )

    def test_compute_apoapsis_radius_molniya_orbit(self):
        """
        Test against the published classic Molniya orbit: a semimajor axis
        of approximately 26,600 km and eccentricity of approximately 0.74
        yields an apogee altitude of approximately 40,000 km.

        See: https://en.wikipedia.org/wiki/Molniya_orbit
        """
        apogee_altitude = (
            compute_apoapsis_radius(26600000, 0.74) - constants.EARTH_MEAN_RADIUS
        )
        self.assertAlmostEqual(apogee_altitude, 40000000, delta=500000)

    def test_compute_j2_raan_rate_sun_synchronous(self):
        """
        Test against the sun-synchronous condition: at Landsat/Terra/Aqua's
        published altitude (705 km) and inclination (98.2 degrees), the
        RAAN precession rate should match the Earth's mean motion around
        the Sun (360 degrees per 365.2422-day year) -- this is exactly why
        that inclination was chosen for those missions.
        """
        semimajor_axis = constants.EARTH_MEAN_RADIUS + 705000
        sun_synchronous_rate = 360 / 365.2422 / 86400
        self.assertAlmostEqual(
            compute_j2_raan_rate(semimajor_axis, 98.2, 0),
            sun_synchronous_rate,
            delta=1e-7,
        )

    def test_compute_j2_raan_rate_zero_for_polar_orbit(self):
        """
        Test that a polar orbit (inclination = 90 degrees) has zero RAAN
        precession rate.
        """
        semimajor_axis = constants.EARTH_MEAN_RADIUS + 705000
        self.assertAlmostEqual(
            compute_j2_raan_rate(semimajor_axis, 90, 0), 0.0, delta=1e-15
        )

    def test_compute_j2_raan_rate_sign_by_inclination(self):
        """
        Test that the RAAN precession rate is negative (westward
        regression) for prograde orbits and positive for retrograde
        orbits.
        """
        semimajor_axis = constants.EARTH_MEAN_RADIUS + 705000
        self.assertLess(compute_j2_raan_rate(semimajor_axis, 51.6, 0), 0.0)
        self.assertGreater(compute_j2_raan_rate(semimajor_axis, 98.2, 0), 0.0)

    def test_compute_j2_raan_rate_increases_with_eccentricity(self):
        """
        Test that the magnitude of the RAAN precession rate increases
        monotonically with eccentricity, for a fixed semimajor axis and
        inclination.
        """
        semimajor_axis = constants.EARTH_MEAN_RADIUS + 705000
        magnitudes = [
            abs(compute_j2_raan_rate(semimajor_axis, 51.6, eccentricity))
            for eccentricity in (0, 0.1, 0.3, 0.5)
        ]
        self.assertEqual(magnitudes, sorted(magnitudes))

    def test_compute_j2_aop_rate_zero_at_critical_inclination(self):
        """
        Test against the well-known critical inclination (~63.43 degrees,
        arccos(1/sqrt(5))) at which the argument of periapsis rate is
        exactly zero -- the reason Molniya-type orbits use this
        inclination to keep their periapsis location "frozen".
        """
        critical_inclination = np.degrees(np.arccos(1 / np.sqrt(5)))
        semimajor_axis = 26600000
        self.assertAlmostEqual(
            compute_j2_aop_rate(semimajor_axis, critical_inclination, 0.74),
            0.0,
            delta=1e-15,
        )

    def test_compute_j2_aop_rate_sign_by_inclination(self):
        """
        Test that the argument of periapsis rate is positive below the
        critical inclination (~63.43 degrees) and negative above it.
        """
        semimajor_axis = constants.EARTH_MEAN_RADIUS + 705000
        self.assertGreater(compute_j2_aop_rate(semimajor_axis, 45, 0), 0.0)
        self.assertLess(compute_j2_aop_rate(semimajor_axis, 90, 0), 0.0)

    def test_compute_j2_mean_motion_rate_nonzero_at_critical_inclination(self):
        """
        Test that, unlike compute_j2_aop_rate, the correction to mean
        anomaly's rate of advance is NOT zero at the critical inclination
        (~63.43 degrees). This is the key distinction that makes this
        correction (not the argument of periapsis rate) the correct one
        for computing a J2-corrected orbital period for Molniya/Tundra
        orbits, which always use the critical inclination.
        """
        critical_inclination = np.degrees(np.arccos(1 / np.sqrt(5)))
        semimajor_axis = 26600000
        self.assertNotAlmostEqual(
            compute_j2_mean_motion_rate(semimajor_axis, critical_inclination, 0.74),
            0.0,
            delta=1e-8,
        )

    def test_compute_j2_mean_motion_rate_zero_crossing(self):
        """
        Test that the correction vanishes at its own zero-crossing
        inclination (~54.74 degrees, arccos(1/sqrt(3)), where
        3*cos^2(inclination) - 1 = 0) -- a different inclination than
        compute_j2_aop_rate's zero crossing, confirming this is a
        genuinely distinct correction term.
        """
        zero_crossing_inclination = np.degrees(np.arccos(1 / np.sqrt(3)))
        semimajor_axis = constants.EARTH_MEAN_RADIUS + 705000
        self.assertAlmostEqual(
            compute_j2_mean_motion_rate(semimajor_axis, zero_crossing_inclination, 0),
            0.0,
            delta=1e-9,
        )

    def test_compute_j2_mean_motion_rate_sign_by_inclination(self):
        """
        Test that the correction is positive below the ~54.74 degree
        zero-crossing inclination, and negative above it.
        """
        semimajor_axis = constants.EARTH_MEAN_RADIUS + 705000
        self.assertGreater(compute_j2_mean_motion_rate(semimajor_axis, 0, 0), 0.0)
        self.assertLess(compute_j2_mean_motion_rate(semimajor_axis, 90, 0), 0.0)

    def test_mean_anomaly_to_true_anomaly(self):
        """
        Test that the mean anomaly can be converted to true anomaly.
        """
        self.assertAlmostEqual(
            mean_anomaly_to_true_anomaly(78.940629, 0.0001492), 78.95065818, delta=0.01
        )

    def test_true_anomaly_to_mean_anomaly(self):
        """
        Test that the true anomaly can be converted to mean anomaly.
        """
        self.assertAlmostEqual(
            true_anomaly_to_mean_anomaly(78.95065818, 0.0001492), 78.940629, delta=0.01
        )
