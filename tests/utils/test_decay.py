"""
Unit tests for the tatc.utils.decay module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest
from datetime import datetime, timezone

import numpy as np

from tatc.constants import EARTH_MEAN_RADIUS, EARTH_MU, EARTH_ROTATION_RATE
from tatc.utils import (
    SolarCycle,
    compute_atmospheric_density,
    compute_decay_rates,
    compute_drag_delta_v,
    compute_orbit_decay,
    compute_orbit_lifetime,
    interpolate_density,
)
from tatc.utils.atmosphere import YEAR_S


def get_gauss_rates(a, e, inclination, ballistic_coefficient, activity, count=200001):
    """
    Gets the rates of decay of an orbit by integrating Gauss's equations of
    the semimajor axis and eccentricity under a tangential drag over a
    revolution in true anomaly (trapezoidal rule), independently of the
    orbit-averaged formulas of `compute_decay_rates`.
    """
    anomaly = np.linspace(0, 2 * np.pi, count)
    p = a * (1 - e**2)
    radius = p / (1 + e * np.cos(anomaly))
    speed = np.sqrt(EARTH_MU * (2 / radius - 1 / a))
    factor = (
        1 - radius * EARTH_ROTATION_RATE * np.cos(np.radians(inclination)) / speed
    ) ** 2
    density = interpolate_density(radius - EARTH_MEAN_RADIUS, activity)
    drag = -density * factor * speed**2 / (2 * ballistic_coefficient)
    dt = radius**2 / np.sqrt(EARTH_MU * p)  # time per radian of true anomaly
    da = 2 * a**2 * speed * drag / EARTH_MU * dt
    de = 2 * (e + np.cos(anomaly)) * drag / speed * dt
    period = 2 * np.pi * np.sqrt(a**3 / EARTH_MU)
    step = anomaly[1] - anomaly[0]
    return (
        np.sum((da[1:] + da[:-1]) / 2) * step / period,
        np.sum((de[1:] + de[:-1]) / 2) * step / period,
    )


class TestDecay(unittest.TestCase):
    """
    Unit tests for the tatc.utils.decay module.
    """

    def test_rates_of_circular_orbits(self):
        """
        Test that a circular orbit shrinks at rho * F * sqrt(mu * a) / B.
        """
        a = EARTH_MEAN_RADIUS + 500e3
        density = compute_atmospheric_density(500e3)
        rate_a, rate_e = compute_decay_rates(a, 0, 90, 100)
        self.assertAlmostEqual(rate_a / (-density * np.sqrt(EARTH_MU * a) / 100), 1)
        self.assertEqual(rate_e, 0)
        # an atmosphere rotating with a prograde orbit drags it less
        self.assertGreater(compute_decay_rates(a, 0, 0, 100)[0], rate_a)
        self.assertLess(compute_decay_rates(a, 0, 180, 100)[0], rate_a)

    def test_rates_match_gauss_equations(self):
        """
        Test that the orbit-averaged rates of eccentric orbits match Gauss's
        equations integrated over a revolution.
        """
        for perigee, apogee, inclination in (
            (400e3, 410e3, 51.6),
            (300e3, 800e3, 28.5),
            (250e3, 1500e3, 98.0),
            (200e3, 35786e3, 7.0),
        ):
            with self.subTest(perigee=perigee, apogee=apogee):
                a = EARTH_MEAN_RADIUS + (perigee + apogee) / 2
                e = (apogee - perigee) / (2 * a)
                rate_a, rate_e = compute_decay_rates(a, e, inclination, 80, 0.7)
                gauss_a, gauss_e = get_gauss_rates(a, e, inclination, 80, 0.7)
                self.assertAlmostEqual(rate_a / gauss_a, 1, delta=2e-3)
                self.assertAlmostEqual(rate_e / gauss_e, 1, delta=2e-3)

    def test_rates_outside_the_atmosphere(self):
        """
        Test that orbits above the atmosphere do not decay and orbits below
        100 km have no rates.
        """
        self.assertEqual(
            compute_decay_rates(EARTH_MEAN_RADIUS + 1600e3, 0, 0, 100), (0.0, 0.0)
        )
        self.assertTrue(
            all(
                map(
                    math.isnan, compute_decay_rates(EARTH_MEAN_RADIUS + 90e3, 0, 0, 100)
                )
            )
        )

    def test_lifetime_of_circular_orbits(self):
        """
        Test that the lifetime of a circular orbit is the integral of the
        time to fall each meter, B / (rho * sqrt(mu * a)), down to 100 km.
        """
        altitudes = np.linspace(100e3, 450e3, 70001)
        density = compute_atmospheric_density(altitudes, "high")
        seconds_per_meter = 50 / (
            density * np.sqrt(EARTH_MU * (EARTH_MEAN_RADIUS + altitudes))
        )
        expected = (
            np.sum((seconds_per_meter[1:] + seconds_per_meter[:-1]) / 2)
            * (altitudes[1] - altitudes[0])
            / YEAR_S
        )
        lifetime = compute_orbit_lifetime(
            EARTH_MEAN_RADIUS + 450e3, 0, 90, 50, activity="high"
        )
        self.assertAlmostEqual(lifetime / expected, 1, delta=0.01)

    def test_lifetime_grows_with_altitude_and_ballistic_coefficient(self):
        """
        Test that lifetimes grow with altitude and the ballistic coefficient,
        and shrink with the activity.
        """
        lifetimes = [
            compute_orbit_lifetime(EARTH_MEAN_RADIUS + h, 0, 51.6, 100)
            for h in (300e3, 400e3, 500e3)
        ]
        self.assertTrue(np.all(np.diff(lifetimes) > 0))
        self.assertAlmostEqual(
            compute_orbit_lifetime(EARTH_MEAN_RADIUS + 400e3, 0, 51.6, 200)
            / lifetimes[1],
            2,
            delta=0.02,
        )
        self.assertLess(
            compute_orbit_lifetime(
                EARTH_MEAN_RADIUS + 400e3, 0, 51.6, 100, activity="high"
            ),
            lifetimes[1],
        )

    def test_lifetime_beyond_horizon(self):
        """
        Test that orbits that outlast the horizon (or the atmosphere) have
        infinite lifetimes.
        """
        self.assertEqual(
            compute_orbit_lifetime(EARTH_MEAN_RADIUS + 800e3, 0, 98, 100, horizon=50),
            math.inf,
        )
        self.assertEqual(
            compute_orbit_lifetime(EARTH_MEAN_RADIUS + 2000e3, 0, 98, 100), math.inf
        )

    def test_eccentric_orbits_circularize(self):
        """
        Test that an eccentric orbit circularizes as it decays: its apogee
        falls faster than its perigee, down to a perigee of 100 km.
        """
        a = EARTH_MEAN_RADIUS + 875e3
        decay = compute_orbit_decay(a, 1250e3 / (2 * a), 28.5, 100)
        self.assertTrue(np.all(np.diff(decay.eccentricity) <= 0))
        self.assertTrue(np.all(np.diff(decay.apogee_altitude) < 0))
        self.assertTrue(np.all(np.diff(decay.time) > 0))
        self.assertAlmostEqual(decay.perigee_altitude.iloc[-1], 100e3, delta=1)
        self.assertLess(decay.eccentricity.iloc[-1], 0.001)
        self.assertLess(
            decay.apogee_altitude.iloc[0] - decay.apogee_altitude.iloc[-1],
            decay.apogee_altitude.iloc[0],
        )

    def test_solar_cycle(self):
        """
        Test that the decay under a solar cycle requires a start, and falls
        between those at low and high activity.
        """
        a = EARTH_MEAN_RADIUS + 500e3
        with self.assertRaises(ValueError):
            compute_orbit_lifetime(a, 0, 98, 100, activity=SolarCycle())
        cycle = compute_orbit_lifetime(
            a,
            0,
            98,
            100,
            activity=SolarCycle(),
            start=datetime(2026, 1, 1, tzinfo=timezone.utc),
        )
        self.assertLess(compute_orbit_lifetime(a, 0, 98, 100, activity="high"), cycle)
        self.assertLess(cycle, compute_orbit_lifetime(a, 0, 98, 100, activity="low"))

    def test_invalid_inputs(self):
        """
        Test that perigees below 100 km and ballistic coefficients that are
        not positive raise errors.
        """
        with self.assertRaises(ValueError):
            compute_orbit_decay(EARTH_MEAN_RADIUS + 90e3, 0, 0, 100)
        with self.assertRaises(ValueError):
            compute_orbit_decay(EARTH_MEAN_RADIUS + 400e3, 0, 0, 0)
        with self.assertRaises(ValueError):
            compute_drag_delta_v(400e3, 0, 100, -1)

    def test_drag_delta_v(self):
        """
        Test that the velocity to hold an orbit is its drag acceleration
        times the duration, averaged over a solar cycle.
        """
        altitude, a = 450e3, EARTH_MEAN_RADIUS + 450e3
        speed = np.sqrt(EARTH_MU / a)
        density = compute_atmospheric_density(altitude, "moderate")
        expected = density * speed**2 / (2 * 100) * 3 * YEAR_S
        self.assertAlmostEqual(compute_drag_delta_v(altitude, 90, 100, 3) / expected, 1)
        self.assertEqual(compute_drag_delta_v(altitude, 90, 100, 0), 0)
        self.assertEqual(compute_drag_delta_v(1600e3, 90, 100, 3), 0)
        # over a whole cycle, the mean of the density over the cycle
        cycle = SolarCycle(period=11)
        start = cycle.minimum
        indices = cycle.get_activity_after(np.linspace(0, 11, 100001))
        mean = np.mean(interpolate_density(altitude, indices))
        self.assertAlmostEqual(
            compute_drag_delta_v(altitude, 90, 100, 11, activity=cycle, start=start)
            / (mean * speed**2 / (2 * 100) * 11 * YEAR_S),
            1,
            delta=1e-3,
        )
