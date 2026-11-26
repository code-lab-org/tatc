"""
Unit tests for the velocity budgets of the tatc.sizing module.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest
from datetime import datetime, timedelta, timezone

from pydantic import ValidationError

from tatc.constants import EARTH_MEAN_RADIUS, EARTH_MU, STANDARD_GRAVITY
from tatc.schemas import CircularOrbit, KeplerianOrbit
from tatc.sizing import (
    DeltaVBudget,
    DisposalRule,
    compute_delta_v_budget,
    compute_disposal_delta_v,
    compute_propellant_mass,
)
from tatc.utils import (
    SolarCycle,
    compute_drag_delta_v,
    compute_hohmann_delta_v,
    compute_orbit_lifetime,
)
from tatc.utils.atmosphere import YEAR_S

EPOCH = datetime(2027, 1, 1, tzinfo=timezone.utc)


class TestPropellant(unittest.TestCase):
    """
    Tests the propellant of velocities and the transfers between orbits.
    """

    def test_rocket_equation(self):
        """
        Tests that the propellant doubles the mass at a velocity of the
        specific impulse times standard gravity times ln(2).
        """
        delta_v = 220 * STANDARD_GRAVITY * math.log(2)
        self.assertAlmostEqual(compute_propellant_mass(150, delta_v, 220), 150)
        self.assertEqual(compute_propellant_mass(150, 0, 220), 0)
        with self.assertRaises(ValueError):
            compute_propellant_mass(150, 10, 0)
        with self.assertRaises(ValueError):
            compute_propellant_mass(-1, 10, 220)

    def test_hohmann_transfer(self):
        """
        Tests that a Hohmann transfer from low Earth orbit (300 km) to
        geostationary orbit takes about 3.9 km/s, either way, and none
        between equal orbits.
        """
        low, geostationary = EARTH_MEAN_RADIUS + 300e3, 42164e3
        self.assertAlmostEqual(
            compute_hohmann_delta_v(low, geostationary), 3.9e3, delta=30
        )
        self.assertAlmostEqual(
            compute_hohmann_delta_v(low, geostationary),
            compute_hohmann_delta_v(geostationary, low),
        )
        self.assertEqual(compute_hohmann_delta_v(low, low), 0)


class TestDisposal(unittest.TestCase):
    """
    Tests the velocity to dispose of spacecraft.
    """

    def test_natural_disposal(self):
        """
        Tests that orbits that drag brings down within the limit need no
        velocity, and that a stricter limit needs some.
        """
        self.assertEqual(compute_disposal_delta_v(500e3, 98, 100), 0)
        self.assertGreater(
            compute_disposal_delta_v(500e3, 98, 100, rule=DisposalRule(limit=2)), 0
        )

    def test_lowered_perigee(self):
        """
        Tests that the perigee lowered by the velocity brings the orbit down
        within the limit, and that a perigee 2 km higher does not.
        """
        rule = DisposalRule(limit=10)
        altitude = 700e3
        delta_v = compute_disposal_delta_v(altitude, 98, 100, rule=rule)
        radius = EARTH_MEAN_RADIUS + altitude
        speed = math.sqrt(EARTH_MU / radius) - delta_v
        # the perigee of the ellipse whose apogee speed is reduced by delta_v
        a = 1 / (2 / radius - speed**2 / EARTH_MU)
        perigee_radius = 2 * a - radius

        def lifetime(perigee_radius):
            return compute_orbit_lifetime(
                (radius + perigee_radius) / 2,
                (radius - perigee_radius) / (radius + perigee_radius),
                98,
                100,
                horizon=20,
            )

        self.assertLessEqual(lifetime(perigee_radius), rule.limit)
        self.assertGreater(lifetime(perigee_radius + 2e3), rule.limit)

    def test_controlled_reentry(self):
        """
        Tests that a controlled reentry lowers the perigee to its altitude,
        and that orbits above the protected region are not disposed of.
        """
        rule = DisposalRule(perigee_altitude=50e3)
        radius = EARTH_MEAN_RADIUS + 700e3
        perigee_radius = EARTH_MEAN_RADIUS + 50e3
        expected = math.sqrt(EARTH_MU / radius) - math.sqrt(
            EARTH_MU * (2 / radius - 2 / (radius + perigee_radius))
        )
        self.assertAlmostEqual(
            compute_disposal_delta_v(700e3, 98, 100, rule=rule), expected
        )
        self.assertEqual(compute_disposal_delta_v(2500e3, 98, 100, rule=rule), 0)
        self.assertEqual(compute_disposal_delta_v(2500e3, 98, 100), 0)
        with self.assertRaises(ValidationError):
            DisposalRule(limit=0)


class TestDeltaVBudget(unittest.TestCase):
    """
    Tests the velocity budgets of orbits.
    """

    def test_parts(self):
        """
        Tests the parts of a budget and its total.
        """
        orbit = CircularOrbit(altitude=450e3, inclination=51.6, epoch=EPOCH)
        budget = compute_delta_v_budget(
            orbit, 100, 5, injection_altitude=300e3, other=10, margin=0.1
        )
        self.assertAlmostEqual(
            budget.raising,
            compute_hohmann_delta_v(
                EARTH_MEAN_RADIUS + 300e3, EARTH_MEAN_RADIUS + 450e3
            ),
        )
        self.assertAlmostEqual(budget.drag, compute_drag_delta_v(450e3, 51.6, 100, 5))
        self.assertEqual(budget.disposal, 0)
        self.assertEqual(budget.other, 10)
        self.assertAlmostEqual(
            budget.margin, 0.1 * (budget.raising + budget.drag + budget.other)
        )
        self.assertAlmostEqual(budget.total, 1.1 * (budget.raising + budget.drag + 10))
        self.assertIn("total", budget.model_dump())

    def test_disposal_after_life(self):
        """
        Tests that the disposal is that at the end of the life, under a solar
        cycle, and none without a rule.
        """
        cycle = SolarCycle()
        orbit = CircularOrbit(altitude=750e3, inclination=98, epoch=EPOCH)
        budget = compute_delta_v_budget(orbit, 100, 3, activity=cycle)
        self.assertAlmostEqual(
            budget.disposal,
            compute_disposal_delta_v(
                750e3,
                98,
                100,
                activity=cycle,
                start=EPOCH + timedelta(seconds=3 * YEAR_S),
            ),
        )
        self.assertGreater(budget.disposal, 0)
        self.assertEqual(
            compute_delta_v_budget(orbit, 100, 3, disposal=None).disposal, 0
        )

    def test_invalid_orbits(self):
        """
        Tests that budgets require near-circular orbits and non-negative
        values.
        """
        orbit = KeplerianOrbit(
            semimajor_axis=EARTH_MEAN_RADIUS + 1000e3,
            eccentricity=0.05,
            inclination=60,
            right_ascension_ascending_node=0,
            perigee_argument=0,
            true_anomaly=0,
            epoch=EPOCH,
        )
        with self.assertRaises(ValueError):
            compute_delta_v_budget(orbit, 100, 5)
        circular = CircularOrbit(altitude=500e3, inclination=60, epoch=EPOCH)
        with self.assertRaises(ValueError):
            compute_delta_v_budget(circular, 100, -1)
        with self.assertRaises(ValidationError):
            DeltaVBudget(drag=-1)
