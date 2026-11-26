"""
Unit tests for the sizing of spacecraft of the tatc.sizing module.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

from pydantic import ValidationError

from tatc.schemas import (
    CircularOrbit,
    Instrument,
    InstrumentResources,
    Satellite,
    SpacecraftResources,
    WalkerConstellation,
)
from tatc.sizing import (
    DisposalRule,
    SpacecraftSizing,
    compute_delta_v_budget,
    compute_propellant_mass,
    size_spacecraft,
)
from tatc.utils import compute_orbit_lifetime

EPOCH = datetime(2027, 1, 1, tzinfo=timezone.utc)
PAYLOAD = Instrument(
    name="imager", field_of_regard=30, resources=InstrumentResources(mass=40, power=60)
)
SIZING = SpacecraftSizing(
    payload_mass_fraction=0.25,
    payload_power_fraction=0.4,
    design_life=5,
    drag_area=1,
    specific_impulse=220,
)


def get_satellite(altitude, resources=None, instruments=(PAYLOAD,)):
    """
    Gets a satellite in a circular orbit at an altitude (meters).
    """
    return Satellite(
        name="satellite",
        orbit=CircularOrbit(altitude=altitude, inclination=97, epoch=EPOCH),
        instruments=list(instruments),
        resources=resources,
    )


class TestSizeSpacecraft(unittest.TestCase):
    """
    Tests the sizing of spacecraft.
    """

    def test_sized_by_payload(self):
        """
        Tests that the dry mass and power are sized by the payload's
        fractions, and the launch mass by the propellant of the budget.
        """
        satellite = get_satellite(500e3)
        sized = size_spacecraft(satellite, SIZING)
        self.assertAlmostEqual(sized.resources.dry_mass, 160)
        self.assertAlmostEqual(sized.resources.power, 150)
        self.assertAlmostEqual(sized.power_margin, 0.6)
        self.assertAlmostEqual(sized.ballistic_coefficient, 160 / 2.2)
        budget = compute_delta_v_budget(satellite.orbit, 160 / 2.2, 5)
        self.assertEqual(sized.budget, budget)
        self.assertAlmostEqual(
            sized.propellant_mass, compute_propellant_mass(160, budget.total, 220)
        )
        self.assertAlmostEqual(sized.resources.launch_mass, 160 + sized.propellant_mass)
        self.assertEqual(sized.resources.design_life, 5)
        self.assertTrue(sized.feasible)

    def test_resources_take_precedence(self):
        """
        Tests that the spacecraft's own resources take precedence over the
        sizing, and an explicit design life over both.
        """
        resources = SpacecraftResources(
            dry_mass=300, launch_mass=400, power=500, design_life=7, drag_area=2
        )
        sized = size_spacecraft(get_satellite(500e3, resources), SIZING)
        self.assertEqual(sized.resources.dry_mass, 300)
        self.assertEqual(sized.resources.launch_mass, 400)
        self.assertEqual(sized.resources.power, 500)
        self.assertEqual(sized.resources.design_life, 7)
        self.assertAlmostEqual(sized.ballistic_coefficient, 300 / (2.2 * 2))
        self.assertEqual(
            size_spacecraft(
                get_satellite(500e3, resources), SIZING, design_life=3
            ).resources.design_life,
            3,
        )

    def test_missing_values(self):
        """
        Tests that sizing without a dry mass, design life, or drag raises
        errors.
        """
        with self.assertRaises(ValueError):
            size_spacecraft(get_satellite(500e3), SpacecraftSizing(design_life=5))
        with self.assertRaises(ValueError):
            size_spacecraft(
                get_satellite(
                    500e3, instruments=[Instrument(name="x", field_of_regard=30)]
                ),
                SIZING,
            )
        with self.assertRaises(ValueError):
            size_spacecraft(
                get_satellite(500e3), SIZING.model_copy(update={"design_life": None})
            )
        with self.assertRaises(ValueError):
            size_spacecraft(
                get_satellite(500e3), SIZING.model_copy(update={"drag_area": None})
            )
        with self.assertRaises(ValidationError):
            SpacecraftSizing(drag_area=1, ballistic_coefficient=50)

    def test_ballistic_coefficient(self):
        """
        Tests that a ballistic coefficient may be given in place of a drag
        area.
        """
        sizing = SIZING.model_copy(
            update={"drag_area": None, "ballistic_coefficient": 50.0}
        )
        self.assertEqual(
            size_spacecraft(get_satellite(500e3), sizing).ballistic_coefficient, 50
        )

    def test_altitude_window_without_propulsion(self):
        """
        Tests that without propulsion an orbit must last the design life and
        come down within the disposal limit after it.
        """
        sizing = SIZING.model_copy(update={"specific_impulse": None})
        low = size_spacecraft(get_satellite(400e3), sizing)
        lifetime = compute_orbit_lifetime(
            get_satellite(400e3).orbit.get_semimajor_axis(),
            0,
            97,
            low.ballistic_coefficient,
        )
        self.assertAlmostEqual(low.lifetime, lifetime)
        self.assertAlmostEqual(low.violations["decay"], (5 - lifetime) / 5)
        self.assertIsNone(low.budget)
        self.assertEqual(low.propellant_mass, 0)
        self.assertEqual(low.resources.launch_mass, low.resources.dry_mass)
        self.assertTrue(size_spacecraft(get_satellite(600e3), sizing).feasible)
        high = size_spacecraft(get_satellite(700e3), sizing)
        self.assertLess(high.lifetime, 200)
        self.assertAlmostEqual(
            high.violations["disposal"], (high.lifetime - 5 - 25) / 25
        )
        # beyond the horizon of lifetimes, the excess of the horizon
        self.assertAlmostEqual(
            size_spacecraft(get_satellite(900e3), sizing).violations["disposal"],
            (200 - 5 - 25) / 25,
        )
        controlled = sizing.model_copy(
            update={"disposal": DisposalRule(perigee_altitude=50e3)}
        )
        self.assertEqual(
            size_spacecraft(get_satellite(600e3), controlled).violations,
            {"disposal": 1.0},
        )
        self.assertTrue(
            size_spacecraft(
                get_satellite(800e3), sizing.model_copy(update={"disposal": None})
            ).feasible
        )

    def test_propellant_and_power_violations(self):
        """
        Tests that a launch mass that does not hold the propellant of the
        budget, and a payload that uses more than the power, are violations.
        """
        resources = SpacecraftResources(dry_mass=160, launch_mass=161, power=50)
        sized = size_spacecraft(get_satellite(400e3, resources), SIZING)
        self.assertGreater(sized.propellant_mass, 1)
        self.assertAlmostEqual(
            sized.violations["propellant"],
            (sized.propellant_mass - 1) / sized.propellant_mass,
        )
        self.assertAlmostEqual(sized.violations["power"], 60 / 50 - 1)
        self.assertFalse(sized.feasible)

    def test_constellation(self):
        """
        Tests that a constellation is sized as each of its members.
        """
        constellation = WalkerConstellation(
            name="walker",
            orbit=CircularOrbit(altitude=500e3, inclination=60, epoch=EPOCH),
            number_satellites=4,
            number_planes=2,
            instruments=[PAYLOAD],
        )
        member = constellation.generate_members()[1]
        self.assertEqual(
            size_spacecraft(constellation, SIZING), size_spacecraft(member, SIZING)
        )
