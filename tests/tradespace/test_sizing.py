"""
Unit tests for the resources and sizing of the satellites of mission plans.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import pandas as pd

from tatc.schemas import Instrument, InstrumentResources, Satellite, SpacecraftResources
from tatc.sizing import SpacecraftSizing, size_spacecraft
from tatc.tradespace import (
    Categorical,
    DeploymentSchedule,
    DesignSpace,
    Infeasible,
    MissionPlan,
    build_constellation,
    walker_space,
)
from tatc.tradespace.mission import Capability
from tatc.tradespace.plans import SIZING_COLUMNS

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


def _design(altitude=500e3, planes=2, per_plane=2):
    """The only design of a design space of a Walker constellation."""
    (design,) = DesignSpace(
        layout=walker_space(
            inclination=97.0, planes=planes, per_plane=per_plane, relative_spacing=0
        ),
        altitude=Categorical([altitude]),
    ).enumerate()
    return design


def _plan(design, plan, payload=(PAYLOAD,)):
    """The mission of a design by a plan, with its built satellites."""
    satellites = build_constellation(
        design.layout,
        payload if callable(payload) else list(payload),
        epoch=EPOCH,
        altitude=design.altitude,
    )
    capability = Capability(satellites, pd.DataFrame({"satellite": []}), len)
    return plan.plan_mission(design, design.layout, capability), satellites


def _with_resources(**resources):
    """A payload function that builds satellites with resources."""

    def build(name, orbit):
        return Satellite(
            name=name,
            orbit=orbit,
            instruments=[PAYLOAD],
            resources=SpacecraftResources(**resources),
        )

    return build


class TestResourceDefaults(unittest.TestCase):
    """
    Tests the masses and design lives that plans take from the resources of
    satellites.
    """

    def test_masses_and_lives(self):
        """
        Tests that a plan without masses or design lives takes the satellites'
        launch masses and design lives.
        """
        design = _design()
        schedule = DeploymentSchedule(years=10, launch_mass_capacity=400)
        mission, _ = _plan(
            design,
            MissionPlan(schedule),
            _with_resources(launch_mass=150, dry_mass=120, design_life=4),
        )
        self.assertEqual(mission.deployment.masses, {None: 150})
        explicit, _ = _plan(
            design,
            MissionPlan(
                DeploymentSchedule(years=10, launch_mass_capacity=400, design_life=4),
                masses=150,
            ),
        )
        self.assertEqual(mission.deployment.launches, explicit.deployment.launches)
        self.assertGreater(len(mission.deployment.launches), 2)
        self.assertIsNone(mission.sizing)

    def test_explicit_values_take_precedence(self):
        """
        Tests that a plan's masses and a schedule's design life take
        precedence over the satellites' resources.
        """
        design = _design()
        mission, _ = _plan(
            design,
            MissionPlan(DeploymentSchedule(years=10, design_life=20), masses=90),
            _with_resources(launch_mass=150, design_life=4),
        )
        self.assertEqual(mission.deployment.masses, {None: 90})
        self.assertEqual(mission.deployment.schedule.design_life, 20)
        # no replacements within the mission: a launch for each plane
        self.assertEqual(len(mission.deployment.launches), 2)

    def test_without_resources(self):
        """
        Tests that satellites without resources leave the plan unchanged.
        """
        mission, _ = _plan(_design(), MissionPlan(DeploymentSchedule()))
        self.assertEqual(mission.deployment.masses, {None: None})

    def test_disagreeing_lives(self):
        """
        Tests that satellites of a type with different design lives raise an
        error.
        """

        def build(name, orbit):
            life = 3 if name.endswith("1") else 5
            return Satellite(
                name=name,
                orbit=orbit,
                resources=SpacecraftResources(design_life=life),
            )

        with self.assertRaises(ValueError):
            _plan(_design(), MissionPlan(DeploymentSchedule(years=10)), build)


class TestSizing(unittest.TestCase):
    """
    Tests the sizing of the satellites of mission plans.
    """

    def test_sized_columns(self):
        """
        Tests that a sized plan takes the sized launch masses and design
        lives, and records the columns of the sizing.
        """
        design = _design()
        mission, satellites = _plan(
            design, MissionPlan(DeploymentSchedule(), sizing=SIZING)
        )
        sized = size_spacecraft(satellites[0], SIZING)
        self.assertEqual(mission.deployment.masses, {None: sized.resources.launch_mass})
        self.assertEqual(mission.deployment.schedule.design_life, {None: 5})
        self.assertEqual(tuple(mission.sizing), SIZING_COLUMNS)
        self.assertAlmostEqual(
            mission.sizing["launch mass"], sized.resources.launch_mass
        )
        self.assertAlmostEqual(mission.sizing["delta-v"], sized.budget.total)
        self.assertAlmostEqual(mission.sizing["power margin"], 0.6)
        record = mission._record()  # pylint: disable=protected-access
        for column in SIZING_COLUMNS:
            self.assertIn(column, record)

    def test_schedule_design_life(self):
        """
        Tests that the schedule's design life sizes the satellites.
        """
        design = _design()
        mission, satellites = _plan(
            design, MissionPlan(DeploymentSchedule(design_life=2), sizing=SIZING)
        )
        sized = size_spacecraft(satellites[0], SIZING, design_life=2)
        self.assertAlmostEqual(mission.sizing["delta-v"], sized.budget.total)

    def test_lower_orbits_are_heavier(self):
        """
        Tests that lower orbits need heavier satellites, so that launches
        limited in mass carry fewer of them.
        """
        schedule = DeploymentSchedule(launch_mass_capacity=340)
        plan = MissionPlan(schedule, sizing=SIZING)
        low, _ = _plan(_design(350e3, 1, 2), plan)
        high, _ = _plan(_design(550e3, 1, 2), plan)
        self.assertGreater(low.sizing["launch mass"], high.sizing["launch mass"])
        self.assertGreater(low.sizing["delta-v"], high.sizing["delta-v"])
        self.assertEqual(len(high.deployment.launches), 1)
        self.assertEqual(len(low.deployment.launches), 2)

    def test_infeasible_designs(self):
        """
        Tests that designs whose satellites violate their sizing, or exceed
        the mass capacity of a launch, are infeasible.
        """
        drifting = SIZING.model_copy(update={"specific_impulse": None})
        with self.assertRaises(Infeasible) as raised:
            _plan(_design(350e3), MissionPlan(DeploymentSchedule(), sizing=drifting))
        self.assertGreater(raised.exception.violation, 0)
        self.assertIn("decay", str(raised.exception))
        with self.assertRaises(Infeasible):
            _plan(
                _design(),
                MissionPlan(
                    DeploymentSchedule(launch_mass_capacity=100), sizing=SIZING
                ),
            )

    def test_sizing_of_the_design(self):
        """
        Tests that the sizing may be a function of the design.
        """
        design = _design()
        mission, _ = _plan(
            design,
            MissionPlan(
                DeploymentSchedule(),
                sizing=lambda d: SIZING.model_copy(update={"design_life": 2.0}),
            ),
        )
        self.assertEqual(mission.deployment.schedule.design_life, {None: 2.0})
        with self.assertRaises(TypeError):
            _plan(design, MissionPlan(DeploymentSchedule(), sizing=lambda d: 5))
