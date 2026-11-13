"""
Unit tests for the deployment of constellations and capability over its
phases.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    DeploymentSchedule,
    DesignSpace,
    Exponential,
    FailureModel,
    Real,
    Requirement,
    RestorationPolicy,
    sso_planes_space,
    walker_space,
)
from tatc.tradespace.mission import Phase, evaluate_phases, summarize_phases


def _walker(planes=3, per_plane=4):
    (design,) = DesignSpace(
        layout=walker_space(
            inclination=55.0, planes=planes, per_plane=per_plane, relative_spacing=0
        )
    ).enumerate()
    return design.layout


def _names(planes, per_plane):
    """The names of the satellites of a Walker layout, in build order."""
    return [f"satellite {k + 1}" for k in range(planes * per_plane)]


class TestPlan(unittest.TestCase):
    """
    Unit tests for planned deployments.
    """

    def test_mass_capacity_by_plane(self):
        """A launch capacity in mass may differ by plane (such as by inclination),
        for the launches of a deployment and those of its restorations."""

        def capacity(index, values):
            del values  # the planes of a Walker constellation have none
            return 100.0 if index == 0 else 200.0

        schedule = DeploymentSchedule(years=5, launch_mass_capacity=capacity)
        deployment = schedule.plan(_walker(2, 4), masses=50.0)
        launches = {
            plane: sum(l.plane == plane for l in deployment.launches)
            for plane in (0, 1)
        }
        self.assertEqual(launches, {0: 2, 1: 1})
        self.assertEqual(deployment.plane_values, [None, None])
        with self.assertRaises(ValueError):
            schedule.launches_of(4, 50.0)
        self.assertEqual(schedule.launches_of(4, 50.0, (1, None)), 1)
        # restorations that wait to share a launch fill it to the plane's capacity
        history = deployment.simulate(
            FailureModel(life=Exponential(2.0)),
            RestorationPolicy(batch_years=1.0),
            seed=0,
        )
        restorations = [l for l in history.launches if l.purpose == "restoration"]
        self.assertTrue(any(l.plane == 1 and l.satellites > 2 for l in restorations))
        self.assertTrue(
            all(l.satellites <= (2 if l.plane == 0 else 4) for l in restorations)
        )

    def test_production_rate(self):
        """A production line deploys planes as their satellites are built."""
        deployment = DeploymentSchedule(development_years=2, production_rate=8).plan(
            _walker(3, 4)
        )
        # 4 satellites a plane at 8 a year: ready every half year
        self.assertEqual([l.time for l in deployment.launches], [2.5, 3.0, 3.5])
        self.assertEqual(deployment.begin, 2.5)
        self.assertEqual(deployment.deployed, 3.5)
        larger = DeploymentSchedule(development_years=2, production_rate=8).plan(
            _walker(3, 6)
        )
        self.assertGreater(larger.deployed, deployment.deployed)

    def test_launch_cadence(self):
        """Launches are at least an interval apart, a plane taking several."""
        schedule = DeploymentSchedule(
            years=6,
            design_life=3,
            development_years=2,
            launch_interval=1,
            launch_capacity=3,
        )
        deployment = schedule.plan(_walker(3, 4))
        first = [l for l in deployment.launches if l.generation == 0]
        self.assertEqual([l.time for l in first], [2, 3, 4, 5, 6, 7])
        self.assertEqual(
            [(l.plane, l.first, l.satellites) for l in first[:2]],
            [(0, 0, 2), (0, 2, 2)],
        )
        # replacements 3 years later, before the end of the mission at year 8
        replacements = [l.time for l in deployment.launches if l.generation > 0]
        self.assertEqual(replacements, [5, 6, 7])
        self.assertEqual(deployment.end, 8)
        self.assertEqual(deployment.units, {None: 12 + 6})

    def test_types_and_order(self):
        """Types are developed over their own durations, planes in an order."""
        schedule = DeploymentSchedule(development_years={"new": 5, "heritage": 2})
        sso = DesignSpace(
            layout=sso_planes_space(crossing=Real(9, 12, step=3), satellites=1, count=2)
        )
        layout = sso.enumerate()[0].layout
        deployment = schedule.plan({"new": _walker(1, 2), "heritage": layout})
        self.assertEqual(
            sorted((l.time, l.kind) for l in deployment.launches),
            [(2.0, "heritage"), (2.0, "heritage"), (5.0, "new")],
        )
        ordered = DeploymentSchedule(
            production_rate=1, order=lambda index, plane: -plane.crossing
        ).plan(layout)
        crossings = [plane.crossing for plane in layout.planes]
        latest = int(np.argmax(crossings))
        self.assertEqual(ordered.launches[0].plane, latest)
        with self.assertRaises(ValueError):
            DeploymentSchedule(development_years={"new": 5}).plan(
                _walker(), kind="heritage"
            )

    def test_design_life_by_type(self):
        """Each type of satellite has its own design life (or none)."""
        schedule = DeploymentSchedule(
            years=6, design_life={"short": 2, "long": None}, development_years=1
        )
        deployment = schedule.plan({"short": _walker(1, 2), "long": _walker(1, 2)})
        replacements = [
            (l.kind, l.time) for l in deployment.launches if l.purpose == "replacement"
        ]
        self.assertEqual(replacements, [("short", 3), ("short", 5)])
        self.assertEqual(deployment.units, {"short": 6, "long": 2})
        self.assertEqual(schedule.design_life_of("short"), 2)
        self.assertIsNone(schedule.design_life_of("long"))
        with self.assertRaises(ValueError):
            DeploymentSchedule(design_life={"short": 2}).plan(_walker(), kind="long")

    def test_sustained(self):
        """A sustained constellation's replacements are produced after the end
        of the mission (from within its development years after the end), but
        not launched in it."""
        terms = {
            "years": 5,
            "design_life": 2,
            "development_years": 3,
            "production_rate": 0.5,
        }
        ending = DeploymentSchedule(**terms).plan(_walker(1, 1))
        sustained = DeploymentSchedule(**terms, sustained=True).plan(_walker(1, 1))
        self.assertEqual(sustained.end, 10)
        self.assertEqual(
            [(b.begin, b.end) for b in ending.production], [(3, 5), (5, 7), (7, 9)]
        )
        self.assertEqual(
            [(b.begin, b.end) for b in sustained.production],
            [(3, 5), (5, 7), (7, 9), (9, 11), (11, 13)],
        )
        self.assertEqual([l.time for l in sustained.launches], [5, 7, 9])
        self.assertEqual(sustained.launches, ending.launches)
        # those in production at the end, not those begun after it
        self.assertEqual(sustained.units, {None: 4})
        self.assertEqual(ending.units, {None: 3})

    def test_mass_capacity(self):
        """A launch carries as many satellites as its mass capacity allows."""
        schedule = DeploymentSchedule(launch_mass_capacity=120)
        deployment = schedule.plan(_walker(2, 4), masses=50)
        self.assertEqual(len(deployment.launches), 4)
        with self.assertRaises(ValueError):
            schedule.plan(_walker(2, 4))
        with self.assertRaises(ValueError):
            schedule.plan(_walker(2, 4), masses=200)

    def test_validation(self):
        """Schedules are validated."""
        for invalid in (
            {"years": 0},
            {"design_life": 0},
            {"design_life": {"a": 2, "b": 0}},
            {"production_rate": 0},
            {"launch_interval": -1},
            {"launch_capacity": 0},
            {"development_years": -1},
        ):
            with self.assertRaises(ValueError):
                DeploymentSchedule(**invalid)


class TestPhases(unittest.TestCase):
    """
    Unit tests for phases of a mission.
    """

    def setUp(self):
        self.schedule = DeploymentSchedule(
            years=4, design_life=2, development_years=1, production_rate=2
        )
        self.deployment = self.schedule.plan(_walker(2, 2))
        self.names = _names(2, 2)

    def test_deployment(self):
        """Phases begin as planes are launched."""
        phases = self.deployment.phases(self.names)
        self.assertEqual([(p.start, p.end) for p in phases], [(2, 3), (3, 6)])
        self.assertEqual(phases[0].satellites, tuple(self.names[:2]))
        self.assertEqual(phases[1].satellites, tuple(self.names))
        with self.assertRaises(ValueError):
            self.deployment.phases(self.names[:3])

    def test_evaluate_and_summarize(self):
        """Phases are evaluated with their satellites' observations, and summarized."""
        observations = pd.DataFrame({"satellite": self.names * 3, "target": range(12)})
        calls = []

        def metric(obs):
            calls.append(len(obs))
            return {"targets": 25.0 * obs["satellite"].nunique()}

        # the first satellite fails at 2.5 and is replaced at 4
        names = tuple(self.names)
        phases = [
            Phase(2, 2.5, names[:2]),
            Phase(2.5, 3, names[1:2]),
            Phase(3, 4, names[1:]),
            Phase(4, 6, names),
        ]
        capability = evaluate_phases(phases, observations, metric)
        self.assertEqual(list(capability.targets), [50, 25, 75, 100])
        self.assertEqual(list(capability.satellites), [2, 1, 3, 4])
        summary = summarize_phases(
            capability,
            requirements=[
                Requirement("targets", ">=", 75),
                Requirement("targets", ">=", 200),
            ],
        )
        expected = (50 * 0.5 + 25 * 0.5 + 75 * 1 + 100 * 2) / 4
        self.assertAlmostEqual(summary["targets mission"], expected)
        self.assertEqual(summary["targets time"], math.inf)  # the last requirement
        first = summarize_phases(
            capability, requirements=[Requirement("targets", ">=", 75)]
        )
        self.assertEqual(first["targets time"], 3)
        # phases of the same satellites are evaluated once
        repeated = evaluate_phases(phases + phases, observations, metric)
        self.assertEqual(len(repeated), 8)
        self.assertEqual(len(calls), 4 + 4)
