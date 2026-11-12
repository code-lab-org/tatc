"""
Unit tests for staged deployments.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from tatc.tradespace import (
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    Exponential,
    FailureModel,
    RestorationPolicy,
    SatelliteCost,
    Stage,
    expected_missions,
    walker_space,
)

# three planes of two satellites
LAYOUT = DesignSpace(layout=walker_space(inclination=60, planes=3, per_plane=2))
SCHEDULE = DeploymentSchedule(
    years=5, design_life=3, development_years=1, production_rate=4, launch_interval=0.5
)
COSTS = CostModel(SatelliteCost(development=10, first_unit=2))


def _layout():
    return LAYOUT.enumerate()[0].layout


def _launches(history, plane):
    return [
        (launch.time, launch.purpose)
        for launch in history.launches
        if launch.plane == plane
    ]


class TestStages(unittest.TestCase):
    """
    Unit tests for stages of deployments.
    """

    def test_planned_delay(self):
        """A stage without a rule is deployed at its first time."""
        deployment = SCHEDULE.plan(_layout(), stages=[Stage(2, 1.0)])
        self.assertEqual({launch.plane for launch in deployment.launches}, {0, 1})
        history = deployment.simulate()
        decided = deployment.begin + 1.0
        self.assertEqual(history.stages, [decided])
        # produced from the decision (after the line's earlier work) and launched
        produced = [b for b in history.production if b.begin == decided]
        self.assertEqual(
            [(b.satellites, b.purpose) for b in produced], [(2, "deployment")]
        )
        launched = _launches(history, 2)
        self.assertEqual(launched[0], (decided + 2 / 4, "deployment"))
        # replaced at the end of its design life, within the mission
        self.assertEqual(launched[1], (decided + 2 / 4 + 3, "replacement"))
        self.assertEqual(history.operating[2, 0][0][0], decided + 2 / 4)
        self.assertEqual(history.metrics()["stages deployed"], 1.0)

    def test_not_deployed(self):
        """A stage that its rule does not deploy is neither produced nor launched."""
        never = SCHEDULE.plan(
            _layout(), stages=[Stage([1, 2], [0.5, 1.0], lambda d: False)]
        )
        history = never.simulate()
        self.assertEqual(history.stages, [None])
        self.assertEqual(_launches(history, 1) + _launches(history, 2), [])
        self.assertNotIn((1, 0), history.operating)
        self.assertEqual(history.metrics()["stages deployed"], 0.0)
        plain = SCHEDULE.plan(_layout())
        self.assertLess(COSTS.estimate(never).total, COSTS.estimate(plain).total)
        self.assertNotIn("stages deployed", plain.simulate().metrics())

    def test_decisions(self):
        """Rules decide at their times, until they deploy, from what is known."""
        seen = []

        def rule(decision):
            seen.append(decision)
            return decision.time >= 1.0

        stages = [Stage(1, 0.25), Stage(2, [0.5, 1.0, 1.5], rule)]
        history = SCHEDULE.plan(_layout(), stages=stages).simulate()
        self.assertEqual([d.time for d in seen], [0.5, 1.0])
        self.assertEqual({d.stage for d in seen}, {1})
        self.assertEqual(seen[0].deployed, (0,))
        # the satellites of the first plane operate (those of the second are
        # still being produced)
        self.assertEqual(seen[0].operating, 2)
        begin = min(launch.time for launch in history.launches)
        self.assertEqual(history.stages, [begin + 0.25, begin + 1.0])

    def test_contingency(self):
        """A rule on failures deploys a stage in the histories that need it."""
        failures = FailureModel(life=Exponential(0.1))
        policy = RestorationPolicy(restore=False)
        stage = Stage(2, [0.5, 1.0, 1.5, 2.0], lambda d: d.operating < 4)
        deployment = SCHEDULE.plan(_layout(), stages=[stage])
        self.assertEqual(deployment.simulate().stages, [None])
        histories = deployment.simulate_many(failures, policy, seed=1, histories=20)
        deployed = [h.stages[0] is not None for h in histories]
        self.assertTrue(any(deployed) and not all(deployed))
        for history in histories:
            failed = [t for t, *_ in history.failures]
            if history.stages[0] is not None:
                self.assertLess(min(failed), history.stages[0])
        # the same history draws the same failures again
        again = deployment.simulate_many(failures, policy, seed=1, histories=20)
        self.assertEqual([h.stages for h in again], [h.stages for h in histories])

    def test_lead_and_launches(self):
        """A lead time delays production, and launch capacity splits launches."""
        schedule = DeploymentSchedule(
            years=5, development_years=1, production_rate=4, launch_capacity=1
        )
        deployment = schedule.plan(
            _layout(), spares=1, stages=[Stage(2, 1.0, lead_years=0.5)]
        )
        history = deployment.simulate()
        # the line is free by the end of the lead time
        decided = history.stages[0]
        launched = [l for l in history.launches if l.plane == 2]
        self.assertEqual([l.satellites for l in launched], [1, 1, 1])
        self.assertEqual(launched[0].time, decided + 0.5 + 3 / 4)
        # the spare of the staged plane waits in orbit
        self.assertTrue(
            any(start >= launched[0].time for start, _, _ in history.parked)
        )

    def test_exact_expectations(self):
        """Exact expectations need a deployment that is not staged."""
        deployment = SCHEDULE.plan(_layout(), stages=[Stage(2, 1.0)])
        with self.assertRaises(ValueError):
            expected_missions(
                deployment, FailureModel(life=Exponential(0.2)), [], None, len
            )

    def test_invalid(self):
        """Stages are of distinct planes of the layout, after some planned ones."""
        for stages in (
            [Stage(1), Stage([1, 2])],
            [Stage(3)],
            [Stage([0, 1, 2])],
        ):
            with self.subTest(stages=stages):
                with self.assertRaises(ValueError):
                    SCHEDULE.plan(_layout(), stages=stages)
        for arguments in (([],), ([1, 1],), ([-1],), (1, -1.0), (1, [])):
            with self.subTest(arguments=arguments):
                with self.assertRaises(ValueError):
                    Stage(*arguments)
        with self.assertRaises(ValueError):
            Stage(1, lead_years=-1)


if __name__ == "__main__":
    unittest.main()
