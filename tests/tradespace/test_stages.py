"""
Unit tests for staged deployments.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import pandas as pd

from tatc.tradespace import (
    Capability,
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    Exponential,
    FailureModel,
    PhaseCache,
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


NAMES = [f"s{i}" for i in range(6)]
# an observation of each satellite, so that the metric counts those operating
OBSERVATIONS = pd.DataFrame({"satellite": NAMES})


def _count(observations):
    """The number of satellites that observe."""
    return {"observing": float(observations.satellite.nunique())}


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

    def test_capability(self):
        """Rules read the capability of the satellites operating at a decision."""
        seen = []

        def rule(decision):
            seen.append((decision.time, decision.operating, decision.capability))
            return decision.capability["observing"] < 4

        deployment = SCHEDULE.plan(_layout(), stages=[Stage(2, [0.0, 0.5], rule)])
        deployment.capability = Capability(NAMES, OBSERVATIONS, _count)
        history = deployment.simulate()
        # the first plane alone at the first launch (the second is launched later)
        self.assertEqual(seen, [(0.0, 2, {"observing": 2.0})])
        self.assertEqual(history.stages, [deployment.begin])
        # the planned history (of the availability) reads it too
        self.assertEqual(history.metrics()["stages deployed"], 1.0)
        self.assertEqual(deployment.planned().stages, [deployment.begin])
        later = SCHEDULE.plan(_layout(), stages=[Stage(2, 0.5, rule)])
        later.capability = deployment.capability
        self.assertEqual(later.simulate().stages, [None])
        # without a capability, a rule cannot read it
        unmeasured = SCHEDULE.plan(_layout(), stages=[Stage(2, 0.0, rule)])
        with self.assertRaisesRegex(ValueError, "capability"):
            unmeasured.simulate()

    def test_capability_of_failures(self):
        """The capability at a decision follows failures, as rules on them do."""
        failures = FailureModel(life=Exponential(0.1))
        policy = RestorationPolicy(restore=False)
        by_count = SCHEDULE.plan(
            _layout(),
            stages=[Stage(2, [0.5, 1.0, 1.5, 2.0], lambda d: d.operating < 4)],
        )
        measured = SCHEDULE.plan(
            _layout(),
            stages=[
                Stage(
                    2,
                    [0.5, 1.0, 1.5, 2.0],
                    lambda d: d.capability["observing"] < 4,
                )
            ],
        )
        cache = PhaseCache(_count)
        measured.capability = Capability(
            NAMES, OBSERVATIONS, cache, context="replicate"
        )
        expected = by_count.simulate_many(failures, policy, seed=1, histories=20)
        histories = measured.simulate_many(failures, policy, seed=1, histories=20)
        self.assertEqual([h.stages for h in histories], [h.stages for h in expected])
        self.assertTrue(any(h.stages[0] is not None for h in histories))
        # each set of satellites read at a decision is evaluated once
        self.assertGreater(cache.hits, 0)

    def test_capability_reads(self):
        """A capability evaluates sets of satellites by name, or by position, only
        when read."""
        calls = []

        def metric(observations):
            calls.append(len(observations))
            return _count(observations)

        capability = Capability(NAMES, OBSERVATIONS, metric)
        self.assertEqual(capability(["s0", "s3"]), {"observing": 2.0})
        planes = SCHEDULE.plan(_layout()).planes
        self.assertEqual(capability.at(planes, [(1, 0), (2, 1)]), {"observing": 2.0})
        self.assertEqual(capability(set()), {"observing": 0.0})
        self.assertEqual(calls, [2, 2, 0])
        # a rule that does not read it does not evaluate it
        deployment = SCHEDULE.plan(_layout(), stages=[Stage(2, 0.0, lambda d: True)])
        deployment.capability = capability
        deployment.simulate()
        self.assertEqual(len(calls), 3)
        with self.assertRaises(ValueError):
            capability(["s9"])
        with self.assertRaises(ValueError):
            Capability(NAMES[:4], OBSERVATIONS, metric).at(planes, [(0, 0)])
        # a cache keys satellites given by name by their names
        cache = PhaseCache(_count)
        cached = Capability(NAMES, OBSERVATIONS, cache)
        cached(["s0", "s1"])
        cached(["s1", "s0"])
        self.assertEqual((cache.hits, cache.misses), (1, 1))

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
