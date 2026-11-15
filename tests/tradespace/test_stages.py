"""
Unit tests for staged deployments.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import dataclasses
import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    Exponential,
    FailureModel,
    RestorationPolicy,
    SatelliteCost,
    Stage,
    walker_space,
)
from tatc.tradespace.expectations import exact_expectations
from tatc.tradespace.mission import Capability, PhaseCache

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
        self.assertEqual(history.decided, [decided])
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
        # the first generation is deployed with the stage
        self.assertEqual(history.deployed, decided + 2 / 4)
        self.assertEqual(deployment.deployed, decided + 2 / 4)
        self.assertLess(SCHEDULE.plan(_layout()).deployed, deployment.deployed)

    def test_sustained(self):
        """A deployed stage of a sustained constellation is replaced after the
        end of the mission as the other planes are: produced, not launched."""
        schedule = dataclasses.replace(SCHEDULE, design_life=3.75, sustained=True)
        history = schedule.plan(_layout(), stages=[Stage(2, 1.0)]).simulate()
        self.assertEqual(history.end, 6.5)
        # launched at 3, its replacement needed at 6.75, after the end
        self.assertEqual(_launches(history, 2), [(3.0, "deployment")])
        after = [(b.begin, b.end) for b in history.production if b.end > history.end]
        self.assertEqual(after, [(6.25, 6.75)])
        # eleven units of 2: five planes launched, and half of the sixth
        self.assertAlmostEqual(COSTS.estimate(history).production, 11 * 2)

    def test_not_deployed(self):
        """A stage that its rule does not deploy is neither produced nor launched."""
        never = SCHEDULE.plan(
            _layout(), stages=[Stage([1, 2], [0.5, 1.0], lambda d: False)]
        )
        history = never.simulate()
        self.assertEqual(history.decided, [None])
        self.assertEqual(_launches(history, 1) + _launches(history, 2), [])
        self.assertNotIn((1, 0), history.operating)
        self.assertEqual(history.metrics()["stages deployed"], 0.0)
        plain = SCHEDULE.plan(_layout())
        self.assertLess(COSTS.estimate(never).total, COSTS.estimate(plain).total)
        self.assertNotIn("stages deployed", plain.simulate().metrics())
        self.assertEqual(never.deployed, history.deployed)

    def test_deployed_after_the_end(self):
        """A stage deployed too late to launch within the mission leaves the
        deployment unfinished: the time it is deployed is infinite."""
        late = SCHEDULE.plan(_layout(), stages=[Stage(2, 4.75)])
        history = late.simulate()
        self.assertEqual(history.decided, [late.begin + 4.75])
        self.assertEqual(_launches(history, 2), [])
        self.assertEqual(history.deployed, float("inf"))
        self.assertEqual(late.deployed, float("inf"))
        # or deployed within it, or not at all
        early = SCHEDULE.plan(_layout(), stages=[Stage(2, 4.0)])
        self.assertEqual(early.deployed, early.begin + 4.0 + 2 / 4)
        never = SCHEDULE.plan(_layout(), stages=[Stage(2, 4.75, lambda d: False)])
        self.assertLess(never.deployed, never.end)

    def test_cost_of_the_option(self):
        """A type of satellite whose planes are all staged costs its
        development, whether or not they are deployed: the cost of the option."""
        plane = (
            DesignSpace(layout=walker_space(inclination=60, planes=1, per_plane=2))
            .enumerate()[0]
            .layout
        )
        model = CostModel(
            {
                "a": SatelliteCost(development=10, first_unit=1),
                "b": SatelliteCost(development=50, first_unit=1),
            }
        )
        never = DeploymentSchedule(years=5, development_years=1).plan(
            {"a": plane, "b": plane}, stages=[Stage(1, 0.0, lambda decision: False)]
        )
        estimate = model.estimate(never)
        self.assertEqual(never.simulate().decided, [None])
        self.assertEqual(estimate.units, {"a": 2})
        self.assertEqual(estimate.development, 10 + 50)

    def test_decisions(self):
        """Rules decide at their times, until they deploy, from what is known."""
        seen = []

        def rule(decision):
            seen.append(decision)
            return decision.since_launch >= 1.0

        stages = [Stage(1, 0.25), Stage(2, [0.5, 1.0, 1.5], rule)]
        history = SCHEDULE.plan(_layout(), stages=stages).simulate()
        self.assertEqual([d.since_launch for d in seen], [0.5, 1.0])
        self.assertEqual({d.stage for d in seen}, {1})
        self.assertEqual(seen[0].stages_deployed, (0,))
        # the satellites of the first plane operate (those of the second are
        # still being produced)
        self.assertEqual(seen[0].satellites_operating, 2)
        begin = min(launch.time for launch in history.launches)
        # times are from the start of development, as every time of a mission
        self.assertEqual([d.time for d in seen], [begin + 0.5, begin + 1.0])
        self.assertEqual(history.decided, [begin + 0.25, begin + 1.0])

    def test_contingency(self):
        """A rule on failures deploys a stage in the histories that need it."""
        failures = FailureModel(life=Exponential(0.1))
        policy = RestorationPolicy(restore=False)
        stage = Stage(2, [0.5, 1.0, 1.5, 2.0], lambda d: d.satellites_operating < 4)
        deployment = SCHEDULE.plan(_layout(), policy=policy, stages=[stage])
        self.assertEqual(deployment.simulate().decided, [None])
        histories = deployment.simulate_many(failures, seed=1, histories=20)
        deployed = [h.decided[0] is not None for h in histories]
        self.assertTrue(any(deployed) and not all(deployed))
        for history in histories:
            failed = [t for t, *_ in history.failures]
            if history.decided[0] is not None:
                self.assertLess(min(failed), history.decided[0])
            # the satellites of a stage count in the availability of the
            # histories that deploy it (and only of those)
            self.assertLessEqual(history.metrics()["availability"], 1.0)
        # the same history draws the same failures again
        again = deployment.simulate_many(failures, seed=1, histories=20)
        self.assertEqual([h.decided for h in again], [h.decided for h in histories])

    def test_capability(self):
        """Rules read the capability of the satellites operating at a decision."""
        seen = []

        def rule(decision):
            seen.append(
                (
                    decision.since_launch,
                    decision.satellites_operating,
                    decision.capability,
                )
            )
            return decision.capability["observing"] < 4

        capability = Capability(NAMES, OBSERVATIONS, _count)
        deployment = SCHEDULE.plan(
            _layout(), stages=[Stage(2, [0.0, 0.5], rule)], capability=capability
        )
        self.assertIs(deployment.capability, capability)
        history = deployment.simulate()
        # the first plane alone at the first launch (the second is launched later)
        self.assertEqual(seen, [(0.0, 2, {"observing": 2.0})])
        self.assertEqual(history.decided, [deployment.begin])
        # the planned history reads it too
        self.assertEqual(history.metrics()["stages deployed"], 1.0)
        self.assertEqual(deployment.planned().decided, [deployment.begin])
        later = SCHEDULE.plan(
            _layout(), stages=[Stage(2, 0.5, rule)], capability=capability
        )
        self.assertEqual(later.simulate().decided, [None])
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
            policy=policy,
            stages=[
                Stage(2, [0.5, 1.0, 1.5, 2.0], lambda d: d.satellites_operating < 4)
            ],
        )
        cache = PhaseCache(_count)
        measured = SCHEDULE.plan(
            _layout(),
            policy=policy,
            stages=[
                Stage(
                    2,
                    [0.5, 1.0, 1.5, 2.0],
                    lambda d: d.capability["observing"] < 4,
                )
            ],
            capability=Capability(NAMES, OBSERVATIONS, cache, context="replicate"),
        )
        expected = by_count.simulate_many(failures, seed=1, histories=20)
        histories = measured.simulate_many(failures, seed=1, histories=20)
        self.assertEqual([h.decided for h in histories], [h.decided for h in expected])
        self.assertTrue(any(h.decided[0] is not None for h in histories))
        # each set of satellites read at a decision is evaluated once
        self.assertGreater(cache.hits, 0)

    def test_capability_reads(self):
        """A capability evaluates sets of satellites by name, or by position, only
        when read."""
        # pylint: disable=protected-access
        calls = []

        def metric(observations):
            calls.append(len(observations))
            return _count(observations)

        capability = Capability(NAMES, OBSERVATIONS, metric)
        self.assertEqual(capability(["s0", "s3"]), {"observing": 2.0})
        planes = SCHEDULE.plan(_layout()).planes
        self.assertEqual(capability._at(planes, [(1, 0), (2, 1)]), {"observing": 2.0})
        self.assertEqual(capability(set()), {"observing": 0.0})
        self.assertEqual(calls, [2, 2, 0])
        # a rule that does not read it does not evaluate it
        deployment = SCHEDULE.plan(
            _layout(), stages=[Stage(2, 0.0, lambda d: True)], capability=capability
        )
        deployment.simulate()
        self.assertEqual(len(calls), 3)
        with self.assertRaises(ValueError):
            capability(["s9"])
        with self.assertRaises(ValueError):
            Capability(NAMES[:4], OBSERVATIONS, metric)._at(planes, [(0, 0)])
        # a cache keys satellites given by name by their names
        cache = PhaseCache(_count)
        cached = Capability(NAMES, OBSERVATIONS, cache)
        cached(["s0", "s1"])
        cached(["s1", "s0"])
        self.assertEqual((cache.hits, cache.misses), (1, 1))
        # which the cache gives by their keys (and the context)
        keys = {"s0": "s0", "s1": "s1"}
        self.assertEqual(cache._metrics_of(keys, OBSERVATIONS), {"observing": 2.0})
        self.assertEqual((cache.hits, cache.misses), (2, 1))
        cache._metrics_of(keys, OBSERVATIONS, context="other")
        self.assertEqual((cache.hits, cache.misses), (2, 2))

    def test_lead_and_launches(self):
        """A lead time delays production, and launch capacity splits launches."""
        schedule = DeploymentSchedule(
            years=5, development_years=1, production_rate=4, launch_capacity=1
        )
        deployment = schedule.plan(
            _layout(),
            policy=RestorationPolicy(orbit_spares=1),
            stages=[Stage(2, 1.0, lead_years=0.5)],
        )
        history = deployment.simulate()
        # the line is free by the end of the lead time
        decided = history.decided[0]
        launched = [l for l in history.launches if l.plane == 2]
        self.assertEqual([l.satellites for l in launched], [1, 1, 1])
        self.assertEqual(launched[0].time, decided + 0.5 + 3 / 4)
        # the spare of the staged plane waits in orbit
        self.assertTrue(
            any(start >= launched[0].time for start, _, _ in history.parked)
        )

    def test_exact_expectations(self):
        """Exact expectations need a deployment that is not staged."""
        deployment = SCHEDULE.plan(
            _layout(),
            policy=RestorationPolicy(restore=False),
            stages=[Stage(2, 1.0)],
        )
        with self.assertRaisesRegex(ValueError, "staged"):
            exact_expectations(
                deployment,
                FailureModel(life=Exponential(0.2)),
                Capability([], None, len),
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

    def test_numpy_values(self):
        """Stages take integers and arrays of numpy as planes and times."""
        stage = Stage(np.int64(1), np.float64(0.5))
        self.assertEqual((stage.planes, stage.times), ((1,), (0.5,)))
        stage = Stage(np.array([2, 1]), np.array([1.0, 0.0]))
        self.assertEqual((stage.planes, stage.times), ((2, 1), (0.0, 1.0)))
        self.assertEqual(Stage([1], range(3)).times, (0.0, 1.0, 2.0))
        deployment = SCHEDULE.plan(_layout(), stages=[Stage(np.int64(2), 1.0)])
        self.assertEqual(deployment.simulate().decided, [deployment.begin + 1.0])


if __name__ == "__main__":
    unittest.main()
