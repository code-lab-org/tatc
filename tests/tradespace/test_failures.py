"""
Unit tests for failures of satellites and launches, their restoration, and
simulations of missions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest
from dataclasses import dataclass

import numpy as np

from tatc.tradespace import (
    Bathtub,
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    Exponential,
    FailureModel,
    Integer,
    LaunchCost,
    Lifetime,
    OperationsCost,
    Piecewise,
    RestorationPolicy,
    SatelliteCost,
    Study,
    Uncertain,
    Weibull,
    aggregate,
    configure,
    pareto_front,
    policy_space,
    walker_space,
)


def _walker(planes=1, per_plane=2):
    (design,) = DesignSpace(
        layout=walker_space(
            inclination=55.0, planes=planes, per_plane=per_plane, relative_spacing=0
        )
    ).enumerate()
    return design.layout


def _names(count):
    return [f"satellite {k + 1}" for k in range(count)]


@dataclass(frozen=True)
class _Fixed(Lifetime):
    """A life of a fixed duration (the first satellite only, if `first`)."""

    years: float

    def cumulative_hazard(self, age):
        return np.where(np.asarray(age) >= self.years, np.inf, 0.0)

    def hazard(self, age):
        return np.zeros_like(np.asarray(age, dtype=float))

    def sample(self, rng):
        return self.years


def _empirical(life, ages, n=40_000, seed=1):
    rng = np.random.default_rng(seed)
    lives = np.array([life.sample(rng) for _ in range(n)])
    return np.array([(lives > age).mean() for age in ages])


class TestLifetimes(unittest.TestCase):
    """
    Unit tests for distributions of life.
    """

    ages = np.array([0.5, 2.0, 5.0, 9.0])

    def test_exponential(self):
        """A constant failure rate."""
        life = Exponential(0.2)
        np.testing.assert_allclose(life.survival(self.ages), np.exp(-0.2 * self.ages))
        np.testing.assert_allclose(life.hazard(self.ages), 0.2)
        np.testing.assert_allclose(
            _empirical(life, self.ages), life.survival(self.ages), atol=0.01
        )
        self.assertEqual(Exponential(0).sample(np.random.default_rng(0)), math.inf)
        with self.assertRaises(ValueError):
            Exponential(-1)

    def test_weibull(self):
        """Decreasing and increasing failure rates."""
        infant, wear = Weibull(0.5, 20), Weibull(4, 6)
        for life in (infant, wear):
            np.testing.assert_allclose(
                _empirical(life, self.ages), life.survival(self.ages), atol=0.01
            )
        self.assertGreater(infant.hazard(0.5), infant.hazard(5))
        self.assertLess(wear.hazard(0.5), wear.hazard(5))
        with self.assertRaises(ValueError):
            Weibull(0, 1)

    def test_bathtub(self):
        """Competing causes: hazards add, the earliest fails."""
        parts = (Weibull(0.5, 30), Exponential(0.05), Weibull(5, 7))
        life = Bathtub(*parts)
        np.testing.assert_allclose(
            life.hazard(self.ages), sum(p.hazard(self.ages) for p in parts)
        )
        np.testing.assert_allclose(
            life.survival(self.ages),
            np.prod([p.survival(self.ages) for p in parts], axis=0),
        )
        np.testing.assert_allclose(
            _empirical(life, self.ages), life.survival(self.ages), atol=0.01
        )
        # high early, low in the middle, high late
        hazard = life.hazard(np.array([0.05, 3.0, 9.0]))
        self.assertTrue(hazard[0] > hazard[1] < hazard[2])
        self.assertEqual(Bathtub().sample(np.random.default_rng(0)), math.inf)

    def test_piecewise(self):
        """Failure rates by intervals of age."""
        life = Piecewise([0, 1, 4], [0.3, 0.05, 0.2])
        np.testing.assert_allclose(life.hazard(np.array([0.5, 2, 6])), [0.3, 0.05, 0.2])
        expected = np.exp(
            -np.array([0.15, 0.3 + 0.05, 0.3 + 0.15 + 0.2, 0.3 + 0.15 + 1.0])
        )
        np.testing.assert_allclose(life.survival(self.ages), expected)
        np.testing.assert_allclose(
            _empirical(life, self.ages), life.survival(self.ages), atol=0.01
        )
        never = Piecewise([0, 1], [0.5, 0.0])
        self.assertGreater(
            np.mean(
                [never.sample(np.random.default_rng(s)) == math.inf for s in range(400)]
            ),
            0.5,
        )
        for breaks, rates in (
            ([1, 2], [0.1, 0.1]),
            ([0, 2], [0.1]),
            ([0, 2, 1], [1, 1, 1]),
        ):
            with self.assertRaises(ValueError):
                Piecewise(breaks, rates)

    def test_models(self):
        """Lives of each type of satellite, and failures of launches."""
        model = FailureModel(life={"a": Exponential(0.1), "b": None}, launch=0.05)
        self.assertEqual(model.life_of("a"), Exponential(0.1))
        self.assertIsNone(model.life_of("b"))
        with self.assertRaises(ValueError):
            model.life_of("c")
        with self.assertRaises(ValueError):
            FailureModel(launch=1.5)


class TestSimulation(unittest.TestCase):
    """
    Unit tests for simulated missions.
    """

    def setUp(self):
        # one plane of 2 satellites, launched at year 1, for a mission of 4 years
        self.schedule = DeploymentSchedule(
            years=4, development_years=1, production_rate=4, launch_interval=0.5
        )
        self.layout = _walker(1, 2)
        self.names = _names(2)
        # every satellite fails 1.5 years after it begins to operate
        self.failures = FailureModel(life=_Fixed(1.5))

    def _simulate(self, policy, spares=0, failures=None, seed=0):
        deployment = self.schedule.plan(self.layout, spares=spares)
        return deployment.simulate(failures or self.failures, policy, seed)

    def test_without_failures(self):
        """Without failures, the history is the deployment's."""
        deployment = self.schedule.plan(self.layout)
        history = deployment.simulate()
        self.assertEqual(history.failures, [])
        self.assertEqual(
            history.operating, {(0, 0): [(1.5, 5.5)], (0, 1): [(1.5, 5.5)]}
        )
        model = CostModel(
            SatelliteCost(development=10, first_unit=5),
            launch=LaunchCost(per_launch=3),
            operations=OperationsCost(per_satellite_year=1),
        )
        self.assertAlmostEqual(
            model.estimate(deployment).total, model.estimate(history).total
        )
        # a failure model without a seed draws no failures
        self.assertEqual(deployment.simulate(self.failures, seed=None).failures, [])
        self.assertEqual(history.metrics()["availability"], 1.0)

    def test_degrade(self):
        """Without restoration, failed satellites stay failed."""
        history = self._simulate(RestorationPolicy(restore=False))
        # both fail at 3.0 (1.5 years after their launch at 1.5)
        self.assertEqual(history.operating[0, 0], [(1.5, 3.0)])
        self.assertEqual(len(history.failures), 2)
        self.assertEqual(history.restorations, [])
        self.assertAlmostEqual(history.metrics()["availability"], 1.5 / 4)
        phases = history.phases(self.names)
        self.assertEqual(
            [(p.start, p.end, len(p.satellites)) for p in phases],
            [(1.5, 3.0, 2), (3.0, 5.5, 0)],
        )

    def test_restore_by_production(self):
        """A failed satellite is produced, then launched at the next slot."""
        policy = RestorationPolicy(detection_years=0.1, procurement_years=0.2)
        history = self._simulate(policy)
        # failed at 3.0, ordered at 3.1, produced from 3.3 to 3.55 (the
        # second after it, to 3.8), each launched when ready (slots 0.5 apart)
        restored = sorted(
            (p, failed, round(when, 6), source)
            for p, failed, when, source in [
                (pos, f, r, s) for _, pos, f, r, s in history.restorations[:2]
            ]
        )
        self.assertEqual([r[2] for r in restored], [3.55, 4.05])
        self.assertTrue(all(r[3] == "production" for r in restored))
        purposes = [b.purpose for b in history.production]
        self.assertIn("restoration", purposes)
        self.assertEqual(
            sorted(l.purpose for l in history.launches).count("restoration"),
            sum(1 for l in history.launches if l.purpose == "restoration"),
        )

    def test_restore_from_orbit(self):
        """A spare in orbit takes a failed satellite's place after drifting."""
        policy = RestorationPolicy(detection_years=0.1, drift_years=0.2)
        history = self._simulate(policy, spares=1)
        first = sorted(history.restorations, key=lambda r: r[3])[0]
        self.assertEqual(first[4], "orbit")
        self.assertAlmostEqual(first[3] - first[2], 0.3)
        # the other satellite, failing at the same time, waits for production
        sources = [r[4] for r in history.restorations]
        self.assertIn("production", sources)
        self.assertGreater(history.metrics()["spares used"], 0)

    def test_restore_from_ground(self):
        """A spare on the ground is launched at the next slot."""
        policy = RestorationPolicy(detection_years=0.1, ground_spares=2)
        history = self._simulate(policy)
        sources = [r[4] for r in history.restorations[:2]]
        self.assertEqual(sources, ["ground", "ground"])
        delays = sorted(round(r[3] - r[2], 6) for r in history.restorations[:2])
        # both known at 3.1, and launched together then
        self.assertEqual(delays, [0.1, 0.1])
        restorations = [l for l in history.launches if l.purpose == "restoration"]
        self.assertEqual(restorations[0].satellites, 2)
        spares = [b for b in history.production if b.purpose == "spare"]
        self.assertEqual(sum(b.satellites for b in spares), 2)

    def test_batching(self):
        """Restorations of a plane waiting within a window share a launch."""
        alone = self._simulate(RestorationPolicy(ground_spares=2))
        shared = self._simulate(RestorationPolicy(ground_spares=2, batch_years=0.2))
        count = lambda h: sum(l.purpose == "restoration" for l in h.launches)
        self.assertGreater(count(alone), count(shared))

    def test_launch_failures(self):
        """A failed launch loses its satellites."""
        always = FailureModel(launch=1.0)
        degraded = self._simulate(RestorationPolicy(restore=False), failures=always)
        self.assertEqual(degraded.operating, {})
        self.assertTrue(all(launch.failed for launch in degraded.launches))
        self.assertEqual(degraded.metrics()["availability"], 0)
        restored = self._simulate(RestorationPolicy(), failures=always)
        # restorations are attempted (and fail) until the end of the mission
        self.assertGreater(len(restored.launches), len(degraded.launches))
        self.assertEqual(restored.operating, {})

    def test_planned_replacement(self):
        """A planned replacement restores a degraded satellite."""
        schedule = DeploymentSchedule(years=4, design_life=2, development_years=1)
        history = schedule.plan(self.layout).simulate(
            self.failures, RestorationPolicy(restore=False), seed=0
        )
        # launched at 1, failed at 2.5, replaced at 3, failed at 4.5
        self.assertEqual(history.operating[0, 0], [(1, 2.5), (3, 4.5)])

    def test_common_random_numbers(self):
        """The satellite at the same position fails at the same age."""
        failures = FailureModel(life=Exponential(0.3), launch=0.1)
        first = self._simulate(
            RestorationPolicy(restore=False), failures=failures, seed=4
        )
        other = self._simulate(RestorationPolicy(), failures=failures, seed=4)
        self.assertEqual(first.failures[0], other.failures[0])
        different = self._simulate(RestorationPolicy(), failures=failures, seed=5)
        self.assertNotEqual(first.failures, different.failures)

    def test_costs_of_policies(self):
        """Restorations and spares cost more than degrading."""
        model = CostModel(
            SatelliteCost(development=10, first_unit=5),
            launch=LaunchCost(per_launch=3),
            operations=OperationsCost(per_satellite_year=1),
        )
        failures = FailureModel(life=Exponential(0.3))
        totals = {}
        for name, policy, spares in (
            ("degrade", RestorationPolicy(restore=False), 0),
            ("produce", RestorationPolicy(), 0),
            ("orbit", RestorationPolicy(replenish=True), 1),
        ):
            costs, availability = [], []
            for seed in range(100):
                history = self._simulate(
                    policy, spares=spares, failures=failures, seed=seed
                )
                costs.append(model.estimate(history).total)
                availability.append(history.metrics()["availability"])
            totals[name] = (np.mean(costs), np.mean(availability))
        self.assertLess(totals["degrade"][0], totals["produce"][0])
        self.assertLess(totals["degrade"][1], totals["produce"][1])
        self.assertLess(totals["produce"][0], totals["orbit"][0])

    def test_policy_validation(self):
        """Policies are validated."""
        for invalid in (
            {"detection_years": -1},
            {"drift_years": -1},
            {"procurement_years": -1},
            {"batch_years": -1},
            {"ground_spares": -1},
        ):
            with self.assertRaises(ValueError):
                RestorationPolicy(**invalid)
        self.assertEqual(
            RestorationPolicy(ground_spares={"a": 2}).ground_spares_of("b"), 0
        )


class TestPolicySpace(unittest.TestCase):
    """
    Unit tests for policies as decisions of designs.
    """

    def test_fragment(self):
        """Policies are enumerated, without spares that are never used."""
        space = policy_space(orbit_spares=Integer(0, 1), ground_spares=Integer(0, 1))
        flats = [d.flat() for d in space.enumerate()]
        self.assertEqual(len(flats), 5)
        self.assertNotIn(
            {"restore": False, "orbit_spares": 1, "ground_spares": 0}, flats
        )
        loose = policy_space(
            orbit_spares=Integer(0, 1),
            ground_spares=Integer(0, 1),
            require_restoration=False,
        )
        self.assertEqual(len(loose.enumerate()), 8)
        with self.assertRaises(ValueError):
            policy_space(drift=0.1)

    def test_configure(self):
        """Values of a design configure a policy or a schedule."""
        space = DesignSpace(
            policy=policy_space(
                restore=True, ground_spares=Integer(1, 2), drift_years=0.5
            ),
            launch_interval=Integer(0, 1),
        )
        design = space.enumerate()[-1]
        base = RestorationPolicy(procurement_years=0.3)
        policy = configure(base, design.policy)
        self.assertEqual(
            (policy.restore, policy.ground_spares, policy.drift_years), (True, 2, 0.5)
        )
        self.assertEqual(policy.procurement_years, 0.3)
        schedule = configure(DeploymentSchedule(years=5), design.values, years=6)
        self.assertEqual((schedule.launch_interval, schedule.years), (1, 6))
        with self.assertRaises(ValueError):
            configure(base, design.policy, unknown=1)
        with self.assertRaises(TypeError):
            configure(RestorationPolicy, design.policy)

    def test_trade_of_policies(self):
        """Policies trade cost against availability on a front."""
        schedule = DeploymentSchedule(
            years=5, development_years=1, production_rate=4, launch_capacity=3
        )
        failures = FailureModel(life=Exponential(0.3), launch=0.05)
        model = CostModel(
            SatelliteCost(development=10, first_unit=5),
            launch=LaunchCost(per_launch=3),
            operations=OperationsCost(per_satellite_year=1),
        )
        space = DesignSpace(
            layout=walker_space(
                inclination=55.0, planes=Integer(1, 2), per_plane=2, relative_spacing=0
            ),
            policy=policy_space(
                orbit_spares=Integer(0, 1), ground_spares=Integer(0, 1)
            ),
            failure_seed=Uncertain.seed(),
        )

        def evaluate(d):
            deployment = schedule.plan(d.layout, spares=d.policy.orbit_spares)
            history = deployment.simulate(
                failures, configure(RestorationPolicy(), d.policy), d.failure_seed
            )
            return {
                "cost": model.estimate(history).total,
                "availability": history.metrics()["availability"],
            }

        designs = space.enumerate()
        self.assertEqual(len(designs), 2 * 5)
        results = Study(space, lambda d: d, evaluate, replicates=30).run(designs)
        summary = aggregate(results, ["cost", "availability"])
        front = summary[pareto_front(summary, {"cost": "min", "availability": "max"})]
        policies = {
            (
                row["policy.restore"],
                row["policy.orbit_spares"],
                row["policy.ground_spares"],
            )
            for _, row in results.drop_duplicates("design")
            .set_index("design")
            .loc[front.design]
            .iterrows()
        }
        # degrading is the cheapest, and a policy that restores is the most available
        self.assertIn((False, 0, 0), policies)
        self.assertGreater(len(policies), 1)
