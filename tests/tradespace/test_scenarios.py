"""
Unit tests for scenarios of studies and the robustness of designs across them.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Criteria,
    DesignSpace,
    Fixed,
    Integer,
    Requirement,
    ResultStore,
    Study,
    Uncertain,
    aggregate,
    attribute_variance,
    robustness,
)

SPACE = DesignSpace(x=Integer(1, 3), level=Fixed(1.0), noise=Uncertain.uniform(0, 1))
SCENARIOS = {"low": {}, "high": {"level": 3.0}, "pinned": {"noise": 0.5}}


def _evaluate(d):
    return {"y": d.x * d.level + d.noise, "cost": float(d.x)}


class TestStudyScenarios(unittest.TestCase):
    """
    Unit tests for studies in scenarios.
    """

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / "results.sqlite"
        self.calls = []

    def tearDown(self):
        self.directory.cleanup()

    def _study(self, scenarios=None, **options):
        def evaluate(d):
            self.calls.append((d.x, d.replicate))
            return _evaluate(d)

        return Study(
            SPACE, lambda d: d, evaluate, replicates=2, scenarios=scenarios, **options
        )

    def test_results(self):
        """Each design is evaluated in each scenario in each replicate."""
        results = self._study(SCENARIOS).run()
        self.assertEqual(len(results), 3 * 3 * 2)
        self.assertEqual(list(results.columns[:3]), ["design", "replicate", "scenario"])
        first = results[results.design == results.design.iloc[0]]
        self.assertEqual(
            list(zip(first.scenario, first.replicate)),
            [(s, r) for s in SCENARIOS for r in range(2)],
        )
        by = results.set_index(["design", "replicate", "scenario"])
        low, high = by.xs("low", level="scenario"), by.xs("high", level="scenario")
        # the same uncertain parameters in every scenario
        np.testing.assert_array_equal(low.noise, high.noise)
        np.testing.assert_allclose(high.y - low.y, 2 * low.cost)
        pinned = by.xs("pinned", level="scenario")
        self.assertTrue((pinned.noise == 0.5).all())
        # judged in each scenario
        summary = aggregate(results, Criteria({"y": "max"}, by="scenario"))
        self.assertEqual(len(summary), 3 * 3)
        with self.assertRaises(ValueError):
            aggregate(results, Criteria({"y": "max"}))

    def test_without_scenarios(self):
        """A study without scenarios has no column of them."""
        results = self._study().run()
        self.assertNotIn("scenario", results)
        self.assertEqual(len(results), 3 * 2)

    def test_parallel(self):
        """Studies in scenarios run in parallel as in sequence."""
        serial = self._study(SCENARIOS).run()
        with ThreadPoolExecutor(2) as executor:
            parallel = self._study(SCENARIOS, executor=executor).run()
        columns = [c for c in serial.columns if c != "seconds"]
        pd.testing.assert_frame_equal(serial[columns], parallel[columns])

    def test_store(self):
        """The evaluations of each scenario are stored by its values."""
        store = ResultStore(self.path, "scenarios")
        plain = self._study(store=store).run()
        self.assertEqual(len(self.calls), 6)
        # a scenario without values shares the evaluations without scenarios
        self.calls.clear()
        both = self._study({"low": {}, "high": {"level": 3.0}}, store=store).run()
        self.assertEqual(len(self.calls), 6)
        low = both[both.scenario == "low"].drop(columns="scenario")
        columns = [c for c in plain.columns if c != "seconds"]
        pd.testing.assert_frame_equal(
            low[columns].reset_index(drop=True), plain[columns]
        )
        # renamed scenarios are taken from the store, changed ones evaluated
        self.calls.clear()
        self._study({"base": {}, "raised": {"level": 3.0}}, store=store).run()
        self.assertEqual(self.calls, [])
        self._study({"raised": {"level": 4.0}}, store=store).run()
        self.assertEqual(len(self.calls), 6)
        # the evaluations of scenarios are stored without their names
        self.calls.clear()
        again = self._study(store=store).run()
        self.assertEqual(self.calls, [])
        self.assertNotIn("scenario", again)

    def test_attribution(self):
        """The variance of designs in each scenario is attributed to their
        uncertain parameters."""
        store = ResultStore(self.path, "scenarios")
        study = self._study(SCENARIOS, store=store)
        a = attribute_variance(
            study, SPACE.enumerate(), ["y"], replicates=20, by=["scenario"]
        )
        shares = a.shares.set_index(["design", "scenario"])["total"]
        key = SPACE.enumerate()[0].key
        # y varies by its noise alone, but for the scenario that pins it
        self.assertAlmostEqual(shares[key, "low"], 1.0, delta=0.5)
        self.assertAlmostEqual(shares[key, "high"], 1.0, delta=0.5)
        self.assertTrue(np.isnan(shares[key, "pinned"]))
        self.calls.clear()
        attribute_variance(
            study, SPACE.enumerate(), ["y"], replicates=20, by=["scenario"]
        )
        self.assertEqual(self.calls, [])

    def test_invalid(self):
        """Scenarios set fixed values and uncertain parameters alone."""
        for scenarios in ({"a": {"x": 2}}, {"a": {"other": 1}}):
            with self.subTest(scenarios=scenarios):
                with self.assertRaises(ValueError):
                    self._study(scenarios)
        for scenarios in ({"a": 1}, {1: {}}, [("a", {})]):
            with self.subTest(scenarios=scenarios):
                with self.assertRaises(TypeError):
                    self._study(scenarios)
        named = DesignSpace(x=Integer(1, 2), scenario=Fixed("a"))
        with self.assertRaises(ValueError):
            Study(named, lambda d: d, _evaluate, scenarios={"a": {}})


def _results(rows):
    """Results of one replicate of designs in scenarios, from (design,
    scenario, perf, cost) rows (`None` for infeasible)."""
    records = [
        {
            "design": design,
            "replicate": 0,
            "scenario": scenario,
            "feasible": perf is not None,
            "perf": math.nan if perf is None else perf,
            "cost": math.nan if cost is None else cost,
        }
        for design, scenario, perf, cost in rows
    ]
    return pd.DataFrame(records)


# in s1, B dominates C (A and B are efficient); in s2, B dominates A and C
RESULTS = _results(
    [
        ("A", "s1", 10, 5),
        ("B", "s1", 8, 3),
        ("C", "s1", 7, 4),
        ("A", "s2", 6, 5),
        ("B", "s2", 9, 3),
        ("C", "s2", 9, 4),
    ]
)
OBJECTIVES = {"perf": "max", "cost": "min"}


class TestRobustness(unittest.TestCase):
    """
    Unit tests for the robustness of designs across scenarios.
    """

    def _robustness(self, results=RESULTS, criteria=OBJECTIVES, **options):
        return robustness(results, criteria, **options).set_index("design")

    def test_trace_and_regret(self):
        """The share of scenarios in which designs are efficient, and their
        regret: the improvement in every objective to match the front."""
        r = self._robustness()
        self.assertEqual(
            list(r.columns), ["pareto trace", "max regret", "s1 regret", "s2 regret"]
        )
        np.testing.assert_allclose(r["pareto trace"], [0.5, 1.0, 0.0])
        # C in s1: 1 short of B in both objectives (3 short of A in perf)
        np.testing.assert_allclose(r["s1 regret"], [0, 0, 1])
        # A in s2: 3 short of B in perf (2 in cost); C: 1 short of B in cost
        np.testing.assert_allclose(r["s2 regret"], [3, 0, 1])
        np.testing.assert_allclose(r["max regret"], [3, 0, 1])

    def test_scale(self):
        """Regret is in units of the scales of the objectives."""
        r = self._robustness(scale={"perf": 2})
        np.testing.assert_allclose(r["s1 regret"], [0, 0, 1])
        np.testing.assert_allclose(r["s2 regret"], [2, 0, 1])
        with self.assertRaises(ValueError):
            self._robustness(scale={"other": 1})
        with self.assertRaises(ValueError):
            self._robustness(scale={"perf": 0})

    def test_tolerance(self):
        """With tolerances, designs near the front are efficient, and regret
        is in tolerances by default."""
        criteria = Criteria(OBJECTIVES, tolerance={"cost": 2})
        r = self._robustness(criteria=criteria)
        np.testing.assert_allclose(r["pareto trace"], [0.5, 1.0, 1.0])
        # A in s2: 3 short of B and C in perf
        np.testing.assert_allclose(r["s2 regret"], [3, 0, 0])

    def test_requirements(self):
        """Designs that miss the requirements, or are infeasible, have
        infinite regret."""
        criteria = Criteria(OBJECTIVES, requirements=[Requirement("perf", ">=", 7)])
        r = self._robustness(criteria=criteria)
        self.assertEqual(r.loc["A", "s2 regret"], math.inf)
        self.assertEqual(r.loc["A", "max regret"], math.inf)
        infeasible = pd.concat([RESULTS, _results([("D", "s1", None, None)])])
        infeasible = pd.concat(
            [infeasible, _results([("D", "s2", 9, 2)])], ignore_index=True
        )
        r = self._robustness(results=infeasible)
        self.assertEqual(r.loc["D", "s1 regret"], math.inf)
        self.assertEqual(r.loc["D", "s2 regret"], 0)
        self.assertEqual(r.loc["D", "pareto trace"], 0.5)

    def test_invalid(self):
        """Robustness needs scenarios, and criteria by them alone."""
        with self.assertRaises(ValueError):
            robustness(RESULTS.drop(columns="scenario"), OBJECTIVES)
        with self.assertRaises(ValueError):
            robustness(RESULTS, Criteria(OBJECTIVES, by="other"))
        r = robustness(RESULTS, Criteria(OBJECTIVES, by="scenario"))
        self.assertEqual(len(r), 3)

    def test_study(self):
        """The robustness of the designs of a study in scenarios."""
        results = Study(SPACE, lambda d: d, _evaluate, scenarios=SCENARIOS).run()
        r = robustness(results, {"y": "max", "cost": "min"})
        # every design trades y against cost, in every scenario
        np.testing.assert_allclose(r["pareto trace"], 1.0)
        self.assertEqual(len(r), 3)


if __name__ == "__main__":
    unittest.main()
