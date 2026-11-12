"""
Unit tests for scenarios of studies and the robustness of designs across them.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import math
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Across,
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
    optimize_optuna,
    optimize_pymoo,
    pareto_front,
    pareto_probability,
    robustness,
)
from tatc.tradespace.optimization import Evaluations

HAS_OPTUNA = importlib.util.find_spec("optuna") is not None
HAS_PYMOO = importlib.util.find_spec("pymoo") is not None

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
            study, ["y"], SPACE.enumerate(), replicates=20, by=["scenario"]
        )
        shares = a.shares.set_index(["design", "scenario"])["total"]
        key = SPACE.enumerate()[0].key
        # y varies by its noise alone, but for the scenario that pins it
        self.assertAlmostEqual(shares[key, "low"], 1.0, delta=0.5)
        self.assertAlmostEqual(shares[key, "high"], 1.0, delta=0.5)
        self.assertTrue(np.isnan(shares[key, "pinned"]))
        self.calls.clear()
        attribute_variance(
            study, ["y"], SPACE.enumerate(), replicates=20, by=["scenario"]
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
            list(r.columns), ["pareto_trace", "max_regret", "s1 regret", "s2 regret"]
        )
        np.testing.assert_allclose(r["pareto_trace"], [0.5, 1.0, 0.0])
        # C in s1: 1 short of B in both objectives (3 short of A in perf)
        np.testing.assert_allclose(r["s1 regret"], [0, 0, 1])
        # A in s2: 3 short of B in perf (2 in cost); C: 1 short of B in cost
        np.testing.assert_allclose(r["s2 regret"], [3, 0, 1])
        np.testing.assert_allclose(r["max_regret"], [3, 0, 1])

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
        np.testing.assert_allclose(r["pareto_trace"], [0.5, 1.0, 1.0])
        # A in s2: 3 short of B and C in perf
        np.testing.assert_allclose(r["s2 regret"], [3, 0, 0])

    def test_requirements(self):
        """Designs that miss the requirements, or are infeasible, have
        infinite regret."""
        criteria = Criteria(OBJECTIVES, requirements=[Requirement("perf", ">=", 7)])
        r = self._robustness(criteria=criteria)
        self.assertEqual(r.loc["A", "s2 regret"], math.inf)
        self.assertEqual(r.loc["A", "max_regret"], math.inf)
        infeasible = pd.concat([RESULTS, _results([("D", "s1", None, None)])])
        infeasible = pd.concat(
            [infeasible, _results([("D", "s2", 9, 2)])], ignore_index=True
        )
        r = self._robustness(results=infeasible)
        self.assertEqual(r.loc["D", "s1 regret"], math.inf)
        self.assertEqual(r.loc["D", "s2 regret"], 0)
        self.assertEqual(r.loc["D", "pareto_trace"], 0.5)

    def test_invalid(self):
        """Robustness needs scenarios, and criteria by them alone."""
        with self.assertRaises(ValueError):
            robustness(RESULTS.drop(columns="scenario"), OBJECTIVES)
        with self.assertRaises(ValueError):
            robustness(RESULTS, Criteria(OBJECTIVES, by="other"))
        r = robustness(RESULTS, Criteria(OBJECTIVES, by="scenario"))
        self.assertEqual(len(r), 3)
        with self.assertRaises(ValueError):
            robustness(RESULTS, Criteria(OBJECTIVES, across=Across()))

    def test_study(self):
        """The robustness of the designs of a study in scenarios."""
        results = Study(SPACE, lambda d: d, _evaluate, scenarios=SCENARIOS).run()
        r = robustness(results, {"y": "max", "cost": "min"})
        # every design trades y against cost, in every scenario
        np.testing.assert_allclose(r["pareto_trace"], 1.0)
        self.assertEqual(len(r), 3)


class TestAcross(unittest.TestCase):
    """
    Unit tests for criteria that combine statistics across scenarios.
    """

    def _summary(self, criteria, results=RESULTS):
        summary = aggregate(results, criteria)
        summary["pareto"] = pareto_front(summary, criteria)
        return summary.set_index("design")

    def test_worst(self):
        """Each objective is its worst over the scenarios, and designs are
        judged once."""
        summary = self._summary(Criteria(OBJECTIVES, across=Across()))
        self.assertEqual(list(summary.index), ["A", "B", "C"])
        self.assertNotIn("scenario", summary.columns)
        np.testing.assert_allclose(summary.perf, [6, 8, 7])
        np.testing.assert_allclose(summary.cost, [5, 3, 4])
        np.testing.assert_array_equal(summary.replicates, [1, 1, 1])
        self.assertEqual(list(summary.pareto), [False, True, False])

    def test_mean(self):
        """Each objective is its mean over the scenarios, weighted by their
        weights (normalized over the scenarios of the results)."""
        summary = self._summary(Criteria(OBJECTIVES, across=Across(rule="mean")))
        np.testing.assert_allclose(summary.perf, [8, 8.5, 8])
        weighted = Across(rule="mean", weights={"s1": 3, "s2": 1, "s3": 2})
        summary = self._summary(Criteria(OBJECTIVES, across=weighted))
        np.testing.assert_allclose(summary.perf, [9, 8.25, 7.5])
        self.assertEqual(list(summary.pareto), [True, True, False])
        unweighted = Across("scenario", "mean", {"s1": 1})
        with self.assertRaises(ValueError):
            aggregate(RESULTS, Criteria(OBJECTIVES, across=unweighted))

    def test_requirements(self):
        """A requirement is met by the combined statistic: in every scenario,
        by the worst; a statistic of another name may be combined otherwise
        than the objectives."""
        required = [Requirement("perf", ">=", 7)]
        criteria = Criteria(OBJECTIVES, requirements=required, across=Across())
        summary = self._summary(criteria)
        self.assertEqual(list(summary.meets_requirements), [False, True, True])
        criteria = Criteria(
            OBJECTIVES,
            statistics={"perf floor": "perf"},
            requirements=[Requirement("perf floor", ">=", 7)],
            across=Across(rule="mean", rules={"perf floor": "worst"}),
        )
        summary = self._summary(criteria)
        np.testing.assert_allclose(summary.perf, [8, 8.5, 8])
        np.testing.assert_allclose(summary["perf floor"], [6, 8, 7])
        self.assertEqual(list(summary.meets_requirements), [False, True, True])
        self.assertAlmostEqual(summary.loc["A", "violation"], 1 / 7)
        # the worst of a requirement to stay at most a value is its largest
        below = Criteria(
            {"perf": "max"},
            statistics={"spent": "cost"},
            requirements=[Requirement("spent", "<=", 4)],
            across=Across(),
        )
        summary = self._summary(below)
        self.assertEqual(list(summary.meets_requirements), [False, True, True])

    def test_incomplete(self):
        """Designs without feasible results in every scenario are not judged."""
        more = _results(
            [("D", "s1", 20, 1), ("E", "s1", 20, 1), ("E", "s2", None, None)]
        )
        results = pd.concat([RESULTS, more], ignore_index=True)
        summary = self._summary(Criteria(OBJECTIVES, across=Across()), results)
        self.assertEqual(list(summary.index), ["A", "B", "C"])

    def test_agrees_with_scenarios(self):
        """The worst and the mean of the statistics in a study's scenarios."""
        results = Study(
            SPACE, lambda d: d, _evaluate, replicates=3, scenarios=SCENARIOS
        ).run()
        objectives = {"y": "max", "cost": "min"}
        by_scenario = aggregate(results, Criteria(objectives, by="scenario"))
        y = by_scenario.pivot(index="design", columns="scenario", values="y")
        for rule, expected in (("worst", y.min(axis=1)), ("mean", y.mean(axis=1))):
            criteria = Criteria(objectives, across=Across(rule=rule))
            summary = aggregate(results, criteria).set_index("design")
            np.testing.assert_allclose(summary.y, expected.loc[summary.index])
        probability = pareto_probability(
            results, Criteria(objectives, across=Across()), n_resamples=20
        ).set_index("design")
        np.testing.assert_allclose(probability.y, y.min(axis=1).loc[probability.index])
        np.testing.assert_array_equal(probability.replicates, 3)
        self.assertTrue(probability.pareto.all())

    def test_resamples_share_replicates(self):
        """A bootstrap draws each replicate in every scenario at once."""
        rows = [
            {
                "design": "A",
                "replicate": r,
                "scenario": s,
                "feasible": True,
                "perf": v,
            }
            for r in range(6)
            for s, v in (("s1", float(r)), ("s2", 10.0 - r))
        ]
        results = pd.DataFrame(rows)
        mean = Criteria({"perf": "max"}, across=Across(rule="mean"))
        probability = pareto_probability(results, mean, n_resamples=50)
        self.assertAlmostEqual(probability["perf"].iloc[0], 5.0)
        self.assertAlmostEqual(probability["perf se"].iloc[0], 0.0)
        worst = Criteria({"perf": "max"}, across=Across())
        probability = pareto_probability(results, worst, n_resamples=50)
        self.assertGreater(probability["perf se"].iloc[0], 0.0)

    def test_invalid(self):
        """Rules, weights, and the statistics combined by their worst are
        checked."""
        with self.assertRaises(ValueError):
            Across(rule="best")
        with self.assertRaises(ValueError):
            Across(rules={"perf": "best"})
        with self.assertRaises(ValueError):
            Across(rule="mean", weights={"s1": 0})
        with self.assertRaises(TypeError):
            Across(column="")
        with self.assertRaises(TypeError):
            Criteria(OBJECTIVES, across="scenario")
        with self.assertRaises(ValueError):
            Criteria(OBJECTIVES, by="scenario", across=Across())
        with self.assertRaises(ValueError):
            Criteria(OBJECTIVES, across=Across(rules={"other": "mean"}))
        # a statistic without a direction has no worst, unless combined by its mean
        with self.assertRaises(ValueError):
            Criteria(OBJECTIVES, statistics={"extra": "perf"}, across=Across())
        Criteria(
            OBJECTIVES,
            statistics={"extra": "perf"},
            across=Across(rules={"extra": "mean"}),
        )
        # an objective and a requirement that disagree on its worst
        with self.assertRaises(ValueError):
            Criteria(
                OBJECTIVES,
                requirements=[Requirement("perf", "<=", 9)],
                across=Across(),
            )
        with self.assertRaises(ValueError):
            aggregate(
                RESULTS.drop(columns="scenario"), Criteria(OBJECTIVES, across=Across())
            )


def _scenario_study():
    return Study(SPACE, lambda d: d, _evaluate, replicates=2, scenarios=SCENARIOS)


class TestSearchAcross(unittest.TestCase):
    """
    Unit tests for searches across scenarios.
    """

    def test_evaluations(self):
        """A search judges each design across the scenarios, as aggregate does."""
        criteria = Criteria({"y": "max", "cost": "min"}, across=Across())
        evaluations = Evaluations(_scenario_study(), criteria)
        designs = SPACE.enumerate()
        values = evaluations.evaluate(designs)
        summary = aggregate(evaluations.results, criteria).set_index("design")
        for design, (y, cost) in zip(designs, values):
            self.assertAlmostEqual(y, summary.loc[design.key, "y"])
            self.assertAlmostEqual(cost, summary.loc[design.key, "cost"])
        self.assertEqual(len(evaluations.results), 3 * 2 * len(SCENARIOS))

    def test_needs_across(self):
        """A search of a study with several scenarios needs criteria across them."""
        evaluations = Evaluations(_scenario_study(), {"y": "max", "cost": "min"})
        with self.assertRaisesRegex(ValueError, "across"):
            evaluations.evaluate(SPACE.enumerate()[:1])
        # a study with one scenario has one record per evaluation
        one = Study(SPACE, lambda d: d, _evaluate, scenarios={"high": {"level": 3.0}})
        Evaluations(one, {"y": "max", "cost": "min"}).evaluate(SPACE.enumerate())

    @unittest.skipUnless(HAS_PYMOO, "requires pymoo")
    def test_pymoo(self):
        """pymoo searches across the scenarios."""
        # pylint: disable-next=import-outside-toplevel
        from pymoo.algorithms.moo.nsga2 import NSGA2

        criteria = Criteria({"y": "max", "cost": "min"}, across=Across(rule="mean"))
        search = optimize_pymoo(
            _scenario_study(), criteria, NSGA2(pop_size=4), ("n_gen", 2), seed=1
        )
        self.assertEqual(set(search.results.scenario), set(SCENARIOS))
        self.assertGreater(len(search.designs), 0)

    @unittest.skipUnless(HAS_OPTUNA, "requires optuna")
    def test_optuna(self):
        """Optuna searches across the scenarios."""
        criteria = Criteria({"y": "max", "cost": "min"}, across=Across())
        search = optimize_optuna(_scenario_study(), criteria, n_trials=6)
        self.assertEqual(set(search.results.scenario), set(SCENARIOS))


if __name__ == "__main__":
    unittest.main()
