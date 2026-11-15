"""
Unit tests for summaries of the results of studies: the statistics of each
design over its replicates (`aggregate` and `summarize`), its judgment by
criteria, and the probability that it is efficient (`pareto_probability`).

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
import warnings

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Across,
    Criteria,
    Probability,
    Quantile,
    Requirement,
    aggregate,
    pareto_probability,
    summarize,
)


def _results():
    """Results of three designs in four replicates."""
    values = {"a": [10, 12, 14, 16], "b": [5, 30, 5, 30], "c": [1, 1, 1, 1]}
    return pd.DataFrame(
        [
            {
                "design": key,
                "replicate": r,
                "feasible": True,
                "cost": {"a": 2, "b": 2, "c": 1}[key],
                "value": v,
                "size": {"a": "small", "b": "small", "c": "large"}[key],
            }
            for key, row in values.items()
            for r, v in enumerate(row)
        ]
    )


def _partly_infeasible(scenario=None, shares=None):
    """Results of three designs in four replicates: a is feasible in every
    replicate, b in two, and c in none (whose records lack metrics)."""
    values = {"a": [1, 2, 3, 4], "b": [10, None, 30, None], "c": [None] * 4}
    values.update(shares or {})
    return pd.DataFrame(
        [
            {
                "design": key,
                "replicate": r,
                **({} if scenario is None else {"scenario": scenario}),
                "feasible": v is not None,
                "cost": (
                    {"a": 2.0, "b": 1.0, "c": 1.0}[key] if v is not None else np.nan
                ),
                "value": np.nan if v is None else float(v),
            }
            for key, row in values.items()
            for r, v in enumerate(row)
        ]
    )


def _replicates(values, cost=1.0):
    """Results of designs (rows) in replicates (columns), with a cost."""
    return pd.DataFrame(
        [
            {
                "design": f"d{i}",
                "replicate": r,
                "cost": cost,
                "value": v,
                "feasible": True,
            }
            for i, row in enumerate(values)
            for r, v in enumerate(row)
        ]
    )


def _scenarios(values):
    """Results of a design in scenarios (by name), from its values in each
    replicate (in order)."""
    return pd.DataFrame(
        [
            {
                "design": "d",
                "replicate": r,
                "scenario": scenario,
                "feasible": True,
                "v": float(v),
            }
            for scenario, row in values.items()
            for r, v in enumerate(row)
        ]
    )


class TestAggregate(unittest.TestCase):
    """
    Unit tests for aggregating results.
    """

    def test_statistics(self):
        """Statistics and design values are aggregated, without violations."""
        frame = aggregate(
            _results(),
            statistics={
                "value": "value",
                "value p25": Quantile("value", 0.25),
                "cost": "cost",
            },
        ).set_index("design")
        self.assertEqual(frame.loc["a", "value"], 13)
        self.assertAlmostEqual(frame.loc["b", "value p25"], 5)
        self.assertEqual(frame.loc["c", "size"], "large")
        self.assertTrue((frame.replicates == 4).all())
        self.assertNotIn("violation", frame.columns)
        self.assertEqual(list(aggregate(_results(), ["value"]).value), [13, 17.5, 1])
        # a mapping of statistics by name is not one of objectives
        with self.assertRaisesRegex(ValueError, "statistics"):
            aggregate(_results(), {"value": "value"})
        with self.assertRaises(ValueError):
            aggregate(_results())

    def test_objectives(self):
        """A mapping of objectives are criteria of their means, with others."""
        frame = aggregate(
            _results(),
            {"value": "max"},
            statistics={"value p25": Quantile("value", 0.25)},
        ).set_index("design")
        self.assertEqual(frame.loc["a", "value"], 13)
        self.assertAlmostEqual(frame.loc["b", "value p25"], 5)
        self.assertTrue(frame.meets_requirements.all())
        # statistics that the criteria name otherwise
        with self.assertRaisesRegex(ValueError, "otherwise"):
            aggregate(
                _results(), {"value": "max"}, statistics={"value": Quantile("value", 1)}
            )
        # a statistic the criteria name the same way is not a clash
        same = aggregate(_results(), {"value": "max"}, statistics={"value": "value"})
        self.assertEqual(list(same.value), [13, 17.5, 1])
        # objectives grouped by a column
        self.assertEqual(len(aggregate(_results(), {"value": "max"}, by="size")), 3)

    def test_criteria(self):
        """Criteria add the violations of their requirements, relative to them."""
        criteria = Criteria(
            {"cost": "min", "value": "max"},
            statistics={"value p25": Quantile("value", 0.25)},
            requirements=[
                Requirement("value p25", ">=", 5),
                Requirement("cost", "<=", 1.5),
            ],
        )
        frame = aggregate(_results(), criteria).set_index("design")
        self.assertEqual(
            list(frame.columns[-8:]),
            [
                "cost",
                "value",
                "value p25",
                "replicates",
                "feasible",
                "violation",
                "meets_requirements",
                "pareto",
            ],
        )
        self.assertAlmostEqual(frame.loc["a", "violation"], 0.5 / 1.5)
        self.assertAlmostEqual(frame.loc["c", "violation"], 4 / 5)
        self.assertEqual(list(frame.meets_requirements), [False, False, False])
        # none meets the requirements, so that none is efficient
        self.assertEqual(list(frame.pareto), [False, False, False])
        with self.assertRaises(ValueError):
            aggregate(_results(), criteria, by=["size"])

    def test_requirement_on_metric(self):
        """A requirement on a metric is on its mean."""
        criteria = Criteria(
            {"value": "max"}, requirements=[Requirement("cost", "<=", 1)]
        )
        frame = aggregate(_results(), criteria)
        self.assertIn("cost", frame.columns)
        self.assertEqual(list(frame.meets_requirements), [False, False, True])
        # a requirement on a metric that the results do not have
        missing = Criteria({"value": "max"}, requirements=[Requirement("w", ">=", 1)])
        with self.assertRaises(ValueError):
            aggregate(_results(), missing)

    def test_statistics_of_names(self):
        """Statistics given by names are the means of the metrics, and other
        statistics are named by convention."""
        summary = aggregate(_results(), {"value": "max"}, statistics=["cost"])
        self.assertEqual(list(summary.cost), [2.0, 2.0, 1.0])
        self.assertEqual(
            list(aggregate(_results(), statistics=[Quantile("value", 0.25)]).columns)[
                -3:
            ],
            ["value p25", "replicates", "feasible"],
        )

    def test_names_of_columns(self):
        """A statistic of the name of a column of the results is its mean."""
        results = _results().assign(**{"P(value >= 10)": 0.5})
        statistic = Probability("value", ">=", 10)
        criteria = Criteria({"P(value >= 10)": "max"}, statistics=[statistic])
        for function in (aggregate, summarize, pareto_probability):
            with self.subTest(function=function.__name__):
                with self.assertRaisesRegex(ValueError, "P\\(value >= 10\\)"):
                    function(results, criteria)
        # the mean of the column, of its name
        summary = aggregate(results, {"P(value >= 10)": "max"})
        np.testing.assert_allclose(summary["P(value >= 10)"], 0.5)
        renamed = Criteria({"reliable": "max"}, statistics={"reliable": statistic})
        np.testing.assert_allclose(
            aggregate(results, renamed).reliable, [1.0, 0.5, 0.0]
        )
        # nor do the columns that summaries add
        with self.assertRaisesRegex(ValueError, "value se"):
            summarize(_results().assign(**{"value se": 1.0}), ["value"])
        with self.assertRaisesRegex(ValueError, "value se"):
            pareto_probability(_results().assign(**{"value se": 1.0}), {"value": "max"})

    def test_records_of_groups(self):
        """Several records of an evaluation need the columns that separate them."""
        results = pd.concat([_results().assign(sensor=s) for s in ("x", "y")])
        with self.assertRaisesRegex(ValueError, r"by.*\['sensor'\]"):
            aggregate(results, ["value"])
        with self.assertRaisesRegex(ValueError, "sensor"):
            pareto_probability(results, {"value": "max"})
        frame = aggregate(results, Criteria({"value": "max"}, by="sensor"))
        self.assertEqual(len(frame), 6)
        self.assertEqual(len(summarize(results, ["value"], by="sensor")), 6)


class TestSummarize(unittest.TestCase):
    """
    Unit tests for summaries with the spread of the means of metrics.
    """

    def test_summarize_criteria(self):
        """Summaries of criteria add the spread of their means of metrics."""
        criteria = Criteria(
            {"cost": "min", "value": "max"},
            statistics={"value p25": Quantile("value", 0.25)},
        )
        summary = summarize(_results(), criteria)
        for column in ("value std", "value se", "pareto"):
            self.assertIn(column, summary.columns)
        self.assertNotIn("value p25 se", summary.columns)
        grouped = summarize(_results(), {"value": "max"}, by="size")
        self.assertEqual(len(grouped), 3)
        with self.assertRaises(ValueError):
            summarize(_results(), Criteria({"value": "max"}), by="size")

    def test_summarize_statistics(self):
        """Summaries take other statistics, with the spread of their means."""
        results = _results().assign(cost=lambda f: f.cost + 0.1 * f.replicate)
        summary = summarize(results, {"value": "max"}, statistics=["cost"])
        for column in ("value", "value se", "cost", "cost std", "cost se", "pareto"):
            self.assertIn(column, summary.columns)
        # each mean followed by its spread, without criteria
        self.assertEqual(
            list(summarize(results, ["value", "cost"]).columns[-8:]),
            [
                "value",
                "value std",
                "value se",
                "cost",
                "cost std",
                "cost se",
                "replicates",
                "feasible",
            ],
        )

    def test_constant_spreads(self):
        """The spreads of metrics that are the same in every replicate of each
        design are left out."""
        summary = summarize(_results(), ["value", "cost"])
        self.assertNotIn("cost std", summary.columns)
        self.assertNotIn("cost se", summary.columns)
        spread = summary.set_index("design")[["value std", "value se"]]
        np.testing.assert_allclose(spread.loc["c"].to_numpy(dtype=float), 0.0)
        # of designs of one replicate, as of summaries without designs
        single = summarize(_results()[_results().replicate == 0], ["value"])
        self.assertNotIn("value se", single.columns)
        empty = summarize(_results().iloc[:0], ["value"])
        self.assertIn("value se", empty.columns)
        # as are the standard errors of the bootstrap
        probability = pareto_probability(_results(), {"cost": "min", "value": "max"})
        self.assertNotIn("cost se", probability.columns)
        self.assertIn("value se", probability.columns)

    def test_spread_of_worst(self):
        """Across groups, the spread of the worst is that of its group, as
        is the standard error of the bootstrap."""
        results = _scenarios({"a": [10] * 8, "b": [5, 9] * 4})
        worst = Criteria({"v": "max"}, across=Across())
        summary = summarize(results, worst).iloc[0]
        b = results[results.scenario == "b"].v
        self.assertEqual(summary.v, 7.0)
        self.assertAlmostEqual(summary["v std"], b.std())
        self.assertAlmostEqual(summary["v se"], b.std() / np.sqrt(8))
        bootstrap = pareto_probability(results, worst, n_resamples=2000)
        # the bootstrap's standard deviation, of n rather than n - 1 replicates
        self.assertAlmostEqual(
            bootstrap["v se"].iloc[0], b.std(ddof=0) / np.sqrt(8), delta=0.05
        )

    def test_spread_of_mean(self):
        """Across groups, the spread of a weighted mean is that of the
        weighted sum of each replicate's values, of groups that share
        replicates, and that of independent groups otherwise."""
        rng = np.random.default_rng(2)
        x = rng.normal(0, 1, 30)
        weights = {"s1": 3, "s2": 1}
        mean = Criteria({"v": "max"}, across=Across(rule="mean", weights=weights))
        # the same replicates, of correlated values
        paired = _scenarios({"s1": x, "s2": 2 * x + rng.normal(0, 0.1, 30)})
        summary = summarize(paired, mean).iloc[0]
        wide = paired.pivot(index="replicate", columns="scenario", values="v")
        combined = 0.75 * wide.s1 + 0.25 * wide.s2
        self.assertAlmostEqual(summary.v, combined.mean())
        self.assertAlmostEqual(summary["v std"], combined.std())
        self.assertAlmostEqual(summary["v se"], combined.std() / np.sqrt(30))
        bootstrap = pareto_probability(paired, mean, n_resamples=2000)
        self.assertAlmostEqual(
            bootstrap["v se"].iloc[0] / summary["v se"], 1.0, delta=0.1
        )
        # replicates that no other group shares: independent groups
        independent = paired.assign(
            replicate=paired.replicate + 100 * (paired.scenario == "s2")
        )
        summary = summarize(independent, mean).iloc[0]
        errors = wide.std() / np.sqrt(30)
        self.assertAlmostEqual(
            summary["v se"],
            np.sqrt(0.75**2 * errors.s1**2 + 0.25**2 * errors.s2**2),
        )
        self.assertAlmostEqual(
            summary["v std"],
            np.sqrt(0.75**2 * wide.s1.var() + 0.25**2 * wide.s2.var()),
        )
        # what the groups share cancels
        cancel = Criteria({"v": "max"}, across=Across(rule="mean"))
        steady = _scenarios({"s1": range(6), "s2": [10 - r for r in range(6)]})
        summary = summarize(steady, cancel)
        self.assertEqual(summary.v.iloc[0], 5.0)
        self.assertNotIn("v se", summary.columns)
        self.assertNotIn("v se", pareto_probability(steady, cancel, n_resamples=50))


class TestFeasibility(unittest.TestCase):
    """
    Unit tests for the judgment of designs by their feasible replicates.
    """

    def test_feasible_share(self):
        """Designs are judged by their feasible replicates, of which
        `replicates` is the number and `feasible` the share; designs without a
        feasible replicate are left out."""
        summary = aggregate(_partly_infeasible(), {"cost": "min", "value": "max"})
        summary = summary.set_index("design")
        self.assertEqual(list(summary.index), ["a", "b"])
        self.assertEqual(summary.loc["b", "value"], 20)
        self.assertEqual(list(summary.replicates), [4, 2])
        self.assertEqual(list(summary.feasible), [1.0, 0.5])
        # b, cheaper and better in its feasible replicates, dominates a
        self.assertEqual(list(summary.pareto), [False, True])
        self.assertEqual(
            list(summary.columns[-5:]),
            ["replicates", "feasible", "violation", "meets_requirements", "pareto"],
        )

    def test_requirement(self):
        """Requirements on `feasible` judge the share of feasible replicates
        over every replicate of a design."""
        every = Criteria(
            {"cost": "min", "value": "max"},
            requirements=[Requirement("feasible", ">=", 1)],
        )
        summary = aggregate(_partly_infeasible(), every).set_index("design")
        self.assertEqual(list(summary.feasible), [1.0, 0.5])
        self.assertEqual(list(summary.meets_requirements), [True, False])
        self.assertAlmostEqual(summary.loc["b", "violation"], 0.5)
        self.assertEqual(list(summary.pareto), [True, False])
        half = Criteria(
            {"cost": "min", "value": "max"},
            requirements=[Requirement("feasible", ">=", 0.5)],
        )
        summary = aggregate(_partly_infeasible(), half)
        self.assertTrue(summary.meets_requirements.all())
        # its spread is not summarized
        summary = summarize(_partly_infeasible(), every)
        self.assertNotIn("feasible se", summary.columns)
        self.assertEqual(list(summary.feasible), [1.0, 0.5])

    def test_groups(self):
        """Infeasible records without a value of a group (such as a sensor)
        count in every group of the design."""
        results = _partly_infeasible()
        feasible = results[results.feasible]
        records = pd.concat(
            [
                feasible.assign(sensor="x"),
                feasible.assign(sensor="y", value=feasible.value + 1),
                results[~results.feasible],
            ],
            ignore_index=True,
        )
        summary = aggregate(records, ["value"], by="sensor")
        summary = summary.set_index(["sensor", "design"])
        self.assertEqual(summary.loc[("y", "b"), "value"], 21)
        self.assertEqual(summary.loc[("x", "b"), "feasible"], 0.5)
        self.assertEqual(summary.loc[("y", "b"), "feasible"], 0.5)
        self.assertEqual(summary.loc[("y", "a"), "feasible"], 1.0)

    def test_across(self):
        """Across groups, `feasible` is over the replicates of every group
        together, and designs need feasible replicates in every group."""
        results = pd.concat(
            [
                _partly_infeasible("s1", {"b": [10, 20, 30, 40]}),
                _partly_infeasible("s2"),
            ],
            ignore_index=True,
        )
        for rule in ("worst", "mean"):
            with self.subTest(rule=rule):
                criteria = Criteria(
                    {"cost": "min", "value": "max"},
                    requirements=[Requirement("feasible", ">=", 0.7)],
                    across=Across(rule=rule),
                )
                summary = aggregate(results, criteria).set_index("design")
                self.assertEqual(list(summary.index), ["a", "b"])
                self.assertEqual(summary.loc["b", "feasible"], 0.75)
                self.assertEqual(summary.loc["b", "replicates"], 2)
                self.assertTrue(summary.meets_requirements.all())

    def test_no_feasible_design(self):
        """Results without a feasible record have empty summaries."""
        results = _partly_infeasible()
        infeasible = results[results.design == "c"].drop(columns=["cost", "value"])
        criteria = Criteria({"cost": "min", "value": "max"}, by="size")
        for summary in (
            aggregate(infeasible, {"cost": "min", "value": "max"}),
            aggregate(infeasible.assign(size=1), criteria),
            aggregate(infeasible.iloc[:0], {"cost": "min", "value": "max"}),
        ):
            self.assertTrue(summary.empty)
            self.assertEqual(
                list(summary.columns[-5:]),
                ["replicates", "feasible", "violation", "meets_requirements", "pareto"],
            )
            # of the types of their columns, which select rows
            self.assertEqual(summary.pareto.dtype, bool)
            self.assertTrue(summary[summary.pareto].empty)
        self.assertTrue(aggregate(infeasible, ["value"]).empty)
        across = Criteria({"cost": "min", "value": "max"}, across=Across(rule="mean"))
        scenarios = pd.concat(
            [_partly_infeasible("s1"), _partly_infeasible("s2")], ignore_index=True
        )
        self.assertEqual(len(aggregate(scenarios, across)), 2)
        self.assertTrue(aggregate(scenarios.iloc[:0], across).empty)


class TestParetoProbability(unittest.TestCase):
    """
    Unit tests for the probability of Pareto efficiency.
    """

    def test_settled(self):
        """Designs separated in every replicate are efficient or not for sure."""
        results = _replicates([[1, 2, 3], [5, 6, 7]])
        probability = pareto_probability(results, {"cost": "min", "value": "max"})
        self.assertEqual(list(probability.pareto_probability), [0, 1])
        self.assertEqual(list(probability.pareto), [False, True])
        self.assertEqual(list(probability.value), [2, 6])

    def test_common_random_numbers(self):
        """Designs that differ by a constant in every replicate are settled."""
        rng = np.random.default_rng(0)
        common = rng.normal(0, 5, 8)
        results = _replicates([common, common + 0.1])
        probability = pareto_probability(results, {"cost": "min", "value": "max"})
        self.assertEqual(list(probability.pareto_probability), [0, 1])

    def test_unsettled(self):
        """Designs whose order changes between replicates are uncertain."""
        # equal means: d0 varies between replicates, d1 does not
        results = _replicates([[0, 2] * 4, [1] * 8])
        probability = pareto_probability(
            results, {"cost": "min", "value": "max"}, n_resamples=300
        )
        p = probability.pareto_probability
        self.assertTrue(((p > 0.05) & (p < 0.95)).all())

    def test_groups(self):
        """Fronts of groups are found separately."""
        results = pd.concat(
            [
                _replicates([[1, 2], [3, 4]]).assign(sensor="a"),
                _replicates([[3, 4], [1, 2]]).assign(sensor="b"),
            ]
        )
        probability = pareto_probability(
            results, Criteria({"cost": "min", "value": "max"}, by=["sensor"])
        )
        winners = probability[probability.pareto_probability == 1]
        self.assertEqual(
            sorted(zip(winners.sensor, winners.design)), [("a", "d1"), ("b", "d0")]
        )

    def test_design_values(self):
        """Probabilities keep the values of the designs, as summaries do."""
        results = _replicates([[1, 2], [3, 4]])
        results["size"] = results["design"].map({"d0": "small", "d1": "large"})
        results["seed"] = results["replicate"]  # an uncertain parameter
        probability = pareto_probability(results, {"cost": "min", "value": "max"})
        self.assertEqual(list(probability["size"]), ["small", "large"])
        self.assertNotIn("seed", probability.columns)
        bare = pareto_probability(
            results, {"cost": "min", "value": "max"}, design_columns=[]
        )
        self.assertNotIn("size", bare.columns)

    def test_quantile_objective(self):
        """A pessimistic objective prefers the design with less spread."""
        objectives = {"cost": "min", "value": "max"}
        by_mean = pareto_probability(_results(), objectives).set_index("design")
        self.assertTrue(by_mean.loc["b", "pareto"])  # b has the larger mean
        pessimistic = pareto_probability(
            _results(),
            Criteria(
                {"cost": "min", "value q": "max"},
                statistics={"value q": Quantile("value", 0.25)},
            ),
        ).set_index("design")
        self.assertTrue(pessimistic.loc["a", "pareto"])
        self.assertFalse(pessimistic.loc["b", "pareto"])

    def test_requirements(self):
        """Designs that miss a requirement are not efficient."""
        probability = pareto_probability(
            _results(),
            Criteria(
                {"cost": "min", "value": "max"},
                requirements=[Requirement("reliable", ">=", 0.75)],
                statistics={"reliable": Probability("value", ">=", 10)},
            ),
        ).set_index("design")
        self.assertEqual(list(probability.meets_requirements), [True, False, False])
        # a is efficient unless b meets the requirement (and dominates it)
        self.assertEqual(probability.loc["a", "requirements_probability"], 1)
        self.assertGreater(probability.loc["a", "pareto_probability"], 0.5)
        self.assertLess(probability.loc["a", "pareto_probability"], 1)
        self.assertEqual(probability.loc["c", "pareto_probability"], 0)
        # b meets the threshold in half of its replicates: it meets the
        # requirement in some resamples, and is then efficient
        self.assertGreater(probability.loc["b", "requirements_probability"], 0.1)
        self.assertLess(probability.loc["b", "requirements_probability"], 0.9)
        self.assertEqual(probability.loc["c", "requirements_probability"], 0)

    def test_unequal_replicates(self):
        """Designs may have different replicates."""
        results = _results()
        results = results[~((results.design == "c") & (results.replicate > 0))]
        probability = pareto_probability(results, {"cost": "min", "value": "max"})
        self.assertEqual(list(probability.replicates), [4, 4, 1])
        self.assertEqual(
            probability.set_index("design").loc["c", "pareto_probability"], 1
        )

    def test_standard_errors(self):
        """The bootstrap gives the standard error of each objective."""
        rng = np.random.default_rng(0)
        results = pd.DataFrame(
            [
                {
                    "design": design,
                    "replicate": r,
                    "feasible": True,
                    "x": float(i),
                    "y": i + rng.normal(0, 2),
                }
                for i, design in enumerate("ab")
                for r in range(40)
            ]
        )
        probability = pareto_probability(
            results, {"x": "min", "y": "max"}, n_resamples=2000
        )
        # x is the same in every replicate of each design
        self.assertNotIn("x se", probability.columns)
        expected = results.groupby("design").y.std(ddof=0) / np.sqrt(40)
        np.testing.assert_allclose(
            probability.set_index("design")["y se"], expected, rtol=0.15
        )

    def test_resamples_share_replicates(self):
        """A bootstrap draws each replicate in every scenario at once."""
        results = _scenarios({"s1": range(6), "s2": [10 - r for r in range(6)]})
        mean = Criteria({"v": "max"}, across=Across(rule="mean"))
        probability = pareto_probability(results, mean, n_resamples=50)
        self.assertAlmostEqual(probability["v"].iloc[0], 5.0)
        # the mean of every resample is 5
        self.assertNotIn("v se", probability.columns)
        worst = Criteria({"v": "max"}, across=Across())
        probability = pareto_probability(results, worst, n_resamples=50)
        self.assertGreater(probability["v se"].iloc[0], 0.0)


class TestFeasibleProbability(unittest.TestCase):
    """
    Unit tests for the probability of efficiency of designs infeasible in
    some replicates.
    """

    @staticmethod
    def _results(feasible):
        """Two designs in eight replicates: d1, better, is feasible only where
        given."""
        return pd.DataFrame(
            [
                {
                    "design": design,
                    "replicate": r,
                    "feasible": True if design == "d0" else feasible[r],
                    "cost": 1.0,
                    "value": 1.0 if design == "d0" else 2.0,
                }
                for design in ("d0", "d1")
                for r in range(len(feasible))
            ]
        )

    def test_feasible_requirement(self):
        """A requirement on the share of feasible replicates is uncertain in
        the bootstrap."""
        results = self._results([True, False] * 4)
        criteria = Criteria(
            {"cost": "min", "value": "max"},
            requirements=[Requirement("feasible", ">=", 0.5)],
        )
        probability = pareto_probability(results, criteria, n_resamples=200)
        probability = probability.set_index("design")
        self.assertEqual(probability.loc["d1", "feasible"], 0.5)
        self.assertEqual(probability.loc["d1", "replicates"], 4)
        met = probability.loc["d1", "requirements_probability"]
        self.assertGreater(met, 0.2)
        self.assertLess(met, 0.9)
        self.assertEqual(probability.loc["d0", "requirements_probability"], 1.0)
        # without the requirement, d1 is judged by its feasible replicates,
        # but for the few resamples that draw none of them (1 in 256)
        plain = pareto_probability(results, {"cost": "min", "value": "max"})
        plain = plain.set_index("design")
        self.assertGreater(plain.loc["d1", "pareto_probability"], 0.98)
        self.assertAlmostEqual(
            plain.loc["d0", "pareto_probability"],
            1 - plain.loc["d1", "pareto_probability"],
        )

    def test_infeasible_draws(self):
        """A design whose drawn replicates are all infeasible is not judged
        by its undrawn replicates."""
        rows = []
        for r in range(5):
            rows.append(
                {
                    "design": "rare",
                    "replicate": r,
                    "feasible": r == 0,
                    "value": 10.0 if r == 0 else np.nan,
                }
            )
            rows.append(
                {"design": "steady", "replicate": r, "feasible": True, "value": 9.0}
            )
        results = pd.DataFrame(rows)
        probability = pareto_probability(
            results, {"value": "max"}, n_resamples=4000, seed=1
        ).set_index("design")
        # rare is efficient in the resamples that draw its replicate 0
        drawn = 1 - 0.8**5
        self.assertAlmostEqual(
            probability.loc["rare", "pareto_probability"], drawn, delta=0.03
        )
        self.assertAlmostEqual(
            probability.loc["steady", "pareto_probability"], 1 - drawn, delta=0.03
        )
        # nor does it meet the requirements without a feasible drawn replicate
        self.assertAlmostEqual(
            probability.loc["rare", "requirements_probability"], drawn, delta=0.03
        )
        self.assertEqual(probability.loc["steady", "requirements_probability"], 1.0)

    def test_no_feasible_design(self):
        """Results without a feasible record have an empty summary."""
        results = self._results([False] * 4)
        results = results[results.design == "d1"].drop(columns=["cost", "value"])
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            probability = pareto_probability(results, {"cost": "min", "value": "max"})
        self.assertTrue(probability.empty)
        self.assertEqual(
            list(probability.columns[-6:]),
            [
                "feasible",
                "violation",
                "meets_requirements",
                "pareto",
                "pareto_probability",
                "requirements_probability",
            ],
        )
        self.assertIn("value se", probability.columns)
