"""
Unit tests for statistics over replicates and requirements on them.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
import warnings

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Criteria,
    DesignSpace,
    Integer,
    Maximum,
    Mean,
    Minimum,
    Probability,
    Quantile,
    Requirement,
    StandardDeviation,
    StandardError,
    Study,
    Uncertain,
    aggregate,
    pareto_front,
    pareto_probability,
    screen,
    summarize,
)

VALUES = np.array(
    [[1.0, 2.0, 3.0, np.nan], [4.0, np.nan, np.nan, np.nan], [np.nan] * 4]
)


class TestStatistics(unittest.TestCase):
    """
    Unit tests for statistics.
    """

    def test_reductions(self):
        """Statistics ignore missing replicates."""
        np.testing.assert_allclose(Mean("x").reduce(VALUES), [2, 4, np.nan])
        np.testing.assert_allclose(Minimum("x").reduce(VALUES), [1, 4, np.nan])
        np.testing.assert_allclose(Maximum("x").reduce(VALUES), [3, 4, np.nan])
        np.testing.assert_allclose(Quantile("x", 0.25).reduce(VALUES), [1.5, 4, np.nan])
        np.testing.assert_allclose(
            StandardDeviation("x").reduce(VALUES), [1, np.nan, np.nan]
        )
        np.testing.assert_allclose(
            StandardError("x").reduce(VALUES), [1 / np.sqrt(3), np.nan, np.nan]
        )
        np.testing.assert_allclose(
            Probability("x", ">=", 2).reduce(VALUES), [2 / 3, 1, np.nan]
        )
        np.testing.assert_allclose(
            Probability("x", "<", 2).reduce(VALUES), [1 / 3, 0, np.nan]
        )

    def test_quantiles_with_missing_values(self):
        """Quantiles are numpy's, of the values of each design."""
        rng = np.random.default_rng(5)
        values = rng.normal(size=(200, 12))
        values[rng.random(values.shape) < 0.4] = np.nan
        values[0] = np.nan
        values[1, 1:] = np.nan
        for q in (0, 0.1, 0.5, 0.7, 1):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                expected = np.nanquantile(values, q, axis=1)
            np.testing.assert_allclose(Quantile("x", q).reduce(values), expected)

    def test_invalid(self):
        """Quantiles and operators are checked."""
        with self.assertRaises(ValueError):
            Quantile("x", 1.5)
        with self.assertRaises(ValueError):
            Probability("x", "==", 1)
        with self.assertRaises(ValueError):
            Requirement("x", "!=", 1)
        with self.assertRaises(ValueError):
            Requirement("x", ">=", 1, scale=0)

    def test_equality(self):
        """Statistics of the same kind and parameters are equal."""
        self.assertEqual(Quantile("x", 0.2), Quantile("x", 0.2))
        self.assertNotEqual(Quantile("x", 0.2), Quantile("x", 0.3))
        self.assertNotEqual(Mean("x"), Maximum("x"))
        self.assertEqual(len({Mean("x"), Mean("x"), Mean("y")}), 2)


class TestRequirement(unittest.TestCase):
    """
    Unit tests for requirements.
    """

    def test_violation(self):
        """Violations are the gaps to the required values, relative to them."""
        statistic = np.array([5.0, 10.0, 12.0, np.nan])
        np.testing.assert_allclose(
            Requirement("x", ">=", 10).violation(statistic), [0.5, 0, 0, np.nan]
        )
        np.testing.assert_allclose(
            Requirement("x", "<=", 10).violation(statistic), [0, 0, 0.2, np.nan]
        )
        strict = Requirement("x", ">", 10).violation(statistic)
        self.assertGreater(strict[1], 0)
        self.assertEqual(strict[2], 0)

    def test_scale(self):
        """Violations are relative to a scale: given, or 1 for a value of 0."""
        statistic = np.array([5.0, 12.0])
        np.testing.assert_allclose(
            Requirement("x", ">=", 10, scale=1).violation(statistic), [5, 0]
        )
        np.testing.assert_allclose(
            Requirement("x", "<=", 0).violation(statistic), [5, 12]
        )
        np.testing.assert_allclose(
            Requirement("x", "<=", -10).violation(statistic), [1.5, 2.2]
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


class TestAggregate(unittest.TestCase):
    """
    Unit tests for aggregating results.
    """

    def test_statistics(self):
        """Statistics and design values are aggregated, without violations."""
        frame = aggregate(
            _results(),
            {"value": "value", "value p25": Quantile("value", 0.25), "cost": "cost"},
        ).set_index("design")
        self.assertEqual(frame.loc["a", "value"], 13)
        self.assertAlmostEqual(frame.loc["b", "value p25"], 5)
        self.assertEqual(frame.loc["c", "size"], "large")
        self.assertTrue((frame.replicates == 4).all())
        self.assertNotIn("violation", frame.columns)
        self.assertEqual(list(aggregate(_results(), ["value"]).value), [13, 17.5, 1])

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
            list(frame.columns[-6:-2]), ["cost", "value", "value p25", "replicates"]
        )
        self.assertAlmostEqual(frame.loc["a", "violation"], 0.5 / 1.5)
        self.assertAlmostEqual(frame.loc["c", "violation"], 4 / 5)
        self.assertEqual(list(frame.meets_requirements), [False, False, False])
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


class TestCriteria(unittest.TestCase):
    """
    Unit tests for criteria.
    """

    def test_names(self):
        """Names are those of objectives, requirements, and other statistics."""
        criteria = Criteria(
            {"a": "min", "b": "max"},
            statistics={"c": Mean("x"), "b": "y"},
            requirements=[Requirement("c", ">=", 1), Requirement("a", "<=", 2)],
            by="sensor",
        )
        self.assertEqual(criteria.names, ["a", "b", "c"])
        self.assertEqual(criteria.metrics, ["a", "x", "y"])
        self.assertEqual(criteria.statistic("a"), Mean("a"))
        self.assertEqual(criteria.statistic("b"), Mean("y"))
        self.assertEqual(criteria.by, ("sensor",))

    def test_invalid(self):
        """Objectives, statistics, requirements, and tolerances are checked."""
        with self.assertRaises(ValueError):
            Criteria({})
        with self.assertRaises(ValueError):
            Criteria({"a": "lowest"})
        with self.assertRaises(ValueError):
            Criteria({"a": "min"}, tolerance={"b": 1})
        with self.assertRaises(ValueError):
            Criteria({"a": "min"}, tolerance={"a": "5"})
        with self.assertRaises(ValueError):
            Criteria({"a": "min"}, tolerance={"a": -1})
        with self.assertRaises(TypeError):
            Criteria({"a": "min"}, statistics={"b": 1})
        with self.assertRaises(TypeError):
            Criteria({"a": "min"}, requirements=[("a", "<=", 1)])
        with self.assertRaises(TypeError):
            pareto_front(pd.DataFrame({"a": [1]}), ["a"])

    def test_shorthand(self):
        """A mapping of objectives is the criteria of their means."""
        frame = pd.DataFrame({"cost": [1, 2, 3], "value": [1, 3, 2]})
        objectives = {"cost": "min", "value": "max"}
        pd.testing.assert_series_equal(
            pareto_front(frame, objectives), pareto_front(frame, Criteria(objectives))
        )

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


class TestParetoProbabilityRequirements(unittest.TestCase):
    """
    Unit tests for the probability of efficiency with statistics and requirements.
    """

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


class TestScreenRequirements(unittest.TestCase):
    """
    Unit tests for screening with requirements.
    """

    def test_misses_ranked_last(self):
        """Designs that miss requirements are ranked after those that meet them."""
        space = DesignSpace(a=Integer(1, 6), noise=Uncertain.uniform(0, 1))
        study = Study(
            space,
            lambda d: d,
            lambda d: {"cost": d.a, "value": d.a**0.5 + d.noise},
            replicates=2,
        )
        screening = screen(
            space.enumerate(),
            [study, study],
            Criteria(
                {"cost": "min", "value": "max"},
                requirements=[Requirement("cost", "<=", 3)],
            ),
            keep=4,
        )
        first = screening.summaries[0].set_index("design")
        meets = first[first.meets_requirements]
        misses = first[~first.meets_requirements]
        self.assertGreater(misses.pareto_rank.min(), meets.pareto_rank.max())
        # the nearest miss (cost 4) is kept to fill four designs
        self.assertEqual(sorted(first[first.kept].cost), [1, 2, 3, 4])
