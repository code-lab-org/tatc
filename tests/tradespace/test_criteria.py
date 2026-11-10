"""
Unit tests for statistics over replicates and requirements on them.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    DesignSpace,
    Integer,
    Maximum,
    Mean,
    Minimum,
    Probability,
    Quantile,
    Requirement,
    StandardDeviation,
    Study,
    Uncertain,
    aggregate,
    pareto_probability,
    screen,
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
            Probability("x", ">=", 2).reduce(VALUES), [2 / 3, 1, np.nan]
        )
        np.testing.assert_allclose(
            Probability("x", "<", 2).reduce(VALUES), [1 / 3, 0, np.nan]
        )

    def test_invalid(self):
        """Quantiles and operators are checked."""
        with self.assertRaises(ValueError):
            Quantile("x", 1.5)
        with self.assertRaises(ValueError):
            Probability("x", "==", 1)
        with self.assertRaises(ValueError):
            Requirement("x", "!=", 1)


class TestRequirement(unittest.TestCase):
    """
    Unit tests for requirements.
    """

    def test_violation(self):
        """Violations are the gaps to the required values."""
        statistic = np.array([5.0, 10.0, 12.0, np.nan])
        np.testing.assert_allclose(
            Requirement("x", ">=", 10).violation(statistic), [5, 0, 0, np.nan]
        )
        np.testing.assert_allclose(
            Requirement("x", "<=", 10).violation(statistic), [0, 0, 2, np.nan]
        )
        strict = Requirement("x", ">", 10).violation(statistic)
        self.assertGreater(strict[1], 0)
        self.assertEqual(strict[2], 0)


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

    def test_statistics_and_requirements(self):
        """Statistics, violations, and design values are aggregated."""
        frame = aggregate(
            _results(),
            {"value": "value", "value p25": Quantile("value", 0.25), "cost": "cost"},
            [Requirement("value p25", ">=", 5), Requirement("cost", "<=", 1.5)],
        ).set_index("design")
        self.assertEqual(frame.loc["a", "value"], 13)
        self.assertAlmostEqual(frame.loc["b", "value p25"], 5)
        self.assertEqual(frame.loc["a", "violation"], 0.5)
        self.assertEqual(frame.loc["c", "violation"], 4)
        self.assertEqual(list(frame.meets_requirements), [False, False, False])
        self.assertEqual(frame.loc["c", "size"], "large")
        self.assertTrue((frame.replicates == 4).all())

    def test_metric_names(self):
        """Metrics given as names are their means, and requirements add theirs."""
        frame = aggregate(_results(), ["value"], [Requirement("cost", "<=", 1)])
        self.assertIn("cost", frame.columns)
        self.assertEqual(list(frame.meets_requirements), [False, False, True])
        # a requirement on a metric that the results do not have
        with self.assertRaises(ValueError):
            aggregate(_results(), {"v": Mean("value")}, [Requirement("w", ">=", 1)])


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
            {"cost": "min", "value q": "max"},
            statistics={"value q": Quantile("value", 0.25)},
        ).set_index("design")
        self.assertTrue(pessimistic.loc["a", "pareto"])
        self.assertFalse(pessimistic.loc["b", "pareto"])

    def test_requirements(self):
        """Designs that miss a requirement are not efficient."""
        probability = pareto_probability(
            _results(),
            {"cost": "min", "value": "max"},
            requirements=[Requirement("reliable", ">=", 0.75)],
            statistics={"reliable": Probability("value", ">=", 10)},
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
            {"cost": "min", "value": "max"},
            keep=4,
            requirements=[Requirement("cost", "<=", 3)],
        )
        first = screening.summaries[0].set_index("design")
        meets = first[first.meets_requirements]
        misses = first[~first.meets_requirements]
        self.assertGreater(misses.pareto_rank.min(), meets.pareto_rank.max())
        # the nearest miss (cost 4) is kept to fill four designs
        self.assertEqual(sorted(first[first.kept].cost), [1, 2, 3, 4])
