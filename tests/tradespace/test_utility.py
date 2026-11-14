"""
Unit tests for utilities of the metrics of designs.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Criteria,
    ExponentialUtility,
    LinearUtility,
    Mean,
    PiecewiseUtility,
    Probability,
    Quantile,
    Requirement,
    Utility,
    aggregate,
    pareto_front,
)


class TestSingleAttribute(unittest.TestCase):
    """
    Unit tests for single-attribute utilities.
    """

    def test_linear(self):
        """Linear utilities are 0 at the worst value and 1 at the best."""
        u = LinearUtility(0, 100)
        np.testing.assert_allclose(
            u([-10, 0, 50, 100, 120, np.nan]), [0, 0, 0.5, 1, 1, np.nan]
        )
        # a metric to minimize
        self.assertAlmostEqual(float(LinearUtility(10, 2)(6)), 0.5)
        with self.assertRaises(ValueError):
            LinearUtility(1, 1)

    def test_exponential(self):
        """Exponential utilities are concave if averse, convex if seeking."""
        np.testing.assert_allclose(ExponentialUtility(0, 10, 2)([0, 10]), [0, 1])
        averse = float(ExponentialUtility(0, 10, 2)(5))
        self.assertAlmostEqual(averse, (1 - math.exp(-1)) / (1 - math.exp(-2)))
        self.assertGreater(averse, 0.5)
        self.assertLess(float(ExponentialUtility(0, 10, -2)(5)), 0.5)
        np.testing.assert_allclose(
            ExponentialUtility(0, 10)([2, 7]), LinearUtility(0, 10)([2, 7])
        )
        self.assertTrue(np.isnan(ExponentialUtility(0, 10, 2)(np.nan)))

    def test_piecewise(self):
        """Piecewise utilities interpolate between points."""
        u = PiecewiseUtility([(50, 0.8), (0, 0), (100, 1)])
        np.testing.assert_allclose(u([-5, 25, 50, 75, 200]), [0, 0.4, 0.8, 0.9, 1])
        decreasing = PiecewiseUtility([(0, 1), (10, 0)])
        self.assertAlmostEqual(float(decreasing(2.5)), 0.75)
        self.assertTrue(np.isnan(u(np.nan)))
        for points in ([(0, 0)], [(0, 0), (0, 1)], [(0, 0), (1, 2)]):
            with self.subTest(points=points):
                with self.assertRaises(ValueError):
                    PiecewiseUtility(points)


class TestUtility(unittest.TestCase):
    """
    Unit tests for multi-attribute utilities.
    """

    def setUp(self):
        self.attributes = {"a": LinearUtility(0, 1), "b": LinearUtility(0, 1)}

    def test_additive(self):
        """Additive utilities are weighted sums, by normalized weights."""
        u = Utility(self.attributes, {"a": 3, "b": 1})
        np.testing.assert_allclose(
            u.evaluate(
                {"a": np.array([1.0, 0.0, 0.5]), "b": np.array([0.0, 1.0, 0.5])}
            ),
            [0.75, 0.25, 0.5],
        )
        self.assertEqual(u.scaling, 0.0)
        self.assertEqual(u.columns, ("a", "b"))

    def test_multiplicative(self):
        """Multiplicative utilities solve their scaling constant from the
        weights: complements if they sum to less than 1, substitutes if more."""
        complements = Utility(self.attributes, {"a": 0.5, "b": 0.3}, "multiplicative")
        # 1 + K = (1 + 0.5 K)(1 + 0.3 K)
        self.assertAlmostEqual(complements.scaling, 4 / 3, places=9)
        substitutes = Utility(self.attributes, {"a": 0.7, "b": 0.6}, "multiplicative")
        self.assertAlmostEqual(substitutes.scaling, -0.3 / 0.42, places=9)
        for u in (complements, substitutes):
            values = {
                "a": np.array([0.0, 1.0, 1.0, 0.0]),
                "b": np.array([0.0, 1.0, 0.0, 1.0]),
            }
            weights = u.weights
            np.testing.assert_allclose(
                u.evaluate(values), [0, 1, weights["a"], weights["b"]], atol=1e-9
            )
        # complements reward doing well in both, substitutes in either
        half = {"a": np.array([0.5]), "b": np.array([0.5])}
        self.assertLess(complements.evaluate(half)[0], substitutes.evaluate(half)[0])
        # weights that sum to 1 are additive
        additive = Utility(self.attributes, {"a": 0.4, "b": 0.6}, "multiplicative")
        self.assertEqual(additive.scaling, 0.0)

    def test_frame(self):
        """Utilities of the rows of a data frame."""
        u = Utility(self.attributes, {"a": 1, "b": 1})
        frame = pd.DataFrame({"a": [1.0, 0.0], "b": [1.0, 0.5]}, index=["x", "y"])
        pd.testing.assert_series_equal(
            u(frame), pd.Series([1.0, 0.25], index=["x", "y"], name="utility")
        )
        with self.assertRaises(ValueError):
            u(frame.drop(columns="b"))

    def test_invalid(self):
        """Utilities need a positive weight of each attribute."""
        for weights, form in (
            ({"a": 1}, "additive"),
            ({"a": 1, "b": 1, "c": 1}, "additive"),
            ({"a": 0, "b": 1}, "additive"),
            ({"a": 1, "b": 1}, "other"),
            ({"a": 1.5, "b": 0.5}, "multiplicative"),
        ):
            with self.subTest(weights=weights, form=form):
                with self.assertRaises(ValueError):
                    Utility(self.attributes, weights, form)
        with self.assertRaises(ValueError):
            Utility({}, {})


def _results():
    """Two designs in two replicates: x varies between replicates, y does not
    (with the same mean of a)."""
    return pd.DataFrame(
        {
            "design": ["x", "x", "y", "y"],
            "replicate": [0, 1, 0, 1],
            "feasible": True,
            "a": [0.0, 1.0, 0.5, 0.5],
            "cost": [1.0, 1.0, 2.0, 2.0],
        }
    )


class TestStatistics(unittest.TestCase):
    """
    Unit tests for statistics of utilities.
    """

    def setUp(self):
        self.utility = Utility({"a": ExponentialUtility(0, 1, 2)}, {"a": 1})

    def test_expected_utility(self):
        """The mean of a utility is its expected value over the replicates,
        which penalizes the spread of a risk-averse utility."""
        criteria = Criteria({"utility": "max"}, statistics={"utility": self.utility})
        self.assertEqual(criteria.metrics, ["a"])
        summary = aggregate(_results(), criteria).set_index("design")
        expected = (0 + 1) / 2
        self.assertAlmostEqual(summary.loc["x", "utility"], expected)
        self.assertAlmostEqual(
            summary.loc["y", "utility"], float(ExponentialUtility(0, 1, 2)(0.5))
        )
        # the utility of the means is the same for both
        means = aggregate(_results(), ["a"]).set_index("design")
        np.testing.assert_allclose(self.utility(means), self.utility(means).iloc[0])

    def test_statistics(self):
        """Utilities are metrics of any statistic."""
        statistics = {
            "worst": Quantile(self.utility, 0),
            "reliable": Probability(self.utility, ">=", 0.5),
            "expected": Mean(self.utility),
        }
        summary = aggregate(_results(), statistics=statistics).set_index("design")
        self.assertEqual(summary.loc["x", "worst"], 0.0)
        self.assertEqual(summary.loc["x", "reliable"], 0.5)
        self.assertEqual(summary.loc["y", "reliable"], 1.0)
        self.assertEqual(Mean(self.utility), Mean(self.utility))
        self.assertEqual(len({Mean(self.utility), Mean(self.utility)}), 1)

    def test_utility_against_cost(self):
        """Designs trade the utility of their metrics against cost."""
        criteria = Criteria(
            {"utility": "max", "cost": "min"},
            statistics={"utility": self.utility},
            requirements=[Requirement("worst", ">=", 0.1)],
        )
        # a requirement on a statistic that the criteria do not name
        with self.assertRaises(ValueError):
            aggregate(_results(), criteria)
        criteria = Criteria(
            {"utility": "max", "cost": "min"},
            statistics={"utility": self.utility, "worst": Quantile(self.utility, 0)},
            requirements=[Requirement("worst", ">=", 0.1)],
        )
        summary = aggregate(_results(), criteria)
        front = summary[pareto_front(summary, criteria)]
        # x is cheaper, but misses the requirement on its worst replicate
        self.assertEqual(list(front.design), ["y"])


if __name__ == "__main__":
    unittest.main()
