"""
Unit tests for statistics of the metrics of designs over their replicates.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
import warnings

import numpy as np

from tatc.tradespace import (
    Maximum,
    Mean,
    Minimum,
    Probability,
    Quantile,
    StandardDeviation,
    StandardError,
    Statistic,
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

    def test_equality(self):
        """Statistics of the same kind and parameters are equal."""
        self.assertEqual(Quantile("x", 0.2), Quantile("x", 0.2))
        self.assertNotEqual(Quantile("x", 0.2), Quantile("x", 0.3))
        self.assertNotEqual(Mean("x"), Maximum("x"))
        self.assertEqual(len({Mean("x"), Mean("x"), Mean("y")}), 2)
        self.assertIsInstance(Mean("x"), Statistic)


class TestNames(unittest.TestCase):
    """
    Unit tests for the conventional names of statistics.
    """

    def test_names(self):
        """Statistics are named by convention."""
        self.assertEqual(Mean("cost").name, "cost")
        self.assertEqual(StandardDeviation("cost").name, "cost std")
        self.assertEqual(StandardError("cost").name, "cost se")
        self.assertEqual(Minimum("cost").name, "cost min")
        self.assertEqual(Maximum("cost").name, "cost max")
        self.assertEqual(Quantile("cost", 0.7).name, "cost p70")
        self.assertEqual(Probability("cost", "<=", 150).name, "P(cost <= 150)")

    def test_columns(self):
        """A statistic needs the column of its metric."""
        self.assertEqual(Quantile("cost", 0.7).columns, ("cost",))
        self.assertEqual(Quantile("cost", 0.7).metric, "cost")


if __name__ == "__main__":
    unittest.main()
