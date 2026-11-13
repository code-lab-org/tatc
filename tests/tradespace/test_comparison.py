"""
Unit tests for comparisons of designs by the replicates of a metric.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import exceedance


def _results(values, feasible=None, **columns):
    """Results of designs (by key) in replicates, from a mapping of keys to
    the metric of each replicate."""
    rows = []
    for design, metrics in values.items():
        for replicate, value in enumerate(metrics):
            rows.append({"design": design, "replicate": replicate, "cost": value})
    frame = pd.DataFrame(rows)
    frame["feasible"] = True if feasible is None else feasible
    for name, value in columns.items():
        frame[name] = value
    return frame


class TestExceedance(unittest.TestCase):
    """
    Unit tests for the probability that the metric of a design exceeds that
    of another.
    """

    def test_pairs(self):
        """Shares of the replicates in which a design's metric is greater, ties
        counting as half."""
        results = _results(
            {"a": [1.0, 2.0, 3.0, 4.0], "b": [2.0, 2.0, 1.0, 5.0], "c": [9.0] * 4}
        )
        compared = exceedance(results, "cost").set_index(["design", "other"])
        # designs in the order of their means: a (2.5), b (2.5), c (9)
        self.assertEqual(
            list(compared.index.get_level_values(0).unique()), ["a", "b", "c"]
        )
        self.assertAlmostEqual(compared.loc[("a", "b"), "exceedance"], 1.5 / 4)
        self.assertAlmostEqual(compared.loc[("b", "a"), "exceedance"], 2.5 / 4)
        self.assertEqual(compared.loc[("a", "a"), "exceedance"], 0.5)
        self.assertEqual(compared.loc[("c", "a"), "exceedance"], 1.0)
        self.assertAlmostEqual(compared.loc[("a", "b"), "cost difference"], 0.0)
        self.assertAlmostEqual(compared.loc[("c", "a"), "cost difference"], 6.5)
        self.assertTrue((compared["replicates"] == 4).all())
        self.assertEqual(len(compared), 9)

    def test_order_and_infeasible(self):
        """Given designs keep their order, and only replicates in which both
        designs are feasible count."""
        results = _results(
            {"a": [1.0, 5.0, 1.0], "b": [2.0, 2.0, 2.0]},
            feasible=[True, False, True, True, True, True],
        )
        compared = exceedance(results, "cost", ["b", "a"])
        self.assertEqual(list(compared["design"]), ["b", "b", "a", "a"])
        pair = compared[(compared.design == "b") & (compared.other == "a")].iloc[0]
        self.assertEqual(pair["replicates"], 2)
        self.assertEqual(pair["exceedance"], 1.0)
        missing = _results({"a": [np.nan, np.nan], "b": [1.0, 2.0]})
        pair = exceedance(missing, "cost").set_index(["design", "other"])
        self.assertTrue(np.isnan(pair.loc[("a", "b"), "exceedance"]))
        self.assertEqual(pair.loc[("a", "b"), "replicates"], 0)

    def test_groups(self):
        """Designs are compared in each group of `by`, which several records of
        an evaluation need."""
        hot = _results({"a": [1.0, 2.0], "b": [3.0, 4.0]}, scenario="hot")
        cold = _results({"a": [5.0, 6.0], "b": [3.0, 4.0]}, scenario="cold")
        results = pd.concat([hot, cold], ignore_index=True)
        with self.assertRaises(ValueError):
            exceedance(results, "cost")
        compared = exceedance(results, "cost", by="scenario").set_index(
            ["scenario", "design", "other"]
        )
        self.assertEqual(compared.loc[("hot", "a", "b"), "exceedance"], 0.0)
        self.assertEqual(compared.loc[("cold", "a", "b"), "exceedance"], 1.0)

    def test_infinite_and_infeasible(self):
        """Infinite values are ordered (equal ones tie), designs infeasible in
        every replicate come last, and results without feasible ones give
        pairs of no replicates."""
        results = _results({"a": [np.inf, 2.0], "b": [np.inf, 1.0]})
        compared = exceedance(results, "cost").set_index(["design", "other"])
        self.assertEqual(compared.loc[("a", "a"), "exceedance"], 0.5)
        self.assertEqual(compared.loc[("a", "b"), "exceedance"], 0.75)
        self.assertEqual(compared.loc[("b", "a"), "exceedance"], 0.25)
        results = _results(
            {"a": [3.0, 4.0], "b": [1.0, 2.0]}, feasible=[False, False, True, True]
        )
        compared = exceedance(results, "cost")
        self.assertEqual(list(compared["design"].unique()), ["b", "a"])
        self.assertTrue(compared[compared.design == "a"]["exceedance"].isna().all())
        results["feasible"] = False
        compared = exceedance(results, "cost")
        self.assertEqual(len(compared), 4)
        self.assertTrue((compared["replicates"] == 0).all())
        empty = exceedance(results.iloc[:0], "cost")
        self.assertTrue(empty.empty)
        self.assertIn("exceedance", empty.columns)

    def test_invalid(self):
        """The metric and the designs must be in the results."""
        results = _results({"a": [1.0], "b": [2.0]})
        with self.assertRaises(ValueError):
            exceedance(results, "mass")
        with self.assertRaises(ValueError):
            exceedance(results, "cost", ["a", "z"])
        with self.assertRaises(TypeError):
            exceedance(results, "cost", [1, 2])


if __name__ == "__main__":
    unittest.main()
