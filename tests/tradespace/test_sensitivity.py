"""
Unit tests for sensitivity analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import main_effects, sensitivity


class TestSensitivity(unittest.TestCase):
    """
    Unit tests for sensitivity analysis.
    """

    def setUp(self):
        rng = np.random.default_rng(0)
        a = np.repeat([1, 2, 3, 4], 50)
        b = np.tile(["x", "y"], 100)
        c = rng.uniform(0, 1, 200)
        branch = np.where(np.arange(200) < 100, 1.0, np.nan)
        self.frame = pd.DataFrame(
            {
                "a": a,
                "b": b,
                "c": c,
                "branch": branch,
                # a drives the metric, b a little, c not at all
                "metric": 10 * a
                + np.where(b == "x", 1.0, 0.0)
                + rng.normal(0, 0.1, 200),
            }
        )

    def test_shares(self):
        """Variables that drive a metric explain most of its variance."""
        result = sensitivity(self.frame, ["a", "b", "c", "branch"], ["metric"])
        self.assertEqual(result.index[0], "a")
        self.assertGreater(result.loc["a", "metric"], 0.95)
        self.assertLess(result.loc["b", "metric"], 0.05)
        self.assertLess(result.loc["c", "metric"], 0.1)
        self.assertEqual(result.loc["branch", "active"], 0.5)
        self.assertEqual(result.loc["a", "active"], 1.0)
        # a variable with one level where active explains nothing within it
        self.assertAlmostEqual(result.loc["branch", "metric"], 0)
        # the correlation ratios are no smaller than the adjusted shares
        ratios = sensitivity(self.frame, ["a", "b", "c"], ["metric"], adjusted=False)
        self.assertTrue(
            (ratios.metric >= result.loc[ratios.index, "metric"] - 1e-12).all()
        )

    def test_main_effects(self):
        """Main effects are the means of the metric at each level."""
        effects = main_effects(self.frame, "a", "metric")
        np.testing.assert_allclose(effects["mean"], [10.5, 20.5, 30.5, 40.5], atol=0.1)
        self.assertEqual(list(effects["count"]), [50] * 4)
        binned = main_effects(self.frame, "c", "metric", max_levels=4)
        self.assertEqual(len(binned), 4)

    def test_results_raise(self):
        """The results of a study (a row per replicate) are not summaries."""
        results = pd.DataFrame(
            {"design": ["x", "x", "y", "y"], "a": [1, 1, 2, 2], "metric": [1, 2, 3, 4]}
        )
        with self.assertRaisesRegex(ValueError, "aggregate"):
            sensitivity(results, ["a"], ["metric"])
        with self.assertRaisesRegex(ValueError, "aggregate"):
            main_effects(results, "a", "metric")

    def test_irrelevant_variable(self):
        """The adjusted share of a variable that does not matter is near 0."""
        rng = np.random.default_rng(1)
        adjusted, ratios = [], []
        for _ in range(100):
            frame = pd.DataFrame(
                {"v": rng.integers(0, 10, 20), "m": rng.normal(size=20)}
            )
            adjusted.append(sensitivity(frame, ["v"], ["m"]).loc["v", "m"])
            ratios.append(
                sensitivity(frame, ["v"], ["m"], adjusted=False).loc["v", "m"]
            )
        # (k - 1) / (n - 1) = 9 / 19 by chance alone
        self.assertGreater(np.mean(ratios), 0.35)
        self.assertLess(np.mean(adjusted), 0.1)
        # no more designs than levels
        frame = pd.DataFrame({"v": [1, 2, 3], "m": [1.0, 2.0, 4.0]})
        self.assertTrue(np.isnan(sensitivity(frame, ["v"], ["m"]).loc["v", "m"]))
