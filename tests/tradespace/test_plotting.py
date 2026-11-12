"""
Unit tests for plots of the results of trade studies.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import pareto_front, plot_designs, plot_means, plot_progress

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    plt = None


@unittest.skipIf(plt is None, "requires matplotlib")
class TestPlots(unittest.TestCase):
    """
    Unit tests for plots.
    """

    def setUp(self):
        self.summary = pd.DataFrame(
            {
                "design": list("abcdef"),
                "satellites": [1, 2, 2, 3, 4, 4],
                "planes": [1, 1, 2, 1, 2, 1],
                "sufficient": [10.0, 30.0, 20.0, 50.0, 60.0, 40.0],
            }
        )

    def tearDown(self):
        plt.close("all")

    def test_designs(self):
        """Designs are colored by group, and the efficient ones ringed and labeled."""
        objectives = {"satellites": "min", "sufficient": "max"}
        self.summary["pareto"] = pareto_front(self.summary, objectives)
        ax = plot_designs(
            self.summary,
            "satellites",
            "sufficient",
            by="planes",
            front="pareto",
            labels="design",
            names=lambda planes: f"{planes} planes",
            jitter=0.2,
        )
        # a collection of each group, and of the efficient designs
        self.assertEqual(len(ax.collections), 3)
        efficient = int(self.summary.pareto.sum())
        self.assertEqual(len(ax.collections[-1].get_offsets()), efficient)
        self.assertEqual(
            sorted(t.get_text() for t in ax.texts),
            sorted(self.summary.design[self.summary.pareto]),
        )
        legend = [t.get_text() for t in ax.get_legend().get_texts()]
        self.assertEqual(legend, ["1 planes", "2 planes", "Pareto-efficient"])
        # jitter shifts only the horizontal positions
        offsets = np.concatenate([c.get_offsets() for c in ax.collections[:2]])
        self.assertTrue(np.all(np.abs(offsets[:, 0] - np.round(offsets[:, 0])) <= 0.2))
        # criteria find the front
        ax = plot_designs(self.summary, "satellites", "sufficient", front=objectives)
        self.assertEqual(len(ax.collections[-1].get_offsets()), efficient)
        self.assertEqual(ax.get_xlabel(), "satellites")

    def test_given_axes(self):
        """Plots draw on given axes."""
        _, ax = plt.subplots()
        self.assertIs(plot_designs(self.summary, "satellites", "sufficient", ax=ax), ax)
        self.assertIsNone(ax.get_legend())

    def test_progress(self):
        """Searches are lines, labeled at their ends."""
        progress = {
            "one": pd.Series([1.0, 2.0, 3.0], index=[10, 20, 30]),
            "two": pd.Series([2.0, 2.5], index=[10, 20]),
        }
        ax = plot_progress(progress)
        self.assertEqual(len(ax.lines), 2)
        self.assertEqual([t.get_text() for t in ax.texts], ["one", "two"])
        self.assertEqual(ax.get_xlim()[0], 0)

    def test_means(self):
        """Means over replicates, with bands, by group (feasible records only)."""
        results = pd.DataFrame(
            {
                "satellites": [1, 1, 2, 2, 1, 2],
                "sensor": ["a", "a", "a", "a", "b", "b"],
                "feasible": [True, True, True, False, True, True],
                "resolved": [1.0, 3.0, 4.0, np.nan, 5.0, 7.0],
            }
        )
        ax = plot_means(results, "satellites", "resolved", by="sensor")
        self.assertEqual(len(ax.lines), 2)
        self.assertEqual(len(ax.collections), 2)  # the bands
        np.testing.assert_allclose(ax.lines[0].get_ydata(), [2.0, 4.0])
        single = plot_means(results, "satellites", "resolved")
        self.assertEqual(len(single.lines), 1)
