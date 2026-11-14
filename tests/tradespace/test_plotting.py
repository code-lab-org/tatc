"""
Unit tests for plots of the results of trade studies.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    exceedance,
    pareto_front,
    plot_designs,
    plot_exceedance,
    plot_means,
    plot_progress,
    plot_shares,
    use_style,
)

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

    def test_exceedance(self):
        """A matrix of the probabilities, in the order of the designs, labeled."""
        results = pd.DataFrame(
            {
                "design": ["a", "a", "b", "b", "c", "c"],
                "replicate": [0, 1, 0, 1, 0, 1],
                "cost": [1.0, 2.0, 1.5, 1.5, 3.0, 4.0],
                "scenario": "nominal",
            }
        )
        compared = exceedance(results, "cost")
        ax = plot_exceedance(compared, names={"a": "first"})
        matrix = ax.images[0].get_array()
        self.assertEqual(matrix.shape, (3, 3))
        self.assertEqual(matrix[2, 0], 1.0)  # c costs more than a
        self.assertEqual(matrix[0, 1], 0.5)  # a and b, each once
        self.assertEqual(
            [t.get_text() for t in ax.get_yticklabels()], ["first", "b", "c"]
        )
        self.assertEqual(len(ax.figure.axes), 2)  # with a color bar
        twice = pd.concat([compared, compared.assign(scenario="other")])
        with self.assertRaises(ValueError):
            plot_exceedance(twice)


class TestHighlights(unittest.TestCase):
    """
    Unit tests for plots of some designs over others.
    """

    def tearDown(self):
        plt.close("all")

    def test_highlight(self):
        """Highlighted designs are colored over the others, in gray."""
        summary = pd.DataFrame(
            {
                "x": [1.0, 2.0, 3.0, 4.0],
                "y": [1.0, 2.0, 2.5, 4.0],
                "pareto": [True, True, False, True],
            }
        )
        ax = plot_designs(
            summary,
            "x",
            "y",
            highlight=[True, False, True, False],
            label="chosen",
            others="the rest",
            front=summary.pareto.to_numpy(),
            front_label="true front",
            colors={"a": "#123456"},
        )
        labels = [text.get_text() for text in ax.get_legend().get_texts()]
        self.assertEqual(labels, ["the rest", "chosen", "true front"])
        self.assertEqual(len(ax.collections[1].get_offsets()), 2)
        with self.assertRaises(ValueError):
            plot_designs(summary, "x", "y", highlight=[True])


class TestShares(unittest.TestCase):
    """
    Unit tests for plots of shares, groups of values, and the style.
    """

    def tearDown(self):
        plt.close("all")

    def test_shares(self):
        """Shares are grouped by factor or by metric, with their errors."""
        table = pd.DataFrame(
            {"a": [0.6, 0.2], "b": [0.1, 0.5], "a se": [0.05, 0.02], "active": [1, 1]},
            index=pd.Index(["x", "y"], name="parameter"),
        )
        ax = plot_shares(table, labels={"a": "first"})
        self.assertEqual(len(ax.patches), 4)  # two factors, two metrics
        self.assertEqual(
            [text.get_text() for text in ax.get_legend().get_texts()], ["first", "b"]
        )
        flipped = plot_shares(
            table, groups="metrics", names={"x": "ex"}, horizontal=True
        )
        self.assertEqual(
            [label.get_text() for label in flipped.get_yticklabels()], ["a", "b"]
        )
        with self.assertRaises(ValueError):
            plot_shares(table, groups="other")

    def test_groups_of_values(self):
        """Designs are colored by values given for their rows."""
        summary = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]})
        ax = plot_designs(summary, "x", "y", by=[1, 1, 2], names={1: "one", 2: "two"})
        self.assertEqual(
            [text.get_text() for text in ax.get_legend().get_texts()], ["one", "two"]
        )

    def test_style(self):
        """The style of plots applies to the figures of the session."""
        import matplotlib as mpl  # pylint: disable=import-outside-toplevel

        with mpl.rc_context():
            use_style()
            self.assertFalse(mpl.rcParams["axes.spines.top"])
            self.assertTrue(mpl.rcParams["figure.constrained_layout.use"])
