"""
Unit tests for the tables of trade studies: the results of studies and the
summaries of their designs.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import Mean
from tatc.tradespace.tables import (
    check_records,
    check_summary,
    design_values,
    draw_resamples,
    efficient_keys,
    feasible_records,
    groups,
    matrices,
    names,
    value_columns,
)


def _results():
    """Results of two designs in two replicates, one infeasible."""
    return pd.DataFrame(
        {
            "design": ["a", "a", "b", "b"],
            "replicate": [0, 1, 0, 1],
            "feasible": [True, True, True, False],
            "size": [1, 1, 2, 2],
            "seed": [7, 8, 7, 8],
            "value": [1.0, 3.0, 2.0, np.nan],
        }
    )


class TestTables(unittest.TestCase):
    """
    Unit tests for tables of results and summaries.
    """

    def test_names(self):
        """Names are given as a name or a sequence of names."""
        self.assertEqual(names("a"), ["a"])
        self.assertEqual(names(("a", "b")), ["a", "b"])

    def test_groups(self):
        """Groups are the positions of the rows of each value."""
        frame = _results()
        self.assertEqual([key for key, _ in groups(frame, ())], [()])
        found = {key: list(rows) for key, rows in groups(frame, ["size"])}
        self.assertEqual(found, {(1,): [0, 1], (2,): [2, 3]})
        with self.assertRaises(ValueError):
            list(groups(frame, ["sensor"]))

    def test_records(self):
        """Results need one record per evaluation, and a summary one row per
        design."""
        results = _results()
        check_records(results)
        doubled = pd.concat([results, results.assign(sensor="y", value=0.0)])
        with self.assertRaisesRegex(ValueError, r"as `by`, such as \['sensor'\]$"):
            check_records(
                doubled.fillna({"sensor": "x"}), hint="as `by`", metrics=["value"]
            )
        # columns of the evaluation, and metrics, do not separate them
        check_records(doubled, ["design", "replicate", "sensor"])
        with self.assertRaisesRegex(ValueError, "aggregate"):
            check_summary(results)
        check_summary(results, by=["replicate"])

    def test_matrices(self):
        """The values of metrics of designs in replicates."""
        designs, replicates, values = matrices(_results(), ["value"])
        self.assertEqual(list(designs), ["a", "b"])
        self.assertEqual(list(replicates), [0, 1])
        np.testing.assert_array_equal(values["value"], [[1.0, 3.0], [2.0, np.nan]])
        with self.assertRaises(ValueError):
            matrices(_results(), ["cost"])

    def test_values_of_designs(self):
        """The values of designs are the columns of their feasible records that
        do not vary."""
        feasible, infeasible = feasible_records(_results())
        self.assertEqual((len(feasible), len(infeasible)), (3, 1))
        statistics = {"value": Mean("value")}
        self.assertEqual(value_columns(_results(), statistics, (), None), ["size"])
        self.assertEqual(value_columns(_results(), statistics, ("size",), None), [])
        self.assertEqual(
            value_columns(_results(), {}, (), ["size", "value"]), ["size", "value"]
        )
        frame = design_values(_results(), ["b", "a"], ["size"])
        self.assertEqual(frame.to_dict("list"), {"design": ["b", "a"], "size": [2, 1]})

    def test_resamples(self):
        """Resamples draw replicates with replacement."""
        drawn = draw_resamples(5, 3, 0)
        self.assertEqual(drawn.shape, (3, 5))
        self.assertTrue(((drawn >= 0) & (drawn < 5)).all())
        np.testing.assert_array_equal(drawn, draw_resamples(5, 3, 0))
        with self.assertRaises(ValueError):
            draw_resamples(5, 1, 0)

    def test_efficient_keys(self):
        """The efficient designs of a summary, each once, in order."""
        summary = pd.DataFrame(
            {"design": ["b", "a", "b", "c"], "pareto": [True, False, True, True]}
        )
        self.assertEqual(efficient_keys(summary), ["b", "c"])
