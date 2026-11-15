"""
Unit tests for criteria: requirements on statistics, and the objectives,
tolerances, and rules across groups by which designs are judged.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import dataclasses
import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Across,
    Criteria,
    Mean,
    Quantile,
    Requirement,
    pareto_front,
)


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

    def test_statements(self):
        """Requirements tell which values meet them, and state themselves."""
        requirement = Requirement("cost", "<=", 1.5)
        np.testing.assert_array_equal(
            requirement.met([1.0, 2.0, np.nan]), [True, False, False]
        )
        self.assertEqual(str(requirement), "cost <= 1.5")

    def test_invalid(self):
        """Operators and scales are checked."""
        with self.assertRaises(ValueError):
            Requirement("x", "!=", 1)
        with self.assertRaises(ValueError):
            Requirement("x", ">=", 1, scale=0)


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
        self.assertEqual(
            criteria.resolved, {"a": Mean("a"), "b": Mean("y"), "c": Mean("x")}
        )
        self.assertEqual(criteria.by, ("sensor",))

    def test_conventional_names(self):
        """Statistics given without names are named by convention."""
        criteria = Criteria(
            {"value p25": "max"}, statistics=[Quantile("value", 0.25), "cost"]
        )
        self.assertEqual(criteria.names, ["value p25", "cost"])

    def test_invalid(self):
        """Objectives, statistics, requirements, and tolerances are checked."""
        with self.assertRaises(ValueError):
            Criteria({})
        with self.assertRaises(ValueError):
            Criteria({"a": "lowest"})
        with self.assertRaises(ValueError):
            Criteria({"a": "min"}, tolerances={"b": 1})
        with self.assertRaises(ValueError):
            Criteria({"a": "min"}, tolerances={"a": "5"})
        with self.assertRaises(ValueError):
            Criteria({"a": "min"}, tolerances={"a": -1})
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
        criteria = Criteria(objectives)
        self.assertIs(Criteria.of(criteria), criteria)
        self.assertEqual(Criteria.of(objectives), criteria)
        with self.assertRaises(TypeError):
            Criteria.of(["cost"])

    def test_costs_and_violation(self):
        """Objectives are costs (maximized ones negated), and the violations
        of requirements add."""
        criteria = Criteria(
            {"cost": "min", "value": "max"},
            requirements=[Requirement("cost", "<=", 2), Requirement("value", ">=", 4)],
        )
        np.testing.assert_array_equal(criteria.signs, [1.0, -1.0])
        values = {"cost": np.array([1.0, 3.0]), "value": np.array([5.0, np.nan])}
        np.testing.assert_array_equal(
            criteria.costs(values), [[1.0, -5.0], [3.0, np.nan]]
        )
        np.testing.assert_allclose(criteria.violation(values, 2), [0.0, np.nan])
        np.testing.assert_allclose(
            criteria.violation({"cost": [3.0], "value": [2.0]}, 1), [0.5 + 0.5]
        )
        self.assertEqual(criteria.costs({"cost": [], "value": []}).shape, (0, 2))
        np.testing.assert_array_equal(Criteria({"a": "min"}).violation({}, 3), 0.0)

    def test_tolerance_array(self):
        """Tolerances are values in the units of objectives, or shares of
        them, in the order of the objectives (zero if not given)."""
        criteria = Criteria(
            {"a": "min", "b": "max", "c": "min"}, tolerances={"c": "5%", "a": 2}
        )
        np.testing.assert_array_equal(
            criteria.tolerance_array, [[2.0, 0.0, 0.0], [0.0, 0.0, 0.05]]
        )
        # derived from the tolerances, which compare criteria
        changed = dataclasses.replace(criteria, tolerances={"b": 1})
        np.testing.assert_array_equal(
            changed.tolerance_array, [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]
        )
        self.assertEqual(criteria, dataclasses.replace(criteria))
        self.assertNotEqual(criteria, changed)

    def test_combiners(self):
        """Statistics across groups are combined by their worst, or their
        mean."""
        self.assertEqual(Criteria({"a": "min"}).combiners, {})
        criteria = Criteria(
            {"a": "min", "b": "max"},
            statistics={"c": "x", "d": "y"},
            requirements=[Requirement("c", ">=", 1)],
            across=Across(rules={"b": "mean"}),
        )
        self.assertEqual(
            criteria.combiners, {"a": "max", "b": "mean", "c": "min", "d": "mean"}
        )
        with self.assertRaisesRegex(ValueError, "no worst"):
            Criteria(
                {"a": "min"}, statistics={"d": "y"}, across=Across(rules={"d": "worst"})
            )
        with self.assertRaisesRegex(ValueError, "disagree"):
            Criteria(
                {"a": "min"}, requirements=[Requirement("a", ">=", 1)], across=Across()
            )
