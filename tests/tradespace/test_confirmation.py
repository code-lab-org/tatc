"""
Unit tests for the confirmation of selected designs.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np

from tatc.tradespace import (
    DesignSpace,
    Integer,
    Probability,
    Requirement,
    Study,
    Uncertain,
    confirm,
    pareto_front,
    summarize,
)

OBJECTIVES = {"cost": "min", "value": "max"}


def _study(replicates=3):
    """Designs of equal true values (but for the cost), with independent noise."""
    space = DesignSpace(
        cost=Integer(1, 2), variant=Integer(0, 49), noise=Uncertain.normal(0, 1)
    )
    return Study(
        space,
        lambda d: d,
        lambda d: {"cost": d.cost, "value": d.cost + d.noise},
        replicates=replicates,
        common_random_numbers=False,
    )


class TestConfirm(unittest.TestCase):
    """
    Unit tests for confirmation.
    """

    def setUp(self):
        self.study = _study()
        self.designs = self.study.space.enumerate()
        self.results = self.study.run(self.designs)
        summary = summarize(self.results, list(OBJECTIVES), design_columns=[])
        # the best looking designs of each cost
        best = summary.sort_values("value", ascending=False).groupby(
            summary.design.map({d.key: d.cost for d in self.designs})
        )
        self.selected_keys = set(best.head(5).design)
        self.selected = [d for d in self.designs if d.key in self.selected_keys]

    def test_removes_optimism(self):
        """New replicates remove the optimism of selection."""
        confirmation = confirm(
            self.study, self.selected, self.results, OBJECTIVES, replicates=20
        )
        self.assertEqual(set(confirmation.summary.design), self.selected_keys)
        self.assertTrue((confirmation.summary.replicates == 20).all())
        # the new replicates follow every replicate of the results
        self.assertEqual(set(confirmation.results.replicate), set(range(3, 23)))
        optimism = confirmation.optimism(OBJECTIVES)
        self.assertAlmostEqual(optimism["cost"], 0)
        # the 5 best of 50 designs of 3 replicates look better by about a
        # standard deviation of their means (0.58), which the new ones remove
        self.assertGreater(optimism["value"], 0.6)
        # the confirmed values are unbiased (their true values are the costs)
        errors = confirmation.summary.value - confirmation.summary.design.map(
            {d.key: d.cost for d in self.designs}
        )
        self.assertLess(abs(errors.mean()), 0.15)

    def test_front(self):
        """The front is found among the confirmed designs by the new replicates."""
        confirmation = confirm(self.study, self.selected, self.results, OBJECTIVES)
        summary = confirmation.summary
        np.testing.assert_array_equal(
            summary.pareto, pareto_front(summary, OBJECTIVES).to_numpy()
        )

    def test_requirements(self):
        """Requirements are judged by the new replicates."""
        confirmation = confirm(
            self.study,
            self.selected,
            self.results,
            OBJECTIVES,
            statistics={"reliability": Probability("value", ">=", 1.5)},
            requirements=[Requirement("reliability", ">=", 0.5)],
        )
        summary = confirmation.summary
        self.assertIn("reliability selected", summary.columns)
        np.testing.assert_array_equal(
            summary.meets_requirements, summary.reliability >= 0.5
        )
        self.assertFalse((summary.pareto & ~summary.meets_requirements).any())
        with self.assertRaises(ValueError):
            confirm(self.study, self.selected, self.results, OBJECTIVES, replicates=0)
