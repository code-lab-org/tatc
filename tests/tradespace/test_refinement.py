"""
Unit tests for adaptive replication.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
import warnings

import numpy as np

from tatc.tradespace import (
    DesignSpace,
    Integer,
    Quantile,
    Study,
    Uncertain,
    pareto_probability,
    refine,
)

OBJECTIVES = {"cost": "min", "value": "max"}


def _study(replicates=5):
    """Designs whose values are close to others' within noise for some costs."""
    space = DesignSpace(
        cost=Integer(1, 8),
        variant=Integer(0, 2),
        noise=Uncertain.normal(0, 1),
    )

    def evaluate(d):
        # variants 0 and 1 differ by a unit of noise; variant 2 is clearly worse
        value = np.sqrt(d.cost) * 3 + {0: 0.0, 1: 1.0, 2: -3.0}[d.variant]
        return {"cost": d.cost, "value": value + d.noise}

    # independent noise between designs, so that close designs are uncertain
    return Study(
        space, lambda d: d, evaluate, replicates=replicates, common_random_numbers=False
    )


class TestRefine(unittest.TestCase):
    """
    Unit tests for adaptive replication.
    """

    def setUp(self):
        self.study = _study()
        self.designs = self.study.space.enumerate()
        self.results = self.study.run(self.designs)

    def test_refines_uncertain_designs(self):
        """Replicates are added only to uncertain designs, which become fewer."""
        refinement = refine(
            self.study,
            self.designs,
            self.results,
            OBJECTIVES,
            step=5,
            max_rounds=6,
            n_resamples=200,
        )
        history = refinement.history
        self.assertGreater(history.refined.iloc[0], 0)
        replicates = refinement.results.groupby("design").replicate.nunique()
        # the clearly worse variant is never refined
        worse = [d.key for d in self.designs if d.variant == 2]
        self.assertTrue((replicates[worse] == 5).all())
        self.assertGreater(replicates.max(), 5)
        before = pareto_probability(self.results, OBJECTIVES, n_resamples=200)
        uncertain_before = (
            (before.pareto_probability > 0.05) & (before.pareto_probability < 0.95)
        ).sum()
        self.assertLess(len(refinement.uncertain), uncertain_before)
        # far fewer evaluations than uniform replication to the same depth
        uniform = replicates.max() * len(self.designs)
        self.assertLess(len(refinement.results), 0.6 * uniform)

    def test_replicates_continue(self):
        """Added replicates follow each design's own."""
        refinement = refine(
            self.study, self.designs, self.results, OBJECTIVES, step=2, max_rounds=1
        )
        added = refinement.results.iloc[len(self.results) :]
        self.assertTrue(set(added.replicate) <= {5, 6})
        self.assertFalse(refinement.results.duplicated(["design", "replicate"]).any())

    def test_limits(self):
        """Refinement stops at the most replicates of a design."""
        refinement = refine(
            self.study,
            self.designs,
            self.results,
            OBJECTIVES,
            step=3,
            max_replicates=8,
            max_rounds=5,
        )
        self.assertLessEqual(
            refinement.results.groupby("design").replicate.nunique().max(), 8
        )
        self.assertEqual(refinement.history.refined.iloc[-1], 0)
        with self.assertRaises(ValueError):
            refine(
                self.study, self.designs, self.results, OBJECTIVES, bounds=(0.9, 0.1)
            )

    def test_statistics(self):
        """Refinement uses the statistics of the objectives."""
        refinement = refine(
            self.study,
            self.designs,
            self.results,
            {"cost": "min", "value p25": "max"},
            statistics={"value p25": Quantile("value", 0.25)},
            max_rounds=1,
        )
        self.assertIn("value p25", refinement.probability.columns)

    def test_tolerance(self):
        """A tolerance settles designs that are within it of the front."""
        strict = pareto_probability(self.results, OBJECTIVES, n_resamples=300)
        tolerant = pareto_probability(
            self.results, OBJECTIVES, n_resamples=300, tolerance={"value": 2.0}
        )
        # designs efficient for sure, of the better variant (within 2 of each other)
        self.assertGreater(
            (tolerant.pareto_probability >= 0.95).sum(),
            (strict.pareto_probability >= 0.95).sum(),
        )
        refinement = refine(
            self.study,
            self.designs,
            self.results,
            OBJECTIVES,
            max_rounds=1,
            tolerance={"value": 2.0},
        )
        settled = refinement.probability.pareto_probability >= 0.95
        self.assertGreaterEqual(
            settled.sum(), (tolerant.pareto_probability >= 0.95).sum()
        )

    def test_min_replicates(self):
        """Designs with few replicates are refined whatever their probability."""
        space = DesignSpace(
            x=Integer(0, 9), a=Uncertain.normal(0, 1), b=Uncertain.normal(0, 1)
        )
        # equal designs, of which two replicates find some surely inefficient
        study = Study(
            space,
            lambda d: d,
            lambda d: {"a": d.a, "b": d.b},
            replicates=2,
            common_random_numbers=False,
        )
        designs = space.enumerate()
        results = study.run(designs)
        objectives = {"a": "min", "b": "min"}
        before = pareto_probability(results, objectives, n_resamples=200)
        self.assertTrue((before.pareto_probability == 0).any())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            refinement = refine(
                study, designs, results, objectives, step=3, max_rounds=1
            )
        replicates = refinement.results.groupby("design").replicate.nunique()
        self.assertTrue((replicates >= 5).all())
        with self.assertRaises(ValueError):
            refine(study, designs, results, objectives, max_replicates=3)

    def test_partial_step(self):
        """The last step of a design adds the replicates left to its limit."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            refinement = refine(
                self.study,
                self.designs,
                self.results,
                OBJECTIVES,
                step=4,
                max_replicates=7,
                max_rounds=3,
            )
        replicates = refinement.results.groupby("design").replicate.nunique()
        self.assertEqual(replicates.max(), 7)
        self.assertTrue(set(replicates) <= {5, 7})

    def test_warns_without_tolerance(self):
        """Designs that remain uncertain without a tolerance give a warning."""
        with self.assertWarns(UserWarning):
            refine(self.study, self.designs, self.results, OBJECTIVES, max_rounds=1)
