"""
Unit tests for adaptive replication.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
import warnings

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Across,
    Criteria,
    DesignSpace,
    Fixed,
    Infeasible,
    Integer,
    Quantile,
    Study,
    Uncertain,
    pareto_probability,
    refine,
)
from tatc.tradespace.refinement import _precise

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

    def test_designs_of_results(self):
        """Without designs, every design of the results may be refined; designs
        may also be given by their keys."""
        given = refine(self.study, OBJECTIVES, self.results, self.designs)
        rebuilt = refine(self.study, OBJECTIVES, self.results)
        pd.testing.assert_frame_equal(
            given.results.drop(columns="seconds"),
            rebuilt.results.drop(columns="seconds"),
        )
        self.assertEqual(
            [d.key for d in rebuilt.candidates], [d.key for d in self.designs]
        )
        # the designs selected for the next step: the efficient designs
        summary = rebuilt.summary
        self.assertEqual(
            [d.key for d in rebuilt.designs], list(summary.design[summary.pareto])
        )
        self.assertEqual(rebuilt.criteria, given.criteria)
        # the evaluations added, over every round
        self.assertEqual(rebuilt.evaluations, len(rebuilt.results) - len(self.results))
        self.assertGreater(rebuilt.evaluations, 0)
        keys = [d.key for d in self.designs[:4]]
        some = refine(self.study, OBJECTIVES, self.results, keys)
        added = some.results.iloc[len(self.results) :]
        self.assertTrue(set(added.design) <= set(keys))
        self.assertIn("pareto_probability", some.summary.columns)
        self.assertIn("variant", some.summary.columns)  # the values of designs

    def test_refines_uncertain_designs(self):
        """Replicates are added only to uncertain designs, which become fewer."""
        refinement = refine(
            self.study,
            OBJECTIVES,
            self.results,
            self.designs,
            step=5,
            max_rounds=6,
            n_resamples=200,
        )
        rounds = refinement.rounds
        self.assertGreater(rounds.refined.iloc[0], 0)
        replicates = refinement.results.groupby("design").replicate.nunique()
        # the clearly worse variant is never refined
        worse = [d.key for d in self.designs if d.variant == 2]
        self.assertTrue((replicates[worse] == 5).all())
        self.assertGreater(replicates.max(), 5)
        before = pareto_probability(self.results, OBJECTIVES, n_resamples=200)
        uncertain_before = (
            (before.pareto_probability > 0.05) & (before.pareto_probability < 0.95)
        ).sum()
        summary = refinement.summary
        unsettled = summary[summary.uncertain & ~summary.indifferent]
        self.assertLess(len(unsettled), uncertain_before)
        # far fewer evaluations than uniform replication to the same depth
        uniform = replicates.max() * len(self.designs)
        self.assertLess(len(refinement.results), 0.6 * uniform)

    def test_replicates_continue(self):
        """Added replicates follow each design's own."""
        refinement = refine(
            self.study, OBJECTIVES, self.results, self.designs, step=2, max_rounds=1
        )
        added = refinement.results.iloc[len(self.results) :]
        self.assertTrue(set(added.replicate) <= {5, 6})
        self.assertFalse(refinement.results.duplicated(["design", "replicate"]).any())

    def test_limits(self):
        """Refinement stops at the most replicates of a design."""
        refinement = refine(
            self.study,
            OBJECTIVES,
            self.results,
            self.designs,
            step=3,
            max_replicates=8,
            max_rounds=5,
        )
        self.assertLessEqual(
            refinement.results.groupby("design").replicate.nunique().max(), 8
        )
        self.assertEqual(refinement.rounds.refined.iloc[-1], 0)
        with self.assertRaises(ValueError):
            refine(
                self.study, OBJECTIVES, self.results, self.designs, bounds=(0.9, 0.1)
            )

    def test_statistics(self):
        """Refinement uses the statistics of the objectives."""
        refinement = refine(
            self.study,
            Criteria(
                {"cost": "min", "value p25": "max"},
                statistics={"value p25": Quantile("value", 0.25)},
            ),
            self.results,
            self.designs,
            max_rounds=1,
        )
        self.assertIn("value p25", refinement.summary.columns)

    def test_tolerance(self):
        """A tolerance settles designs that are within it of the front."""
        strict = pareto_probability(self.results, OBJECTIVES, n_resamples=300)
        criteria = Criteria(OBJECTIVES, tolerances={"value": 2.0})
        tolerant = pareto_probability(self.results, criteria, n_resamples=300)
        # designs efficient for sure, of the better variant (within 2 of each other)
        self.assertGreater(
            (tolerant.pareto_probability >= 0.95).sum(),
            (strict.pareto_probability >= 0.95).sum(),
        )
        refinement = refine(
            self.study, criteria, self.results, self.designs, max_rounds=1
        )
        settled = refinement.summary.pareto_probability >= 0.95
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
                study, objectives, results, designs, step=3, max_rounds=1
            )
        replicates = refinement.results.groupby("design").replicate.nunique()
        self.assertTrue((replicates >= 5).all())
        with self.assertRaises(ValueError):
            refine(study, objectives, results, designs, max_replicates=3)

    def test_partial_step(self):
        """The last step of a design adds the replicates left to its limit."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            refinement = refine(
                self.study,
                OBJECTIVES,
                self.results,
                self.designs,
                step=4,
                max_replicates=7,
                max_rounds=3,
            )
        replicates = refinement.results.groupby("design").replicate.nunique()
        self.assertEqual(replicates.max(), 7)
        self.assertTrue(set(replicates) <= {5, 7})

    def test_warns_without_tolerance(self):
        """Designs that remain uncertain without a tolerance give a warning,
        at the call of the refinement."""
        with self.assertWarns(UserWarning) as context:
            refine(self.study, OBJECTIVES, self.results, self.designs, max_rounds=1)
        self.assertEqual(context.filename, __file__)

    def test_precision(self):
        """A design about a tolerance from the front stops once precise enough."""
        space = DesignSpace(variant=Integer(0, 1), noise=Uncertain.normal(0, 1))
        # variant 1 is one tolerance below variant 0, with noise: whether it is
        # within the tolerance of the front stays uncertain (with this seed)
        study = Study(
            space,
            lambda d: d,
            lambda d: {"cost": 1, "value": 10.0 - d.variant * (2.0 - d.noise)},
            replicates=5,
            seed=2,
        )
        designs = space.enumerate()
        results = study.run(designs)
        criteria = Criteria(OBJECTIVES, tolerances={"value": 2.0})

        def refined(precision):
            refinement = refine(
                study,
                criteria,
                results,
                designs,
                step=5,
                max_replicates=40,
                precision=precision,
            )
            replicates = refinement.results.groupby("design").replicate.nunique()
            return refinement, replicates[designs[1].key]

        _, unlimited = refined(None)
        self.assertEqual(unlimited, 40)
        # a standard error below a tenth of the tolerance (0.2) takes about 25
        refinement, precise = refined(0.1)
        self.assertGreaterEqual(precise, 20)
        self.assertLessEqual(precise, 30)
        summary = refinement.summary
        self.assertEqual(summary.design[summary.indifferent].tolist(), [designs[1].key])
        self.assertFalse((summary.uncertain & ~summary.indifferent).any())
        self.assertEqual(refinement.rounds.indifferent.iloc[-1], 1)
        with self.assertRaises(ValueError):
            refine(study, criteria, results, designs, precision=0)

    def test_budget(self):
        """A budget limits the evaluations, refining the least settled designs first."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            refinement = refine(
                self.study, OBJECTIVES, self.results, self.designs, step=5, budget=12
            )
            nothing = refine(
                self.study, OBJECTIVES, self.results, self.designs, budget=0
            )
        self.assertEqual(refinement.rounds.evaluations.sum(), 12)
        replicates = refinement.results.groupby("design").replicate.nunique()
        self.assertEqual(sorted(replicates[replicates > 5]), [7, 10, 10])
        # the designs whose probability was nearest 0.5
        before = pareto_probability(self.results, OBJECTIVES, n_resamples=500, seed=0)
        before = before.sort_values("design")
        distance = (before.pareto_probability - 0.5).abs()
        nearest = before.design[distance.sort_values(kind="stable").index[:3]]
        self.assertEqual(set(replicates[replicates > 5].index), set(nearest))
        self.assertEqual(len(nothing.results), len(self.results))

    def test_budget_of_scenarios(self):
        """A replicate of a study with scenarios spends an evaluation of the
        budget in each scenario, as the refinement counts them."""
        space = DesignSpace(
            cost=Integer(1, 6),
            variant=Integer(0, 2),
            shift=Fixed(0.0),
            noise=Uncertain.normal(0, 1),
        )

        def evaluate(d):
            value = np.sqrt(d.cost) * 3 + {0: 0.0, 1: 1.0, 2: -3.0}[d.variant]
            return {"cost": d.cost, "value": value + d.shift + d.noise}

        study = Study(
            space,
            lambda d: d,
            evaluate,
            replicates=3,
            common_random_numbers=False,
            scenarios={"a": {}, "b": {"shift": 0.5}, "c": {"shift": -0.5}},
        )
        results = study.run(space.enumerate())
        criteria = Criteria(OBJECTIVES, across=Across(rule="mean"))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            refinement = refine(
                study, criteria, results, budget=13, step=2, max_rounds=3
            )
        added = len(refinement.results) - len(results)
        self.assertEqual(refinement.evaluations, added)
        # four replicates of the budget's 13 evaluations, in each of 3 scenarios
        self.assertEqual(added, 12)

    def test_precise(self):
        """Designs are precise if each varying objective is within its tolerance."""
        probability = pd.DataFrame(
            {
                "cost": [1.0, 1.0, 1.0],
                "cost se": [0.0, 0.0, 0.1],
                "value": [100.0, 100.0, 100.0],
                "value se": [1.0, 6.0, 1.0],
            }
        )
        # a tolerance of 10% of 100: precise below a standard error of 5
        criteria = Criteria(OBJECTIVES, tolerances={"value": "10%"})
        np.testing.assert_array_equal(
            _precise(probability, criteria, 0.5), [True, False, False]
        )
        np.testing.assert_array_equal(
            _precise(probability, criteria, None), [False, False, False]
        )
        np.testing.assert_array_equal(
            _precise(probability, Criteria(OBJECTIVES), 0.5), [False, False, False]
        )
        # an objective without a standard error does not vary
        np.testing.assert_array_equal(
            _precise(probability.drop(columns="cost se"), criteria, 0.5),
            [True, False, True],
        )


class TestRefineInfeasible(unittest.TestCase):
    """
    Unit tests for the refinement of results without a feasible design.
    """

    def test_no_feasible_design(self):
        """Nothing is refined, and nothing is selected."""

        def reject(design):
            raise Infeasible("never")

        study = Study(_study().space, reject, lambda built: built, replicates=2)
        with self.assertWarnsRegex(UserWarning, "feasible"):
            results = study.run(study.space.enumerate()[:4])
        refinement = refine(study, OBJECTIVES, results, max_rounds=2, n_resamples=10)
        self.assertTrue(refinement.summary.empty)
        self.assertEqual(refinement.designs, [])
        self.assertEqual(len(refinement.candidates), 4)
        self.assertEqual(refinement.evaluations, 0)
        self.assertEqual(len(refinement.results), len(results))
