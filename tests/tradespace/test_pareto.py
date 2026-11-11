"""
Unit tests for Pareto efficiency.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import itertools
import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    hypervolume,
    pareto_front,
    pareto_probability,
    pareto_ranks,
)


def _brute_hypervolume(costs: np.ndarray, reference: np.ndarray) -> float:
    """Hypervolume (to minimize) by inclusion-exclusion over subsets."""
    total = 0.0
    for size in range(1, len(costs) + 1):
        for subset in itertools.combinations(range(len(costs)), size):
            corner = costs[list(subset)].max(axis=0)
            total += (-1) ** (size + 1) * np.prod(np.clip(reference - corner, 0, None))
    return total


class TestPareto(unittest.TestCase):
    """
    Unit tests for Pareto efficiency.
    """

    def setUp(self):
        self.frame = pd.DataFrame(
            {
                "cost": [1, 2, 3, 2, 4, 1],
                "value": [1, 3, 4, 2, 4, 1],
            },
            index=list("abcdef"),
        )
        self.objectives = {"cost": "min", "value": "max"}

    def test_front(self):
        """The front is the non-dominated rows, including ties."""
        front = pareto_front(self.frame, self.objectives)
        self.assertEqual(list(front[front].index), ["a", "b", "c", "f"])
        flipped = pareto_front(self.frame, {"cost": "max", "value": "max"})
        self.assertEqual(list(flipped[flipped].index), ["e"])

    def test_ranks(self):
        """Ranks peel successive fronts."""
        ranks = pareto_ranks(self.frame, self.objectives)
        self.assertEqual(list(ranks), [0, 0, 0, 1, 1, 0])

    def test_hypervolume_2d(self):
        """The two-dimensional hypervolume is the area of the dominated staircase."""
        # costs (cost, -value) to (5, 0): a (1,-1), b (2,-3), c (3,-4)
        area = hypervolume(self.frame, self.objectives, {"cost": 5, "value": 0})
        self.assertAlmostEqual(area, 4 * 1 + 3 * 2 + 2 * 1)
        self.assertEqual(hypervolume(self.frame, self.objectives, [0, 10]), 0.0)
        self.assertEqual(hypervolume(self.frame.iloc[:0], self.objectives, [5, 0]), 0)

    def test_hypervolume_matches_brute_force(self):
        """Hypervolumes in three and four dimensions match inclusion-exclusion."""
        rng = np.random.default_rng(4)
        for dimensions in (3, 4):
            costs = rng.random((9, dimensions))
            frame = pd.DataFrame(costs, columns=[f"f{i}" for i in range(dimensions)])
            objectives = {column: "min" for column in frame.columns}
            reference = np.full(dimensions, 1.1)
            self.assertAlmostEqual(
                hypervolume(frame, objectives, reference),
                _brute_hypervolume(costs, reference),
            )

    def test_invalid(self):
        """Directions and references are checked."""
        with self.assertRaises(ValueError):
            pareto_front(self.frame, {"cost": "lowest"})
        with self.assertRaises(ValueError):
            hypervolume(self.frame, self.objectives, [1])

    def test_tolerance(self):
        """Rows within the tolerances of the front are efficient."""
        frame = pd.DataFrame(
            {"cost": [1.0, 1.05, 1.5, 2.0], "value": [1.0, 1.02, 0.99, 0.5]}
        )
        objectives = {"cost": "min", "value": "max"}
        self.assertEqual(
            list(pareto_front(frame, objectives)), [True, True, False, False]
        )
        # c is worse than a and b by 0.01 and 0.03 in value, within the
        # tolerance; d is not
        loose = pareto_front(frame, objectives, tolerance={"value": 0.05})
        self.assertEqual(list(loose), [True, True, True, False])
        ranks = pareto_ranks(frame, objectives, tolerance={"value": 0.05})
        self.assertEqual(list(ranks), [0, 0, 0, 1])
        # a larger tolerance keeps d (worse than a by 0.5 in value)
        self.assertTrue(pareto_front(frame, objectives, tolerance={"value": 0.6}).all())
        with self.assertRaises(ValueError):
            pareto_front(frame, objectives, tolerance={"other": 1})


class TestRelativeTolerance(unittest.TestCase):
    """
    Unit tests for tolerances proportional to objectives.
    """

    def test_relative(self):
        """A percentage tolerates differences in proportion to the values."""
        frame = pd.DataFrame({"cost": [1, 1, 2, 2], "value": [100.0, 96.0, 4.0, 3.0]})
        objectives = {"cost": "min", "value": "max"}
        # 5% of 96 (4.8) tolerates 96 against 100; 5% of 3 (0.15) does not
        # tolerate 3 against 4
        front = pareto_front(frame, objectives, tolerance={"value": "5%"})
        self.assertEqual(list(front), [True, True, False, False])
        # an absolute tolerance of 1 does not tolerate 96 against 100
        absolute = pareto_front(frame, objectives, tolerance={"value": 1.0})
        self.assertEqual(list(absolute), [True, False, False, False])
        with self.assertRaises(ValueError):
            pareto_front(frame, objectives, tolerance={"value": "5"})


class TestParetoProbability(unittest.TestCase):
    """
    Unit tests for the probability of Pareto efficiency.
    """

    @staticmethod
    def _results(values):
        """Results of designs (rows) in replicates (columns), with a cost of 1."""
        return pd.DataFrame(
            [
                {
                    "design": f"d{i}",
                    "replicate": r,
                    "cost": 1.0,
                    "value": v,
                    "feasible": True,
                }
                for i, row in enumerate(values)
                for r, v in enumerate(row)
            ]
        )

    def test_settled(self):
        """Designs separated in every replicate are efficient or not for sure."""
        results = self._results([[1, 2, 3], [5, 6, 7]])
        probability = pareto_probability(results, {"cost": "min", "value": "max"})
        self.assertEqual(list(probability.pareto_probability), [0, 1])
        self.assertEqual(list(probability.pareto), [False, True])
        self.assertEqual(list(probability.value), [2, 6])

    def test_common_random_numbers(self):
        """Designs that differ by a constant in every replicate are settled."""
        rng = np.random.default_rng(0)
        common = rng.normal(0, 5, 8)
        results = self._results([common, common + 0.1])
        probability = pareto_probability(results, {"cost": "min", "value": "max"})
        self.assertEqual(list(probability.pareto_probability), [0, 1])

    def test_unsettled(self):
        """Designs whose order changes between replicates are uncertain."""
        # equal means: d0 varies between replicates, d1 does not
        results = self._results([[0, 2] * 4, [1] * 8])
        probability = pareto_probability(
            results, {"cost": "min", "value": "max"}, n_resamples=300
        )
        p = probability.pareto_probability
        self.assertTrue(((p > 0.05) & (p < 0.95)).all())

    def test_groups(self):
        """Fronts of groups are found separately."""
        results = pd.concat(
            [
                self._results([[1, 2], [3, 4]]).assign(sensor="a"),
                self._results([[3, 4], [1, 2]]).assign(sensor="b"),
            ]
        )
        probability = pareto_probability(
            results, {"cost": "min", "value": "max"}, by=["sensor"]
        )
        winners = probability[probability.pareto_probability == 1]
        self.assertEqual(
            sorted(zip(winners.sensor, winners.design)), [("a", "d1"), ("b", "d0")]
        )


class TestManyPoints(unittest.TestCase):
    """
    Unit tests for the dominance of many points.
    """

    def test_sweep(self):
        """A sweep finds the same dominated points as comparisons of every pair."""
        # pylint: disable=import-outside-toplevel,protected-access
        from tatc.tradespace import pareto

        rng = np.random.default_rng(3)
        for costs in (
            rng.random((500, 3)),
            rng.integers(0, 4, (500, 3)).astype(float),
            np.vstack([rng.random((100, 2))] * 3),
        ):
            np.testing.assert_array_equal(
                pareto._dominated_sweep(costs), pareto._dominated(costs)
            )
        many = rng.random((pareto._SWEEP_POINTS + 10, 3))
        frame = pd.DataFrame(many, columns=["a", "b", "c"])
        front = pareto_front(frame, {"a": "min", "b": "min", "c": "min"})
        np.testing.assert_array_equal(~front.to_numpy(), pareto._dominated_sweep(many))
