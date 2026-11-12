"""
Unit tests for Pareto efficiency.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import itertools
import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Criteria,
    Requirement,
    hypervolume,
    hypervolume_progress,
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
        loose = Criteria(objectives, tolerance={"value": 0.05})
        self.assertEqual(list(pareto_front(frame, loose)), [True, True, True, False])
        self.assertEqual(list(pareto_ranks(frame, loose)), [0, 0, 0, 1])
        # a larger tolerance keeps d (worse than a by 0.5 in value)
        looser = Criteria(objectives, tolerance={"value": 0.6})
        self.assertTrue(pareto_front(frame, looser).all())

    def test_requirements_and_groups(self):
        """Rows that miss requirements are not efficient, and groups are separate."""
        frame = self.frame.assign(
            reliability=[1.0, 1.0, 0.2, 1.0, 1.0, np.nan], sensor=list("xxxxxy")
        )
        criteria = Criteria(
            self.objectives, requirements=[Requirement("reliability", ">=", 0.5)]
        )
        front = pareto_front(frame, criteria)
        # c (efficient by its objectives) and f (unknown) miss the
        # requirement, and e is efficient without c
        self.assertEqual(list(front[front].index), ["a", "b", "e"])
        ranks = pareto_ranks(frame, criteria)
        self.assertEqual(list(ranks[list("abde")]), [0, 0, 1, 0])
        # those that miss it rank after, in order of their violations
        self.assertEqual(list(ranks[["c", "f"]]), [2, 3])
        area = hypervolume(frame, criteria, {"cost": 5, "value": 0})
        self.assertAlmostEqual(area, 4 * 1 + 3 * 2 + 1 * 1)
        # f is alone in its group
        grouped = Criteria(self.objectives, by="sensor")
        front = pareto_front(frame, grouped)
        self.assertEqual(list(front[front].index), ["a", "b", "c", "f"])
        self.assertEqual(pareto_ranks(frame, grouped)["f"], 0)
        with self.assertRaises(ValueError):
            hypervolume(frame, grouped, [5, 0])
        with self.assertRaises(ValueError):
            pareto_front(self.frame, criteria)  # no column of the requirement


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
        front = pareto_front(frame, Criteria(objectives, tolerance={"value": "5%"}))
        self.assertEqual(list(front), [True, True, False, False])
        # an absolute tolerance of 1 does not tolerate 96 against 100
        absolute = pareto_front(frame, Criteria(objectives, tolerance={"value": 1.0}))
        self.assertEqual(list(absolute), [True, False, False, False])


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
            results, Criteria({"cost": "min", "value": "max"}, by=["sensor"])
        )
        winners = probability[probability.pareto_probability == 1]
        self.assertEqual(
            sorted(zip(winners.sensor, winners.design)), [("a", "d1"), ("b", "d0")]
        )


class TestManyPoints(unittest.TestCase):
    """
    Unit tests for the dominance of many points.
    """

    @staticmethod
    def _every_pair(costs, tolerances=None):
        """Dominance by comparing every pair of points."""
        absolute, relative = (
            np.zeros((2, costs.shape[1])) if tolerances is None else tolerances
        )
        dominated = np.zeros(len(costs), dtype=bool)
        for i, point in enumerate(costs):
            penalized = costs + absolute + relative * np.abs(point)
            no_worse = np.all(penalized <= point, axis=1)
            better = np.any(penalized < point, axis=1)
            dominated[i] = np.any(no_worse & better)
        return dominated

    def test_matches_every_pair(self):
        """Dominance, with or without tolerances, is that of every pair."""
        # pylint: disable=import-outside-toplevel,protected-access
        from tatc.tradespace import pareto

        rng = np.random.default_rng(3)
        sets = (
            rng.random((700, 3)),
            rng.integers(0, 4, (700, 3)).astype(float),
            np.vstack([rng.random((150, 2))] * 3),
            np.column_stack([rng.integers(1, 11, 900), rng.normal(0, 30, 900)]),
        )
        tolerances = (
            None,
            np.array([[0.0, 0.0, 0.1], [0.0, 0.0, 0.0]]),
            np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.1]]),
        )
        for costs in sets:
            for tolerance in tolerances:
                if tolerance is not None:
                    tolerance = tolerance[:, : costs.shape[1]]
                np.testing.assert_array_equal(
                    pareto._dominated(costs, tolerance),
                    self._every_pair(costs, tolerance),
                )

    def test_many_points(self):
        """Fronts of many points are found quickly."""
        rng = np.random.default_rng(4)
        frame = pd.DataFrame(rng.random((20_000, 3)), columns=["a", "b", "c"])
        front = pareto_front(frame, {"a": "min", "b": "min", "c": "max"})
        self.assertGreater(front.sum(), 0)
        tolerant = pareto_front(
            frame, Criteria({"a": "min", "b": "min", "c": "max"}, tolerance={"c": "5%"})
        )
        self.assertTrue((tolerant >= front).all())


class TestHypervolumeProgress(unittest.TestCase):
    """
    Unit tests for the hypervolume of searches as they proceed.
    """

    def test_progress(self):
        """The hypervolume of the first designs evaluated, in steps."""
        results = pd.DataFrame(
            {
                "design": ["a", "a", "b", "c", "d"],
                "replicate": [0, 1, 0, 0, 0],
                "feasible": [True, True, True, True, False],
                "x": [1.0, 3.0, 4.0, 1.0, np.nan],
                "y": [4.0, 6.0, 8.0, 2.0, np.nan],
            }
        )
        objectives = {"x": "min", "y": "max"}
        # a (2, 5) and b (4, 8) dominate 58 of the box to (10, 0), and c (1,
        # 2) adds 2; the infeasible design d adds nothing
        progress = hypervolume_progress(results, objectives, (10, 0), step=2)
        self.assertEqual(progress.index.tolist(), [2, 4])
        np.testing.assert_allclose(progress.values, [58.0, 60.0])
        last = hypervolume_progress(results, objectives, (10, 0), step=3)
        self.assertEqual(last.index.tolist(), [3, 4])
        np.testing.assert_allclose(last.values, [60.0, 60.0])
        with self.assertRaises(ValueError):
            hypervolume_progress(results, objectives, (10, 0), step=0)


class TestStandardErrors(unittest.TestCase):
    """
    Unit tests for the standard errors of the objectives of designs.
    """

    def test_standard_errors(self):
        """The bootstrap gives the standard error of each objective."""
        rng = np.random.default_rng(0)
        results = pd.DataFrame(
            [
                {
                    "design": design,
                    "replicate": r,
                    "feasible": True,
                    "x": float(i),
                    "y": i + rng.normal(0, 2),
                }
                for i, design in enumerate("ab")
                for r in range(40)
            ]
        )
        probability = pareto_probability(
            results, {"x": "min", "y": "max"}, n_resamples=2000
        )
        np.testing.assert_allclose(probability["x se"], 0.0)
        expected = results.groupby("design").y.std(ddof=0) / np.sqrt(40)
        np.testing.assert_allclose(
            probability.set_index("design")["y se"], expected, rtol=0.15
        )
