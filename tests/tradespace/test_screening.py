"""
Unit tests for screening with analyses of increasing fidelity.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
import warnings

import numpy as np

from tatc.tradespace import (
    Criteria,
    DesignSpace,
    Integer,
    Requirement,
    Study,
    Uncertain,
    confirm,
    pareto_front,
    refine,
    screen,
    summarize,
)

OBJECTIVES = {"cost": "min", "value": "max"}


def _space():
    return DesignSpace(
        a=Integer(1, 10),
        b=Integer(1, 10),
        noise=Uncertain.normal(0, 1),
        constraints=[lambda d: d.a + d.b <= 16],
    )


def _metrics(design, noise_scale):
    cost = design.a + design.b
    value = np.sqrt(design.a) * np.log1p(design.b) + noise_scale * design.noise
    return {"cost": cost, "value": value}


def _study(noise_scale, replicates):
    # independent noise of each design, as of a cheap analysis
    return Study(
        _space(),
        lambda d: _metrics(d, noise_scale),
        lambda metrics: metrics,
        replicates=replicates,
        common_random_numbers=False,
    )


class TestScreen(unittest.TestCase):
    """
    Unit tests for screening.
    """

    def setUp(self):
        self.designs = _space().enumerate()
        exact = summarize(_study(0, 1).run(self.designs), statistics=list(OBJECTIVES))
        self.front = set(exact.design[pareto_front(exact, OBJECTIVES)])

    def test_keeps_whole_fronts(self):
        """Each stage keeps whole fronts until it keeps enough designs."""
        screening = screen(
            [_study(0.3, 1), _study(0, 1)], OBJECTIVES, self.designs, keep=20
        )
        first = screening.summaries[screening.summaries.stage == 0]
        self.assertGreaterEqual(first.kept.sum(), 20)
        last_rank = first[first.kept].pareto_rank.max()
        self.assertTrue((first[first.pareto_rank <= last_rank].kept).all())
        self.assertFalse((first[first.pareto_rank > last_rank].kept).any())
        self.assertEqual(len(screening.designs), first.kept.sum())
        self.assertEqual(set(screening.history.stage), {0, 1})
        self.assertEqual(len(screening.summary), len(screening.designs))
        # the results of the last stage, which other steps continue
        last = screening.history[screening.history.stage == 1]
        self.assertEqual(len(screening.results), len(last))
        self.assertNotIn("stage", screening.results.columns)
        self.assertEqual(screening.criteria.objectives, OBJECTIVES)
        # the summaries of the stages are those of summarize, with spreads
        noisy = screen(
            [_study(0.3, 2), _study(0.3, 2)], OBJECTIVES, self.designs, keep=20
        )
        for column in ("value std", "value se"):
            self.assertIn(column, noisy.summary.columns)
            self.assertIn(column, noisy.summaries.columns)
        self.assertNotIn("value se", screening.summary.columns)  # without noise

    def test_finds_front(self):
        """A tolerance of the noise of the first stage keeps the efficient designs."""
        stages = [_study(0.3, 1), _study(0, 1)]
        # about two standard deviations of the noise of the first stage
        cheap = Criteria(OBJECTIVES, tolerances={"value": 0.6})
        screening = screen(stages, [cheap, OBJECTIVES], self.designs, keep=0.1)
        final = screening.summary
        found = set(final.design[pareto_front(final, OBJECTIVES)])
        self.assertEqual(found, self.front)
        self.assertLess(len(screening.designs), 0.6 * len(self.designs))
        # without a tolerance, the noise loses efficient designs
        screening = screen(stages, OBJECTIVES, self.designs, keep=0.1)
        final = screening.summary
        found = set(final.design[pareto_front(final, OBJECTIVES)])
        self.assertLess(len(found & self.front), len(self.front))

    def test_three_stages(self):
        """Values of keep and criteria can be given for each stage."""
        screening = screen(
            [_study(0.5, 1), _study(0.1, 2), _study(0, 1)],
            [Criteria(OBJECTIVES, tolerances={"value": 1.0}), OBJECTIVES, OBJECTIVES],
            self.designs,
            keep=[0.5, 10],
        )
        counts = screening.summaries.groupby("stage").size().tolist()
        self.assertEqual(counts[0], len(self.designs))
        self.assertGreaterEqual(counts[1], 0.5 * counts[0])
        self.assertGreaterEqual(counts[2], 10)
        with self.assertRaises(ValueError):
            screen([_study(0, 1)] * 3, OBJECTIVES, self.designs, keep=[0.5])
        with self.assertRaises(ValueError):
            screen([_study(0, 1)] * 3, [OBJECTIVES] * 2, self.designs)
        with self.assertRaises(ValueError):
            screen([_study(0, 1)], Criteria(OBJECTIVES, by="sensor"), self.designs)

    def test_next_steps(self):
        """The results and designs of a screening continue in other steps."""
        full = _study(0, 2)
        screening = screen([_study(0.3, 1), full], OBJECTIVES, self.designs, keep=10)
        confirmation = confirm(full, OBJECTIVES, screening.results, screening.designs)
        self.assertEqual(
            {d.key for d in confirmation.designs}, {d.key for d in screening.designs}
        )
        confirm(full, OBJECTIVES, screening.results)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # designs still uncertain
            refinement = refine(
                full,
                OBJECTIVES,
                screening.results,
                screening.designs,
                max_rounds=1,
                n_resamples=20,
            )
        self.assertEqual(len(refinement.candidates), len(screening.designs))

    def test_keep(self):
        """Designs are kept by a number of at least 1, or a share of at most 1."""
        stages = [_study(0, 1), _study(0, 1)]
        designs = self.designs[:20]
        everything = screen(stages, OBJECTIVES, designs, keep=1.0)
        self.assertEqual(len(everything.designs), len(designs))
        for keep in (0, -1, 0.0, 1.5, -0.5, [0.5, 0]):
            with self.subTest(keep=keep):
                stages_of = (
                    stages if not isinstance(keep, list) else stages + stages[:1]
                )
                with self.assertRaises(ValueError):
                    screen(stages_of, OBJECTIVES, designs, keep=keep)
        for keep in (True, "10", None):
            with self.subTest(keep=keep):
                with self.assertRaises(TypeError):
                    screen(stages, OBJECTIVES, designs, keep=keep)

    def test_criteria_tolerance(self):
        """Every stage uses the tolerance of the criteria unless given its own."""
        stages = [_study(0.3, 1), _study(0, 1)]
        criteria = Criteria(OBJECTIVES, tolerances={"value": 0.6})
        tolerant = screen(stages, criteria, self.designs, keep=0.1)
        overridden = screen(stages, [OBJECTIVES, criteria], self.designs, keep=0.1)
        self.assertGreater(len(tolerant.designs), len(overridden.designs))
        # the last stage ranks with the criteria's tolerance
        final = tolerant.summary
        self.assertTrue(
            (final.pareto_rank == 0).sum() >= pareto_front(final, OBJECTIVES).sum()
        )


class TestScreenRequirements(unittest.TestCase):
    """
    Unit tests for screening with requirements.
    """

    def test_misses_ranked_last(self):
        """Designs that miss requirements are ranked after those that meet them."""
        space = DesignSpace(a=Integer(1, 6), noise=Uncertain.uniform(0, 1))
        study = Study(
            space,
            lambda d: d,
            lambda d: {"cost": d.a, "value": d.a**0.5 + d.noise},
            replicates=2,
        )
        screening = screen(
            [study, study],
            Criteria(
                {"cost": "min", "value": "max"},
                requirements=[Requirement("cost", "<=", 3)],
            ),
            space.enumerate(),
            keep=4,
        )
        summaries = screening.summaries
        first = summaries[summaries.stage == 0].set_index("design")
        meets = first[first.meets_requirements]
        misses = first[~first.meets_requirements]
        self.assertGreater(misses.pareto_rank.min(), meets.pareto_rank.max())
        # the nearest miss (cost 4) is kept to fill four designs
        self.assertEqual(sorted(first[first.kept].cost), [1, 2, 3, 4])
