"""
Unit tests for screening with analyses of increasing fidelity.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np

from tatc.tradespace import (
    DesignSpace,
    Integer,
    Study,
    Uncertain,
    pareto_front,
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
        exact = summarize(_study(0, 1).run(self.designs), list(OBJECTIVES))
        self.front = set(exact.design[pareto_front(exact, OBJECTIVES)])

    def test_keeps_whole_fronts(self):
        """Each stage keeps whole fronts until it keeps enough designs."""
        screening = screen(
            self.designs, [_study(0.3, 1), _study(0, 1)], OBJECTIVES, keep=20
        )
        first = screening.summaries[0]
        self.assertGreaterEqual(first.kept.sum(), 20)
        last_rank = first[first.kept].pareto_rank.max()
        self.assertTrue((first[first.pareto_rank <= last_rank].kept).all())
        self.assertFalse((first[first.pareto_rank > last_rank].kept).any())
        self.assertEqual(len(screening.designs), first.kept.sum())
        self.assertEqual(set(screening.results.stage), {0, 1})
        self.assertEqual(len(screening.final), len(screening.designs))

    def test_finds_front(self):
        """A tolerance of the noise of the first stage keeps the efficient designs."""
        stages = [_study(0.3, 1), _study(0, 1)]
        # about two standard deviations of the noise of the first stage
        screening = screen(
            self.designs, stages, OBJECTIVES, keep=0.1, tolerance={"value": 0.6}
        )
        final = screening.final
        found = set(final.design[pareto_front(final, OBJECTIVES)])
        self.assertEqual(found, self.front)
        self.assertLess(len(screening.designs), 0.6 * len(self.designs))
        # without a tolerance, the noise loses efficient designs
        screening = screen(self.designs, stages, OBJECTIVES, keep=0.1)
        final = screening.final
        found = set(final.design[pareto_front(final, OBJECTIVES)])
        self.assertLess(len(found & self.front), len(self.front))

    def test_three_stages(self):
        """Values of keep and tolerance can be given for each transition."""
        screening = screen(
            self.designs,
            [_study(0.5, 1), _study(0.1, 2), _study(0, 1)],
            OBJECTIVES,
            keep=[0.5, 10],
            tolerance=[{"value": 1.0}, None],
        )
        counts = [len(summary) for summary in screening.summaries]
        self.assertEqual(counts[0], len(self.designs))
        self.assertGreaterEqual(counts[1], 0.5 * counts[0])
        self.assertGreaterEqual(counts[2], 10)
        with self.assertRaises(ValueError):
            screen(self.designs, [_study(0, 1)] * 3, OBJECTIVES, keep=[0.5])
