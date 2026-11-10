"""
Unit tests for trade studies.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from tatc.schemas import CircularOrbit, Satellite, WalkerConstellation
from tatc.tradespace import DesignSpace, Integer, Real, Study, Uncertain, summarize


class TestStudy(unittest.TestCase):
    """
    Unit tests for trade studies.
    """

    def setUp(self):
        self.space = DesignSpace(
            planes=Integer(1, 3),
            per_plane=Integer(1, 3),
            noise=Uncertain.uniform(0, 1),
        )

    def test_run_records(self):
        """Every design is evaluated in every replicate, with its values and metrics."""
        study = Study(
            self.space,
            build=lambda d: (d.planes, d.per_plane, d.noise),
            evaluate=lambda built: {
                "satellites": built[0] * built[1],
                "score": built[0] + built[2],
            },
            replicates=4,
        )
        results = study.run()
        self.assertEqual(len(results), 9 * 4)
        self.assertTrue(results.feasible.all())
        self.assertTrue(
            (results.satellites == results.planes * results.per_plane).all()
        )
        # common random numbers: the same noise in the same replicate of every design
        self.assertEqual(results.groupby("replicate").noise.nunique().max(), 1)
        self.assertEqual(results.groupby("design").noise.nunique().min(), 4)

    def test_independent_draws(self):
        """Without common random numbers, designs draw independently."""
        study = Study(
            self.space,
            build=lambda d: d.noise,
            evaluate=lambda noise: {"noise_out": noise},
            replicates=2,
            common_random_numbers=False,
        )
        results = study.run()
        self.assertEqual(results.groupby("replicate").noise.nunique().min(), 9)

    def test_frame_outcomes_and_duplicates(self):
        """Data frame outcomes give a record per row, and repeated designs run once."""
        designs = self.space.enumerate()[:2] * 3
        study = Study(
            self.space,
            build=lambda d: d.planes,
            evaluate=lambda planes: pd.DataFrame(
                {"sensor": ["a", "b"], "value": [planes, 2 * planes]}
            ).set_index("sensor"),
        )
        results = study.run(designs)
        self.assertEqual(len(results), 2 * 2)
        self.assertEqual(list(results.sensor), ["a", "b", "a", "b"])

    def test_infeasible_designs(self):
        """Designs rejected by validation are recorded as infeasible."""
        epoch = datetime(2026, 1, 1, tzinfo=timezone.utc)
        space = DesignSpace(
            altitude=Real(500e3, 600e3, step=100e3),
            planes=Integer(1, 3),
            satellites=Integer(1, 2),
        )

        def build(d):
            return WalkerConstellation(
                name="walker",
                orbit=CircularOrbit(altitude=d.altitude, inclination=60, epoch=epoch),
                number_planes=d.planes,
                number_satellites=d.satellites,
            ).generate_members()

        study = Study(space, build, lambda members: {"members": len(members)})
        results = study.run()
        self.assertEqual(len(results), 12)
        infeasible = results[~results.feasible]
        # more planes than satellites
        self.assertEqual(len(infeasible), 2 * 3)
        self.assertTrue(infeasible.error.str.contains("number planes").all())
        feasible = results[results.feasible]
        self.assertTrue((feasible.members == feasible.satellites).all())
        self.assertIsInstance(build(space.enumerate()[0])[0], Satellite)

    def test_progress(self):
        """Progress is reported as evaluations complete."""
        calls = []
        Study(
            self.space,
            build=lambda d: None,
            evaluate=lambda built: {},
            replicates=2,
        ).run(progress=lambda done, total: calls.append((done, total)))
        self.assertEqual(calls[0], (1, 18))
        self.assertEqual(calls[-1], (18, 18))

    def _noise_study(self, **kwargs):
        return Study(
            self.space,
            build=lambda d: (d.planes, d.per_plane, d.noise),
            evaluate=lambda built: {"score": built[0] * built[1] + built[2]},
            replicates=3,
            **kwargs,
        )

    def test_parallel_matches_sequential(self):
        """Evaluations in worker processes give the results of a sequential run."""
        sequential = self._noise_study().run()
        parallel = self._noise_study(n_jobs=2).run()
        columns = [c for c in sequential.columns if c != "seconds"]
        pd.testing.assert_frame_equal(sequential[columns], parallel[columns])

    def test_executor(self):
        """Evaluations run on a given executor, in the order of the designs."""
        with ThreadPoolExecutor(max_workers=3) as executor:
            threaded = self._noise_study(executor=executor).run()
        sequential = self._noise_study().run()
        pd.testing.assert_series_equal(threaded.score, sequential.score)
        self.assertEqual(self._noise_study(executor=executor).workers, 3)

    def test_summarize(self):
        """Summaries give the mean and spread of metrics over replicates."""
        results = self._noise_study().run()
        summary = summarize(results, ["score"])
        self.assertEqual(len(summary), 9)
        self.assertIn("planes", summary.columns)
        self.assertNotIn("noise", summary.columns)
        self.assertTrue((summary.replicates == 3).all())
        first = results[results.design == summary.design[0]]
        self.assertAlmostEqual(summary["score"][0], first.score.mean())
        self.assertAlmostEqual(summary["score std"][0], first.score.std())
        self.assertAlmostEqual(summary["score sem"][0], first.score.std() / np.sqrt(3))
        # a metric may also be a value of the design
        overlap = summarize(results, ["score", "planes"], design_columns=["planes"])
        self.assertEqual(list(overlap.columns).count("planes"), 1)

    def test_summarize_by(self):
        """Summaries group the records of an evaluation."""
        study = Study(
            self.space,
            build=lambda d: d.planes,
            evaluate=lambda planes: pd.DataFrame(
                {"sensor": ["a", "b"], "value": [planes, 2 * planes]}
            ).set_index("sensor"),
            replicates=2,
        )
        summary = summarize(study.run(), ["value"], by=["sensor"])
        self.assertEqual(len(summary), 18)
        b = summary[(summary.sensor == "b") & (summary.planes == 3)]
        self.assertTrue((b.value == 6).all())
