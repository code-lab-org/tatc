"""
Unit tests for trade studies.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import tempfile
import unittest
import warnings
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, time, timezone

import numpy as np
import pandas as pd

from tatc.schemas import (
    CircularOrbit,
    Satellite,
    SunSynchronousOrbit,
    WalkerConstellation,
)
from tatc.tradespace import (
    DesignSpace,
    Infeasible,
    Integer,
    Real,
    ResultStore,
    Study,
    Uncertain,
    summarize,
)


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

    def test_design_argument(self):
        """An evaluation with a parameter named design is given the realized design."""
        build = lambda d: d.planes  # pylint: disable=unnecessary-lambda-assignment
        given = Study(
            self.space,
            build,
            lambda planes, design: {"score": planes + design.noise},
            replicates=2,
        ).run()
        np.testing.assert_allclose(given.score, given.planes + given.noise)

        def keyword(planes, *, design):
            return {"satellites": planes * design.per_plane}

        results = Study(self.space, build, keyword).run()
        self.assertTrue(
            (results.satellites == results.planes * results.per_plane).all()
        )
        # functions without the parameter, or without signatures, are given the
        # built object only
        Study(self.space, build, lambda planes: {"score": planes}).run()
        results = Study(self.space, lambda d: {"value": d.planes}, dict).run()
        self.assertEqual(list(results.value), list(results.planes))

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
        self.assertAlmostEqual(summary["score se"][0], first.score.std() / np.sqrt(3))
        self.assertNotIn("score sem", summary.columns)
        # a metric may also be a value of the design
        overlap = summarize(results, ["score", "planes"], design_columns=["planes"])
        self.assertEqual(list(overlap.columns).count("planes"), 1)

    def test_replace_and_frames(self):
        """Variants of a study replace arguments; designs may be data frames."""
        study = Study(
            self.space,
            lambda d: d,
            lambda d: {"score": d.planes + d.noise},
            replicates=2,
            seed=3,
        )
        more = study.replace(replicates=4)
        self.assertEqual((more.replicates, more.seed), (4, 3))
        self.assertIs(more.evaluate, study.evaluate)
        self.assertEqual(study.replicates, 2)
        results = more.run()
        self.assertEqual(len(results), 9 * 4)
        with self.assertRaisesRegex(TypeError, "no arguments"):
            study.replace(replicate=3)
        # designs given by a data frame of their values (such as a summary)
        summary = summarize(results, ["score"])
        best = summary.nlargest(2, "score")
        again = study.run(best)
        self.assertEqual(list(again.design.unique()), list(best.design))
        with self.assertRaisesRegex(ValueError, "not given"):
            study.run(list(best.design))

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


class TestErrors(unittest.TestCase):
    """
    Unit tests for errors of evaluations.
    """

    def setUp(self):
        self.space = DesignSpace(a=Integer(1, 4))
        self.fail = True

    def _evaluate(self, d):
        if self.fail and d.a == 3:
            raise RuntimeError("diverged")
        return {"value": d.a}

    def test_raise(self):
        """Errors other than rejections of designs are raised by default."""
        with self.assertRaises(RuntimeError):
            Study(self.space, lambda d: d, self._evaluate).run()
        with self.assertRaises(ValueError):
            Study(self.space, lambda d: d, self._evaluate, errors="ignore")

    def test_record(self):
        """Errors are recorded as infeasible evaluations, which are not stored."""
        with tempfile.TemporaryDirectory() as directory:
            store = ResultStore(Path(directory) / "results.sqlite", "errors")
            study = Study(
                self.space, lambda d: d, self._evaluate, errors="record", store=store
            )
            results = study.run()
            failed = results[~results.feasible]
            self.assertEqual(list(failed.a), [3])
            self.assertEqual(failed.error.iloc[0], "RuntimeError: diverged")
            self.assertFalse([c for c in results.columns if c.startswith("_")])
            self.assertEqual(len(store), 3)
            # a later run evaluates the failed evaluation again
            self.fail = False
            again = study.run()
            self.assertTrue(again.feasible.all())
            self.assertEqual(len(store), 4)


class TestFeasibility(unittest.TestCase):
    """
    Unit tests for designs that studies reject, and mistakes in their functions.
    """

    def setUp(self):
        self.space = DesignSpace(a=Integer(1, 4))
        self.mistaken = True
        self.built = []

    @staticmethod
    def _store(directory: str) -> ResultStore:
        return ResultStore(Path(directory) / "results.sqlite", "feasibility")

    def test_rejected(self):
        """Designs that either function rejects with Infeasible are stored as infeasible."""

        def build(d):
            if d.a == 1:
                raise Infeasible("too small")
            return d.a

        def evaluate(a):
            if a == 2:
                raise Infeasible("too slow")
            return {"value": a}

        with tempfile.TemporaryDirectory() as directory:
            store = self._store(directory)
            results = Study(self.space, build, evaluate, store=store).run()
            self.assertEqual(list(results.feasible), [False, False, True, True])
            self.assertEqual(list(results.error[:2]), ["too small", "too slow"])
            self.assertEqual(len(store), 4)

    def _build_orbit(self, d):
        self.built.append(d.a)
        return SunSynchronousOrbit(
            altitude=d.a * 2000e3,
            equator_crossing_time=time(10, 30),
            epoch=datetime(2026, 1, 1, tzinfo=timezone.utc),
        )

    def test_validation_not_stored(self):
        """Designs that fail validation are infeasible, and evaluated again later."""
        with tempfile.TemporaryDirectory() as directory:
            store = self._store(directory)
            study = Study(
                self.space,
                self._build_orbit,
                lambda orbit: {"inclination": orbit.get_inclination()},
                store=store,
            )
            results = study.run()
            # too high to be sun-synchronous
            self.assertEqual(list(results.feasible), [True, True, False, False])
            self.assertTrue(results.error[2:].str.contains("altitude").all())
            self.assertEqual(len(store), 2)
            self.built.clear()
            again = study.run()
            self.assertEqual(self.built, [3, 4])
            self.assertEqual(list(again.feasible), list(results.feasible))

    def _build_mistaken(self, d):
        if self.mistaken:
            # pylint: disable-next=unbalanced-tuple-unpacking,unused-variable
            first, second = [d.a]
        return d.a

    def test_mistakes_not_infeasible(self):
        """Other errors, such as a ValueError of a mistake, fail rather than reject designs."""

        def evaluate(a):
            return {"value": a}

        with self.assertRaisesRegex(ValueError, "not enough values to unpack"):
            Study(self.space, self._build_mistaken, evaluate).run()
        with tempfile.TemporaryDirectory() as directory:
            store = self._store(directory)
            study = Study(
                self.space,
                self._build_mistaken,
                evaluate,
                errors="record",
                store=store,
            )
            with self.assertWarnsRegex(UserWarning, "none of the 4 designs"):
                failed = study.run()
            self.assertTrue(failed.error.str.startswith("ValueError: not enough").all())
            self.assertEqual(len(store), 0)
            # a fix that the study's fingerprint cannot see
            self.mistaken = False
            self.assertTrue(study.run().feasible.all())
            self.assertEqual(len(store), 4)

    def test_warning_without_feasible(self):
        """A run of several designs, none feasible, warns with the most common error."""

        def build(d):
            raise Infeasible("never" if d.a > 1 else "once")

        study = Study(self.space, build, lambda built: {})
        with self.assertWarnsRegex(UserWarning, r"\(3 of 4 evaluations\): never"):
            study.run()
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            # a single design, or any feasible design
            study.run(self.space.enumerate()[:1])
            Study(self.space, lambda d: d, lambda d: {}).run()


class TestMetricNames(unittest.TestCase):
    """
    Unit tests for the names of the metrics of evaluations.
    """

    def test_reserved_names(self):
        """Metrics named as status columns, or replacing values of designs, raise errors."""
        space = DesignSpace(a=Integer(1, 2), noise=Uncertain.uniform(0, 1))
        for metrics in ({"replicate": 0}, {"error": 0.1}, {"a": 9}, {"noise": 0.5}):
            for errors in ("raise", "record"):
                with self.subTest(metrics=metrics, errors=errors):
                    study = Study(
                        space, lambda d: d, lambda d, m=metrics: m, errors=errors
                    )
                    with self.assertRaisesRegex(ValueError, "metrics named"):
                        study.run()

    def test_restated_values(self):
        """A metric may restate a value of the design, such as a variable that is an objective."""
        space = DesignSpace(a=Integer(1, 2), noise=Uncertain.uniform(0, 1))
        study = Study(
            space,
            lambda d: d,
            lambda d: {"a": float(d.a), "noise": d.noise, "value": 2 * d.a},
            replicates=2,
        )
        results = study.run()
        self.assertEqual(list(results.a), [1, 1, 2, 2])
        self.assertEqual(list(results.value), [2, 2, 4, 4])

    def test_reserved_index(self):
        """Named indices of data frame outcomes are metrics too."""
        space = DesignSpace(a=Integer(1, 2))
        frame = pd.DataFrame({"value": [1]}, index=pd.Index([0], name="a"))
        with self.assertRaisesRegex(ValueError, r"\['a'\]"):
            Study(space, lambda d: d, lambda d: frame).run()
        # other names are kept
        renamed = frame.rename_axis("sensor")
        results = Study(space, lambda d: d, lambda d: renamed).run()
        self.assertEqual(list(results.a), [1, 2])
        self.assertEqual(list(results.sensor), [0, 0])
