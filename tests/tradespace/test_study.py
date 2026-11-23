"""
Unit tests for trade studies.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import math
import tempfile
import unittest
from types import SimpleNamespace
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
    Derived,
    DesignSpace,
    Fixed,
    Infeasible,
    Integer,
    Real,
    ResultStore,
    Study,
    Uncertain,
    summarize,
)

HAS_JOBLIB = importlib.util.find_spec("joblib") is not None


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
        # a first parameter named design is given the built object (as here,
        # where it is the design)
        results = Study(
            self.space, lambda d: d, lambda design: {"value": design.planes}
        ).run()
        self.assertTrue((results.value == results.planes).all())
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

    @unittest.skipUnless(HAS_JOBLIB, "requires joblib")
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
        summary = summarize(results, statistics=["score"])
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
        overlap = summarize(results, statistics=["score", "planes"])
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
        summary = summarize(results, statistics=["score"])
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
        summary = summarize(study.run(), statistics=["value"], by=["sensor"])
        self.assertEqual(len(summary), 18)
        b = summary[(summary.sensor == "b") & (summary.planes == 3)]
        self.assertTrue((b.value == 6).all())


class TestErrors(unittest.TestCase):
    """
    Unit tests for errors of evaluations.
    """

    def setUp(self):
        self.space = DesignSpace(a=Integer(1, 4))
        # the state of the evaluations, which their fingerprint does not see
        # (an object of an installed class is described by its type)
        self.state = SimpleNamespace(fail=True)

    def _evaluate(self, d):
        if self.state.fail and d.a == 3:
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
            store = ResultStore(Path(directory) / "results.sqlite", namespace="errors")
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
            self.state.fail = False
            again = study.run()
            self.assertTrue(again.feasible.all())
            self.assertEqual(len(store), 4)

    def test_realization(self):
        """Errors of the realization of a design (such as of a value derived
        from a drawn uncertain parameter) are recorded too, with the values of
        the design."""
        space = DesignSpace(
            n=Integer(1, 2),
            u=Uncertain.uniform(-1, 1),
            root=Derived(lambda d: math.sqrt(d.u)),
        )
        with tempfile.TemporaryDirectory() as directory:
            store = ResultStore(Path(directory) / "results.sqlite", namespace="draws")
            study = Study(
                space,
                lambda d: d.n,
                lambda n: {"m": n},
                replicates=8,
                errors="record",
                store=store,
            )
            results = study.run()
            failed = results[~results.feasible]
            self.assertTrue(len(failed) > 0)
            self.assertTrue(failed.error.str.startswith("ValueError: ").all())
            self.assertEqual(sorted(failed.n.unique()), [1, 2])
            self.assertEqual(len(store), int(results.feasible.sum()))
        with self.assertRaises(ValueError):
            study.replace(errors="raise", store=None).run()


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
        return ResultStore(Path(directory) / "results.sqlite", namespace="feasibility")

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
        # also `scenario` (and the marker of unstored records) without scenarios
        for metrics in (
            {"replicate": 0},
            {"error": 0.1},
            {"a": 9},
            {"noise": 0.5},
            {"scenario": "mine"},
            {"_unstored": False},
        ):
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


class TestDesignsOfOtherSpaces(unittest.TestCase):
    """
    Unit tests for designs of other design spaces given to a study.
    """

    def test_evaluated_in_study_space(self):
        """Designs of another design space are evaluated as the designs of the
        study's design space, also with a store and new uncertain parameters."""
        space = DesignSpace(n=Integer(1, 2), altitude=Fixed(500))
        designs = space.enumerate()
        higher = space.extend(altitude=Fixed(700))
        with tempfile.TemporaryDirectory() as directory:
            store = ResultStore(Path(directory) / "results.sqlite")
            study = Study(
                higher, lambda d: d.altitude * d.n, lambda x: {"score": x}, store=store
            )
            self.assertEqual(list(study.run(designs).score), [700, 1400])
            self.assertEqual(list(study.run().score), [700, 1400])
        self.assertEqual(study.realize(designs[0], 0).altitude, 700)
        risky = space.extend(growth=Uncertain.uniform(1, 2))
        results = Study(risky, lambda d: d.n * d.growth, lambda x: {"s": x}).run(
            designs
        )
        self.assertTrue(results.growth.between(1, 2).all())

    def test_other_decisions(self):
        """Designs with decisions that the study's design space does not have
        (such as of a design space that extends it) raise, rather than being
        evaluated without them."""
        base = DesignSpace(a=Integer(1, 2))
        spared = base.extend(spares=Integer(0, 1)).design({"a": 1, "spares": 1})
        study = Study(base, lambda d: d.a, lambda a: {"m": a})
        with self.assertRaisesRegex(ValueError, "another design space.*spares"):
            study.run([spared])
        with self.assertRaisesRegex(ValueError, "another design space.*spares"):
            study.realize(spared)


class TestStudyScenarios(unittest.TestCase):
    """
    Unit tests for the constraints of design spaces in the scenarios of a
    study, and the names of scenarios of realized designs.
    """

    def test_infeasible_in_scenario(self):
        """A design that violates a constraint in a scenario is infeasible in it."""
        space = DesignSpace(
            n=Integer(1, 4), limit=Fixed(4), constraints=[lambda d: d.n - d.limit]
        )
        study = Study(
            space,
            lambda d: d.n,
            lambda x: {"score": x},
            scenarios={"strict": {"limit": 2}, "loose": {}},
        )
        results = study.run()
        feasible = results.pivot(index="n", columns="scenario", values="feasible")
        self.assertEqual(list(feasible.loose), [True] * 4)
        self.assertEqual(list(feasible.strict), [True, True, False, False])
        errors = results[~results.feasible].error
        self.assertTrue(errors.str.contains("constraint 0").all())
        # with the values of the design, from which designs are rebuilt
        rebuilt = space.designs_of(results)
        self.assertEqual(list(rebuilt.map(lambda d: d.n)), list(results.n))
        with self.assertRaises(Infeasible):
            study.realize(space.design({"n": 4}), 0, "strict")

    def test_scenario_of_design(self):
        """An evaluation is given the design realized in its scenario."""
        space = DesignSpace(n=Integer(1, 2), level=Fixed(1.0))
        study = Study(
            space,
            lambda d: d,
            lambda d, design: {"named": design.scenario, "level": design.level},
            scenarios={"high": {"level": 3.0}, "low": {}},
        )
        results = study.run()
        self.assertEqual(list(results.named), list(results.scenario))
        self.assertEqual(set(results[results.scenario == "high"].level), {3.0})
        self.assertEqual(
            study.realize(space.enumerate()[0], 0, "high").scenario, "high"
        )
        self.assertIsNone(
            Study(space, lambda d: d, dict).realize(space.enumerate()[0]).scenario
        )


class TestArguments(unittest.TestCase):
    """
    Unit tests for the arguments of studies.
    """

    def test_seed(self):
        """A seed is a non-negative integer, checked by the study."""
        space = DesignSpace(n=Integer(1, 2), u=Uncertain.uniform(0, 1))
        with self.assertRaisesRegex(ValueError, "seed"):
            Study(space, lambda d: d.u, lambda x: {"s": x}, seed=-1)

    def test_keywords(self):
        """The settings of a study, and the options of a run, are keywords."""
        space = DesignSpace(n=Integer(1, 2))
        with self.assertRaises(TypeError):
            Study(space, lambda d: d, dict, 3)  # pylint: disable=too-many-function-args
        study = Study(space, lambda d: d.n, lambda n: {"m": n})
        with self.assertRaises(TypeError):
            study.run(None, print)  # pylint: disable=too-many-function-args

    def test_indices(self):
        """Replicate indices and the offsets of shifts are checked before any
        evaluation."""
        calls = []
        space = DesignSpace(n=Integer(1, 2), u=Uncertain.uniform(0, 1))
        study = Study(space, lambda d: d.n, lambda n: calls.append(n) or {"m": n})
        for options in (
            {"replicate_indices": [0, -1]},
            {"replicate_indices": [0.5]},
            {"shift": {"u": -1}},
            {"shift": {"u": 1.5}},
        ):
            with self.subTest(options=options):
                with self.assertRaisesRegex(ValueError, "non-negative integer"):
                    study.run(**options)
        self.assertEqual(calls, [])
        shifted = study.run(replicate_indices=[2], shift={"u": 1})
        plain = study.run(replicate_indices=[3])
        self.assertEqual(list(shifted.u), list(plain.u))


class TestEmptyResults(unittest.TestCase):
    """
    Unit tests for evaluations without records and runs without evaluations.
    """

    def setUp(self):
        self.space = DesignSpace(n=Integer(0, 2))

    @staticmethod
    def _evaluate(n):
        # a record per sensor, of which a design may have none
        frame = pd.DataFrame({"sensor": list(range(n)), "score": [1.0] * n})
        return frame.set_index("sensor")

    def test_empty_frame(self):
        """An evaluation that gives an empty data frame raises an error."""
        for errors in ("raise", "record"):
            with self.subTest(errors=errors):
                study = Study(self.space, lambda d: d.n, self._evaluate, errors=errors)
                with self.assertRaisesRegex(ValueError, "empty data frame"):
                    study.run()

    def test_no_designs(self):
        """A run of no designs gives the status columns of results."""
        results = Study(self.space, lambda d: d.n, self._evaluate).run([])
        self.assertEqual(
            list(results.columns),
            ["design", "replicate", "feasible", "error", "seconds"],
        )
        self.assertEqual(len(results), 0)
        scenarios = Study(
            self.space, lambda d: d.n, self._evaluate, scenarios={"a": {}}
        ).run([])
        self.assertEqual(
            list(scenarios.columns)[:3], ["design", "replicate", "scenario"]
        )


class TestWarningLocation(unittest.TestCase):
    """
    Unit tests for the location of the warnings of studies.
    """

    def test_no_feasible_designs(self):
        """The warning of no feasible designs points at the run."""

        def build(d):
            raise Infeasible("never")

        study = Study(DesignSpace(a=Integer(1, 2)), build, dict)
        with self.assertWarnsRegex(UserWarning, "none of the 2 designs") as context:
            study.run()
        self.assertEqual(Path(context.filename).resolve(), Path(__file__).resolve())


class TestDescriptions(unittest.TestCase):
    """
    Unit tests for the descriptions of designs in the results of studies.
    """

    def setUp(self):
        self.space = DesignSpace(
            planes=Integer(1, 2),
            per_plane=Integer(1, 2),
            noise=Uncertain.uniform(0, 1),
        )
        self.calls = []

        def describe(d):
            self.calls.append(d.key)
            return f"{d.planes} x {d.per_plane}"

        self.describe = describe

    def _study(self, **options):
        return Study(
            self.space,
            lambda d: d,
            lambda d: {"satellites": d.planes * d.per_plane, "score": d.noise},
            replicates=3,
            describe=self.describe,
            **options,
        )

    def test_column(self):
        """A description of each design follows its values, given once per design
        (of the design, not its replicates), and summaries keep it."""
        results = self._study().run()
        columns = list(results.columns)
        self.assertEqual(columns.index("description"), columns.index("feasible") - 1)
        self.assertEqual(columns.index("description"), columns.index("noise") + 1)
        labels = [f"{p} x {s}" for p in (1, 2) for s in (1, 2)]
        self.assertEqual(list(results.description), [x for x in labels for _ in "abc"])
        self.assertEqual(len(self.calls), 4)
        summary = summarize(results, {"satellites": "min", "score": "max"})
        self.assertEqual(list(summary.description), labels)
        # without a function, no descriptions
        self.assertNotIn("description", self._study().replace(describe=None).run())
        with self.assertRaises(TypeError):
            self._study().replace(describe="planes")

    def test_not_stored(self):
        """Descriptions are not stored, nor part of the fingerprint, so that
        another function needs no new evaluations."""
        with tempfile.TemporaryDirectory() as directory:
            store = ResultStore(Path(directory) / "results.sqlite")
            study = self._study(store=store)
            first = study.run()
            other = study.replace(describe=lambda d: f"{d.planes} planes")
            self.assertEqual(other.fingerprint, study.fingerprint)
            progress = []
            again = other.run(progress=lambda done, total: progress.append(done))
            self.assertEqual(progress, [])  # every evaluation from the store
            self.assertEqual(list(again.description.unique()), ["1 planes", "2 planes"])
            pd.testing.assert_frame_equal(
                again.drop(columns="description"), first.drop(columns="description")
            )
            records = store._get(  # pylint: disable=protected-access
                [(key, 0) for key in first.design.unique()], study.fingerprint
            )
            self.assertEqual(len(records), 4)
            self.assertTrue(
                all("description" not in r for rs in records.values() for r in rs)
            )

    def test_other_spaces_and_scenarios(self):
        """Designs of another design space are described as designs of the
        study's, in each scenario."""
        extended = self.space.extend(limit=Fixed(1))
        study = self._study().replace(
            space=extended, scenarios={"a": {}, "b": {"limit": 2}}
        )
        results = study.run(self.space.enumerate()[:2])
        self.assertEqual(list(results.description.unique()), ["1 x 1", "1 x 2"])
        self.assertEqual(len(self.calls), 2)
        self.assertEqual(list(results.columns[:3]), ["design", "replicate", "scenario"])

    def test_names(self):
        """Metrics and values may not be named as the descriptions."""
        study = Study(
            self.space,
            lambda d: d,
            lambda d: {"description": "a metric"},
            describe=self.describe,
        )
        with self.assertRaisesRegex(ValueError, "description"):
            study.run()
        # but may without descriptions
        self.assertIn("description", study.replace(describe=None).run())
        with self.assertRaisesRegex(ValueError, "named 'description'"):
            self._study().replace(space=self.space.extend(description=Fixed("x")))

    def test_no_designs(self):
        """A run of no designs gives the column of descriptions."""
        results = self._study().run([])
        self.assertEqual(
            list(results.columns),
            ["design", "replicate", "description", "feasible", "error", "seconds"],
        )
