"""
Unit tests for searches of design spaces with optimizers.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import unittest

import numpy as np
import pandas as pd

from tatc.tradespace import (
    Categorical,
    Choice,
    DesignSpace,
    Evaluations,
    Infeasible,
    Integer,
    Minimum,
    Real,
    Repeat,
    Requirement,
    Study,
    Uncertain,
    hypervolume,
    make_pymoo_problem,
    optimize_optuna,
    optimize_pymoo,
    pareto_front,
    summarize,
)

HAS_OPTUNA = importlib.util.find_spec("optuna") is not None
HAS_PYMOO = importlib.util.find_spec("pymoo") is not None

OBJECTIVES = {"cost": "min", "value": "max"}
REFERENCE = {"cost": 30, "value": 0}


def _space():
    """Constellations of planes, or a single satellite of a quality."""
    return DesignSpace(
        layout=Choice(
            planes=DesignSpace(
                planes=Repeat(
                    DesignSpace(satellites=Integer(1, 4), size=Real(0.5, 2, step=0.5)),
                    count=Integer(1, 3),
                    ordered=False,
                )
            ),
            single=DesignSpace(quality=Categorical([1, 2, 3])),
        ),
        noise=Uncertain.uniform(-0.01, 0.01),
        # at most 8 satellites, violated by the satellites in excess
        constraints=[
            lambda d: (
                sum(p.satellites for p in d.layout.planes) - 8
                if d.layout.choice == "planes"
                else True
            )
        ],
    )


def _build(design):
    if design.layout.choice == "single":
        return {"cost": 2.0 * design.layout.quality, "value": design.layout.quality}
    if any(p.size == 2 and p.satellites == 4 for p in design.layout.planes):
        raise ValueError("too heavy")
    cost = sum(p.satellites * p.size + 1 for p in design.layout.planes)
    value = sum(np.sqrt(p.satellites) * p.size for p in design.layout.planes)
    return {"cost": cost, "value": value + design.noise}


def _study(**kwargs):
    return Study(_space(), _build, lambda metrics: metrics, replicates=2, **kwargs)


class TestEncoding(unittest.TestCase):
    """
    Unit tests for encoding designs in the unit hypercube.
    """

    def test_unit_paths(self):
        """Paths cover every decision, with repeated items up to the largest count."""
        paths = _space().unit_paths()
        self.assertEqual(
            paths[:3],
            [
                "layout",
                "layout.planes.planes.count",
                "layout.planes.planes[0].satellites",
            ],
        )
        self.assertEqual(len(paths), 1 + 1 + 3 * 2 + 1)
        self.assertEqual(paths[-1], "layout.single.quality")

    def test_round_trip(self):
        """Every design decodes from its encoding."""
        space = _space()
        for design in space.enumerate():
            self.assertEqual(space.decode(space.encode(design)).key, design.key)

    def test_decode_random(self):
        """Random points decode to designs of the space (or infeasible ones)."""
        space = _space()
        keys = {design.key for design in space.enumerate()}
        rng = np.random.default_rng(0)
        decoded = 0
        for _ in range(200):
            try:
                design = space.decode(rng.random(len(space.unit_paths())))
            except Infeasible:
                continue
            decoded += 1
            self.assertIn(design.key, keys)
        self.assertGreater(decoded, 100)

    def test_dynamic_count_bound(self):
        """Counts need fixed upper bounds to be encoded."""
        space = DesignSpace(
            n=Integer(1, 3),
            items=Repeat(Integer(0, 1), count=Integer(1, lambda d: d.n)),
        )
        with self.assertRaises(ValueError):
            space.unit_paths()


class TestEvaluations(unittest.TestCase):
    """
    Unit tests for evaluations of searches.
    """

    def test_memo(self):
        """Designs are evaluated once, and infeasible designs have no objectives."""
        evaluations = Evaluations(_study(), OBJECTIVES)
        designs = _space().enumerate()[:5]
        first = evaluations.evaluate(designs + designs)
        self.assertEqual(first[:5], first[5:])
        evaluations.evaluate(designs)
        self.assertEqual(len(evaluations.results), 5 * 2)

    def test_one_record_per_evaluation(self):
        """Objectives need one record per evaluation."""
        study = Study(
            _space(),
            lambda d: None,
            lambda built: pd.DataFrame({"cost": [1, 2], "value": [1, 2]}),
        )
        with self.assertRaises(ValueError):
            Evaluations(study, OBJECTIVES).evaluate(_space().enumerate()[:1])


def _exact_hypervolume():
    summary = summarize(_study().run(), list(OBJECTIVES))
    return summary, hypervolume(summary, OBJECTIVES, REFERENCE)


@unittest.skipUnless(HAS_OPTUNA, "requires optuna")
class TestOptuna(unittest.TestCase):
    """
    Unit tests for searches with Optuna.
    """

    def test_suggest(self):
        """Suggestions follow the kinds of the variables."""
        import optuna  # pylint: disable=import-outside-toplevel

        trial = optuna.trial.FixedTrial(
            {
                "layout": "planes",
                "layout.planes.planes.count": 2,
                "layout.planes.planes[0].satellites": 3,
                "layout.planes.planes[0].size": 1.5,
                "layout.planes.planes[1].satellites": 1,
                "layout.planes.planes[1].size": 0.5,
            }
        )
        design = _space().suggest(trial)
        self.assertEqual(
            [(p.satellites, p.size) for p in design.layout.planes], [(1, 0.5), (3, 1.5)]
        )

    def test_optimize(self):
        """A search finds most of the hypervolume of the true front."""
        import optuna  # pylint: disable=import-outside-toplevel

        from tatc.tradespace import (
            optimize_optuna,
        )  # pylint: disable=import-outside-toplevel

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        summary, exact = _exact_hypervolume()
        optuna_study, results = optimize_optuna(
            _study(),
            OBJECTIVES,
            n_trials=150,
            sampler=optuna.samplers.NSGAIISampler(seed=1, population_size=20),
        )
        self.assertEqual(len(optuna_study.trials), 150)
        found = summarize(results, list(OBJECTIVES))
        self.assertLessEqual(found.design.nunique(), 150)
        self.assertGreater(hypervolume(found, OBJECTIVES, REFERENCE), 0.9 * exact)
        self.assertTrue(pareto_front(summary, OBJECTIVES).any())


@unittest.skipUnless(HAS_PYMOO, "requires pymoo")
class TestPymoo(unittest.TestCase):
    """
    Unit tests for searches with pymoo.
    """

    def test_optimize(self):
        """NSGA-II on the encoded space finds most of the hypervolume of the true front."""
        from tatc.tradespace import (
            optimize_pymoo,
        )  # pylint: disable=import-outside-toplevel

        _, exact = _exact_hypervolume()
        result, results = optimize_pymoo(
            _study(), OBJECTIVES, termination=("n_gen", 15), seed=1
        )
        self.assertEqual(result.F.shape[1], 2)
        found = summarize(results, list(OBJECTIVES))
        self.assertGreater(hypervolume(found, OBJECTIVES, REFERENCE), 0.9 * exact)
        self.assertFalse(results.feasible.all())  # some designs were too heavy


@unittest.skipUnless(HAS_OPTUNA, "requires optuna")
class TestOptunaConstraints(unittest.TestCase):
    """
    Unit tests for constraints in searches with Optuna.
    """

    def setUp(self):
        import optuna  # pylint: disable=import-outside-toplevel

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        self.optuna = optuna

    def test_violations(self):
        """Trials record the violations of their designs' constraints."""
        optuna_study, _ = optimize_optuna(
            _study(),
            OBJECTIVES,
            n_trials=80,
            sampler=self.optuna.samplers.RandomSampler(seed=3),
        )
        states = {trial.state for trial in optuna_study.trials}
        self.assertEqual(states, {self.optuna.trial.TrialState.COMPLETE})
        violations = [trial.constraints["feasibility"] for trial in optuna_study.trials]
        self.assertTrue(any(v == 0 for v in violations))
        excess = [v for v in violations if v > 0]
        self.assertTrue(excess)
        # violations of the space count satellites in excess; the study's are 1
        for trial in optuna_study.trials:
            if "error" in trial.user_attrs:
                planes = [
                    trial.params[f"layout.planes.planes[{k}].satellites"]
                    for k in range(trial.params["layout.planes.planes.count"])
                ]
                self.assertEqual(trial.constraints["feasibility"], sum(planes) - 8)

    def test_pruning(self):
        """Without constraints, infeasible trials are pruned."""
        optuna_study, _ = optimize_optuna(
            _study(),
            OBJECTIVES,
            n_trials=60,
            sampler=self.optuna.samplers.RandomSampler(seed=3),
            constraints=False,
        )
        states = {trial.state for trial in optuna_study.trials}
        self.assertIn(self.optuna.trial.TrialState.PRUNED, states)
        self.assertTrue(all(not trial.constraints for trial in optuna_study.trials))


@unittest.skipUnless(HAS_PYMOO, "requires pymoo")
class TestPymooConstraints(unittest.TestCase):
    """
    Unit tests for constraints of pymoo problems.
    """

    def test_violations(self):
        """Violations of the space are the constraint values of the problem."""
        from pymoo.core.evaluator import (
            Evaluator,
        )  # pylint: disable=import-outside-toplevel
        from pymoo.core.population import (
            Population,
        )  # pylint: disable=import-outside-toplevel

        space = _space()
        problem = make_pymoo_problem(_study(), OBJECTIVES)
        light = space.design({"layout": ["single", {"quality": 2}]})
        # three planes of four satellites: the first branch, the largest count
        # and number of satellites (the last of their levels)
        paths = space.unit_paths()
        heavy = np.full(len(paths), 0.1)
        for i, path in enumerate(paths):
            if path.endswith("count") or path.endswith("satellites"):
                heavy[i] = 0.99
        with self.assertRaises(Infeasible) as context:
            space.decode(heavy)
        self.assertEqual(context.exception.violation, 4)
        x = np.array([space.encode(light), heavy])
        population = Evaluator().eval(problem, Population.new(X=x))
        self.assertEqual(list(population.get("G")[:, 0]), [0, 4])


@unittest.skipUnless(HAS_OPTUNA and HAS_PYMOO, "requires optuna and pymoo")
class TestRequirements(unittest.TestCase):
    """
    Unit tests for requirements in searches.
    """

    REQUIREMENTS = [Requirement("value p0", ">=", 2)]
    STATISTICS = {"value p0": Minimum("value")}

    def test_evaluations(self):
        """Evaluations give the violations of requirements on statistics."""
        evaluations = Evaluations(
            _study(), OBJECTIVES, self.STATISTICS, self.REQUIREMENTS
        )
        designs = [d for d in _space().enumerate() if d.layout.choice == "single"]
        values = evaluations.evaluate(designs)
        for design, value in zip(designs, values):
            violation = evaluations.violations[design.key][0]
            # the least value over replicates is the quality less the noise
            self.assertAlmostEqual(violation, max(2 - value[1], 0), delta=0.02)

    def test_optuna(self):
        """Trials record the violation of each requirement."""
        import optuna  # pylint: disable=import-outside-toplevel

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        optuna_study, _ = optimize_optuna(
            _study(),
            OBJECTIVES,
            n_trials=40,
            sampler=optuna.samplers.RandomSampler(seed=1),
            statistics=self.STATISTICS,
            requirements=self.REQUIREMENTS,
        )
        keys = {tuple(sorted(trial.constraints)) for trial in optuna_study.trials}
        self.assertEqual(keys, {("feasibility", "value p0 >= 2")})
        violations = [t.constraints["value p0 >= 2"] for t in optuna_study.trials]
        self.assertTrue(any(v > 0 for v in violations))
        with self.assertRaises(ValueError):
            optimize_optuna(
                _study(),
                OBJECTIVES,
                n_trials=1,
                constraints=False,
                requirements=self.REQUIREMENTS,
                statistics=self.STATISTICS,
            )

    def test_pymoo(self):
        """Problems have a constraint for each requirement."""
        problem = make_pymoo_problem(
            _study(), OBJECTIVES, self.STATISTICS, self.REQUIREMENTS
        )
        self.assertEqual(problem.n_ieq_constr, 2)
        result, _ = optimize_pymoo(
            _study(),
            OBJECTIVES,
            termination=("n_gen", 3),
            seed=1,
            statistics=self.STATISTICS,
            requirements=self.REQUIREMENTS,
        )
        self.assertEqual(result.pop.get("G").shape[1], 2)


@unittest.skipUnless(HAS_OPTUNA, "optuna is not installed")
class TestOptunaSearch(unittest.TestCase):
    """
    Unit tests for searches with Optuna: dynamic options, batches, and budgets.
    """

    def test_dynamic_options(self):
        """Categorical options that depend on earlier decisions are searched."""
        import optuna  # pylint: disable=import-outside-toplevel

        space = DesignSpace(
            kind=Categorical(["a", "b"]),
            option=Categorical(lambda d: ["x", "y"] if d.kind == "a" else ["z"]),
        )
        study = Study(space, lambda d: d, lambda d: {"f": float(d.option == "x")})
        optuna_study, results = optimize_optuna(
            study,
            {"f": "max"},
            n_trials=20,
            sampler=optuna.samplers.RandomSampler(seed=1),
        )
        self.assertEqual(set(results.option), {"x", "y", "z"})
        self.assertEqual(optuna_study.best_value, 1.0)

    def test_default_sampler(self):
        """The default sampler lies about pending trials in batches."""
        study = Study(DesignSpace(x=Integer(0, 9)), lambda d: d, lambda d: {"f": d.x})
        batched, _ = optimize_optuna(study, {"f": "min"}, n_trials=4, batch_size=2)
        single, _ = optimize_optuna(study, {"f": "min"}, n_trials=2, batch_size=1)
        # pylint: disable=protected-access
        self.assertTrue(batched.sampler._constant_liar)
        self.assertFalse(single.sampler._constant_liar)

    def test_max_designs(self):
        """The search stops once it has evaluated the most designs."""
        study = Study(DesignSpace(x=Integer(0, 2)), lambda d: d, lambda d: {"f": d.x})
        optuna_study, results = optimize_optuna(
            study, {"f": "min"}, n_trials=50, max_designs=2, batch_size=1
        )
        self.assertEqual(results.design.nunique(), 2)
        self.assertLessEqual(len(optuna_study.trials), 50)

    def test_partially_infeasible(self):
        """A design infeasible in some replicates violates by their share."""
        space = DesignSpace(x=Integer(0, 1), noise=Uncertain.uniform(0, 1))

        def build(d):
            if d.x == 1 and d.noise < 0.5:
                raise ValueError("unlucky")
            return d

        study = Study(space, build, lambda d: {"f": d.x}, replicates=8)
        evaluations = Evaluations(study, {"f": "max"})
        values = evaluations.evaluate(space.enumerate())
        self.assertIsNotNone(values[0])
        self.assertIsNone(values[1])
        share = evaluations.infeasibility[space.enumerate()[1].key]
        self.assertGreater(share, 0)
        self.assertLess(share, 1)
        optuna_study, _ = optimize_optuna(study, {"f": "max"}, n_trials=6)
        for trial in optuna_study.trials:
            if trial.params.get("x") == 1:
                self.assertAlmostEqual(trial.constraints["feasibility"], share)

    def test_records_checked_first(self):
        """A study of several records per evaluation fails after one design."""
        calls = []

        def evaluate(d):
            calls.append(d.x)
            return pd.DataFrame({"f": [d.x, d.x]}, index=pd.Index([0, 1], name="k"))

        study = Study(DesignSpace(x=Integer(0, 9)), lambda d: d, evaluate)
        with self.assertRaises(ValueError):
            Evaluations(study, {"f": "min"}).evaluate(study.space.enumerate())
        self.assertEqual(len(calls), 1)
