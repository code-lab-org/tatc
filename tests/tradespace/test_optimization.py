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
    Criteria,
    DesignSpace,
    Infeasible,
    Integer,
    Minimum,
    Real,
    Repeat,
    Requirement,
    Study,
    Uncertain,
    hypervolume,
    optimize_optuna,
    optimize_pymoo,
    pareto_front,
    summarize,
)
from tatc.tradespace.optimization import (
    Evaluations,
    _constraint_names,
    make_pymoo_problem,
)
from tatc.tradespace.values import Scope
from tatc.tradespace.walks import _Context, _walk

HAS_OPTUNA = importlib.util.find_spec("optuna") is not None
HAS_PYMOO = importlib.util.find_spec("pymoo") is not None

OBJECTIVES = {"cost": "min", "value": "max"}
REFERENCE = {"cost": 30, "value": 0}


def _encode(space, design):
    """The center of each decided level's part of [0, 1], and 0.5 for inactive paths."""
    paths = space._unit_paths()
    index = {path: i for i, path in enumerate(paths)}
    unit = np.full(len(index), 0.5)

    def record(path, domain, value):
        unit[index[path]] = (domain.levels.index(value) + 0.5) / len(domain.levels)
        return value

    _walk(space, "", Scope(), _Context(normalize=record), design.decisions)
    return unit


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
        raise Infeasible("too heavy")
    cost = sum(p.satellites * p.size + 1 for p in design.layout.planes)
    value = sum(np.sqrt(p.satellites) * p.size for p in design.layout.planes)
    return {"cost": cost, "value": value + design.noise}


def _study(**kwargs):
    return Study(_space(), _build, lambda metrics: metrics, replicates=2, **kwargs)


def _partly_infeasible():
    """A study of three designs, of which the first is feasible in every
    replicate, the second in some, and the third in none, and its designs."""
    space = DesignSpace(x=Integer(0, 2), noise=Uncertain.uniform(0, 1))

    def build(d):
        if d.x == 2 or (d.x == 1 and d.noise < 0.5):
            raise Infeasible("unlucky")
        return d

    study = Study(space, build, lambda d: {"f": d.x}, replicates=8)
    return study, space.enumerate()


class TestEncoding(unittest.TestCase):
    """
    Unit tests for encoding designs in the unit hypercube.
    """

    def test_unit_paths(self):
        """Paths cover every decision, with repeated items up to the largest count."""
        paths = _space()._unit_paths()
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
            self.assertEqual(space._decode(_encode(space, design)).key, design.key)

    def test_decode_random(self):
        """Random points decode to designs of the space (or infeasible ones)."""
        space = _space()
        keys = {design.key for design in space.enumerate()}
        rng = np.random.default_rng(0)
        decoded = 0
        for _ in range(200):
            try:
                design = space._decode(rng.random(len(space._unit_paths())))
            except Infeasible:
                continue
            decoded += 1
            self.assertIn(design.key, keys)
        self.assertGreater(decoded, 100)

    def test_dynamic_count_bound(self):
        """Counts need fixed upper bounds to be encoded."""
        space = DesignSpace(
            n=Integer(1, 3),
            parts=Repeat(Integer(0, 1), count=Integer(1, lambda d: d.n)),
        )
        with self.assertRaises(ValueError):
            space._unit_paths()


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

    def test_feasible_replicates(self):
        """A design infeasible in some replicates is judged by the others, as
        summaries judge it, unless a requirement on its share of feasible
        replicates rejects it; a design without a feasible replicate is
        infeasible."""
        study, designs = _partly_infeasible()
        results = study.run(designs)
        share = results[results.design == designs[1].key].feasible.mean()
        self.assertGreater(share, 0)
        self.assertLess(share, 1)
        evaluations = Evaluations(study, {"f": "max"})
        values = evaluations.evaluate(designs)
        self.assertEqual(values[:2], [(0.0,), (1.0,)])
        self.assertIsNone(values[2])
        self.assertEqual(
            [evaluations.infeasibility[d.key] for d in designs], [0.0, 0.0, 1.0]
        )
        # the summary of a search is that of its results (see summarize),
        # with the spreads of means
        search = evaluations.search(None)
        pd.testing.assert_frame_equal(search.summary, summarize(results, {"f": "max"}))
        self.assertEqual([d.key for d in search.designs], [designs[1].key])
        self.assertEqual(len(search.evaluated), 3)
        # a requirement on the share of feasible replicates
        reliable = Criteria(
            {"f": "max"}, requirements=[Requirement("feasible", ">=", 1)]
        )
        evaluations = Evaluations(study, reliable)
        evaluations.evaluate(designs)
        self.assertEqual(evaluations.violations[designs[0].key], (0.0,))
        self.assertAlmostEqual(evaluations.violations[designs[1].key][0], 1 - share)
        # a search of designs without a feasible replicate finds none
        none = Evaluations(study, reliable)
        none.evaluate(designs[2:])
        self.assertTrue(none.search(None).summary.empty)
        self.assertEqual(none.search(None).designs, [])

    def test_one_record_per_evaluation(self):
        """Objectives need one record per evaluation."""
        study = Study(
            _space(),
            lambda d: None,
            lambda built: pd.DataFrame({"cost": [1, 2], "value": [1, 2]}),
        )
        with self.assertRaises(ValueError):
            Evaluations(study, OBJECTIVES).evaluate(_space().enumerate()[:1])

    def test_constraint_names(self):
        """Requirements are constraints of distinct names: their statements,
        their values in full where these are alike, or their order."""
        requirements = [
            Requirement("a", ">=", 2.0000001),
            Requirement("a", ">=", 2.0000002),
            Requirement("b", "<=", 1),
            Requirement("b", "<=", 1, scale=2),
        ]
        self.assertEqual(
            _constraint_names(requirements),
            ["a >= 2.0000001", "a >= 2.0000002", "b <= 1", "b <= 1 (2)"],
        )


def _exact_hypervolume():
    summary = summarize(_study().run(), statistics=list(OBJECTIVES))
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
        design = _space()._suggest(trial)
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
        search = optimize_optuna(
            _study(),
            OBJECTIVES,
            n_trials=150,
            sampler=optuna.samplers.NSGAIISampler(seed=1, population_size=20),
        )
        optuna_study, results = search.state, search.results
        self.assertEqual(len(optuna_study.trials), 150)
        found = summarize(results, statistics=list(OBJECTIVES))
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
        search = optimize_pymoo(_study(), OBJECTIVES, termination=("n_gen", 15), seed=1)
        result, results = search.state, search.results
        self.assertEqual(result.F.shape[1], 2)
        self.assertEqual({d.key for d in search.evaluated}, set(results.design))
        # the designs of a search are the efficient designs it found
        self.assertEqual(
            [d.key for d in search.designs],
            list(search.summary.design[search.summary.pareto]),
        )
        self.assertEqual(search.criteria.objectives, OBJECTIVES)
        found = summarize(results, statistics=list(OBJECTIVES))
        self.assertGreater(hypervolume(found, OBJECTIVES, REFERENCE), 0.9 * exact)
        self.assertFalse(results.feasible.all())  # some designs were too heavy
        # the summary of the feasible designs, and their front
        summary = search.summary
        self.assertEqual(set(summary.design), set(results.design[results.feasible]))
        self.assertIn("value se", summary.columns)  # with the spreads of means
        np.testing.assert_array_equal(
            summary.pareto, pareto_front(summary, OBJECTIVES).to_numpy()
        )

    def test_duplicates(self):
        """Proposals of the same design, or of designs evaluated before, are
        duplicates; proposals that violate a constraint are not."""
        from pymoo.core.population import (  # pylint: disable=import-outside-toplevel
            Population,
        )

        study = _study()
        problem = make_pymoo_problem(study, OBJECTIVES)
        single = study.space.design({"layout": ["single", {"quality": 2}]})
        x = _encode(study.space, single)
        # the same level of each variable, and other values of the branch not chosen
        same = x.copy()
        same[study.space._unit_paths().index("layout.single.quality")] += 0.1
        other = x.copy()
        other[study.space._unit_paths().index("layout.planes.planes.count")] = 0.9
        too_many = _encode(
            study.space,
            study.space.design(
                {"layout": ["planes", {"planes": [{"satellites": 4, "size": 0.5}] * 2}]}
            ),
        )
        too_many[study.space._unit_paths().index("layout.planes.planes.count")] = 0.99
        for k in range(3):
            too_many[
                study.space._unit_paths().index(f"layout.planes.planes[{k}].satellites")
            ] = 0.99
        with self.assertRaises(Infeasible):
            study.space._decode(too_many)
        pop = Population.new("X", np.array([x, same, other, too_many, too_many]))
        kept = problem.duplicates.do(pop)
        self.assertEqual(len(kept), 3)  # one of the same design, and both infeasible
        # designs evaluated before are duplicates
        problem.evaluations.evaluate([single])
        self.assertEqual(len(problem.duplicates.do(pop)), 2)

    def test_distinct(self):
        """Distinct searches evaluate more designs in as many generations."""
        distinct = optimize_pymoo(
            _study(), OBJECTIVES, termination=("n_gen", 5), seed=1
        )
        plain = optimize_pymoo(
            _study(), OBJECTIVES, termination=("n_gen", 5), seed=1, distinct=False
        )
        self.assertGreater(len(distinct.evaluated), len(plain.evaluated))
        with self.assertRaises(ValueError):
            optimize_pymoo(_study(), OBJECTIVES, algorithm=object(), seed=1)

    def test_algorithm_unchanged(self):
        """A distinct search does not change the algorithm it is given, which
        another search may take."""
        # pylint: disable-next=import-outside-toplevel
        from pymoo.algorithms.moo.nsga2 import NSGA2

        algorithm = NSGA2(pop_size=8)
        own = algorithm.mating.eliminate_duplicates
        optimize_pymoo(
            _study(), OBJECTIVES, algorithm=algorithm, termination=("n_gen", 2), seed=1
        )
        self.assertIs(algorithm.mating.eliminate_duplicates, own)
        self.assertIs(algorithm.initialization.eliminate_duplicates, own)
        plain = optimize_pymoo(
            _study(),
            OBJECTIVES,
            algorithm=algorithm,
            termination=("n_gen", 2),
            seed=1,
            distinct=False,
        )
        self.assertGreater(len(plain.evaluated), 0)

    def test_max_designs(self):
        """A search stops after the generation in which it has evaluated the
        most designs, or by its termination."""
        from pymoo.algorithms.moo.nsga2 import (  # pylint: disable=import-outside-toplevel
            NSGA2,
        )

        search = optimize_pymoo(
            _study(),
            OBJECTIVES,
            algorithm=NSGA2(pop_size=10),
            termination=("n_gen", 50),
            seed=1,
            max_designs=25,
        )
        evaluated = search.results.design.nunique()
        self.assertGreaterEqual(evaluated, 25)
        self.assertLess(evaluated, 25 + 10)
        self.assertLess(search.state.algorithm.n_gen, 50)
        # the termination, if first
        short = optimize_pymoo(
            _study(),
            OBJECTIVES,
            algorithm=NSGA2(pop_size=10),
            termination=("n_gen", 2),
            seed=1,
            max_designs=1000,
        )
        self.assertLessEqual(short.results.design.nunique(), 20)
        with self.assertRaises(ValueError):
            optimize_pymoo(_study(), OBJECTIVES, max_designs=0)


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
        optuna_study = optimize_optuna(
            _study(),
            OBJECTIVES,
            n_trials=80,
            sampler=self.optuna.samplers.RandomSampler(seed=3),
        ).state
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
        optuna_study = optimize_optuna(
            _study(),
            OBJECTIVES,
            n_trials=60,
            sampler=self.optuna.samplers.RandomSampler(seed=3),
            constraints=False,
        ).state
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
        paths = space._unit_paths()
        heavy = np.full(len(paths), 0.1)
        for i, path in enumerate(paths):
            if path.endswith("count") or path.endswith("satellites"):
                heavy[i] = 0.99
        with self.assertRaises(Infeasible) as context:
            space._decode(heavy)
        self.assertEqual(context.exception.violation, 4)
        x = np.array([_encode(space, light), heavy])
        population = Evaluator().eval(problem, Population.new(X=x))
        self.assertEqual(list(population.get("G")[:, 0]), [0, 4])


@unittest.skipUnless(HAS_OPTUNA and HAS_PYMOO, "requires optuna and pymoo")
class TestRequirements(unittest.TestCase):
    """
    Unit tests for requirements in searches.
    """

    CRITERIA = Criteria(
        OBJECTIVES,
        statistics={"value p0": Minimum("value")},
        requirements=[Requirement("value p0", ">=", 2)],
    )

    def test_evaluations(self):
        """Evaluations give the violations of requirements on statistics."""
        evaluations = Evaluations(_study(), self.CRITERIA)
        designs = [d for d in _space().enumerate() if d.layout.choice == "single"]
        values = evaluations.evaluate(designs)
        for design, value in zip(designs, values):
            violation = evaluations.violations[design.key][0]
            # the least value over replicates is the quality less the noise,
            # and violations are relative to the required value
            self.assertAlmostEqual(violation, max(2 - value[1], 0) / 2, delta=0.01)
        with self.assertRaises(ValueError):
            Evaluations(_study(), Criteria(OBJECTIVES, by="sensor"))

    def test_optuna(self):
        """Trials record the violation of each requirement."""
        import optuna  # pylint: disable=import-outside-toplevel

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        optuna_study = optimize_optuna(
            _study(),
            self.CRITERIA,
            n_trials=40,
            sampler=optuna.samplers.RandomSampler(seed=1),
        ).state
        keys = {tuple(sorted(trial.constraints)) for trial in optuna_study.trials}
        self.assertEqual(keys, {("feasibility", "value p0 >= 2")})
        violations = [t.constraints["value p0 >= 2"] for t in optuna_study.trials]
        self.assertTrue(any(v > 0 for v in violations))
        with self.assertRaises(ValueError):
            optimize_optuna(_study(), self.CRITERIA, n_trials=1, constraints=False)

    def test_optuna_names(self):
        """Requirements that state themselves alike are distinct constraints."""
        import optuna  # pylint: disable=import-outside-toplevel

        optuna.logging.set_verbosity(optuna.logging.WARNING)
        criteria = Criteria(
            OBJECTIVES,
            statistics={"value p0": Minimum("value")},
            requirements=[
                Requirement("value p0", ">=", 2.0000001),
                Requirement("value p0", ">=", 2.5000001),
                Requirement("value p0", ">=", 2.5000002),
            ],
        )
        optuna_study = optimize_optuna(
            _study(), criteria, n_trials=10, sampler=optuna.samplers.RandomSampler(1)
        ).state
        names = ("value p0 >= 2", "value p0 >= 2.5000001", "value p0 >= 2.5000002")
        keys = {tuple(sorted(trial.constraints)) for trial in optuna_study.trials}
        self.assertEqual(keys, {("feasibility", *names)})

    def test_pymoo(self):
        """Problems have a constraint for each requirement."""
        problem = make_pymoo_problem(_study(), self.CRITERIA)
        self.assertEqual(problem.n_ieq_constr, 2)
        result = optimize_pymoo(
            _study(), self.CRITERIA, termination=("n_gen", 3), seed=1
        ).state
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
        search = optimize_optuna(
            study,
            {"f": "max"},
            n_trials=20,
            sampler=optuna.samplers.RandomSampler(seed=1),
        )
        optuna_study, results = search.state, search.results
        self.assertEqual(set(results.option), {"x", "y", "z"})
        self.assertEqual(optuna_study.best_value, 1.0)

    def test_default_sampler(self):
        """The default sampler lies about pending trials in batches."""
        study = Study(DesignSpace(x=Integer(0, 9)), lambda d: d, lambda d: {"f": d.x})
        batched = optimize_optuna(study, {"f": "min"}, n_trials=4, batch_size=2).state
        single = optimize_optuna(study, {"f": "min"}, n_trials=2, batch_size=1).state
        # pylint: disable=protected-access
        self.assertTrue(batched.sampler._constant_liar)
        self.assertFalse(single.sampler._constant_liar)

    def test_max_designs(self):
        """The search stops once it has evaluated the most designs."""
        study = Study(DesignSpace(x=Integer(0, 2)), lambda d: d, lambda d: {"f": d.x})
        search = optimize_optuna(
            study, {"f": "min"}, n_trials=50, max_designs=2, batch_size=1
        )
        optuna_study, results = search.state, search.results
        self.assertEqual(results.design.nunique(), 2)
        # the designs evaluated by the search, in the order of their evaluation
        self.assertEqual(
            [design.key for design in search.evaluated],
            results.design.drop_duplicates().tolist(),
        )
        self.assertLessEqual(len(optuna_study.trials), 50)
        with self.assertRaises(ValueError):
            optimize_optuna(study, {"f": "min"}, n_trials=5, max_designs=0)

    def test_partially_infeasible(self):
        """A design without a feasible replicate violates the constraint of
        feasibility; one infeasible in some replicates does not, but violates
        a requirement on its share of feasible replicates."""
        study, designs = _partly_infeasible()
        results = study.run(designs)
        share = results[results.design == designs[1].key].feasible.mean()
        reliable = Criteria(
            {"f": "max"}, requirements=[Requirement("feasible", ">=", 1)]
        )
        optuna_study = optimize_optuna(study, reliable, n_trials=12).state
        for trial in optuna_study.trials:
            x = trial.params.get("x")
            feasibility = trial.constraints["feasibility"]
            self.assertEqual(feasibility, 1.0 if x == 2 else 0.0)
            if x == 1:
                self.assertAlmostEqual(trial.constraints["feasible >= 1"], 1 - share)

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
