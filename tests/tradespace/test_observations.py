"""
Unit tests for studies of the observations of constellations.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

# pylint: disable=protected-access

import dataclasses
import importlib.util
import math
import pickle
import tempfile
import unittest
import warnings
from datetime import datetime, timedelta, timezone
from functools import partial
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

from tatc.analysis import collect_observations, evaluate_architecture
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import (
    ActiveWindow,
    Architecture,
    DataFlow,
    DataGeneration,
    GroundStation,
    Instrument,
    Satellite,
    SpacecraftCommunications,
    StationCommunications,
    SunSynchronousOrbit,
    Transmitter,
)
from tatc.tradespace import (
    Categorical,
    Component,
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    ExponentialLifetime,
    FailureModel,
    Failures,
    Fixed,
    Integer,
    LaunchCost,
    Mean,
    MissionPlan,
    Observations,
    ObservationStudy,
    Probability,
    Quantile,
    Real,
    Requirement,
    RestorationPolicy,
    SatelliteCache,
    SatelliteCost,
    Stage,
    StandardDeviation,
    Study,
    Targets,
    Triangular,
    Uncertain,
    build_constellation,
    configure,
    count_planes,
    count_satellites,
    policy_space,
    risk_space,
    walker_space,
)
from tatc.tradespace.expectations import exact_expectations
from tatc.tradespace.fingerprint import get_fingerprint
from tatc.tradespace.histories import summarize_histories
from tatc.tradespace.mission import summarize_phases

HAS_JOBLIB = importlib.util.find_spec("joblib") is not None

START = datetime(2026, 7, 1, tzinfo=timezone.utc)
HOURS = timedelta(hours=3)
N_POINTS = 20
SPACE = DesignSpace(
    altitude=Real(600e3, 700e3, step=100e3),
    layout=walker_space(inclination=60, planes=Integer(1, 2), per_plane=Integer(1, 2)),
    growth=Uncertain.uniform(1.0, 1.5),
    points_seed=Uncertain.seed(),
    failure_seed=Uncertain.seed(),
)
TARGETS = Targets(
    generate_points_fibonacci_lattice, count=N_POINTS, rotate=True, seed="points_seed"
)
SCHEDULE = DeploymentSchedule(years=3, production_rate=4, launch_interval=0.5)
COSTS = CostModel(SatelliteCost(development=10, first_unit=2))
PLAN = MissionPlan(SCHEDULE, COSTS, growth={"development": "growth"})
FAILURES = FailureModel(life=ExponentialLifetime(0.2))
DEGRADE = dataclasses.replace(PLAN, policy=RestorationPolicy(restore=False))


def _build(d):
    imager = Instrument(name="imager", field_of_regard=60)
    return build_constellation(d.layout, imager, epoch=START, altitude=d.altitude)


def _metrics(observations):
    """The share (%) of the points observed."""
    if len(observations) == 0:
        return {"observed": 0.0}
    return {"observed": 100 * observations.target_hash.nunique() / N_POINTS}


def _observations(directory=None):
    return Observations(TARGETS, START, HOURS, directory=directory)


def _grid(count, seed=None):
    """Targets at the equator, of a seed (which a replicate draws)."""
    del seed
    return pd.DataFrame({"latitude": [0.0] * count, "longitude": range(count)})


class TestTargets(unittest.TestCase):
    """
    Unit tests for the targets of observations.
    """

    def test_generated(self):
        """Targets generated from a seed of the design differ by replicate."""
        design = SPACE.enumerate()[0]
        first, second = design.realize(0), design.realize(1)
        targets = TARGETS._of(first)
        pd.testing.assert_frame_equal(
            targets,
            generate_points_fibonacci_lattice(
                count=N_POINTS, rotate=True, seed=first.points_seed
            ),
        )
        self.assertIs(TARGETS._of(first), targets)  # kept in memory
        self.assertFalse(TARGETS._of(second).equals(targets))
        self.assertEqual(pickle.loads(pickle.dumps(TARGETS))._memory, {})
        with self.assertRaises(ValueError):
            TARGETS._of()

    def test_fixed(self):
        """Fixed targets are the same in every replicate."""
        points = generate_points_fibonacci_lattice(count=5, rotate=True, seed=1)
        targets = Targets(points)
        self.assertIs(targets._of(), points)
        self.assertIsNone(targets._seed_of(None))
        same = Targets(generate_points_fibonacci_lattice, count=5, rotate=True)
        self.assertTrue(same._of().equals(same._of()))
        with self.assertRaises(ValueError):
            Targets(points, seed="points_seed")
        with self.assertRaises(TypeError):
            Targets([1, 2])

    def test_fixed_seed(self):
        """An integer seed generates the same targets in every replicate."""
        design = SPACE.enumerate()[0]
        first, second = design.realize(0), design.realize(1)
        targets = Targets(
            generate_points_fibonacci_lattice, count=N_POINTS, rotate=True, seed=7
        )
        expected = generate_points_fibonacci_lattice(
            count=N_POINTS, rotate=True, seed=7
        )
        pd.testing.assert_frame_equal(targets._of(first), expected)
        self.assertIs(targets._of(second), targets._of(first))
        self.assertEqual(targets._seed_of(None), 7)
        pd.testing.assert_frame_equal(targets._of(), expected)
        # observations of a fixed seed need no design
        observations = Observations(targets, START, HOURS)
        satellites = _build(first)
        pd.testing.assert_frame_equal(
            observations(satellites), observations(satellites, second)
        )

    def test_invalid_seeds(self):
        """A seed is a reference, an integer, or None."""
        for seed in (7.5, True, [1], timedelta(hours=1)):
            with self.assertRaisesRegex(TypeError, "seed of targets"):
                Targets(generate_points_fibonacci_lattice, count=5, seed=seed)
        Targets(generate_points_fibonacci_lattice, count=5, seed=np.int64(3))
        Targets(generate_points_fibonacci_lattice, count=5, seed=lambda d: 3)


class TestObservations(unittest.TestCase):
    """
    Unit tests for the observations of targets.
    """

    def test_matches_analysis(self):
        """The observations are those of the analysis of the targets."""
        d = SPACE.enumerate()[3].realize(0)
        satellites = _build(d)
        observations = _observations()(satellites, d)
        direct = collect_observations(TARGETS._of(d), satellites, START, START + HOURS)
        columns = ["target_hash", "satellite", "epoch"]
        pd.testing.assert_frame_equal(
            observations.sort_values(columns, ignore_index=True)[direct.columns],
            direct.sort_values(columns, ignore_index=True),
        )
        self.assertEqual(_observations()._context_of(d), (d.points_seed, START, HOURS))

    def test_disk_cache(self):
        """Observations on disk are kept by the seed of the targets and the period."""
        designs = [d.realize(r) for d in SPACE.enumerate()[:2] for r in (0, 1)]
        with tempfile.TemporaryDirectory() as directory:
            observations = _observations(directory)
            for d in designs:
                observations(_build(d), d)
            namespaces = sorted(
                (p / "namespace.txt").read_text() for p in Path(directory).iterdir()
            )
            self.assertEqual(len(namespaces), 2)  # the targets of two replicates
            self.assertTrue(all("targets of seed" in n for n in namespaces))
            again = _observations(directory)
            again(_build(designs[0]), designs[0])
            self.assertEqual(again._cache_of(designs[0]).misses, 0)

    def test_targets_generated_when_needed(self):
        """Caches on disk of generated targets are keyed by their recipe and
        seed, so that opening (and pruning) them generates no targets, nor do
        results on disk."""
        designs = [d.realize(r) for d in SPACE.enumerate()[:1] for r in (0, 1)]
        original = Targets._of_seed  # pylint: disable=protected-access
        with (
            tempfile.TemporaryDirectory() as directory,
            mock.patch.object(
                Targets, "_of_seed", autospec=True, side_effect=original
            ) as generated,
        ):
            first = _observations(directory)
            for d in designs:
                first(_build(d)[:1], d)
            self.assertEqual(generated.call_count, len(designs))
            generated.reset_mock()
            again = _observations(directory)
            self.assertEqual(again.prune(), 0)
            again(_build(designs[0])[:1], designs[0])  # from the results on disk
            generated.assert_not_called()
            self.assertEqual(
                sorted(SatelliteCache._contexts(directory).values()),
                sorted((d.points_seed, START) for d in designs),
            )
        # the fingerprint of a cache follows the recipe of the targets (their
        # function and its arguments), not the targets kept in memory
        d = designs[0]

        def fingerprint(targets):
            return get_fingerprint(
                Observations(targets, START, HOURS)._cache_of(d).function
            )

        fresh = Targets(_grid, 5, seed="points_seed")
        recipes = [fresh, Targets(_grid, 6, seed="points_seed"), TARGETS]
        prints = [fingerprint(targets) for targets in recipes]
        self.assertEqual(len(set(prints)), 3)
        fresh._of(d)
        self.assertEqual(fingerprint(fresh), prints[0])

    def test_directory_not_in_fingerprint(self):
        """The directory of a cache does not change the fingerprint of a study."""
        prints = {
            ObservationStudy(
                SPACE, _build, _observations(directory), _metrics
            ).fingerprint
            for directory in (None, "a", "b")
        }
        self.assertEqual(len(prints), 1)
        other = Observations(TARGETS, START, 2 * HOURS)
        self.assertNotIn(
            ObservationStudy(SPACE, _build, other, _metrics).fingerprint, prints
        )


class TestObservationStudy(unittest.TestCase):
    """
    Unit tests for studies of observations.
    """

    def test_matches_study(self):
        """A study of observations matches a study of the same parts."""
        designs = SPACE.enumerate()

        def evaluate(satellites, design):
            deployment = SCHEDULE.plan(design.layout)
            observations = collect_observations(
                TARGETS._of(design), satellites, START, START + HOURS
            )
            return {
                "satellites": len(satellites),
                "planes": count_planes(design.layout),
                **COSTS.estimate(deployment, {"development": design.growth}).metrics(),
                **_metrics(observations),
            }

        plain = Study(SPACE, _build, evaluate, replicates=2).run(designs)
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            layout="layout",
            plan=PLAN,
            replicates=2,
        )
        results = study.run(designs)
        columns = [c for c in plain.columns if c != "seconds"]
        pd.testing.assert_frame_equal(results[columns], plain[columns])
        # the parts of a design, for other analyses
        d = designs[3]
        self.assertEqual(len(study.satellites_of(d, 1)), count_satellites(d.layout))
        self.assertEqual(len(study.targets_of(d, 1)), N_POINTS)
        self.assertGreater(len(study.observations_of(d, 1)), 0)

    def test_phases_and_failures(self):
        """Studies summarize phases, or the histories of failures."""
        designs = SPACE.enumerate()[:2]
        phased = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            layout="layout",
            plan=PLAN,
            phases=[Requirement("observed", ">=", 50)],
        ).run(designs)
        self.assertIn("observed mission", phased)
        self.assertIn("observed time", phased)
        # the years by which the plan deploys every satellite
        planned = [SCHEDULE.plan(d.layout).deployed for d in designs for _ in range(1)]
        self.assertEqual(list(phased.deployed), planned)
        failures = Failures(FAILURES, histories=4, seed="failure_seed")
        degraded = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            layout="layout",
            plan=DEGRADE,
            failures=failures,
        ).run(designs)
        self.assertTrue({"observed mission", "availability", "cost"} <= set(degraded))
        # the policy of each design: by default, its values `policy`
        trade = DesignSpace(
            **SPACE.children, policy=policy_space(orbit_spares=Integer(0, 1))
        )
        restored = ObservationStudy(
            trade, _build, _observations(), _metrics, plan=PLAN, failures=failures
        ).run(trade.enumerate()[:4])
        self.assertTrue(restored.feasible.all())
        spares = restored.groupby("policy.orbit_spares")["cost production"].mean()
        self.assertGreater(spares[1], spares[0])

    def test_scenarios(self):
        """Metrics given the design read the values of scenarios, also of the
        phases of missions (cached by scenario)."""

        def observed_at(observations, design):
            """The share (%) of the points observed at least a number of times."""
            if len(observations) == 0:
                return {"observed": 0.0}
            counts = observations.groupby("target_hash").size()
            return {"observed": 100 * (counts >= design.times).sum() / N_POINTS}

        def space(times):
            return DesignSpace(**SPACE.children, times=times)

        def study(times, **options):
            return ObservationStudy(
                space(times),
                _build,
                _observations(),
                observed_at,
                layout="layout",
                plan=DEGRADE,
                failures=Failures(FAILURES, histories=4, seed="failure_seed"),
                replicates=2,
                **options,
            )

        scenarios = {"once": {}, "twice": {"times": 2}}
        both = study(1, scenarios=scenarios).run(space(1).enumerate())
        columns = ["observed", "observed mission", "availability", "cost"]
        once = both[both.scenario == "once"].reset_index(drop=True)
        twice = both[both.scenario == "twice"].reset_index(drop=True)
        self.assertTrue((once.observed >= twice.observed).all())
        self.assertTrue((once["observed mission"] > twice["observed mission"]).any())
        # as studies with the values of the scenarios
        for times, expected in ((1, once), (2, twice)):
            alone = study(times).run(space(times).enumerate())
            pd.testing.assert_frame_equal(
                alone[columns], expected[columns], check_names=False
            )

    def test_epoch(self):
        """The period may start at an epoch of the design, such as one drawn for
        each replicate, whose observations are those of the period at that epoch."""
        later = START + timedelta(days=100)
        space = DesignSpace(**SPACE.children, epoch=Uncertain.choice([START, later]))
        observations = Observations(TARGETS, "epoch", HOURS)
        design = space.enumerate()[3]
        realized = {}
        for replicate in range(10):
            d = design.realize(replicate)
            realized.setdefault(d.epoch, d)
        self.assertEqual(set(realized), {START, later})
        for epoch, d in realized.items():
            self.assertEqual(observations.start_of(d), epoch)
            self.assertEqual(observations.end_of(d), epoch + HOURS)
            self.assertEqual(observations._context_of(d)[1], epoch)
            satellites = _build(d)
            fixed = Observations(TARGETS, epoch, HOURS)(satellites, d)
            pd.testing.assert_frame_equal(observations(satellites, d), fixed)
        satellites = _build(realized[START])
        self.assertFalse(
            observations(satellites, realized[START]).equals(
                observations(satellites, realized[later])
            )
        )
        with self.assertRaises(ValueError):
            observations.start_of(None)
        self.assertEqual(Observations(TARGETS, START, HOURS).end_of(), START + HOURS)
        # a study with an epoch drawn for each replicate
        results = ObservationStudy(
            space, _build, observations, _metrics, layout="layout", replicates=4
        ).run(space.enumerate()[:2])
        self.assertEqual(set(results.epoch), {START, later})

    def test_sizing(self):
        """Plans size the schedule and costs of each design by its values, and
        limit launches by the masses of the plan."""
        designs = [d for d in SPACE.enumerate() if d.layout.per_plane == 2]

        def cost_model(d):
            """Satellites that cost more to develop higher (an illustrative relationship)."""
            development = 10 * d.altitude / 600e3
            return CostModel(
                SatelliteCost(development=development, first_unit=2),
                launch=LaunchCost(per_launch=1),
            )

        def schedule(d):
            """Launches that lift two satellites to 600 km, and one to 700 km."""
            capacity = 250.0 if d.altitude < 650e3 else 150.0
            return dataclasses.replace(SCHEDULE, launch_mass_capacity=capacity)

        def study(plan):
            return ObservationStudy(
                SPACE, _build, _observations(), _metrics, layout="layout", plan=plan
            )

        results = study(MissionPlan(schedule, cost_model, masses=100.0)).run(designs)
        by_altitude = results.groupby("altitude")[["cost development", "cost launch"]]
        means = by_altitude.mean()
        self.assertAlmostEqual(
            means.loc[700e3, "cost development"] / means.loc[600e3, "cost development"],
            7 / 6,
        )
        # a launch of each plane at 600 km, and of each satellite at 700 km
        planes = results.groupby("altitude")["planes"].mean()
        self.assertAlmostEqual(means.loc[600e3, "cost launch"], planes[600e3])
        self.assertAlmostEqual(means.loc[700e3, "cost launch"], 2 * planes[700e3])
        # masses of each design (a reference to its values)
        light = study(MissionPlan(schedule, cost_model, masses=lambda d: 50.0))
        light = light.run(designs)
        self.assertAlmostEqual(
            light.groupby("altitude")["cost launch"].mean()[700e3], planes[700e3]
        )

    def test_cost_risk(self):
        """Plans draw the growth of the components of each design from its
        scores, and resolve growth by type of satellite."""
        growth = Triangular(0.9, 1.5)
        model = CostModel(
            SatelliteCost(
                components=[Component("imager", development=10, growth=growth)]
            ),
            correlation=0.5,
        )
        space = SPACE.extend(cost_risk=risk_space(model))
        designs = [d for d in space.enumerate() if d.layout.per_plane == 2]

        def study(plan):
            return ObservationStudy(
                space,
                _build,
                _observations(),
                _metrics,
                layout="layout",
                plan=plan,
                replicates=2,
            )

        results = study(MissionPlan(SCHEDULE, model, scores="cost_risk")).run(designs)
        score = math.sqrt(0.5) * (
            results["cost_risk.common"] + results["cost_risk.imager"]
        )
        probability = 0.5 * (1 + score.map(lambda z: math.erf(z / math.sqrt(2))))
        expected = 10 * probability.map(growth.ppf)
        pd.testing.assert_series_equal(
            results["cost development"], expected, check_names=False
        )
        # designs of a replicate share the growth (common random numbers)
        self.assertEqual(
            results.groupby("replicate")["cost development"].nunique().max(), 1
        )
        # without scores, the expected growth
        nominal = study(MissionPlan(SCHEDULE, model)).run(designs)
        self.assertTrue((nominal["cost development"] == 10 * growth.mean()).all())
        # growth by type of satellite resolves the values of the design
        plan = MissionPlan(
            SCHEDULE, model, growth={"development": {"a": "growth", "b": np.int64(2)}}
        )
        design = designs[0].realize(0)
        self.assertEqual(
            plan._growth_of(design), {"development": {"a": design.growth, "b": 2}}
        )
        # scores of a group nested in another
        nested = SPACE.extend(risk=DesignSpace(cost_risk=risk_space(model)))
        design = nested.enumerate()[0].realize(3)
        plan = MissionPlan(SCHEDULE, model, scores="risk.cost_risk")
        self.assertEqual(plan._cost_model_of(design), model.draw(design.risk.cost_risk))

    def test_stages(self):
        """Plans stage the deployments of designs that decide to."""
        space = DesignSpace(**SPACE.children, staged=Categorical([False, True]))
        designs = [d for d in space.enumerate() if d.layout.planes == 2]

        def stages(design):
            """A second plane, never deployed, for a staged design."""
            return [Stage(1, 1.0, lambda decision: False)] if design.staged else []

        plan = MissionPlan(
            SCHEDULE, COSTS, growth={"development": "growth"}, stages=stages
        )
        options = {"layout": "layout", "plan": plan, "replicates": 2}
        phased = ObservationStudy(
            space,
            _build,
            _observations(),
            _metrics,
            phases=[Requirement("observed", ">=", 50)],
            **options,
        ).run(designs)
        means = phased.groupby("staged")[["observed mission", "cost"]].mean()
        self.assertLess(
            means.loc[True, "observed mission"], means.loc[False, "observed mission"]
        )
        self.assertLess(means.loc[True, "cost"], means.loc[False, "cost"])
        # the stages deployed, and when the first generation is deployed
        staged = phased.groupby("staged")[["stages deployed", "deployed"]].mean()
        self.assertEqual(staged.loc[True, "stages deployed"], 0.0)
        self.assertTrue(np.isnan(staged.loc[False, "stages deployed"]))  # no stages
        self.assertLess(staged.loc[True, "deployed"], staged.loc[False, "deployed"])
        # failures without restoration: histories are sampled rather than exact
        failures = Failures(FAILURES, seed="failure_seed", histories=4)
        options["plan"] = dataclasses.replace(
            plan, policy=RestorationPolicy(restore=False)
        )
        degraded = ObservationStudy(
            space, _build, _observations(), _metrics, failures=failures, **options
        ).run(designs)
        self.assertTrue(degraded.feasible.all())
        self.assertTrue(any("stages deployed" in column for column in degraded))

    def test_capability_rules(self):
        """Rules of stages read the capability of the satellites operating at a
        decision, by the study's metrics, in phases and in histories."""
        stagings = ["delayed", "measured", "never", "never measured"]
        space = DesignSpace(**SPACE.children, staging=Categorical(stagings))

        def stages(design):
            """A second plane deployed at the first launch, if a rule decides to."""
            rules = {
                "delayed": None,
                "measured": lambda decision: decision.capability["observed"] <= 100,
                "never": lambda decision: False,
                "never measured": lambda decision: decision.capability["observed"] < 0,
            }
            return [Stage(1, 0.0, rules[design.staging])]

        designs = [
            d
            for d in space.enumerate()
            if d.layout.planes == 2 and d.layout.per_plane == 1
        ]
        plan = MissionPlan(
            SCHEDULE, COSTS, stages=stages, policy=RestorationPolicy(restore=False)
        )
        options = {"layout": "layout", "plan": plan, "replicates": 2}
        failures = Failures(FAILURES, histories=4, seed="failure_seed")
        for summary in ({"phases": True}, {"failures": failures}):
            with self.subTest(summary=list(summary)):
                results = ObservationStudy(
                    space, _build, _observations(), _metrics, **options, **summary
                ).run(designs)
                means = results.groupby("staging")[["observed mission", "cost"]].mean()
                pd.testing.assert_series_equal(
                    means.loc["measured"], means.loc["delayed"], check_names=False
                )
                pd.testing.assert_series_equal(
                    means.loc["never measured"], means.loc["never"], check_names=False
                )
                self.assertLess(
                    means.loc["never", "cost"], means.loc["delayed", "cost"]
                )

    def test_metrics_of_uncertain_values(self):
        """The metrics of phases that read uncertain values of the design are
        those of its replicate, and shared by the designs of a replicate."""
        space = SPACE.extend(u=Uncertain.uniform(0, 1))
        fixed = Targets(
            generate_points_fibonacci_lattice(count=N_POINTS, rotate=True, seed=1)
        )

        def metric(observations, design):
            del observations
            return {"m": design.u}

        study = ObservationStudy(
            space,
            _build,
            Observations(fixed, START, HOURS),
            metric,
            plan=PLAN,
            phases=True,
            replicates=3,
        )
        designs = space.enumerate()[:2]
        results = study.run(designs)
        self.assertEqual(results.u.nunique(), 3)
        np.testing.assert_allclose(results["m mission"], results["m"])
        # with common random numbers, the designs of a replicate share a context
        first, second = (study.realize(d, 1) for d in designs)
        evaluation = study.evaluate
        self.assertEqual(evaluation.context_of(first), evaluation.context_of(second))
        other = study.realize(designs[0], 2)
        self.assertNotEqual(evaluation.context_of(first), evaluation.context_of(other))
        # metrics not given the design are cached by the observations alone
        plain = study.replace(metrics=_metrics).evaluate
        self.assertEqual(plain.context_of(first), plain.context_of(other))

    def test_references(self):
        """Every part of a plan may refer to a value of the design, by name or
        by a function of the design."""
        space = SPACE.extend(
            staging=Fixed((Stage(1, 1.0),)), schedule=Fixed(SCHEDULE), mass=Fixed(5.0)
        )
        design = [d for d in space.enumerate() if d.layout.planes == 2][0]
        plan = MissionPlan("schedule", lambda d: COSTS, stages="staging", masses="mass")
        study = ObservationStudy(space, _build, _observations(), _metrics, plan=plan)
        mission = study.mission_of(design)
        self.assertEqual(mission.deployment.stages, (Stage(1, 1.0),))
        self.assertEqual(mission.deployment.masses, {None: 5.0})
        self.assertIs(mission.deployment.schedule, SCHEDULE)
        self.assertEqual(mission.cost_model, COSTS)

    def test_staging_recorded(self):
        """Studies record the stages that rules deploy (without failures), and
        when the first generation is deployed, with them."""
        space = DesignSpace(**SPACE.children, go=Categorical([False, True]))
        designs = [
            d
            for d in space.enumerate()
            if d.layout.planes == 2 and d.layout.per_plane == 1
        ]
        plan = MissionPlan(
            SCHEDULE, COSTS, stages=lambda d: [Stage(1, 1.0, lambda _: d.go)]
        )
        study = ObservationStudy(
            space, _build, _observations(), _metrics, plan=plan, phases=True
        )
        results = study.run(designs)
        means = results.groupby("go")[["stages deployed", "deployed"]].mean()
        self.assertEqual(list(means["stages deployed"]), [0.0, 1.0])
        # the first plane launched at 0.25, the second produced from the
        # decision a year later, and launched at 1.5
        self.assertEqual(list(means["deployed"]), [0.25, 1.5])

    def test_metrics_of_every_satellite_once(self):
        """With a plan, the metrics of every satellite are those of the last
        phase of the mission, evaluated once."""
        calls = []

        def metrics(observations):
            if len(observations) > 0:
                calls.append(frozenset(observations.satellite))
            return _metrics(observations)

        design = SPACE.enumerate()[3]
        study = ObservationStudy(
            SPACE, _build, _observations(), metrics, plan=PLAN, phases=True
        )
        result = study.run([design]).iloc[0]
        self.assertEqual(len(calls), len(set(calls)))
        satellites = study.satellites_of(design)
        observed = study.observations_of(design)
        self.assertEqual(result["observed"], _metrics(observed)["observed"])
        self.assertIn(frozenset(s.name for s in satellites), calls)

    def test_spares_of_the_policy(self):
        """The policy of each design gives its spares, which are produced
        whether or not the study has failures."""
        designs = [d for d in SPACE.enumerate() if d.altitude == 600e3]

        def study(policy):
            plan = dataclasses.replace(PLAN, policy=policy)
            return ObservationStudy(
                SPACE, _build, _observations(), _metrics, plan=plan, phases=True
            )

        plain = study(RestorationPolicy()).run(designs)
        spared = study(RestorationPolicy(orbit_spares=1, ground_spares=1))
        results = spared.run(designs)
        # a spare in orbit in each plane and one on the ground, of 2 each
        np.testing.assert_allclose(
            results["cost production"] - plain["cost production"],
            2 * (plain.planes + 1),
        )
        # (whose production delays the deployment of the planes)
        self.assertTrue(
            (results["observed mission"] <= plain["observed mission"]).all()
        )
        self.assertTrue((results["deployed"] >= plain["deployed"]).all())
        mission = spared.mission_of(designs[-1])
        self.assertEqual({p.orbit_spares for p in mission.deployment.planes}, {1})
        # a policy of each design (or type), from its values
        by_type = study(lambda d: RestorationPolicy(orbit_spares={None: 2}))
        mission = by_type.mission_of(designs[-1])
        self.assertEqual({p.orbit_spares for p in mission.deployment.planes}, {2})
        with self.assertRaises(TypeError):
            study(lambda d: "spares").run(designs[:1])
        # a group of values that configures the default policy, or a base one
        space = SPACE.extend(policy=policy_space(restore=True, orbit_spares=1))
        realized = space.enumerate()[0].realize(0)
        named = dataclasses.replace(PLAN, policy="policy")._policy_of(realized)
        self.assertEqual(named, RestorationPolicy(orbit_spares=1))
        base = RestorationPolicy(drift_years=1.0)
        configured = dataclasses.replace(
            PLAN, policy=lambda d: configure(base, d.policy)
        )._policy_of(realized)
        self.assertEqual(configured, RestorationPolicy(orbit_spares=1, drift_years=1.0))
        # by default, the design's values `policy`, which must configure a
        # policy (others need a function that gives it), or else none
        self.assertEqual(PLAN._policy_of(realized), RestorationPolicy(orbit_spares=1))
        other = SPACE.extend(policy=Categorical(["a", "b"])).enumerate()[0]
        with self.assertRaisesRegex(TypeError, "policy= a function"):
            PLAN._policy_of(other.realize(0))
        unset = SPACE.extend(policy=Fixed(None)).enumerate()[0]
        self.assertEqual(PLAN._policy_of(unset.realize(0)), RestorationPolicy())
        self.assertEqual(PLAN._policy_of(designs[0].realize(0)), RestorationPolicy())

    def test_failure_times(self):
        """Failures give each history the time until it meets requirements on
        its phases, summarized over the histories."""
        failures = Failures(
            FAILURES,
            seed="failure_seed",
            histories=6,
            statistics=[
                Quantile("observed time", 0.5),
                Probability("observed time", "<=", 1),
            ],
        )
        phased = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            plan=PLAN,
            phases=[Requirement("observed", ">=", 50)],
        )
        study = phased.replace(failures=failures)  # the requirements of the study
        design = SPACE.enumerate()[3]
        result = study.run([design]).iloc[0]
        histories = study.mission_of(design).histories()
        self.assertIn("observed time", histories)
        summary = summarize_histories(histories, statistics=failures.statistics)
        for column in (
            "observed time",
            "observed time histories p50",
            "P(observed time <= 1) histories",
        ):
            self.assertEqual(result[column], summary[column])
        self.assertAlmostEqual(
            result["P(observed time <= 1) histories"],
            (histories["observed time"] <= 1).mean(),
        )

    def test_failure_validation(self):
        """Failures are validated as they are made, and quantiles of costs need a
        cost model."""
        for invalid in (
            {"sampling": "antithetic"},
            {"stress": 0.0},
            {"histories": 0},
            {"control_histories": 0},
            {"seed": None},
        ):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    Failures(FAILURES, **{"seed": "failure_seed", **invalid})
        with self.assertRaises(TypeError):
            Failures(FAILURES)  # pylint: disable=missing-kwoa
        with self.assertRaises(TypeError):
            Failures(
                FAILURES, seed="failure_seed", statistics=[Requirement("x", ">", 1)]
            )
        with self.assertRaises(TypeError):
            ObservationStudy(
                SPACE, _build, _observations(), _metrics, plan=PLAN, phases=["x"]
            )
        with self.assertRaisesRegex(ValueError, "cost model"):
            ObservationStudy(
                SPACE,
                _build,
                _observations(),
                _metrics,
                plan=MissionPlan(SCHEDULE),
                failures=Failures(
                    FAILURES, seed="failure_seed", statistics=[Quantile("cost", 0.7)]
                ),
            )
        # quantiles of other columns that the histories do not have
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            plan=PLAN,
            failures=Failures(
                FAILURES, seed="failure_seed", statistics=[Quantile("other", 0.5)]
            ),
        )
        with self.assertRaisesRegex(ValueError, "other"):
            study.run(SPACE.enumerate()[:1])

    def test_replace(self):
        """Variants of an observation study replace its arguments."""
        study = ObservationStudy(SPACE, _build, _observations(), _metrics)
        phased = study.replace(layout="layout", plan=PLAN, phases=True)
        self.assertEqual(phased.phases, ())
        self.assertIs(phased.build, study.build)
        self.assertIs(phased.observations, study.observations)
        results = phased.replace(replicates=2).run(SPACE.enumerate()[:1])
        self.assertEqual(len(results), 2)
        self.assertIn("observed mission", results)
        with self.assertRaises(TypeError):
            study.replace(evaluate=_metrics)

    def test_invalid(self):
        """Failures and phases need a plan, and a plan needs a layout (by
        default, the value `layout` of the design space)."""
        observations = _observations()
        with self.assertRaises(ValueError):
            ObservationStudy(SPACE, _build, observations, _metrics, phases=True)
        # one requirement on the phases of each metric
        with self.assertRaisesRegex(ValueError, r"metrics \['observed'\]"):
            ObservationStudy(
                SPACE,
                _build,
                observations,
                _metrics,
                plan=PLAN,
                phases=[
                    Requirement("observed", ">=", 50),
                    Requirement("observed", ">=", 80),
                ],
            )
        study = ObservationStudy(SPACE, _build, observations, _metrics, plan=PLAN)
        self.assertEqual(study.layout, "layout")
        shells = DesignSpace(
            altitude=SPACE.children["altitude"], shell=SPACE.children["layout"]
        )
        with self.assertRaises(ValueError):
            ObservationStudy(shells, _build, observations, _metrics, plan=PLAN)
        with self.assertRaises(ValueError):
            ObservationStudy(
                SPACE,
                _build,
                observations,
                _metrics,
                failures=Failures(FAILURES, seed="failure_seed"),
            )

    def test_misspelled_references(self):
        """References of the parts of a study that name no value of the
        designs (such as misspelled names) raise when the study is defined,
        each with the names there are; functions are not checked."""
        misspelled = Targets(
            generate_points_fibonacci_lattice, count=N_POINTS, rotate=True, seed="seed"
        )
        cases = {
            "Targets seed 'seed'": (Observations(misspelled, START, HOURS), {}),
            "Observations start 'start'": (Observations(TARGETS, "start", HOURS), {}),
            "MissionPlan kind 'kinds'": (
                _observations(),
                {"plan": dataclasses.replace(PLAN, kind="kinds")},
            ),
            r"MissionPlan growth\['development'\] 'grwth'.*did you mean 'growth'": (
                _observations(),
                {"plan": MissionPlan(SCHEDULE, COSTS, growth={"development": "grwth"})},
            ),
            r"MissionPlan masses\['imager'\] 'mass'": (
                _observations(),
                {"plan": dataclasses.replace(PLAN, masses={"imager": "mass"})},
            ),
            "Failures seed 'failures_seed'.*did you mean 'failure_seed'": (
                _observations(),
                {"plan": PLAN, "failures": Failures(FAILURES, seed="failures_seed")},
            ),
            "ObservationStudy layout 'layout.plane'": (
                _observations(),
                {"layout": "layout.plane"},
            ),
        }
        for message, (observations, options) in cases.items():
            with self.subTest(message), self.assertRaisesRegex(ValueError, message):
                ObservationStudy(SPACE, _build, observations, _metrics, **options)
        # paths, numbers, and functions
        ObservationStudy(
            SPACE,
            _build,
            Observations(TARGETS, lambda d: START, HOURS),
            _metrics,
            plan=dataclasses.replace(
                PLAN,
                kind=lambda d: "imager",
                growth={"development": "growth", "launch": 1.2},
                masses={"imager": 100},
            ),
            layout="layout",
        )
        # of a design space that the study's is replaced by
        study = ObservationStudy(SPACE, _build, _observations(), _metrics)
        with self.assertRaisesRegex(ValueError, "Targets seed 'points_seed'"):
            study.replace(space=DesignSpace(layout=SPACE.children["layout"]))

    def test_requirements_on_phases_of_metrics(self):
        """Requirements on phases are on metrics, not statistics."""
        with self.assertRaisesRegex(ValueError, "not on statistics"):
            ObservationStudy(
                SPACE,
                _build,
                _observations(),
                _metrics,
                plan=PLAN,
                phases=[Requirement(Quantile("observed", 0.1), ">=", 50)],
            )
        ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            plan=PLAN,
            phases=[Requirement(Mean("observed"), ">=", 50)],
        )

    @unittest.skipUnless(HAS_JOBLIB, "requires joblib")
    def test_parallel(self):
        """Studies of observations run in parallel, with a cache on disk."""
        designs = SPACE.enumerate()
        with tempfile.TemporaryDirectory() as directory:
            study = partial(
                ObservationStudy,
                SPACE,
                _build,
                _observations(directory),
                _metrics,
                layout="layout",
                replicates=2,
            )
            serial = study().run(designs)
            parallel = study(n_jobs=2).run(designs)
        columns = [c for c in serial.columns if c != "seconds"]
        pd.testing.assert_frame_equal(serial[columns], parallel[columns])


class TestOneDesign(unittest.TestCase):
    """
    Unit tests for the parts of the evaluations of one design.
    """

    def setUp(self):
        self.study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            layout="layout",
            plan=PLAN,
            phases=True,
            replicates=2,
        )
        self.design = SPACE.enumerate()[3]

    def test_mission(self):
        """The deployment, phases, and costs of a design are its evaluation's."""
        result = self.study.run([self.design]).iloc[1]
        realized = self.study.realize(self.design, 1)
        mission = self.study.mission_of(self.design, 1)
        self.assertIs(mission.design.decisions, realized.decisions)
        self.assertEqual(mission.design.growth, realized.growth)
        deployment = mission.deployment
        self.assertEqual(deployment.planes, SCHEDULE.plan(realized.layout).planes)
        self.assertIs(mission.capability, deployment.capability)
        phases = mission.phases()
        self.assertAlmostEqual(
            summarize_phases(phases)["observed mission"], result["observed mission"]
        )
        self.assertAlmostEqual(mission.estimate().total, result["cost"])
        self.assertEqual(deployment.planned().deployed, result["deployed"])
        self.assertIs(deployment.planned(), deployment.planned())
        self.assertEqual(deployment.policy, RestorationPolicy())
        with self.assertRaisesRegex(ValueError, "failures"):
            mission.histories()
        unplanned = self.study.replace(plan=None, phases=None)
        with self.assertRaisesRegex(ValueError, "plan"):
            unplanned.mission_of(self.design)

    def test_realize(self):
        """A study realizes designs with its seed, and in its scenarios."""
        self.assertEqual(
            self.study.realize(self.design, 1).growth, self.design.realize(1).growth
        )
        self.assertNotEqual(
            self.study.replace(seed=3).realize(self.design, 1).growth,
            self.design.realize(1).growth,
        )
        slow = self.study.replace(scenarios={"slow": {"growth": 1.2}})
        self.assertEqual(slow.realize(self.design, 1, "slow").growth, 1.2)
        self.assertEqual(
            slow.satellites_of(self.design, 1, "slow"),
            slow.satellites_of(self.design, 1),
        )
        with self.assertRaises(KeyError):
            slow.realize(self.design, 1, "fast")

    def test_prune(self):
        """Observations remove their caches' results of the analysis with other
        code, and keep those of other observations."""

        def analysis(source):
            namespace = {"collect_observations": collect_observations}
            exec(source, namespace)  # pylint: disable=exec-used
            return namespace["analysis"]

        plain = analysis(
            "def analysis(targets, satellite, start, end):\n"
            "    return collect_observations(targets, satellite, start, end)\n"
        )
        edited = analysis(
            "def analysis(targets, satellite, start, end):\n"
            "    found = collect_observations(targets, satellite, start, end)\n"
            "    return found\n"
        )
        more = Targets(
            generate_points_fibonacci_lattice,
            count=N_POINTS + 5,
            rotate=True,
            seed="points_seed",
        )
        with tempfile.TemporaryDirectory() as directory:
            first = Observations(TARGETS, START, HOURS, plain, directory=directory)
            others = Observations(more, START, HOURS, plain, directory=directory)
            latest = Observations(TARGETS, START, HOURS, edited, directory=directory)
            realized = self.design.realize(0)
            satellites = _build(realized)
            first(satellites, realized)
            others(satellites, realized)
            with self.assertWarnsRegex(UserWarning, "other code"):
                latest(satellites, realized)
            stale = first._cache_of(realized).folder
            self.assertEqual(len(list(stale.glob("*.pkl"))), len(satellites))
            # a study removes them too, and counts the evaluations it removes
            study = ObservationStudy(SPACE, _build, latest, _metrics)
            self.assertEqual(study.prune(), 0)  # (without a store)
            self.assertEqual(list(stale.glob("*.pkl")), [])
            self.assertEqual(latest.prune(), 0)
            # the results of the other targets (another definition) are kept
            kept = others._cache_of(realized).folder
            self.assertEqual(len(list(kept.glob("*.pkl"))), len(satellites))
            self.assertEqual(Observations(TARGETS, START, HOURS).prune(), 0)

    def test_failure_options(self):
        """Failures sample stressed histories, with controls and spreads."""
        realized = self.design.realize(0)
        default = MissionPlan(SCHEDULE)._policy_of(realized)
        self.assertTrue(default.restore)  # restored by default
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            layout="layout",
            plan=PLAN,
            failures=Failures(
                FAILURES,
                seed="failure_seed",
                histories=6,
                stress=2.0,
                controls=["availability"],
                control_histories=50,
                statistics=[StandardDeviation("observed mission")],
            ),
        )
        result = study.run([self.design]).iloc[0]
        self.assertIn("observed mission histories std", result.index)
        self.assertLess(result["effective histories"], 6)
        self.assertTrue(math.isfinite(result["observed mission"]))
        with self.assertRaisesRegex(ValueError, "controls"):
            study.replace(
                failures=Failures(FAILURES, seed="failure_seed", controls=["other"])
            ).run([self.design])


class TestFailureStudies(unittest.TestCase):
    """
    Unit tests for the summaries of failures of observation studies.
    """

    def test_costs_and_histories(self):
        """The costs of a study of failures are the means of its histories."""
        design = SPACE.enumerate()[3]
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            plan=PLAN,
            failures=Failures(FAILURES, histories=5, seed="failure_seed"),
        )
        result = study.run([design]).iloc[0]
        histories = study.mission_of(design).histories()
        self.assertEqual(len(histories), 5)
        self.assertAlmostEqual(result["cost"], histories["cost"].mean())
        self.assertAlmostEqual(
            result["observed mission"], histories["observed mission"].mean()
        )
        with self.assertRaisesRegex(ValueError, "failures"):
            study.replace(failures=None).mission_of(design).histories()

    def test_exact_or_sampled(self):
        """Without restoration, the means are exact expectations if the sets of
        satellites are few enough, and otherwise those of the histories."""
        design = SPACE.enumerate()[3]
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            plan=DEGRADE,
            failures=Failures(FAILURES, histories=4, seed="failure_seed"),
        )
        mission = study.mission_of(design)
        exact = exact_expectations(mission.deployment, FAILURES, mission.capability)
        sampled = mission.histories()
        result = study.run([design]).iloc[0]
        for column in ("observed mission", "availability"):
            self.assertAlmostEqual(result[column], exact[column])
        bounded = partial(exact_expectations, max_sets=1)
        with mock.patch("tatc.tradespace.plans.exact_expectations", bounded):
            fallback = study.run([design]).iloc[0]
        for column in ("observed mission", "availability"):
            self.assertAlmostEqual(fallback[column], sampled[column].mean())


BASELINE = Architecture(
    satellites=[
        Satellite(
            name="baseline",
            orbit=SunSynchronousOrbit(altitude=700e3, equator_crossing_time="13:30"),
            instruments=[Instrument(name="imager", field_of_regard=60)],
        )
    ],
    ground_stations=[GroundStation(name="svalbard", latitude=78.2, longitude=15.4)],
)


def _latency_metrics(observations):
    """The share of the points observed and the mean latency (minutes)."""
    latency = observations["latency"].mean() if len(observations) > 0 else None
    return {
        **_metrics(observations),
        "latency": None if pd.isna(latency) else latency / timedelta(minutes=1),
    }


class TestArchitectureStudy(unittest.TestCase):
    """
    Unit tests for studies of the observations of architectures.
    """

    def setUp(self):
        self.designs = SPACE.enumerate()[:3]

    def _expected(self, study, design, replicate=0):
        realized = study.realize(design, replicate)
        return evaluate_architecture(
            study.architecture_of(design, replicate),
            TARGETS._of(realized),
            START,
            START + HOURS,
            _latency_metrics,
        )

    def test_baseline(self):
        """A study evaluates the architecture of each design with the baseline."""
        study = ObservationStudy(
            SPACE, _build, _observations(), _latency_metrics, baseline=BASELINE
        )
        results = study.run(self.designs)
        for (_, row), design in zip(results.iterrows(), self.designs):
            architecture = study.architecture_of(design)
            self.assertEqual(architecture.satellites[-1].name, "baseline")
            self.assertEqual(architecture.ground_stations, BASELINE.ground_stations)
            self.assertEqual(row["satellites"], count_satellites(design.layout))
            self.assertEqual(len(study.satellites_of(design)), row["satellites"])
            expected = self._expected(study, design).metrics
            self.assertAlmostEqual(row["observed"], expected["observed"])
            self.assertAlmostEqual(row["latency"], expected["latency"])
        # the same as an architecture that each design builds
        built = ObservationStudy(
            SPACE,
            lambda d: Architecture(satellites=_build(d)).merge(BASELINE),
            _observations(),
            _latency_metrics,
        ).run(self.designs)
        columns = ["observed", "latency"]
        pd.testing.assert_frame_equal(results[columns], built[columns])
        with self.assertRaises(TypeError):
            ObservationStudy(SPACE, _build, _observations(), _metrics, baseline=[])

    def test_phases_keep_baseline(self):
        """The phases of a mission keep the observations of the baseline."""
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            baseline=BASELINE,
            layout="layout",
            plan=PLAN,
            phases=True,
        )
        design = self.designs[2]
        capability = study.mission_of(design).capability
        alone = evaluate_architecture(
            BASELINE,
            TARGETS._of(study.realize(design, 0)),
            START,
            START + HOURS,
            _metrics,
        )
        self.assertGreater(capability([])["observed"], 0)
        self.assertAlmostEqual(capability([])["observed"], alone.metrics["observed"])
        # phases of other baselines are cached apart
        plain = study.replace(baseline=None)
        self.assertEqual(plain.mission_of(design).capability([])["observed"], 0.0)
        with_baseline = study.run([design]).iloc[0]
        without = plain.run([design]).iloc[0]
        self.assertGreaterEqual(
            with_baseline["observed mission"], without["observed mission"]
        )

    def test_evaluation_of(self):
        """The value of each member of a design's architecture."""
        study = ObservationStudy(
            SPACE, _build, _observations(), _latency_metrics, baseline=BASELINE
        )
        design = self.designs[2]
        evaluation = study.evaluation_of(design)
        expected = self._expected(study, design)
        self.assertEqual(evaluation.metrics, expected.metrics)
        values = evaluation.leave_one_out()
        pd.testing.assert_frame_equal(values, expected.leave_one_out())
        self.assertIn("svalbard", values.index)

    def test_windows_and_contacts_cached(self):
        """Members are observed in their windows, and contacts are cached on
        disk, and pruned."""
        baseline = BASELINE.model_copy(
            update={"windows": {"baseline": ActiveWindow(end=START + HOURS / 3)}}
        )
        with tempfile.TemporaryDirectory() as directory:
            observations = Observations(TARGETS, START, HOURS, directory=directory)
            study = ObservationStudy(
                SPACE, _build, observations, _latency_metrics, baseline=baseline
            )
            design = self.designs[0]
            result = study.run([design]).iloc[0]
            self.assertAlmostEqual(
                result["observed"], self._expected(study, design).metrics["observed"]
            )
            seed = study.realize(design, 0).points_seed
            self.assertIn(
                ("observations", seed, START, START + HOURS / 3), observations._caches
            )
            self.assertIn(("contacts", START, START + HOURS), observations._caches)
            contexts = SatelliteCache._contexts(directory)
            self.assertIn(("contacts", START, START + HOURS), contexts.values())
            self.assertIn((seed, START, START + HOURS / 3), contexts.values())
            folder = observations._contacts(START, START + HOURS).folder
            self.assertGreater(len(list(folder.glob("*.pkl"))), 0)
            # a fresh session reads them from disk
            again = Observations(TARGETS, START, HOURS, directory=directory)
            again.evaluate(study.architecture_of(design), study.realize(design, 0))
            self.assertEqual(again._contacts(START, START + HOURS).misses, 0)
            self.assertEqual(observations.prune(), 0)


def _communicating(satellites, rate=150e6):
    """Satellites with a downlink at a rate and data generated at 10 Mbps."""
    return [
        satellite.model_copy(
            update={
                "communications": SpacecraftCommunications(
                    transmitters=[Transmitter(band="X", max_data_rate=rate)],
                    generation=DataGeneration(rate=10e6),
                )
            }
        )
        for satellite in satellites
    ]


class TestDataFlowStudy(unittest.TestCase):
    """
    Unit tests for studies of the data flow of architectures.
    """

    def setUp(self):
        self.design = SPACE.enumerate()[3]
        self.baseline = BASELINE.model_copy(
            update={
                "ground_stations": [
                    BASELINE.ground_stations[0].model_copy(
                        update={"communications": StationCommunications(antennas=1)}
                    )
                ]
            }
        )

    def test_data_flow(self):
        """
        A study evaluates architectures with a data flow model, as
        `evaluate_architecture` does.
        """
        observations = Observations(TARGETS, START, HOURS, data_flow=DataFlow())
        study = ObservationStudy(
            SPACE,
            lambda d: _communicating(_build(d)),
            observations,
            _latency_metrics,
            baseline=self.baseline,
        )
        row = study.run([self.design]).iloc[0]
        realized = study.realize(self.design, 0)
        expected = evaluate_architecture(
            study.architecture_of(self.design),
            TARGETS._of(realized),
            START,
            START + HOURS,
            _latency_metrics,
            data_flow=DataFlow(),
        )
        self.assertAlmostEqual(row["latency"], expected.metrics["latency"])
        evaluation = observations.evaluate(study.architecture_of(self.design), realized)
        self.assertIn("lost", evaluation.observations.columns)
        self.assertIsNotNone(evaluation.flow)

    def test_communications_share_cache(self):
        """
        Architectures that differ only in communications share the cached
        observations and contacts.
        """
        observations = Observations(TARGETS, START, HOURS, data_flow=DataFlow())
        study = ObservationStudy(SPACE, _build, observations, _latency_metrics)
        realized = study.realize(self.design, 0)
        satellites = _build(self.design)
        slow = Architecture(satellites=_communicating(satellites, 10e6)).merge(
            self.baseline
        )
        fast = Architecture(satellites=_communicating(satellites, 300e6)).merge(
            self.baseline
        )
        first = observations.evaluate(slow, realized)
        seed = realized.points_seed
        cache = observations._cache(seed, START)
        contacts = observations._contacts(START, START + HOURS)
        misses = (cache.misses, contacts.misses)
        second = observations.evaluate(fast, realized)
        self.assertEqual((cache.misses, contacts.misses), misses)
        # but the faster downlink delivers more
        self.assertGreater(
            second.flow.summary["delivered"].sum(),
            first.flow.summary["delivered"].sum(),
        )

    def test_phase_context(self):
        """
        With a data flow model, the context of the metrics of phases
        identifies all the planned satellites, which share antennas.
        """
        satellites = _communicating(_build(self.design))
        timing = ObservationStudy(
            SPACE, _build, _observations(), _metrics, layout="layout", plan=PLAN
        )
        flowing = ObservationStudy(
            SPACE,
            _build,
            Observations(TARGETS, START, HOURS, data_flow=DataFlow()),
            _metrics,
            layout="layout",
            plan=PLAN,
        )
        realized = flowing.realize(self.design, 0)
        architecture = Architecture(satellites=satellites)
        # pylint: disable=protected-access
        plain = timing.evaluate.context_of(realized, architecture, satellites)
        context = flowing.evaluate.context_of(realized, architecture, satellites)
        self.assertNotEqual(plain, context)
        fewer = flowing.evaluate.context_of(realized, architecture, satellites[:1])
        self.assertNotEqual(context, fewer)
