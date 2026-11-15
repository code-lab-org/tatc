"""
Unit tests for studies of the observations of constellations.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

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

from tatc.analysis import collect_observations
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import Instrument
from tatc.tradespace import (
    Categorical,
    Component,
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    Exponential,
    FailureModel,
    Failures,
    Integer,
    LaunchCost,
    MissionPlan,
    ObservationStudy,
    Observations,
    Real,
    Requirement,
    RestorationPolicy,
    SatelliteCache,
    SatelliteCost,
    Stage,
    Study,
    Targets,
    Triangular,
    Uncertain,
    build_constellation,
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
from tatc.tradespace import configure

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
FAILURES = FailureModel(life=Exponential(0.2))
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
        targets = TARGETS.of(first)
        pd.testing.assert_frame_equal(
            targets,
            generate_points_fibonacci_lattice(
                count=N_POINTS, rotate=True, seed=first.points_seed
            ),
        )
        self.assertIs(TARGETS.of(first), targets)  # kept in memory
        self.assertFalse(TARGETS.of(second).equals(targets))
        self.assertEqual(pickle.loads(pickle.dumps(TARGETS))._memory, {})
        with self.assertRaises(ValueError):
            TARGETS.of()

    def test_fixed(self):
        """Fixed targets are the same in every replicate."""
        points = generate_points_fibonacci_lattice(count=5, rotate=True, seed=1)
        targets = Targets(points)
        self.assertIs(targets.of(), points)
        self.assertIsNone(targets.seed_of(None))
        same = Targets(generate_points_fibonacci_lattice, count=5, rotate=True)
        self.assertTrue(same.of().equals(same.of()))
        with self.assertRaises(ValueError):
            Targets(points, seed="points_seed")
        with self.assertRaises(TypeError):
            Targets([1, 2])


class TestObservations(unittest.TestCase):
    """
    Unit tests for the observations of targets.
    """

    def test_matches_analysis(self):
        """The observations are those of the analysis of the targets."""
        d = SPACE.enumerate()[3].realize(0)
        satellites = _build(d)
        observations = _observations()(satellites, d)
        direct = collect_observations(TARGETS.of(d), satellites, START, START + HOURS)
        columns = ["target_hash", "satellite", "epoch"]
        pd.testing.assert_frame_equal(
            observations.sort_values(columns, ignore_index=True)[direct.columns],
            direct.sort_values(columns, ignore_index=True),
        )
        self.assertEqual(_observations().context_of(d), (d.points_seed, START, HOURS))

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
            self.assertEqual(again.cache_of(designs[0]).misses, 0)

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
                sorted(SatelliteCache.contexts(directory).values()),
                sorted((d.points_seed, START) for d in designs),
            )
        # the fingerprint of a cache follows the recipe of the targets (their
        # function and its arguments), not the targets kept in memory
        d = designs[0]

        def fingerprint(targets):
            return get_fingerprint(
                Observations(targets, START, HOURS).cache_of(d).function
            )

        fresh = Targets(_grid, 5, seed="points_seed")
        recipes = [fresh, Targets(_grid, 6, seed="points_seed"), TARGETS]
        prints = [fingerprint(targets) for targets in recipes]
        self.assertEqual(len(set(prints)), 3)
        fresh.of(d)
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
                TARGETS.of(design), satellites, START, START + HOURS
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
        # the policy of each design
        trade = DesignSpace(
            **SPACE.children, policy=policy_space(orbit_spares=Integer(0, 1))
        )
        restored = ObservationStudy(
            trade,
            _build,
            _observations(),
            _metrics,
            plan=dataclasses.replace(PLAN, policy="policy"),
            failures=failures,
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
            self.assertEqual(observations.context_of(d)[1], epoch)
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
            _ = observations.end
        with self.assertRaises(ValueError):
            observations.start_of(None)
        self.assertEqual(Observations(TARGETS, START, HOURS).end, START + HOURS)
        # a study with an epoch drawn for each replicate
        results = ObservationStudy(
            space, _build, observations, _metrics, layout="layout", replicates=4
        ).run(space.enumerate()[:2])
        self.assertEqual(set(results.epoch), {START, later})

    def test_sizing(self):
        """Plans size the schedule and costs of each design by its values, and
        limit launches by the masses of the cost model."""
        designs = [d for d in SPACE.enumerate() if d.layout.per_plane == 2]

        def cost_model(d):
            """Satellites that cost more to develop higher (an illustrative relationship)."""
            development = 10 * d.altitude / 600e3
            return CostModel(
                SatelliteCost(development=development, first_unit=2, mass=100),
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

        results = study(MissionPlan(schedule, cost_model)).run(designs)
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
        # given masses take the place of those of the cost model
        light = study(MissionPlan(schedule, cost_model, masses=50.0)).run(designs)
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
            plan.growth_of(design), {"development": {"a": design.growth, "b": 2}}
        )
        # scores of a group nested in another
        nested = SPACE.extend(risk=DesignSpace(cost_risk=risk_space(model)))
        design = nested.enumerate()[0].realize(3)
        plan = MissionPlan(SCHEDULE, model, scores="risk.cost_risk")
        self.assertEqual(plan.cost_model_of(design), model.draw(design.risk.cost_risk))

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
        failures = Failures(FAILURES, histories=4, seed="failure_seed")
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
        named = dataclasses.replace(PLAN, policy="policy").policy_of(realized)
        self.assertEqual(named, RestorationPolicy(orbit_spares=1))
        base = RestorationPolicy(drift_years=1.0)
        configured = dataclasses.replace(
            PLAN, policy=lambda d: configure(base, d.policy)
        ).policy_of(realized)
        self.assertEqual(configured, RestorationPolicy(orbit_spares=1, drift_years=1.0))

    def test_failure_times(self):
        """Failures give each history the time until it meets requirements on
        its phases, summarized over the histories."""
        failures = Failures(
            FAILURES,
            histories=6,
            seed="failure_seed",
            phases=[Requirement("observed", ">=", 50)],
            quantiles={"observed time": 0.5},
            requirements=[Requirement("observed time", "<=", 1)],
        )
        study = ObservationStudy(
            SPACE, _build, _observations(), _metrics, plan=PLAN, failures=failures
        )
        design = SPACE.enumerate()[3]
        result = study.run([design]).iloc[0]
        histories = study.mission_of(design).histories()
        self.assertIn("observed time", histories)
        summary = summarize_histories(
            histories,
            quantiles=failures.quantiles,
            requirements=failures.requirements,
        )
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
            {"launch_stress": -1.0},
            {"histories": 0},
        ):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    Failures(FAILURES, **invalid)
        with self.assertRaises(TypeError):
            Failures(FAILURES, phases=["observed >= 50"])
        with self.assertRaisesRegex(ValueError, "cost model"):
            ObservationStudy(
                SPACE,
                _build,
                _observations(),
                _metrics,
                plan=MissionPlan(SCHEDULE),
                failures=Failures(FAILURES, quantiles={"cost": 0.7}),
            )
        # quantiles of other columns that the histories do not have
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            plan=PLAN,
            failures=Failures(FAILURES, quantiles={"other": 0.5}),
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
                layout="layout",
                plan=PLAN,
                failures=Failures(FAILURES),
                phases=True,
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
        self.assertIs(deployment.capability, mission.capability)
        phases = mission.phases()
        self.assertAlmostEqual(
            summarize_phases(phases)["observed mission"], result["observed mission"]
        )
        self.assertAlmostEqual(mission.estimate().total, result["cost"])
        self.assertEqual(deployment.planned().deployed, result["deployed"])
        self.assertIs(deployment.planned(), deployment.planned())
        self.assertIs(deployment.policy, PLAN.policy)
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
            stale = first.cache_of(realized).folder
            self.assertEqual(len(list(stale.glob("*.pkl"))), len(satellites))
            # a study removes them too, and counts the evaluations it removes
            study = ObservationStudy(SPACE, _build, latest, _metrics)
            self.assertEqual(study.prune(), 0)  # (without a store)
            self.assertEqual(list(stale.glob("*.pkl")), [])
            self.assertEqual(latest.prune(), 0)
            # the results of the other targets (another definition) are kept
            kept = others.cache_of(realized).folder
            self.assertEqual(len(list(kept.glob("*.pkl"))), len(satellites))
            self.assertEqual(Observations(TARGETS, START, HOURS).prune(), 0)

    def test_failure_options(self):
        """Failures sample stressed histories, with controls and spreads."""
        self.assertTrue(MissionPlan(SCHEDULE).policy.restore)  # restored by default
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            layout="layout",
            plan=PLAN,
            failures=Failures(
                FAILURES,
                histories=6,
                seed="failure_seed",
                stress=2.0,
                controls=["availability"],
                control_histories=50,
                spread=True,
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

    def test_exact_without_seed(self):
        """Without a seed, the one history without failures is not replaced by
        exact expectations over failures."""
        design = SPACE.enumerate()[3]
        study = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            plan=DEGRADE,
            failures=Failures(FAILURES),
        )
        result = study.run([design]).iloc[0]
        self.assertEqual(result["availability"], 1.0)
        self.assertEqual(result["failures"], 0.0)
