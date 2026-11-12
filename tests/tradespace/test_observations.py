"""
Unit tests for studies of the observations of constellations.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import dataclasses
import pickle
import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from functools import partial
from pathlib import Path

import pandas as pd

from tatc.analysis import collect_observations
from tatc.generation import generate_points_rotated_lattice
from tatc.schemas import Instrument
from tatc.tradespace import (
    Categorical,
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
    SatelliteCost,
    Stage,
    Study,
    Targets,
    Uncertain,
    build_constellation,
    count_planes,
    count_satellites,
    policy_space,
    value_of,
    walker_space,
)

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
TARGETS = Targets(generate_points_rotated_lattice, N_POINTS, seed="points_seed")
SCHEDULE = DeploymentSchedule(years=3, production_rate=4, launch_interval=0.5)
COSTS = CostModel(SatelliteCost(development=10, first_unit=2))
PLAN = MissionPlan(SCHEDULE, COSTS, growth={"development": "growth"})
FAILURES = FailureModel(life=Exponential(0.2))


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


class TestReferences(unittest.TestCase):
    """
    Unit tests for references to the values of designs.
    """

    def test_value_of(self):
        """References name values, paths of values, or compute them."""
        d = SPACE.enumerate()[0].realize(0)
        self.assertEqual(value_of(d, "altitude"), d.altitude)
        self.assertEqual(value_of(d, "layout.planes"), d.layout.planes)
        self.assertEqual(value_of(d, lambda d: 2 * d.altitude), 2 * d.altitude)
        self.assertIsNone(value_of(d, None))
        with self.assertRaises(KeyError):
            value_of(d, "layout.other")
        with self.assertRaises(TypeError):
            value_of(d, 3)


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
            targets, generate_points_rotated_lattice(N_POINTS, seed=first.points_seed)
        )
        self.assertIs(TARGETS.of(first), targets)  # kept in memory
        self.assertFalse(TARGETS.of(second).equals(targets))
        self.assertEqual(pickle.loads(pickle.dumps(TARGETS))._memory, {})
        with self.assertRaises(ValueError):
            TARGETS.of()

    def test_fixed(self):
        """Fixed targets are the same in every replicate."""
        points = generate_points_rotated_lattice(5, seed=1)
        targets = Targets(points)
        self.assertIs(targets.of(), points)
        self.assertIsNone(targets.seed_of(None))
        same = Targets(generate_points_rotated_lattice, 5, seed=None)
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
        self.assertIn("observed mean", phased)
        self.assertIn("observed time", phased)
        failures = Failures(
            FAILURES, RestorationPolicy(restore=False), histories=4, seed="failure_seed"
        )
        degraded = ObservationStudy(
            SPACE,
            _build,
            _observations(),
            _metrics,
            layout="layout",
            plan=PLAN,
            failures=failures,
        ).run(designs)
        self.assertTrue({"observed mean", "availability", "cost"} <= set(degraded))
        # the policy of each design
        trade = DesignSpace(
            **SPACE.children, policy=policy_space(orbit_spares=Integer(0, 1))
        )
        restored = ObservationStudy(
            trade,
            _build,
            _observations(),
            _metrics,
            layout="layout",
            plan=PLAN,
            failures=Failures(FAILURES, histories=4, seed="failure_seed"),
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
                plan=PLAN,
                failures=Failures(
                    FAILURES,
                    RestorationPolicy(restore=False),
                    histories=4,
                    seed="failure_seed",
                ),
                replicates=2,
                **options,
            )

        scenarios = {"once": {}, "twice": {"times": 2}}
        both = study(1, scenarios=scenarios).run(space(1).enumerate())
        columns = ["observed", "observed mean", "availability", "cost"]
        once = both[both.scenario == "once"].reset_index(drop=True)
        twice = both[both.scenario == "twice"].reset_index(drop=True)
        self.assertTrue((once.observed >= twice.observed).all())
        self.assertTrue((once["observed mean"] > twice["observed mean"]).any())
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
        means = phased.groupby("staged")[["observed mean", "cost"]].mean()
        self.assertLess(
            means.loc[True, "observed mean"], means.loc[False, "observed mean"]
        )
        self.assertLess(means.loc[True, "cost"], means.loc[False, "cost"])
        # failures without restoration: histories are sampled rather than exact
        failures = Failures(
            FAILURES, RestorationPolicy(restore=False), histories=4, seed="failure_seed"
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
        plan = MissionPlan(SCHEDULE, COSTS, stages=stages)
        options = {"layout": "layout", "plan": plan, "replicates": 2}
        failures = Failures(
            FAILURES, RestorationPolicy(restore=False), histories=4, seed="failure_seed"
        )
        for summary in ({"phases": True}, {"failures": failures}):
            with self.subTest(summary=list(summary)):
                results = ObservationStudy(
                    space, _build, _observations(), _metrics, **options, **summary
                ).run(designs)
                means = results.groupby("staging")[["observed mean", "cost"]].mean()
                pd.testing.assert_series_equal(
                    means.loc["measured"], means.loc["delayed"], check_names=False
                )
                pd.testing.assert_series_equal(
                    means.loc["never measured"], means.loc["never"], check_names=False
                )
                self.assertLess(
                    means.loc["never", "cost"], means.loc["delayed", "cost"]
                )

    def test_invalid(self):
        """Failures and phases need a plan, and a plan needs a layout."""
        observations = _observations()
        with self.assertRaises(ValueError):
            ObservationStudy(SPACE, _build, observations, _metrics, phases=True)
        with self.assertRaises(ValueError):
            ObservationStudy(SPACE, _build, observations, _metrics, plan=PLAN)
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
