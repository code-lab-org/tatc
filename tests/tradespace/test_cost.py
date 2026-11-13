"""
Unit tests for life-cycle costs of constellations.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import dataclasses
import math
import unittest

import numpy as np

from tatc.tradespace import (
    Categorical,
    Component,
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    Integer,
    LaunchCost,
    OperationsCost,
    Quantile,
    SatelliteCost,
    Study,
    Triangular,
    Uncertain,
    aggregate,
    count_satellites,
    exceedance,
    planes_space,
    risk_space,
    sso_planes_space,
    walker_space,
)
from tatc.tradespace.constellations import plane_sizes


def _design(**children):
    """The only design of a design space."""
    (design,) = DesignSpace(**children).enumerate()
    return design


def _walker(planes=3, per_plane=4):
    return _design(
        layout=walker_space(
            inclination=55.0, planes=planes, per_plane=per_plane, relative_spacing=0
        )
    ).layout


class TestSatelliteCost(unittest.TestCase):
    """
    Unit tests for costs of satellites.
    """

    def test_learning_curve(self):
        """Each doubling of units costs the learning rate times as much per unit."""
        flat = SatelliteCost(first_unit=10)
        self.assertEqual(flat.production(5), 50)
        learning = SatelliteCost(first_unit=10, learning_rate=0.9)
        self.assertAlmostEqual(learning.production(1), 10)
        self.assertAlmostEqual(learning.production(2), 10 + 9)
        self.assertAlmostEqual(learning.production(4) - learning.production(3), 8.1)
        self.assertEqual(learning.production(0), 0)
        with self.assertRaises(ValueError):
            SatelliteCost(first_unit=10, learning_rate=1.2)
        with self.assertRaises(ValueError):
            SatelliteCost(development=-1)


class TestPlaneSizes(unittest.TestCase):
    """
    Unit tests for the satellites of each plane of layouts.
    """

    def test_layouts(self):
        """Plane sizes of each layout, and of several."""
        walker = _walker(planes=2, per_plane=3)
        sso = _design(
            layout=sso_planes_space(crossing=13.5, satellites=2, count=1)
        ).layout
        inclined = _design(
            layout=planes_space(inclination=30.0, node=0.0, satellites=4, count=2)
        ).layout
        self.assertEqual(plane_sizes(walker), [3, 3])
        self.assertEqual(plane_sizes(sso), [2])
        self.assertEqual(plane_sizes(inclined), [4, 4])
        self.assertEqual(plane_sizes([walker, sso]), [3, 3, 2])


class TestCostModel(unittest.TestCase):
    """
    Unit tests for life-cycle costs of deployments.
    """

    def setUp(self):
        self.satellite = SatelliteCost(
            development=100, first_unit=10, learning_rate=0.9, mass=50
        )

    def test_elements(self):
        """Costs of development, production, launches, and operations."""
        model = CostModel(
            self.satellite,
            launch=LaunchCost(per_launch=30),
            operations=OperationsCost(per_year=5, per_satellite_year=1),
        )
        schedule = DeploymentSchedule(years=6, design_life=3, launch_capacity=3)
        estimate = model.estimate(schedule.plan(_walker(planes=3, per_plane=4)))
        # two generations of 3 planes of 4 satellites, 2 launches each
        self.assertEqual(estimate.units, {"satellites": 24})
        self.assertEqual(estimate.launches, 12)
        self.assertEqual(estimate.development, 100)
        self.assertAlmostEqual(estimate.production, self.satellite.production(24))
        self.assertEqual(estimate.launch, 12 * 30)
        self.assertEqual(estimate.operations, 6 * (5 + 12))
        self.assertAlmostEqual(
            estimate.total, 100 + self.satellite.production(24) + 360 + 102
        )
        # without discounting, the present value is the total
        self.assertAlmostEqual(estimate.present_value, estimate.total)
        metrics = estimate.metrics()
        self.assertEqual(
            list(metrics),
            [
                "cost",
                "cost development",
                "cost production",
                "cost launch",
                "cost operations",
                "cost present",
                "cost peak annual",
            ],
        )
        self.assertEqual(metrics["cost"], estimate.total)

    def test_rideshare(self):
        """Rideshares are priced per satellite and per kg."""
        model = CostModel(
            self.satellite, launch=LaunchCost(per_satellite=2, per_kg=0.1)
        )
        estimate = model.estimate(DeploymentSchedule().plan(_walker(1, 4)))
        self.assertEqual(estimate.launches, 0)
        self.assertAlmostEqual(estimate.launch, 4 * (2 + 0.1 * 50))
        with self.assertRaises(ValueError):
            CostModel(SatelliteCost(), launch=LaunchCost(per_kg=1)).estimate(
                DeploymentSchedule().plan(_walker())
            )

    def test_spares(self):
        """Spares in orbit are built, launched, and operated with each plane."""
        model = CostModel(
            self.satellite,
            launch=LaunchCost(per_launch=30),
            operations=OperationsCost(per_satellite_year=1),
        )
        schedule = DeploymentSchedule(launch_capacity=4)
        estimate = model.estimate(schedule.plan(_walker(2, 4), spares=1))
        self.assertEqual(estimate.units["satellites"], 10)
        self.assertEqual(estimate.launches, 4)
        self.assertEqual(estimate.operations, 10)

    def test_types(self):
        """Types of satellites are developed once each, with their own curves."""
        small = SatelliteCost(development=10, first_unit=1)
        model = CostModel({"large": self.satellite, "small": small})
        sso = _design(
            layout=sso_planes_space(crossing=13.5, satellites=3, count=1)
        ).layout
        mixed = model.estimate(
            DeploymentSchedule().plan({"large": _walker(1, 2), "small": sso})
        )
        self.assertEqual(mixed.units, {"large": 2, "small": 3})
        self.assertEqual(mixed.development, 110)
        self.assertAlmostEqual(mixed.production, self.satellite.production(2) + 3)
        one = model.estimate(DeploymentSchedule().plan(_walker(1, 2), kind="small"))
        self.assertEqual(one.development, 10)
        with self.assertRaises(ValueError):
            model.estimate(DeploymentSchedule().plan(_walker(1, 2)))
        with self.assertRaises(ValueError):
            model.estimate(DeploymentSchedule().plan(_walker(1, 2), kind="medium"))

    def test_growth(self):
        """Growth factors multiply the cost, or its elements."""
        model = CostModel(
            self.satellite, operations=OperationsCost(per_year=5), discount_rate=0.05
        )
        deployment = DeploymentSchedule(years=2, development_years=2).plan(_walker())
        base = model.estimate(deployment)
        grown = model.estimate(deployment, growth=1.5)
        self.assertAlmostEqual(grown.total, 1.5 * base.total)
        partial = model.estimate(deployment, growth={"development": 2.0})
        self.assertAlmostEqual(partial.development, 200)
        self.assertAlmostEqual(partial.operations, base.operations)
        self.assertAlmostEqual(
            partial.present["development"], 2 * base.present["development"]
        )
        with self.assertRaises(ValueError):
            model.estimate(deployment, growth={"insurance": 1.1})
        with self.assertRaises(ValueError):
            CostModel(self.satellite, discount_rate=-1)

    def test_trade_study(self):
        """Costs are constraints of design spaces and metrics of studies."""
        model = CostModel(
            {
                "narrow": SatelliteCost(development=50, first_unit=10),
                "wide": SatelliteCost(development=80, first_unit=16),
            },
            launch=LaunchCost(per_launch=20),
        )
        schedule = DeploymentSchedule()
        budget = 200
        space = DesignSpace(
            imager=Categorical(["narrow", "wide"]),
            layout=walker_space(
                inclination=55.0,
                planes=Integer(1, 3),
                per_plane=Integer(1, 3),
                relative_spacing=0,
            ),
            growth=Uncertain(lambda rng: float(rng.triangular(0.9, 1.0, 1.5))),
            constraints=[
                lambda d: model.estimate(schedule.plan(d.layout, kind=d.imager)).total
                - budget
            ],
        )
        designs = space.enumerate()
        nominal = [
            model.estimate(schedule.plan(d.layout, kind=d.imager)).total
            for d in designs
        ]
        self.assertTrue(all(cost <= budget for cost in nominal))
        self.assertLess(len(designs), 2 * 9)

        def evaluate(d):
            deployment = schedule.plan(d.layout, kind=d.imager)
            cost = model.estimate(deployment, growth=d.growth)
            return {**cost.metrics(), "satellites": count_satellites(d.layout)}

        results = Study(space, lambda d: d, evaluate, replicates=20).run(designs)
        summary = aggregate(
            results, statistics={"cost": "cost", "cost p70": Quantile("cost", 0.7)}
        )
        # the 0.7 quantile of a growth skewed upward exceeds the mean
        self.assertTrue((summary["cost p70"] > summary["cost"]).all())
        ratio = summary["cost"] / np.array(nominal)
        np.testing.assert_allclose(ratio, ratio.iloc[0])  # common random numbers


class TestDatedCosts(unittest.TestCase):
    """
    Unit tests for dated costs and their present value.
    """

    def setUp(self):
        self.terms = {
            "satellites": SatelliteCost(
                development=90, first_unit=10, learning_rate=0.9
            ),
            "launch": LaunchCost(per_launch=30),
            "operations": OperationsCost(per_year=5, per_satellite_year=1),
        }
        self.layout = _walker(planes=3, per_plane=4)

    def _estimate(self, discount_rate=0.1, start=0.0, **schedule):
        terms = {"years": 6, "design_life": 3, "launch_capacity": 3, **schedule}
        model = CostModel(**self.terms, discount_rate=discount_rate, start=start)
        return model.estimate(DeploymentSchedule(**terms).plan(self.layout))

    def test_present_value(self):
        """Costs are discounted from the present to their dates."""
        dated = self._estimate(development_years=3)
        # development of 90, evenly over 3 years, dated at mid-year
        expected = sum(30 / 1.1 ** (year + 0.5) for year in range(3))
        self.assertAlmostEqual(dated.present["development"], expected)
        # launches of the first generation at year 3, the second at year 6
        self.assertAlmostEqual(dated.present["launch"], 180 / 1.1**3 + 180 / 1.1**6)
        later = self._estimate(development_years=3, start=2)
        self.assertAlmostEqual(later.present_value, dated.present_value / 1.1**2)
        self.assertLess(dated.metrics()["cost present"], dated.total)

    def test_production(self):
        """Production is paid while produced, along the learning curve."""
        dated = self._estimate(development_years=2, production_rate=8)
        production = [f for f in dated.flows if f[1] == "production"]
        self.assertAlmostEqual(production[0][0], 2.25)  # middle of [2, 2.5)
        self.assertAlmostEqual(production[0][2], self.terms["satellites"].production(4))

    def test_operations_from_launch(self):
        """Satellites launched later are operated for less of the mission."""
        dated = self._estimate(development_years=2, launch_interval=1)
        launches = [t for t, element, _ in dated.flows if element == "launch"]
        self.assertEqual(launches, [2, 3, 4, 5, 5, 6, 6, 7, 7])
        expected = 5 * 6 + sum(2 * (8 - t) for t in [2, 3, 4, 5, 6, 7])
        self.assertAlmostEqual(dated.operations, expected)

    def test_annual(self):
        """Annual costs add to the total, and their peak is a metric."""
        dated = self._estimate(development_years=2, launch_interval=1)
        annual = dated.annual()
        self.assertAlmostEqual(annual["total"].sum(), dated.total)
        self.assertEqual(list(annual.index), list(range(0, 8)))
        self.assertAlmostEqual(
            dated.metrics()["cost peak annual"], annual["total"].max()
        )


NEW = Triangular(0.827, 2.274)
"""The growth of new components (Maier et al., 2023, Table I)."""
EXISTING = Triangular(0.981, 1.140)
"""The growth of existing components (Maier et al., 2023, Table I)."""


def _normal_probability(score):
    return 0.5 * (1 + math.erf(score / math.sqrt(2)))


class TestTriangular(unittest.TestCase):
    """
    Unit tests for triangular distributions of growth.
    """

    def test_moments(self):
        """The mean and standard deviation of a triangular distribution."""
        self.assertAlmostEqual(NEW.mean(), 1.367, places=3)
        self.assertAlmostEqual(EXISTING.mean(), 1.0403, places=4)
        rng = np.random.default_rng(1)
        draws = rng.triangular(0.827, 1.0, 2.274, 200000)
        self.assertAlmostEqual(NEW.std(), draws.std(), places=2)

    def test_quantiles(self):
        """Quantiles run from the lowest to the highest growth, through the mode."""
        self.assertAlmostEqual(NEW.ppf(0.0), 0.827)
        self.assertAlmostEqual(NEW.ppf(1.0), 2.274)
        self.assertAlmostEqual(NEW.ppf((1 - 0.827) / (2.274 - 0.827)), 1.0)
        draws = [NEW.ppf(q) for q in np.random.default_rng(2).random(100000)]
        self.assertAlmostEqual(np.mean(draws), NEW.mean(), places=2)
        with self.assertRaises(ValueError):
            Triangular(1.2, 1.5)  # mode below the lowest growth
        with self.assertRaises(ValueError):
            Triangular(1.0, 1.0)


class TestComponents(unittest.TestCase):
    """
    Unit tests for the costs of the components of satellites.
    """

    def setUp(self):
        self.satellite = SatelliteCost(
            development=10,
            first_unit=1,
            components=[
                Component("bus", development=20, first_unit=4),
                Component("imager", development=30, first_unit=6, growth=NEW),
            ],
        )

    def test_expected_costs(self):
        """Without drawn growth, components cost their expected growth."""
        expected_development = 10 + 20 + 30 * NEW.mean()
        expected_first_unit = 1 + 4 + 6 * NEW.mean()
        self.assertAlmostEqual(self.satellite.total_development, expected_development)
        self.assertAlmostEqual(self.satellite.total_first_unit, expected_first_unit)
        self.assertAlmostEqual(self.satellite.production(2), 2 * expected_first_unit)
        estimate = CostModel(self.satellite).estimate(
            DeploymentSchedule().plan(_walker(1, 2))
        )
        self.assertAlmostEqual(estimate.development, expected_development)
        self.assertAlmostEqual(estimate.production, 2 * expected_first_unit)
        # a satellite without components costs as before
        self.assertEqual(SatelliteCost(development=5).total_development, 5)

    def test_invalid(self):
        """Components need names, costs that are not negative, and growth."""
        with self.assertRaises(ValueError):
            Component("")
        with self.assertRaises(ValueError):
            Component("bus", development=-1)
        with self.assertRaises(ValueError):
            Component("bus", growth=-0.5)
        with self.assertRaises(TypeError):
            Component("bus", growth="new")
        with self.assertRaises(ValueError):
            SatelliteCost(components=[Component("bus"), Component("bus")])
        with self.assertRaises(ValueError):
            CostModel(self.satellite, correlation=1.5)
        with self.assertRaises(TypeError):
            SatelliteCost(components=["bus"])
        # the name of the common score is not a component's
        with self.assertRaises(ValueError):
            CostModel(SatelliteCost(components=[Component("common", growth=NEW)]))
        # types that share a component share its costs and growth
        other = SatelliteCost(
            components=[Component("bus", development=25, first_unit=4)]
        )
        with self.assertRaises(ValueError):
            CostModel({"a": self.satellite, "b": other})

    def test_score_names(self):
        """Scores are of uncertain components, after the common score if they
        are correlated."""
        wide = SatelliteCost(
            components=[
                Component("bus", development=20, first_unit=4),
                Component("wide imager", development=50, growth=NEW),
            ]
        )
        model = CostModel({"narrow": self.satellite, "wide": wide}, correlation=0.5)
        self.assertEqual(model.score_names(), ["common", "imager", "wide imager"])
        self.assertEqual(CostModel(self.satellite).score_names(), ["imager"])
        self.assertEqual(CostModel(SatelliteCost(), correlation=0.5).score_names(), [])

    def test_draw(self):
        """Growth is the quantile at the probability of the scores of the
        component and the common score, mixed by the correlation."""
        model = CostModel({"narrow": self.satellite}, correlation=0.5)
        drawn = model.draw({"common": 1.0, "imager": -0.5, "other": 3.0})
        imager = drawn.satellites["narrow"].components[1]
        score = math.sqrt(0.5) * 1.0 + math.sqrt(0.5) * -0.5
        self.assertAlmostEqual(imager.growth, NEW.ppf(_normal_probability(score)))
        self.assertFalse(imager.uncertain)
        # components of fixed growth, and the model, are unchanged
        self.assertEqual(drawn.satellites["narrow"].components[0].growth, 1.0)
        self.assertTrue(model.satellites["narrow"].components[1].uncertain)
        self.assertAlmostEqual(
            drawn.satellites["narrow"].total_development, 10 + 20 + 30 * imager.growth
        )
        # without correlation, the component's score alone
        independent = CostModel(
            SatelliteCost(components=[Component("imager", 30, growth=NEW)])
        )
        self.assertEqual(independent.score_names(), ["imager"])
        growth = independent.draw({"imager": 0.0}).satellites.components[0].growth
        self.assertAlmostEqual(growth, NEW.ppf(0.5))
        with self.assertRaises(ValueError):
            model.draw({"imager": 0.0})  # no common score

    def test_correlation(self):
        """Drawn growths are correlated between every pair of components, and a
        component that types share has the same growth in each."""
        bus = Component("bus", 1, growth=NEW)
        model = CostModel(
            {
                "a": SatelliteCost(
                    components=[bus, Component("imager", 1, growth=NEW)]
                ),
                "b": SatelliteCost(
                    components=[bus, Component("sounder", 1, growth=NEW)]
                ),
            },
            correlation=0.5,
        )
        rng = np.random.default_rng(3)
        draws = []
        for _ in range(4000):
            scores = {name: rng.normal() for name in model.score_names()}
            drawn = model.draw(scores).satellites
            draws.append(
                [c.growth for c in drawn["a"].components]
                + [c.growth for c in drawn["b"].components]
            )
        draws = np.array(draws)
        np.testing.assert_array_equal(draws[:, 0], draws[:, 2])
        r = np.corrcoef(draws.T)
        self.assertAlmostEqual(r[0, 1], 0.5, delta=0.06)
        self.assertAlmostEqual(r[1, 3], 0.5, delta=0.06)
        independent = dataclasses.replace(model, correlation=0.0)
        draws = []
        for _ in range(4000):
            scores = {name: rng.normal() for name in independent.score_names()}
            drawn = independent.draw(scores).satellites
            draws.append(
                [c.growth for c in drawn["a"].components + drawn["b"].components]
            )
        r = np.corrcoef(np.array(draws).T)
        self.assertAlmostEqual(r[1, 3], 0.0, delta=0.06)

    def test_shared_development(self):
        """A component that types of a deployment share is developed once, and
        produced with each type."""
        bus = Component("bus", development=20, first_unit=4)
        model = CostModel(
            {
                "a": SatelliteCost(
                    development=5, components=[bus, Component("imager", 30, 6)]
                ),
                "b": SatelliteCost(
                    development=7, components=[bus, Component("sounder", 10, 2)]
                ),
            }
        )
        sso = _design(
            layout=sso_planes_space(crossing=13.5, satellites=3, count=1)
        ).layout
        mixed = model.estimate(
            DeploymentSchedule().plan({"a": _walker(1, 2), "b": sso})
        )
        self.assertAlmostEqual(mixed.development, 5 + 30 + 7 + 10 + 20)
        self.assertAlmostEqual(mixed.production, 2 * (4 + 6) + 3 * (4 + 2))
        alone = model.estimate(DeploymentSchedule().plan(_walker(1, 2), kind="b"))
        self.assertAlmostEqual(alone.development, 7 + 10 + 20)


class TestGrowthByType(unittest.TestCase):
    """
    Unit tests for growth factors of each type of satellite.
    """

    def setUp(self):
        self.model = CostModel(
            {
                "narrow": SatelliteCost(development=50, first_unit=10),
                "wide": SatelliteCost(development=80, first_unit=16),
            },
            launch=LaunchCost(per_launch=20),
            operations=OperationsCost(per_year=1),
        )
        sso = _design(
            layout=sso_planes_space(crossing=13.5, satellites=2, count=1)
        ).layout
        self.deployment = DeploymentSchedule().plan(
            {"narrow": _walker(1, 2), "wide": sso}
        )

    def test_by_type(self):
        """Each type of satellite grows by its own factor."""
        base = self.model.estimate(self.deployment)
        grown = self.model.estimate(
            self.deployment,
            growth={
                "development": {"narrow": 2.0, "wide": 1.0},
                "production": {"narrow": 1.0, "wide": 1.5},
                "launch": {"narrow": 1.0, "wide": 3.0},
                "operations": 2.0,
            },
        )
        self.assertAlmostEqual(grown.development, 2 * 50 + 80)
        self.assertAlmostEqual(grown.production, 20 + 1.5 * 32)
        self.assertAlmostEqual(grown.launch, 20 + 3 * 20)
        self.assertAlmostEqual(grown.operations, 2 * base.operations)

    def test_invalid(self):
        """Growth by type needs every kind of the deployment, and is not of
        operations."""
        with self.assertRaises(ValueError):
            self.model.estimate(self.deployment, growth={"development": {"narrow": 2}})
        with self.assertRaises(ValueError):
            self.model.estimate(
                self.deployment, growth={"operations": {"narrow": 1, "wide": 1}}
            )
        single = CostModel(SatelliteCost(development=10))
        with self.assertRaises(ValueError):
            single.estimate(
                DeploymentSchedule().plan(_walker(1, 2)),
                growth={"development": {"satellites": 2}},
            )

    def test_shared_types(self):
        """Common random numbers correlate the costs of designs of the same type
        (as in Maier et al., 2023): the order of such designs is settled, and
        that of designs of different types is not."""
        space = DesignSpace(
            imager=Categorical(["narrow", "wide"]),
            layout=walker_space(
                inclination=55.0, planes=1, per_plane=Integer(2, 3), relative_spacing=0
            ),
            narrow=Uncertain.normal(1.0, 0.3),
            wide=Uncertain.normal(1.0, 0.3),
        )
        schedule = DeploymentSchedule()

        def evaluate(d):
            deployment = schedule.plan(d.layout, kind=d.imager)
            growth = {"narrow": d.narrow, "wide": d.wide}
            return self.model.estimate(
                deployment, growth={"development": growth, "production": growth}
            ).metrics()

        designs = space.enumerate()

        def compare(common):
            study = Study(
                space,
                lambda d: d,
                evaluate,
                replicates=200,
                common_random_numbers=common,
            )
            compared = exceedance(study.run(designs), "cost", designs)
            return compared.set_index(["design", "other"])["exceedance"]

        shared, independent = compare(True), compare(False)
        narrow2, narrow3 = (d.key for d in designs if d.imager == "narrow")
        wide2 = next(d.key for d in designs if d.imager == "wide")
        # three narrow satellites cost more than two in every replicate
        self.assertEqual(shared[(narrow3, narrow2)], 1.0)
        self.assertLess(independent[(narrow3, narrow2)], 0.99)
        # designs of different types are not ordered by their shared draws
        self.assertLess(shared[(wide2, narrow3)], 1.0)
        self.assertGreater(shared[(wide2, narrow3)], 0.0)


class TestRiskSpace(unittest.TestCase):
    """
    Unit tests for the scores of the cost risk of designs.
    """

    def setUp(self):
        self.model = CostModel(
            {
                30: SatelliteCost(
                    components=[
                        Component("bus", 20, 4, growth=EXISTING),
                        Component("imager 30", 15, 3, growth=NEW),
                    ],
                ),
                60: SatelliteCost(
                    components=[
                        Component("bus", 20, 4, growth=EXISTING),
                        Component("imager 60", 30, 6, growth=NEW),
                    ],
                ),
            },
            correlation=0.5,
        )

    def test_scores(self):
        """A standard normal score of each name, the same in every design of a
        replicate."""
        group = risk_space(self.model)
        self.assertEqual(
            list(group.children), ["common", "bus", "imager 30", "imager 60"]
        )
        space = DesignSpace(
            imager=Categorical([30, 60]),
            layout=walker_space(
                inclination=55.0, planes=1, per_plane=2, relative_spacing=0
            ),
            cost_risk=group,
        )
        self.assertIn("cost_risk.imager 30", space.uncertain_paths())
        narrow, wide = space.enumerate()
        a, b = narrow.realize(4), wide.realize(4)
        self.assertEqual(a.cost_risk["bus"], b.cost_risk["bus"])
        # every design of a replicate draws the same growths, and the bus grows
        # alike in both types
        drawn = self.model.draw(a.cost_risk)
        self.assertEqual(drawn, self.model.draw(b.cost_risk))
        self.assertEqual(
            drawn.satellites[30].components[0], drawn.satellites[60].components[0]
        )
        scores = np.array([narrow.realize(r).cost_risk["bus"] for r in range(500)])
        self.assertAlmostEqual(scores.mean(), 0.0, delta=0.15)
        self.assertAlmostEqual(scores.std(), 1.0, delta=0.15)

    def test_invalid(self):
        """Names are not reserved, nor have characters of paths."""
        reserved = CostModel(SatelliteCost(components=[Component("count", growth=NEW)]))
        with self.assertRaisesRegex(ValueError, "count"):
            risk_space(reserved)
        dotted = CostModel(
            SatelliteCost(components=[Component("imager.v2", growth=NEW)])
        )
        with self.assertRaises(ValueError):
            risk_space(dotted)
