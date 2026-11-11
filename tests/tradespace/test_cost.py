"""
Unit tests for life-cycle costs of constellations.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np

from tatc.tradespace import (
    Categorical,
    CostModel,
    DeploymentSchedule,
    DesignSpace,
    Integer,
    LaunchCost,
    OperationsCost,
    Quantile,
    SatelliteCost,
    Study,
    Uncertain,
    aggregate,
    count_satellites,
    plane_sizes,
    planes_space,
    sso_planes_space,
    walker_space,
)


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
            results, {"cost": "cost", "cost p70": Quantile("cost", 0.7)}
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
