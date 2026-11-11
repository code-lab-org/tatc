"""
Unit tests for life-cycle costs of constellations.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np

from tatc.tradespace import (
    Categorical,
    CostModel,
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
    Unit tests for life-cycle costs.
    """

    def setUp(self):
        self.satellite = SatelliteCost(
            development=100, first_unit=10, learning_rate=0.9, mass=50
        )

    def test_dedicated_launches(self):
        """Planes are launched in as many launches as their capacity requires."""
        model = CostModel(
            self.satellite,
            launch=LaunchCost(per_launch=30, capacity=3),
            operations=OperationsCost(per_year=5, per_satellite_year=1),
            years=6,
            design_life=3,
        )
        self.assertEqual(model.generations, 2)
        estimate = model.estimate(_walker(planes=3, per_plane=4))
        # two generations of 3 planes of 4 satellites, 2 launches each
        self.assertEqual(estimate.units, {"satellites": 24})
        self.assertEqual(estimate.launches, 12)
        self.assertEqual(estimate.development, 100)
        self.assertAlmostEqual(estimate.production, self.satellite.production(24))
        self.assertEqual(estimate.launch, 12 * 30)
        self.assertEqual(estimate.operations, 6 * (5 + 12))
        self.assertAlmostEqual(
            estimate.total,
            100 + self.satellite.production(24) + 360 + 102,
        )
        metrics = estimate.metrics()
        self.assertEqual(
            list(metrics),
            [
                "cost",
                "cost development",
                "cost production",
                "cost launch",
                "cost operations",
            ],
        )
        self.assertEqual(metrics["cost"], estimate.total)

    def test_rideshare(self):
        """Rideshares are priced per satellite and per kg."""
        model = CostModel(
            self.satellite, launch=LaunchCost(per_satellite=2, per_kg=0.1)
        )
        estimate = model.estimate(_walker(planes=1, per_plane=4))
        self.assertEqual(estimate.launches, 0)
        self.assertAlmostEqual(estimate.launch, 4 * (2 + 0.1 * 50))
        with self.assertRaises(ValueError):
            CostModel(SatelliteCost(), launch=LaunchCost(per_kg=1)).estimate(_walker())

    def test_mass_capacity(self):
        """A launch carries as many satellites as its mass capacity allows."""
        model = CostModel(
            self.satellite, launch=LaunchCost(per_launch=30, mass_capacity=120)
        )
        # 2 satellites of 50 kg per launch: 2 launches for each plane of 4
        self.assertEqual(model.estimate(_walker(planes=2, per_plane=4)).launches, 4)
        with self.assertRaises(ValueError):
            CostModel(
                self.satellite, launch=LaunchCost(per_launch=30, mass_capacity=40)
            ).estimate(_walker())

    def test_spares(self):
        """Spares in orbit are built, launched, and operated with each plane."""
        model = CostModel(
            self.satellite,
            launch=LaunchCost(per_launch=30, capacity=4),
            operations=OperationsCost(per_satellite_year=1),
        )
        estimate = model.estimate(_walker(planes=2, per_plane=4), spares=1)
        self.assertEqual(estimate.units["satellites"], 10)
        self.assertEqual(estimate.launches, 4)
        self.assertEqual(estimate.operations, 10)

    def test_types(self):
        """Types of satellites are developed once each, with their own curves."""
        small = SatelliteCost(development=10, first_unit=1)
        model = CostModel({"large": self.satellite, "small": small})
        walker, sso = (
            _walker(planes=1, per_plane=2),
            _design(
                layout=sso_planes_space(crossing=13.5, satellites=3, count=1)
            ).layout,
        )
        mixed = model.estimate({"large": walker, "small": sso})
        self.assertEqual(mixed.units, {"large": 2, "small": 3})
        self.assertEqual(mixed.development, 110)
        self.assertAlmostEqual(mixed.production, self.satellite.production(2) + 3)
        one = model.estimate(walker, kind="small")
        self.assertEqual(one.development, 10)
        with self.assertRaises(ValueError):
            model.estimate(walker)
        with self.assertRaises(ValueError):
            model.estimate(walker, kind="medium")

    def test_growth(self):
        """Growth factors multiply the cost, or its elements."""
        model = CostModel(
            self.satellite, operations=OperationsCost(per_year=5), years=2
        )
        base = model.estimate(_walker())
        grown = model.estimate(_walker(), growth=1.5)
        self.assertAlmostEqual(grown.total, 1.5 * base.total)
        partial = model.estimate(_walker(), growth={"development": 2.0})
        self.assertAlmostEqual(partial.development, 200)
        self.assertAlmostEqual(partial.operations, base.operations)
        with self.assertRaises(ValueError):
            model.estimate(_walker(), growth={"insurance": 1.1})

    def test_trade_study(self):
        """Costs are constraints of design spaces and metrics of studies."""
        model = CostModel(
            {
                "narrow": SatelliteCost(development=50, first_unit=10),
                "wide": SatelliteCost(development=80, first_unit=16),
            },
            launch=LaunchCost(per_launch=20),
        )
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
                lambda d: model.estimate(d.layout, kind=d.imager).total - budget
            ],
        )
        designs = space.enumerate()
        nominal = [model.estimate(d.layout, kind=d.imager).total for d in designs]
        self.assertTrue(all(cost <= budget for cost in nominal))
        self.assertLess(len(designs), 2 * 9)

        def evaluate(d):
            cost = model.estimate(d.layout, kind=d.imager, growth=d.growth)
            return {**cost.metrics(), "satellites": count_satellites(d.layout)}

        results = Study(space, lambda d: d, evaluate, replicates=20).run(designs)
        summary = aggregate(
            results, {"cost": "cost", "cost p70": Quantile("cost", 0.7)}
        )
        # the 0.7 quantile of a growth skewed upward exceeds the mean
        self.assertTrue((summary["cost p70"] > summary["cost"]).all())
        ratio = summary["cost"] / np.array(nominal)
        np.testing.assert_allclose(ratio, ratio.iloc[0])  # common random numbers
