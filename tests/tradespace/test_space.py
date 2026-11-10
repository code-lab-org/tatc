"""
Unit tests for design spaces.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest

from tatc.tradespace import (
    Categorical,
    Choice,
    Derived,
    DesignSpace,
    Infeasible,
    Integer,
    Real,
    Repeat,
    Uncertain,
)


def _planes_space():
    """Planes spread over 12 hours with a random phase, for at most 6 satellites."""
    return DesignSpace(
        per_plane=Integer(1, 6),
        planes=Repeat(
            DesignSpace(
                crossing=Derived(lambda d: (13.5 + 12 * d.index / d.planes.count) % 24),
                phase=Uncertain.uniform(0, 360),
            ),
            count=Integer(1, 6),
        ),
        constraints=[lambda d: d.planes.count * d.per_plane <= 6],
    )


def _hierarchical_space():
    """A choice between Walker constellations and unordered sun-synchronous planes."""
    plane = DesignSpace(
        crossing=Real(0, 18, step=6),
        satellites=Integer(1, 2),
        payload=Categorical(["atms", "cris"]),
        phase=Uncertain.uniform(0, 360),
    )
    return DesignSpace(
        altitude=Real(700e3, 850e3, step=75e3),
        layout=Choice(
            walker=DesignSpace(
                inclination=Real(40, 70, step=15),
                planes=Integer(1, 3),
                per_plane=Integer(1, 2),
                relative_spacing=Integer(0, lambda d: d.planes - 1),
            ),
            sso=DesignSpace(planes=Repeat(plane, count=Integer(1, 2), ordered=False)),
        ),
    )


class TestVariables(unittest.TestCase):
    """
    Unit tests for design variables.
    """

    def test_levels(self):
        """Discrete variables enumerate their levels, inclusive of the upper bound."""
        space = DesignSpace(
            a=Real(0, 1, step=0.1), b=Integer(2, 10, step=4), c=Categorical("xyz")
        )
        designs = space.enumerate()
        self.assertEqual(len(designs), 11 * 3 * 3)
        self.assertEqual(sorted({d.a for d in designs})[-1], 1.0)
        self.assertEqual(sorted({d.b for d in designs}), [2, 6, 10])
        self.assertIsInstance(designs[0].b, int)

    def test_continuous_cannot_enumerate(self):
        """A continuous variable can be sampled but not enumerated."""
        space = DesignSpace(a=Real(1, 100, log=True))
        with self.assertRaises(ValueError):
            space.enumerate()
        values = [d.a for d in space.sample(50, seed=0)]
        self.assertTrue(all(1 <= v <= 100 for v in values))
        self.assertEqual(len(set(values)), 50)

    def test_dependent_bounds(self):
        """Bounds may depend on earlier values."""
        space = DesignSpace(
            planes=Integer(1, 4), spacing=Integer(0, lambda d: d.planes - 1)
        )
        self.assertEqual(space.size(), 1 + 2 + 3 + 4)

    def test_empty_domain_is_infeasible(self):
        """A variable without levels makes a design infeasible."""
        space = DesignSpace(a=Integer(0, 2), b=Integer(1, lambda d: d.a))
        self.assertEqual(
            {(d.a, d.b) for d in space.enumerate()}, {(1, 1), (2, 1), (2, 2)}
        )

    def test_reserved_names(self):
        """Names used by views cannot be children."""
        with self.assertRaises(ValueError):
            DesignSpace(index=Integer(0, 1))


class TestDesignSpace(unittest.TestCase):
    """
    Unit tests for design spaces.
    """

    def test_constraints(self):
        """Constraints remove designs."""
        designs = _planes_space().enumerate()
        self.assertEqual(len(designs), 14)
        pairs = {(d.planes.count, d.per_plane) for d in designs}
        self.assertIn((2, 3), pairs)
        self.assertNotIn((2, 4), pairs)

    def test_derived_and_scope(self):
        """Derived values resolve names in enclosing groups and the item index."""
        design = [d for d in _planes_space().enumerate() if d.planes.count == 4][0]
        self.assertEqual([p.crossing for p in design.planes], [13.5, 16.5, 19.5, 22.5])
        self.assertEqual(design.planes[2].index, 2)
        self.assertEqual(len(design.planes), 4)

    def test_hierarchical_enumeration(self):
        """Choices and repetitions enumerate the designs of each branch."""
        designs = _hierarchical_space().enumerate()
        walker = [d for d in designs if d.layout.choice == "walker"]
        sso = [d for d in designs if d.layout.choice == "sso"]
        # 3 altitudes x 3 inclinations x (1 x 2 + 2 x 2 + 3 x 2) (planes x spacing x per plane)
        self.assertEqual(len(walker), 3 * 3 * 12)
        # 16 plane designs (4 x 2 x 2): one plane, or two unordered planes (16 x 17 / 2)
        self.assertEqual(len(sso), 3 * (16 + 136))
        self.assertEqual(len({d.key for d in designs}), len(designs))

    def test_unordered_items_are_canonical(self):
        """Designs that differ only in the order of unordered items are one design."""
        space = _hierarchical_space()
        first = {"crossing": 0.0, "satellites": 1, "payload": "atms"}
        second = {"crossing": 6.0, "satellites": 2, "payload": "cris"}
        a = space.design(
            {"altitude": 700e3, "layout": ["sso", {"planes": [first, second]}]}
        )
        b = space.design(
            {"altitude": 700e3, "layout": ["sso", {"planes": [second, first]}]}
        )
        self.assertEqual(a.key, b.key)
        self.assertEqual(
            [p.payload for p in a.layout.planes], [p.payload for p in b.layout.planes]
        )
        self.assertEqual(a.decisions, b.decisions)

    def test_flat(self):
        """Flat values are keyed by path, without fixed values."""
        space = DesignSpace(fixed=3, layout=Choice(a=DesignSpace(x=Integer(1, 1))))
        self.assertEqual(space.enumerate()[0].flat(), {"layout": "a", "layout.a.x": 1})
        flat = _planes_space().enumerate()[1].flat()
        self.assertEqual(flat["planes.count"], 2)
        self.assertEqual(flat["planes[1].crossing"], 19.5)
        self.assertTrue(math.isnan(flat["planes[1].phase"]))

    def test_design_from_decisions(self):
        """A design is rebuilt from its decisions, and its constraints are checked."""
        space = _planes_space()
        for design in space.enumerate():
            self.assertEqual(space.design(design.decisions).key, design.key)
        with self.assertRaises(Infeasible):
            space.design({"per_plane": 4, "planes": [{}, {}]})

    def test_sample(self):
        """Sampled designs are feasible, and Latin hypercube samples are stratified."""
        space = _planes_space()
        keys = {d.key for d in space.enumerate()}
        for method in ("lhs", "random"):
            designs = space.sample(30, method=method, seed=3)
            self.assertEqual(len(designs), 30)
            self.assertTrue({d.key for d in designs} <= keys)
        stratified = DesignSpace(a=Real(0, 1)).sample(10, seed=0)
        self.assertEqual(sorted(int(d.a * 10) for d in stratified), list(range(10)))

    def test_constraint_violation(self):
        """Constraints may give the amount of their violation."""
        space = DesignSpace(
            a=Integer(0, 5), constraints=[lambda d: d.a - 3, lambda d: d.a != 1]
        )
        self.assertEqual([d.a for d in space.enumerate()], [0, 2, 3])
        with self.assertRaises(Infeasible) as context:
            space.design({"a": 5})
        self.assertEqual(context.exception.violation, 2)
        with self.assertRaises(Infeasible) as context:
            space.design({"a": 1})
        self.assertEqual(context.exception.violation, 1)

    def test_infeasible_sample(self):
        """Sampling gives up on a space without feasible designs."""
        space = DesignSpace(a=Integer(0, 1), constraints=[lambda d: False])
        with self.assertRaises(Infeasible):
            space.sample(1, max_tries=5)
        self.assertEqual(space.enumerate(), [])


class TestRealize(unittest.TestCase):
    """
    Unit tests for realizing designs.
    """

    def test_common_random_numbers(self):
        """Uncertain parameters at the same path are common to every design."""
        designs = _planes_space().enumerate()
        phases = {d.realize(2, seed=7).planes[0].phase for d in designs}
        self.assertEqual(len(phases), 1)
        other = designs[0].realize(3, seed=7).planes[0].phase
        self.assertNotIn(other, phases)
        independent = {
            d.realize(2, seed=7, common=False).planes[0].phase for d in designs
        }
        self.assertEqual(len(independent), len(designs))

    def test_realize_is_reproducible(self):
        """Realizing a design again gives the same values, and keeps its key."""
        design = _planes_space().enumerate()[5]
        a, b = design.realize(1), design.realize(1)
        self.assertEqual(a.flat(), b.flat())
        self.assertEqual(a.key, design.key)
        self.assertEqual(a.replicate, 1)
        self.assertTrue(all(0 <= p.phase < 360 for p in a.planes))

    def test_uncertain_choice(self):
        """Uncertain choices are drawn from their options."""
        space = DesignSpace(epoch=Uncertain.choice(["jan", "jul"]))
        design = space.enumerate()[0]
        self.assertIsNone(design.epoch)
        self.assertEqual({design.realize(r).epoch for r in range(20)}, {"jan", "jul"})
