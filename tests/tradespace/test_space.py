"""
Unit tests for design spaces.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest

import pandas as pd

from tatc.tradespace import (
    Categorical,
    Choice,
    Derived,
    DesignSpace,
    Fixed,
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
            d.realize(2, seed=7, common_random_numbers=False).planes[0].phase
            for d in designs
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

    def test_scenario(self):
        """A scenario gives fixed values and uncertain parameters other values."""
        space = DesignSpace(
            a=Integer(1, 2),
            level=Fixed(1.0),
            double=Derived(lambda d: 2 * d.level),
            u=Uncertain.uniform(0, 1),
            v=Uncertain.uniform(0, 1),
            parts=Repeat(DesignSpace(k=Fixed(1)), count=2),
        )
        design = space.enumerate()[0]
        plain = design.realize(3, seed=5)
        scenario = {"level": 3.0, "u": 0.25, "parts[*].k": 7}
        changed = design.realize(3, seed=5, overrides=scenario)
        self.assertEqual((changed.level, changed.double, changed.u), (3.0, 6.0, 0.25))
        self.assertEqual([part.k for part in changed.parts], [7, 7])
        # other parameters are drawn as without the scenario
        self.assertEqual(changed.v, plain.v)
        self.assertEqual((plain.level, plain.double), (1.0, 2.0))
        self.assertEqual(design.realize(3, seed=5, overrides={}).flat(), plain.flat())
        for values in ({"a": 2}, {"double": 1.0}, {"w": 1}):
            with self.subTest(values=values):
                with self.assertRaises(ValueError):
                    design.realize(3, overrides=values)

    def test_uncertain_paths(self):
        """Uncertain parameters are listed by path, with [*] for items."""
        self.assertEqual(_planes_space().uncertain_paths(), ["planes[*].phase"])
        space = DesignSpace(
            seed=Uncertain.seed(),
            layout=Choice(
                a=DesignSpace(phase=Uncertain.uniform(0, 1)),
                b=DesignSpace(planes=Repeat(Uncertain.uniform(0, 1), count=2)),
            ),
        )
        self.assertEqual(
            space.uncertain_paths(), ["seed", "layout.a.phase", "layout.b.planes[*]"]
        )
        self.assertEqual(DesignSpace(a=Integer(1, 2)).uncertain_paths(), [])

    def test_shift(self):
        """A shift draws named parameters as in other replicates."""
        space = DesignSpace(
            u=Uncertain.uniform(0, 1),
            v=Uncertain.uniform(0, 1),
            planes=Repeat(DesignSpace(phase=Uncertain.uniform(0, 360)), count=2),
        )
        design = space.enumerate()[0]
        mixed = design.realize(3, seed=5, shift={"u": 10, "planes[*].phase": 1})
        self.assertEqual(mixed.u, design.realize(13, seed=5).u)
        self.assertEqual(mixed.v, design.realize(3, seed=5).v)
        self.assertEqual(
            [p.phase for p in mixed.planes],
            [p.phase for p in design.realize(4, seed=5).planes],
        )
        self.assertEqual(mixed.replicate, 3)
        self.assertEqual(
            design.realize(3, shift={"u": 0}).flat(), design.realize(3).flat()
        )
        independent = design.realize(3, common_random_numbers=False, shift={"u": 10})
        self.assertEqual(
            independent.u, design.realize(13, common_random_numbers=False).u
        )
        with self.assertRaises(ValueError):
            design.realize(3, shift={"w": 1})
        with self.assertRaises(ValueError):
            design.realize(3, shift={"planes[0].phase": 1})


class TestEnumeration(unittest.TestCase):
    """
    Unit tests for the enumeration of unordered repetitions, which walks only
    one order of their items.
    """

    @staticmethod
    def _every_walk(space):
        """Enumerates by walking every order of the items."""
        # pylint: disable=import-outside-toplevel,protected-access
        from tatc.tradespace.walks import _Context, _Odometer

        odometer = _Odometer()
        ctx = _Context(choose=odometer.choose)
        designs = {}
        while True:
            try:
                design = space._design(ctx)
                designs.setdefault(design.key, design)
            except Infeasible:
                pass
            if not odometer.advance():
                return list(designs.values())

    def test_same_designs(self):
        """Walking one order of unordered items gives the same designs."""
        spaces = [
            DesignSpace(
                parts=Repeat(
                    DesignSpace(
                        kind=Choice(
                            a=DesignSpace(x=Integer(0, 2)),
                            b=DesignSpace(y=Categorical(["p", "q"]), z=Integer(5, 12)),
                        )
                    ),
                    count=Integer(0, 3),
                    ordered=False,
                )
            ),
            DesignSpace(
                outer=Repeat(
                    DesignSpace(
                        inner=Repeat(
                            DesignSpace(v=Integer(0, 2)),
                            count=Integer(1, 2),
                            ordered=False,
                        )
                    ),
                    count=Integer(1, 2),
                    ordered=False,
                ),
                constraints=[lambda d: len(d.outer) < 2 or d.outer[0].inner[0].v != 1],
            ),
        ]
        for space in spaces:
            expected = self._every_walk(space)
            designs = space.enumerate()
            self.assertEqual([d.key for d in designs], [d.key for d in expected])
            self.assertEqual([d.flat() for d in designs], [d.flat() for d in expected])
            self.assertEqual(space.size(), len(expected))


class TestDesignsOf(unittest.TestCase):
    """
    Unit tests for designs rebuilt from data frames of their values, and for
    designs of extended design spaces.
    """

    def setUp(self):
        self.space = _hierarchical_space()
        self.designs = self.space.enumerate()
        # values by path, as in the results of a study (with missing values of
        # the branches a design does not choose, and integers as floats)
        self.frame = pd.DataFrame(
            [{"design": d.key, **d.realize(0).flat()} for d in self.designs]
        )

    def test_round_trip(self):
        """Designs rebuilt from their values are the designs."""
        rebuilt = self.space.designs_of(self.frame)
        self.assertEqual([d.key for d in rebuilt], [d.key for d in self.designs])
        self.assertEqual([d.flat() for d in rebuilt][5], self.designs[5].flat())
        # from a summary indexed by design, of some keys, in their order
        keys = [self.designs[3].key, self.designs[1].key]
        some = self.space.designs_of(self.frame.set_index("design"), keys)
        self.assertEqual([d.key for d in some], keys)
        # a design for each row, aligned with the rows (as of results of
        # several rows of a design), each design rebuilt once
        twice = pd.concat([self.frame, self.frame], ignore_index=True)
        aligned = self.space.designs_of(twice)
        self.assertEqual(list(aligned.index), list(twice.index))
        self.assertEqual(list(aligned.map(lambda d: d.key)), list(twice.design))
        self.assertIs(aligned.iloc[0], aligned.iloc[len(self.frame)])
        # without keys, every row's design
        bare = self.frame.drop(columns="design")
        rebuilt = self.space.designs_of(pd.concat([bare, bare]))
        self.assertEqual(len(rebuilt), 2 * len(self.designs))

    def test_errors(self):
        """Rows without the values of decisions, or of other designs, raise."""
        with self.assertRaises(KeyError):
            self.space.designs_of(self.frame, ["unknown"])
        with self.assertRaisesRegex(ValueError, "no value of altitude"):
            self.space.designs_of(self.frame.drop(columns="altitude"))
        other = self.frame.assign(design=self.frame.design.iloc[::-1].to_numpy())
        with self.assertRaisesRegex(ValueError, "another design"):
            self.space.designs_of(other)
        wrong = self.frame.assign(altitude=1.0)
        with self.assertRaisesRegex(ValueError, "not a level"):
            self.space.designs_of(wrong)

    def test_extend(self):
        """Extended design spaces replace or add children, and keep constraints."""
        space = DesignSpace(a=Integer(1, 3), b=2.0, constraints=[lambda d: d.a <= 2])
        extended = space.extend(c=Categorical(["x", "y"]), b=Integer(0, 1))
        self.assertEqual(list(extended.children), ["a", "b", "c"])
        self.assertEqual(len(extended.enumerate()), 2 * 2 * 2)
        self.assertEqual(list(space.children), ["a", "b"])  # unchanged
        design = space.enumerate()[1]
        derived = extended.design(design, b=1, c="y")
        self.assertEqual((derived.a, derived.b, derived.c), (2, 1, "y"))
        # a design of a group's own space gives its decisions
        group = DesignSpace(k=Integer(0, 1))
        nested = space.extend(g=group).design(design, g=group.enumerate()[1])
        self.assertEqual(nested.g.k, 1)
        with self.assertRaisesRegex(ValueError, "no children"):
            extended.design(design, z=1)
        with self.assertRaises(Infeasible):
            extended.design({"a": 3, "b": 0, "c": "x"})


def _layouts_space():
    """A choice of a group of planes (a number) or a repetition of planes."""
    return DesignSpace(
        altitude=Real(500e3, 600e3, step=100e3),
        layout=Choice(
            walker=DesignSpace(planes=Integer(1, 2), per_plane=Fixed(2)),
            sso=DesignSpace(
                planes=Repeat(
                    DesignSpace(crossing=Real(6, 9, step=3)), count=Categorical([1, 2])
                )
            ),
        ),
    )


class TestValuesByPath(unittest.TestCase):
    """
    Unit tests for the values of designs by their paths.
    """

    def setUp(self):
        self.designs = _layouts_space().enumerate()
        self.walker = next(d for d in self.designs if d.layout == "walker")
        self.sso = [d for d in self.designs if d.layout == "sso"][-1]

    def test_paths(self):
        """Paths of the columns of results (and of fixed values) give values."""
        self.assertEqual(self.walker["layout.walker.planes"], self.walker.layout.planes)
        self.assertEqual(self.walker["layout.walker.per_plane"], 2)  # fixed
        self.assertEqual(self.walker["layout.walker"], self.walker.layout.value)
        self.assertEqual(self.sso["layout.sso.planes.count"], 2)
        self.assertEqual(
            self.sso["layout.sso.planes[1].crossing"],
            self.sso.layout.planes[1].crossing,
        )
        for design in (self.walker, self.sso):
            for path, value in design.flat().items():
                self.assertEqual(design[path], value, path)
        with self.assertRaisesRegex(KeyError, "chose 'sso'"):
            self.sso["layout.walker.planes"]  # pylint: disable=pointless-statement
        with self.assertRaisesRegex(KeyError, "2 items"):
            self.sso[
                "layout.sso.planes[2].crossing"
            ]  # pylint: disable=pointless-statement

    def test_choices(self):
        """A choice equals the name of its branch, and steps into it."""
        layout = self.walker.layout
        self.assertTrue(layout == "walker")
        self.assertFalse(layout == "sso")
        self.assertIn("planes", layout)
        self.assertNotIn("crossing", layout)
        self.assertEqual(layout.walker.planes, layout.planes)
        with self.assertRaisesRegex(AttributeError, "chose 'walker'"):
            layout.sso  # pylint: disable=pointless-statement

    def test_items(self):
        """Repetitions give views of their items, also of slices."""
        planes = self.sso.layout.planes
        self.assertEqual(planes[-1].crossing, planes[1].crossing)
        self.assertEqual([plane.crossing for plane in planes[:1]], [planes[0].crossing])
        self.assertEqual(planes, list(planes))

    def test_constraints(self):
        """The callables of a design space read values by paths too."""
        space = DesignSpace(
            layout=Choice(walker=DesignSpace(planes=Integer(1, 3))),
            constraints=[lambda d: d["layout.walker.planes"] <= 2],
        )
        self.assertEqual(len(space.enumerate()), 2)


class TestReservedNames(unittest.TestCase):
    """
    Unit tests for names that the views of values give otherwise.
    """

    def test_reserved(self):
        """Children and branches may not be named as attributes of views, or
        as the status columns of the results of studies."""
        for name in ("values", "items", "keys", "get", "key", "flat", "space"):
            with self.assertRaisesRegex(ValueError, "reserved"):
                DesignSpace(**{name: 1})
        for name in ("design", "feasible", "error", "seconds", "scenario"):
            with self.assertRaisesRegex(ValueError, "reserved"):
                DesignSpace(**{name: Integer(1, 2)})
        for name in ("choice", "value"):
            with self.assertRaisesRegex(ValueError, "reserved"):
                Choice(**{name: DesignSpace()})


class TestCounts(unittest.TestCase):
    """
    Unit tests for the counts of repetitions.
    """

    def test_kinds(self):
        """Counts may be categorical, fixed, or derived."""
        categorical = DesignSpace(
            parts=Repeat(Integer(0, 1), count=Categorical([1, 3]))
        )
        designs = categorical.enumerate()
        self.assertEqual(len(designs), 2 + 8)
        self.assertEqual(
            categorical.unit_paths(),
            ["parts.count", "parts[0]", "parts[1]", "parts[2]"],
        )
        rebuilt = categorical.designs_of(
            pd.DataFrame([{"design": d.key, **d.flat()} for d in designs])
        )
        self.assertEqual([d.key for d in rebuilt], [d.key for d in designs])
        fixed = DesignSpace(parts=Repeat(Integer(0, 1), count=Fixed(2)))
        self.assertEqual(len(fixed.enumerate()), 4)
        derived = DesignSpace(
            n=Integer(1, 2), parts=Repeat(Integer(0, 1), count=Derived(lambda d: d.n))
        )
        self.assertEqual(len(derived.enumerate()), 2 + 4)
        with self.assertRaisesRegex(ValueError, "upper bound"):
            derived.unit_paths()

    def test_invalid(self):
        """Counts of other nodes, or negative, are rejected at once."""
        with self.assertRaisesRegex(TypeError, "count"):
            Repeat(Integer(0, 1), count=Real(1, 2))
        with self.assertRaises(ValueError):
            Repeat(Integer(0, 1), count=-1)


class TestDomains(unittest.TestCase):
    """
    Unit tests for the decisions of designs given by their values.
    """

    def test_checked(self):
        """Decisions outside their variables' domains are rejected."""
        space = DesignSpace(a=Integer(1, 6), b=Categorical(["x", "y"]), c=Real(0, 1))
        self.assertEqual(space.design({"a": 2, "b": "x", "c": 0.5}).a, 2)
        for decisions in ({"a": 99}, {"b": "z"}, {"c": 2.0}):
            with self.assertRaises(ValueError):
                space.design({"a": 2, "b": "x", "c": 0.5, **decisions})


class TestUncertainConstructors(unittest.TestCase):
    """
    Unit tests for the constructors of uncertain parameters.
    """

    def test_triangular(self):
        """A triangular parameter draws in its range, and names it."""
        growth = Uncertain.triangular(0.9, 1.0, 1.6)
        self.assertEqual(repr(growth), "Uncertain(triangular(0.9, 1.0, 1.6))")
        self.assertTrue(math.isnan(growth.nominal))
        space = DesignSpace(growth=growth)
        draws = [space.enumerate()[0].realize(r).growth for r in range(50)]
        self.assertTrue(all(0.9 <= draw <= 1.6 for draw in draws))
        with self.assertRaises(ValueError):
            Uncertain.triangular(1.0, 0.5, 2.0)

    def test_nominal(self):
        """The constructors take nominal values."""
        space = DesignSpace(x=Uncertain.uniform(0, 1, nominal=0.5))
        self.assertEqual(space.enumerate()[0].x, 0.5)
        self.assertEqual(Uncertain.normal(0, 1, nominal=0.0).nominal, 0.0)
        self.assertIn("seed()", repr(Uncertain.seed()))


class TestDecisionsAndViews(unittest.TestCase):
    """
    Unit tests for decisions given by value, and the views of values.
    """

    def test_levels(self):
        """Decisions given by value are the levels of their domains."""
        import numpy as np  # pylint: disable=import-outside-toplevel

        space = DesignSpace(a=Integer(1, 4), b=Real(0.5, 1.5, step=0.5))
        (enumerated,) = [d for d in space.enumerate() if d.a == 2 and d.b == 0.5]
        for a in (2, 2.0, np.int64(2)):
            design = space.design({"a": a, "b": 0.5})
            self.assertEqual(design.key, enumerated.key)
            self.assertIs(type(design.a), int)

    def test_views(self):
        """Designs and scopes tell whether they have values, and choices print
        as the names of their branches."""
        space = DesignSpace(
            a=Integer(1, 2),
            layout=Choice(walker=DesignSpace(planes=Integer(1, 2))),
            constraints=[lambda d: d.get("a") >= 1 and "layout.walker.planes" in d],
        )
        design = space.enumerate()[0]
        self.assertIn("a", design)
        self.assertIn("layout.walker.planes", design)
        self.assertNotIn("other", design)
        self.assertEqual(str(design.layout), "walker")
        with self.assertRaises(TypeError):
            list(design)
        with self.assertRaisesRegex(ValueError, "reserved"):
            Choice(walker=DesignSpace(value=1))

    def test_items(self):
        """Items are given by index within the repetition only."""
        space = DesignSpace(parts=Repeat(Integer(0, 1), count=3))
        parts = space.enumerate()[0].parts
        self.assertEqual(parts[-1], parts[2])
        for index in (3, -4):
            with self.assertRaises(IndexError):
                parts[index]  # pylint: disable=pointless-statement

    def test_counts_given(self):
        """Items given for a count that is not decided must match it."""
        space = DesignSpace(
            n=Integer(1, 2), parts=Repeat(Integer(0, 1), count=Derived(lambda d: d.n))
        )
        with self.assertRaisesRegex(ValueError, "items are given"):
            space.design({"n": 1, "parts": [0, 1]})


class TestBranchNames(unittest.TestCase):
    """
    Unit tests for choices given the names of their branches.
    """

    def test_override(self):
        """A choice given its branch's name keeps the decisions given of it, or
        decides it at its first levels."""
        space = _layouts_space()
        walker = space.design({"altitude": 600e3, "layout": ["walker", {"planes": 2}]})
        sso = space.design(walker, layout="sso")
        self.assertEqual(sso.layout, "sso")
        self.assertEqual(sso.altitude, 600e3)
        self.assertEqual(
            sso.decisions["layout"], ["sso", {"planes": [{"crossing": 6.0}]}]
        )
        again = space.design(walker, layout="walker")
        self.assertEqual(again.key, walker.key)
        # in the decisions too
        first = space.design({"altitude": 500e3, "layout": "walker"})
        self.assertEqual(first.layout.planes, 1)
        with self.assertRaisesRegex(ValueError, "not a level of layout"):
            space.design(walker, layout="other")


class TestUnorderedLevels(unittest.TestCase):
    """
    Unit tests for the order of the items of unordered repetitions.
    """

    def test_level_order(self):
        """Unordered items are in the order of the levels of their decisions
        (not of their text, as 10.5 precedes 9.0)."""
        space = DesignSpace(
            parts=Repeat(
                DesignSpace(size=Real(9.0, 10.5, step=1.5), kind=Categorical("ba")),
                count=2,
                ordered=False,
            )
        )
        design = space.design(
            {"parts": [{"size": 10.5, "kind": "b"}, {"size": 9.0, "kind": "a"}]}
        )
        self.assertEqual(
            design.decisions["parts"],
            [{"size": 9.0, "kind": "a"}, {"size": 10.5, "kind": "b"}],
        )
        # by the first decision that differs, in the order of declaration
        tied = space.design(
            {"parts": [{"size": 9.0, "kind": "a"}, {"size": 9.0, "kind": "b"}]}
        )
        self.assertEqual([p.kind for p in tied.parts], ["b", "a"])
        # enumerations and samples give the same order
        keys = {d.key for d in space.enumerate()}
        self.assertIn(design.key, keys)
        self.assertLessEqual({d.key for d in space.sample(20, seed=1)}, keys)


class TestScenarioOverrides(unittest.TestCase):
    """
    Unit tests for the overrides of scenarios, which are checked against the
    constraints of design spaces.
    """

    def setUp(self):
        self.space = DesignSpace(
            n=Integer(1, 4),
            limit=Fixed(4),
            u=Uncertain.uniform(0, 1),
            constraints=[lambda d: d.n - d.limit],
        )

    def test_constraints(self):
        """A design that violates a constraint with the overrides is infeasible."""
        design = self.space.design({"n": 3})
        self.assertEqual(design.realize(0, overrides={"limit": 3}).limit, 3)
        with self.assertRaisesRegex(Infeasible, "constraint 0"):
            design.realize(0, overrides={"limit": 2})

    def test_nan(self):
        """A constraint that is NaN is not satisfied."""
        space = DesignSpace(
            n=Integer(1, 4),
            u=Uncertain.uniform(0, 1),
            constraints=[lambda d: d.n - 10 * d.u],
        )
        self.assertEqual(space.size(), 0)
        with self.assertRaisesRegex(Infeasible, "NaN") as context:
            space.design({"n": 1})
        self.assertEqual(context.exception.violation, 1)

    def test_scenario_name(self):
        """A realized design gives the name of its scenario."""
        design = self.space.design({"n": 1})
        self.assertIsNone(design.scenario)
        self.assertIsNone(design.realize(0).scenario)
        named = design.realize(0, overrides={"u": 0.5}, scenario="half")
        self.assertEqual((named.scenario, named.u), ("half", 0.5))
        self.assertIn("scenario='half'", repr(named))

    def test_checks(self):
        """Paths of shifts and overrides are checked by the design space."""
        space = _planes_space()
        space.check_uncertain(["planes[*].phase"])
        space.check_overrides({"planes[*].phase": 1.0})
        self.space.check_overrides({"limit": 1, "u": 0.5})
        for paths in (["planes[0].phase"], ["per_plane"], ["planes[*]"], [1]):
            with self.subTest(paths=paths):
                with self.assertRaisesRegex(ValueError, "not uncertain"):
                    space.check_uncertain(paths)
        with self.assertRaisesRegex(ValueError, "not fixed values"):
            self.space.check_overrides({"n": 2})

    def test_seed(self):
        """Seeds are non-negative integers."""
        design = self.space.design({"n": 1})
        for seed in (-1, 1.5, True):
            with self.subTest(seed=seed):
                with self.assertRaisesRegex(ValueError, "seed"):
                    design.realize(0, seed=seed)


class TestKeywords(unittest.TestCase):
    """
    Unit tests for options given by keyword only.
    """

    def test_keywords(self):
        """The options of samples and realizations are keywords."""
        space = _planes_space()
        with self.assertRaises(TypeError):
            space.sample(5, "random")  # pylint: disable=too-many-function-args
        with self.assertRaises(TypeError):
            space.enumerate()[0].realize(0, 1)  # pylint: disable=too-many-function-args


class TestDesignsOfSpaces(unittest.TestCase):
    """
    Unit tests for designs given to the functions of trade studies.
    """

    def test_other_spaces(self):
        """Designs of other design spaces are the designs of their decisions."""
        # pylint: disable=import-outside-toplevel,protected-access
        from tatc.tradespace.space import _designs

        space = DesignSpace(n=Integer(1, 2), altitude=Fixed(500))
        higher = space.extend(altitude=Fixed(700))
        designs = _designs(higher, space.enumerate())
        self.assertTrue(all(d.space is higher for d in designs))
        self.assertEqual([d.altitude for d in designs], [700, 700])
        same = higher.enumerate()
        self.assertTrue(all(a is b for a, b in zip(_designs(higher, same), same)))

    def test_frames(self):
        """Designs rebuilt from frames read the frame's keys once, and a frame
        that lacks their values gives the reason."""
        # pylint: disable=import-outside-toplevel,protected-access
        from tatc.tradespace.space import _designs

        class Counted(pd.DataFrame):
            """A data frame that counts the reads of its keys of designs."""

            reads = 0

            @property
            def _constructor(self):
                return Counted

            def __getitem__(self, key):
                if isinstance(key, str) and key == "design":
                    Counted.reads += 1
                return super().__getitem__(key)

        space = DesignSpace(a=Integer(0, 99), b=Integer(0, 4))
        designs = space.enumerate()
        rows = pd.DataFrame(
            [
                {"design": d.key, "replicate": r, **d.flat()}
                for d in designs
                for r in (0, 1)
            ]
        )
        keys = [d.key for d in designs]
        rebuilt = _designs(space, keys, Counted(rows))
        self.assertEqual([d.key for d in rebuilt], keys)
        self.assertLess(Counted.reads, 10)
        with self.assertRaisesRegex(ValueError, "not given") as context:
            _designs(space, keys[:2], rows.drop(columns="b"))
        self.assertIn("no value of b", str(context.exception.__cause__))


class TestConfigure(unittest.TestCase):
    """
    Unit tests for dataclasses configured from the values of designs.
    """

    def test_design(self):
        """A dataclass is configured from a design, or a view of its values."""
        # pylint: disable=import-outside-toplevel
        import dataclasses

        from tatc.tradespace.values import configure

        @dataclasses.dataclass(frozen=True)
        class Policy:
            """A policy of two fields."""

            spares: int = 0
            years: float = 1.0

        design = DesignSpace(spares=Integer(1, 2), other=Integer(0, 1)).enumerate()[-1]
        self.assertEqual(configure(Policy(), design), Policy(spares=2))
        self.assertEqual(configure(Policy(), design.values), Policy(spares=2))
        self.assertEqual(configure(Policy(), design, years=2.0).years, 2.0)
        with self.assertRaisesRegex(TypeError, "values of a design"):
            configure(Policy(), 3)
