"""
Unit tests for reading designs: references to their values, dataclasses
configured from them, and functions given them.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

# pylint: disable=protected-access

import dataclasses
import unittest

from tatc.tradespace import (
    Choice,
    DesignSpace,
    Fixed,
    Integer,
    Real,
    Repeat,
    Uncertain,
    configure,
    value_of,
)
from tatc.tradespace.references import (
    _check_references,
    _given_design,
    _takes_design,
)

SPACE = DesignSpace(
    altitude=Real(600e3, 700e3, step=100e3),
    layout=DesignSpace(planes=Integer(1, 2), per_plane=Integer(1, 2)),
    growth=Uncertain.uniform(1.0, 1.5),
)


@dataclasses.dataclass(frozen=True)
class _Policy:
    """A policy of two fields."""

    spares: int = 0
    years: float = 1.0


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

    def test_misspelled_names(self):
        """The error of a misspelled name lists the names, and suggests one."""
        d = SPACE.enumerate()[0]
        with self.assertRaisesRegex(KeyError, "the design has no value 'altitud'"):
            value_of(d, "altitud")
        with self.assertRaisesRegex(KeyError, "did you mean 'altitude'"):
            value_of(d, "altitud")
        with self.assertRaisesRegex(KeyError, "did you mean 'planes'"):
            value_of(d, "layout.plane")
        with self.assertRaisesRegex(AttributeError, "did you mean 'growth'"):
            _ = d.grwth
        self.assertNotIn("altitud", d)
        self.assertIn("layout.planes", d)
        self.assertNotIn("layout.plane", d)

    def test_check_references(self):
        """References that are names or paths must name values of the designs
        (in any branch of a choice); functions are not checked."""
        space = DesignSpace(
            seed=Uncertain.seed(),
            layout=Choice(
                walker=DesignSpace(planes=Integer(1, 2)),
                sso=DesignSpace(crossing=Real(6, 18, step=6)),
            ),
            parts=Repeat(DesignSpace(size=Integer(1, 2)), 2),
            fixed=Fixed({"a": 1}),
        )
        valid = [
            "seed",
            "layout",
            "layout.walker",
            "layout.walker.planes",
            "layout.planes",  # of the branch chosen
            "layout.crossing",
            "parts.count",
            "parts[1].size",
            "fixed.a",  # within a value
        ]
        _check_references(space, [(f"reference {r}", r) for r in valid])
        _check_references(space, [("a function", lambda d: d.other), ("none", None)])
        with self.assertRaisesRegex(
            ValueError, r"Targets seed 'sed': the designs have no value 'sed'"
        ):
            _check_references(space, [("Targets seed", "sed")])
        with self.assertRaisesRegex(ValueError, "did you mean 'planes'"):
            _check_references(space, [("kind", "layout.walker.plane")])
        with self.assertRaisesRegex(ValueError, "did you mean 'crossing'"):
            _check_references(space, [("kind", "layout.crosing")])
        with self.assertRaisesRegex(ValueError, "parts has no value 'size'"):
            _check_references(space, [("kind", "parts.size")])
        with self.assertRaisesRegex(ValueError, "parts.0. has no value 'sise'"):
            _check_references(space, [("kind", "parts[0].sise")])
        with self.assertRaisesRegex(ValueError, "not a name or path"):
            _check_references(space, [("kind", "a..b")])
        with self.assertRaisesRegex(ValueError, "'x'.*'y'"):  # every one
            _check_references(space, [("one", "x"), ("another", "y")])


class TestConfigure(unittest.TestCase):
    """
    Unit tests for dataclasses configured from the values of designs.
    """

    def test_design(self):
        """A dataclass is configured from a design, or a view of its values."""
        design = DesignSpace(spares=Integer(1, 2), other=Integer(0, 1)).enumerate()[-1]
        self.assertEqual(configure(_Policy(), design), _Policy(spares=2))
        self.assertEqual(configure(_Policy(), design.values), _Policy(spares=2))
        self.assertEqual(configure(_Policy(), design, years=2.0).years, 2.0)
        self.assertEqual(configure(_Policy(), {"spares": 3}), _Policy(spares=3))
        with self.assertRaisesRegex(TypeError, "values of a design"):
            configure(_Policy(), 3)
        with self.assertRaisesRegex(TypeError, "dataclass instance"):
            configure(_Policy, design)
        with self.assertRaisesRegex(ValueError, "no fields"):
            configure(_Policy(), design, other=1)

    def test_choice(self):
        """A choice configures by the values of its branch, which must be a
        group (or a fixed `None`)."""
        space = DesignSpace(
            policy=Choice(
                spared=DesignSpace(spares=Integer(1, 2)),
                plain=Integer(1, 2),
                none=None,
            )
        )
        spared = space.design({"policy": ["spared", {"spares": 2}]})
        self.assertEqual(configure(_Policy(), spared.policy), _Policy(spares=2))
        none = space.design({"policy": "none"})
        self.assertEqual(configure(_Policy(), none.policy), _Policy())
        plain = space.design({"policy": ["plain", 2]})
        with self.assertRaisesRegex(TypeError, "chose 'plain', whose value is 2"):
            configure(_Policy(), plain.policy)


class TestGivenDesign(unittest.TestCase):
    """
    Unit tests for functions given the design by a keyword.
    """

    def test_takes_design(self):
        """A function takes the design if it has a parameter named `design`
        that may be given by keyword."""

        def positional(x, design):
            return x, design

        def keyword(x, *, design=None):
            return x, design

        def only_position(x, design, /):
            return x, design

        self.assertTrue(_takes_design(positional))
        self.assertTrue(_takes_design(keyword))
        self.assertFalse(_takes_design(only_position))
        # the first parameter is given what the function is of
        self.assertFalse(_takes_design(lambda design: design))
        self.assertFalse(_takes_design(lambda design, other=None: design))
        self.assertTrue(_takes_design(lambda x, design=None: x))
        self.assertFalse(_takes_design(lambda x: x))
        self.assertFalse(_takes_design(lambda x, **kwargs: x))
        self.assertFalse(_takes_design(len))

    def test_given_design(self):
        """A function that takes the design is given it; others are not."""
        design = SPACE.enumerate()[0]
        self.assertEqual(
            _given_design(lambda x, design: (x, design.key), design)(1),
            (1, design.key),
        )

        def alone(x):
            return x

        self.assertIs(_given_design(alone, design), alone)
