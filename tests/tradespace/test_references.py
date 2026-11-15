"""
Unit tests for reading designs: references to their values, dataclasses
configured from them, and functions given them.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import dataclasses
import unittest

from tatc.tradespace import (
    Choice,
    DesignSpace,
    Integer,
    Real,
    Uncertain,
    configure,
    value_of,
)
from tatc.tradespace.references import given_design, takes_design

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

        self.assertTrue(takes_design(positional))
        self.assertTrue(takes_design(keyword))
        self.assertFalse(takes_design(only_position))
        self.assertFalse(takes_design(lambda x: x))
        self.assertFalse(takes_design(lambda x, **kwargs: x))
        self.assertFalse(takes_design(len))

    def test_given_design(self):
        """A function that takes the design is given it; others are not."""
        design = SPACE.enumerate()[0]
        self.assertEqual(
            given_design(lambda x, design: (x, design.key), design)(1),
            (1, design.key),
        )

        def alone(x):
            return x

        self.assertIs(given_design(alone, design), alone)
