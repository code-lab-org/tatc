"""
Unit tests for the SpaceSystem schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from pydantic import ValidationError

from tatc.schemas import Instrument, PointedInstrument
from tatc.schemas.space.base import SpaceSystem


class TestSpaceSystem(unittest.TestCase):
    """
    Unit tests for the SpaceSystem schema.
    """

    def test_good_data(self):
        """
        Test that a SpaceSystem can be created with explicit instruments.
        """
        o = SpaceSystem(
            name="Test System",
            instruments=[Instrument(name="A"), Instrument(name="B")],
        )
        self.assertEqual(o.name, "Test System")
        self.assertEqual(len(o.instruments), 2)
        self.assertEqual(o.instruments[0].name, "A")
        self.assertEqual(o.instruments[1].name, "B")

    def test_name_required(self):
        """
        Test that name is required.
        """
        with self.assertRaises(ValidationError):
            SpaceSystem()

    def test_instruments_default_single_generic_instrument(self):
        """
        Test that omitting instruments defaults to a single generic
        (default-constructed) Instrument.
        """
        o = SpaceSystem(name="Test System")
        self.assertEqual(len(o.instruments), 1)
        self.assertEqual(o.instruments[0], Instrument())

    def test_instruments_default_independent_across_instances(self):
        """
        Test that the default instruments list is not shared (aliased)
        across separate SpaceSystem instances, since Field(default=[...])
        uses a mutable list/model as its default value: mutating one
        instance's default instruments must not affect another's.
        """
        first = SpaceSystem(name="First")
        second = SpaceSystem(name="Second")
        self.assertIsNot(first.instruments, second.instruments)
        self.assertIsNot(first.instruments[0], second.instruments[0])
        first.instruments[0].name = "Renamed"
        self.assertEqual(second.instruments[0].name, "Default")

    def test_instruments_empty_list_invalid(self):
        """
        Test that an empty instruments list is rejected (min_length=1): a
        space system must carry at least one instrument.
        """
        with self.assertRaises(ValidationError):
            SpaceSystem(name="Test System", instruments=[])

    def test_instruments_accepts_pointed_instrument(self):
        """
        Test that instruments accepts any AllInstruments member, including
        PointedInstrument, not just the base Instrument.
        """
        pointed = PointedInstrument(
            name="Pointed",
            cross_track_field_of_view=20.0,
            along_track_field_of_view=10.0,
        )
        o = SpaceSystem(name="Test System", instruments=[pointed])
        self.assertEqual(len(o.instruments), 1)
        self.assertIsInstance(o.instruments[0], PointedInstrument)


if __name__ == "__main__":
    unittest.main()
