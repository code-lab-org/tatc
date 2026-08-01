"""
Unit tests for the BaseConstellation schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from tatc.schemas import Instrument
from tatc.schemas.space.base_constellation import BaseConstellation


class TestBaseConstellation(unittest.TestCase):
    """
    Unit tests for the BaseConstellation schema.
    """
    def test_generate_members_not_implemented_on_base(self):
        """
        Test that generate_members raises NotImplementedError on the bare
        base class, since BaseConstellation has no universal way to derive
        member satellites.
        """
        with self.assertRaises(NotImplementedError):
            BaseConstellation(name="Test Constellation").generate_members()

    def test_inherits_space_system_fields(self):
        """
        Test that BaseConstellation inherits SpaceSystem's fields
        (name, instruments) unchanged.
        """
        o = BaseConstellation(name="Test Constellation")
        self.assertEqual(o.name, "Test Constellation")
        self.assertEqual(len(o.instruments), 1)
        self.assertEqual(o.instruments[0], Instrument())


if __name__ == "__main__":
    unittest.main()
