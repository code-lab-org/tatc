"""
Unit tests for the tatc.schemas.agility module.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np
from pydantic import ValidationError

from tatc.schemas import Agility


class TestAgility(unittest.TestCase):
    """
    Tests the slew times of agility.
    """

    def test_constant_rate(self):
        """
        Tests that slews without an acceleration turn at the maximum rate.
        """
        agility = Agility(max_rate=2, settle_time=5)
        self.assertAlmostEqual(agility.get_slew_time(30), 20)
        self.assertEqual(agility.get_slew_time(0), 0)

    def test_trapezoidal_profile(self):
        """
        Tests the triangular and trapezoidal profiles of rate, which meet at
        the angle at which a slew reaches the maximum rate.
        """
        agility = Agility(max_rate=2, max_acceleration=0.5)
        # rate ** 2 / acceleration = 8 degrees
        np.testing.assert_allclose(
            agility.get_slew_time([2, 8, 30]),
            [2 * np.sqrt(2 / 0.5), 8.0, 30 / 2 + 2 / 0.5],
        )
        times = agility.get_slew_time(np.linspace(0, 40, 401))
        self.assertTrue(np.all(np.diff(times) > 0))
        # slower than turning at the maximum rate throughout
        self.assertGreater(
            agility.get_slew_time(30), Agility(max_rate=2).get_slew_time(30)
        )

    def test_settle_time(self):
        """
        Tests that the settle time is added to every slew.
        """
        settled = Agility(max_rate=2, max_acceleration=0.5, settle_time=5)
        unsettled = Agility(max_rate=2, max_acceleration=0.5)
        self.assertAlmostEqual(
            settled.get_slew_time(30), unsettled.get_slew_time(30) + 5
        )

    def test_validation(self):
        """
        Tests that rates, accelerations, and settle times are validated, and
        that negative angles raise errors.
        """
        with self.assertRaises(ValidationError):
            Agility(max_rate=0)
        with self.assertRaises(ValidationError):
            Agility(max_rate=1, max_acceleration=0)
        with self.assertRaises(ValidationError):
            Agility(max_rate=1, settle_time=-1)
        with self.assertRaises(ValueError):
            Agility(max_rate=1).get_slew_time(-1)
