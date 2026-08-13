"""
Unit tests for the tatc.utils.formatting module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from tatc.utils import zero_pad


class TestFormatting(unittest.TestCase):
    """
    Unit tests for the tatc.utils.formatting module.
    """

    def test_zero_pad_pads_to_max_number_width(self):
        """
        Test that the current number is zero-padded to the digit width of the max number.
        """
        self.assertEqual(zero_pad("Sat", 12, 3), "Sat 03")

    def test_zero_pad_current_equals_max(self):
        """
        Test that no padding is added when the current number already has the
        same digit width as the max number.
        """
        self.assertEqual(zero_pad("Sat", 12, 12), "Sat 12")

    def test_zero_pad_single_digit_max(self):
        """
        Test that no padding is added when the max number is a single digit.
        """
        self.assertEqual(zero_pad("Sat", 5, 3), "Sat 3")

    def test_zero_pad_current_number_exceeds_max_number_width(self):
        """
        Test that a current number with more digits than the max number is
        left at its natural width rather than truncated.
        """
        self.assertEqual(zero_pad("Sat", 9, 15), "Sat 15")

    def test_zero_pad_preserves_numeric_sort_order(self):
        """
        Test that labels generated in increasing numeric order also sort in
        increasing lexicographic (string) order, which is the purpose of
        the zero padding.
        """
        count = 12
        labels = [zero_pad("Sat", count, i) for i in range(1, count + 1)]
        self.assertEqual(sorted(labels), labels)
