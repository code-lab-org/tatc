"""
Unit tests for the catalog's validation metrics.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

import numpy as np

from tatc.catalog.metrics import (
    compute_mask_iou,
    compute_pass_scores,
    compute_product_limit_error,
    summarize,
)


class TestSummarize(unittest.TestCase):
    """
    Unit tests for summarize.
    """

    def test_values(self):
        """
        Test the statistics of a signed sample, whose percentile and maximum are
        of absolute values while its median is signed.
        """
        result = summarize([-4, -1, 0, 1, 2])
        self.assertEqual(result["n"], 5)
        self.assertEqual(result["median"], 0)
        self.assertEqual(result["max"], 4)
        self.assertAlmostEqual(result["p95"], np.percentile([4, 1, 0, 1, 2], 95))

    def test_ignores_not_finite(self):
        """
        Test that values that are not finite are ignored.
        """
        self.assertEqual(summarize([1, np.nan, np.inf])["n"], 1)

    def test_empty_invalid(self):
        """
        Test that a sample without finite values is an error, not a zero.
        """
        with self.assertRaises(ValueError):
            summarize([np.nan])


class TestPassScores(unittest.TestCase):
    """
    Unit tests for compute_pass_scores.
    """

    def test_values(self):
        """
        Test recall and precision.
        """
        result = compute_pass_scores(90, 100, 180)
        self.assertEqual(result["pass_recall"], 90)
        self.assertEqual(result["pass_precision"], 50)

    def test_invalid_counts(self):
        """
        Test that impossible counts are errors.
        """
        for counts in ((1, 0, 1), (1, 1, 0), (2, 1, 2), (2, 2, 1), (-1, 1, 1)):
            with self.assertRaises(ValueError):
                compute_pass_scores(*counts)


class TestMaskIou(unittest.TestCase):
    """
    Unit tests for compute_mask_iou.
    """

    def test_equal_weights(self):
        """
        Test the intersection over union of overlapping masks.
        """
        self.assertAlmostEqual(
            compute_mask_iou([True, True, False], [False, True, True]), 1 / 3
        )

    def test_weights(self):
        """
        Test that cells count by their weights.
        """
        self.assertAlmostEqual(
            compute_mask_iou([True, True, False], [False, True, True], [1, 2, 1]), 0.5
        )

    def test_empty_union_invalid(self):
        """
        Test that empty masks are an error rather than perfect agreement.
        """
        with self.assertRaises(ValueError):
            compute_mask_iou([False], [False])


class TestProductLimitError(unittest.TestCase):
    """
    Unit tests for compute_product_limit_error.
    """

    def test_observed_beyond_bound(self):
        """
        Test an observation beyond the bound.
        """
        self.assertAlmostEqual(
            compute_product_limit_error([66.9, 67.2, 67.3], [True, True, False], 67),
            0.3,
        )

    def test_unobserved_within_bound(self):
        """
        Test a pixel within the bound that is not observed.
        """
        self.assertAlmostEqual(
            compute_product_limit_error([66.5, 66.6, 67.1], [True, False, False], 67),
            0.5,
        )

    def test_one_sided_invalid(self):
        """
        Test that pixels on one side of the boundary do not locate it.
        """
        with self.assertRaises(ValueError):
            compute_product_limit_error([66, 67], [True, True], 67)
        with self.assertRaises(ValueError):
            compute_product_limit_error([66, 67], [False, False], 67)


if __name__ == "__main__":
    unittest.main()
