"""
Unit tests for the tatc.utils.orbital module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

from tatc.utils import mean_anomaly_to_true_anomaly, true_anomaly_to_mean_anomaly


class TestOrbital(unittest.TestCase):
    """
    Unit tests for the tatc.utils.orbital module.
    """
    def test_mean_anomaly_to_true_anomaly(self):
        """
        Test that the mean anomaly can be converted to true anomaly.
        """
        self.assertAlmostEqual(
            mean_anomaly_to_true_anomaly(78.940629, 0.0001492), 78.95065818, delta=0.01
        )

    def test_true_anomaly_to_mean_anomaly(self):
        """
        Test that the true anomaly can be converted to mean anomaly.
        """
        self.assertAlmostEqual(
            true_anomaly_to_mean_anomaly(78.95065818, 0.0001492), 78.940629, delta=0.01
        )
