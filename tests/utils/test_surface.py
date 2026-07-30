"""
Unit tests for the tatc.utils.surface module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""
import unittest

import numpy as np

from tatc import constants
from tatc.utils import compute_number_samples


class TestSurface(unittest.TestCase):
    """
    Unit tests for the tatc.utils.surface module.
    """
    def test_compute_number_samples(self):
        """
        Test that the number of samples can be computed for a given sample distance.
        """
        # rough approximation based on flat sample areas
        sample_distance = 10000
        num_samples = int(
            constants.EARTH_SURFACE_AREA / (np.pi * (sample_distance / 2) ** 2)
        )
        self.assertEqual(compute_number_samples(sample_distance), num_samples)
