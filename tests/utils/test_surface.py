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

    def test_compute_number_samples_matches_flat_approximation_at_small_scale(self):
        """
        Test that, for sample distances much smaller than the Earth's
        radius, the spherical-cap-based sample count closely approximates
        the flat-plane (circular disk area) approximation.
        """
        for sample_distance in (1000, 50000, 200000):
            spherical = compute_number_samples(sample_distance)
            flat = constants.EARTH_SURFACE_AREA / (np.pi * (sample_distance / 2) ** 2)
            self.assertAlmostEqual(spherical / flat, 1.0, delta=1e-4)

    def test_compute_number_samples_exceeds_flat_approximation_at_large_scale(self):
        """
        Test that, at a sample distance large enough for Earth's curvature
        to matter, the spherical-cap-based sample count is greater than
        the flat-plane approximation: a spherical cap covers less area
        than a flat disk of the same angular radius, so more of them are
        needed to cover the same total surface area.
        """
        sample_distance = 5000000
        spherical = compute_number_samples(sample_distance)
        flat = int(
            constants.EARTH_SURFACE_AREA / (np.pi * (sample_distance / 2) ** 2)
        )
        self.assertGreater(spherical, flat)

    def test_compute_number_samples_decreases_with_distance(self):
        """
        Test that the number of samples decreases monotonically as the
        sample distance increases.
        """
        counts = [
            compute_number_samples(d) for d in (1000, 10000, 100000, 1000000)
        ]
        self.assertEqual(counts, sorted(counts, reverse=True))

    def test_compute_number_samples_zero_distance_raises(self):
        """
        Test that a zero sample distance (a meaningless, infinitely dense
        request) raises a clear ValueError rather than returning a
        misleading value. (Previously this reached an unguarded division
        by zero deep in the geometry, which happened to raise
        OverflowError -- via a RuntimeWarning-then-inf-then-int(inf)
        chain -- purely by numeric coincidence, and did not extend to
        negative distances at all; see
        test_compute_number_samples_negative_distance_raises.)
        """
        with self.assertRaises(ValueError):
            compute_number_samples(0)

    def test_compute_number_samples_negative_distance_raises(self):
        """
        Regression test: a negative sample distance must also be
        rejected. Previously this was silently accepted and treated as
        if it were the corresponding positive distance, since the
        geometry's cos(theta/2) term is an even function of distance and
        so cannot distinguish the sign on its own.
        """
        with self.assertRaises(ValueError):
            compute_number_samples(-2000000)
