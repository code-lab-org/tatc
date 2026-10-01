"""
Unit tests for the tatc.utils.radar module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest

from tatc.utils import (
    compute_radar_beam_height,
    compute_radar_ground_range,
    compute_radar_ground_range_bounds,
    compute_radar_slant_range,
    compute_terrain_elevation_angle,
)


class TestRadar(unittest.TestCase):
    """
    Unit tests for the tatc.utils.radar module.
    """

    def test_compute_radar_beam_height_zero_range(self):
        """
        Test that the beam height at zero slant range equals the station height.
        """
        for elevation_angle in (0, 0.5, 10, 45, 90):
            self.assertAlmostEqual(
                compute_radar_beam_height(0, elevation_angle, 370), 370, delta=1e-6
            )

    def test_compute_radar_beam_height_increases_with_range(self):
        """
        Test that beam height increases monotonically with slant range for
        a non-negative elevation angle.
        """
        for elevation_angle in (0, 0.5, 5, 19.5):
            heights = [
                compute_radar_beam_height(r, elevation_angle)
                for r in (0, 10000, 50000, 100000, 230000)
            ]
            self.assertEqual(heights, sorted(heights))

    def test_compute_radar_ground_range_zero_range(self):
        """
        Test that the ground range at zero slant range is zero.
        """
        self.assertEqual(compute_radar_ground_range(0, 0.5), 0)

    def test_compute_radar_ground_range_increases_with_range(self):
        """
        Test that ground range increases monotonically with slant range.
        """
        ground_ranges = [
            compute_radar_ground_range(r, 0.5)
            for r in (0, 10000, 50000, 100000, 230000)
        ]
        self.assertEqual(ground_ranges, sorted(ground_ranges))

    def test_compute_radar_slant_range_inverts_beam_height(self):
        """
        Test that compute_radar_slant_range is the inverse of
        compute_radar_beam_height across a range of elevation angles, slant
        ranges, and station heights.
        """
        for elevation_angle in (0, 0.5, 5, 19.5, 45):
            for slant_range in (1000, 50000, 100000, 230000):
                for station_height in (0, 370):
                    height = compute_radar_beam_height(
                        slant_range, elevation_angle, station_height
                    )
                    self.assertAlmostEqual(
                        compute_radar_slant_range(
                            elevation_angle, height, station_height
                        ),
                        slant_range,
                        delta=1e-3,
                    )

    def test_compute_radar_slant_range_decreases_with_elevation_angle(self):
        """
        Test that, for a fixed target height, the slant range required to
        reach it decreases monotonically as elevation angle increases (a
        steeper beam reaches a given height sooner).
        """
        slant_ranges = [
            compute_radar_slant_range(elevation_angle, 3048)
            for elevation_angle in (0.5, 2, 5, 10, 19.5)
        ]
        self.assertEqual(slant_ranges, sorted(slant_ranges, reverse=True))

    def test_compute_radar_slant_range_nan_below_station_height(self):
        """
        Test that no real (non-negative) slant range exists for a target
        strictly below the station height, for any elevation angle.
        """
        for elevation_angle in (0, 0.5, 19.5, 45, 90):
            self.assertTrue(
                math.isnan(compute_radar_slant_range(elevation_angle, 0, 370))
            )

    def test_compute_radar_slant_range_zero_at_station_height(self):
        """
        Test that a target exactly at the station height is trivially
        reached at a slant range of zero, for any elevation angle (the
        beam's own starting point).
        """
        for elevation_angle in (0, 0.5, 19.5, 45, 90):
            self.assertAlmostEqual(
                compute_radar_slant_range(elevation_angle, 370, 370), 0, delta=1e-6
            )

    def test_compute_radar_ground_range_bounds_target_at_station_height(self):
        """
        Test that a target at or below the station elevation yields a full
        disk (zero inner ground range) out to the ground range reached by
        max_range at the minimum elevation angle.
        """
        bounds = compute_radar_ground_range_bounds(0.5, 19.5, 230000, 0, 0)
        self.assertIsNotNone(bounds)
        self.assertEqual(bounds[0], 0.0)
        self.assertAlmostEqual(
            bounds[1], compute_radar_ground_range(230000, 0.5, 0), delta=1e-6
        )

    def test_compute_radar_ground_range_bounds_annulus(self):
        """
        Test that a target well above the station elevation yields a
        nontrivial annulus (positive inner bound, inner < outer).
        """
        bounds = compute_radar_ground_range_bounds(0.5, 19.5, 230000, 3048, 0)
        self.assertIsNotNone(bounds)
        inner, outer = bounds
        self.assertGreater(inner, 0)
        self.assertLess(inner, outer)

    def test_compute_radar_ground_range_bounds_cone_shrinks_with_max_elevation(self):
        """
        Test that increasing max_elevation_angle (holding everything else
        fixed) shrinks (or holds equal) the inner "cone of silence" ground
        range, since a steeper top tilt reaches a given height sooner.
        """
        inner_ranges = [
            compute_radar_ground_range_bounds(
                0.5, max_elevation_angle, 230000, 3048, 0
            )[0]
            for max_elevation_angle in (5, 10, 19.5, 45)
        ]
        self.assertEqual(inner_ranges, sorted(inner_ranges, reverse=True))

    def test_compute_radar_ground_range_bounds_unreachable_height(self):
        """
        Test that a target too high to be reached by the highest elevation
        angle within max_range yields no coverage (None).
        """
        self.assertIsNone(compute_radar_ground_range_bounds(0.5, 19.5, 1000, 100000, 0))

    def test_compute_radar_ground_range_bounds_degenerate_single_tilt(self):
        """
        Test that a degenerate single scanned elevation angle (min equals
        max) yields no coverage (a zero-width ring) for an elevated target.
        """
        self.assertIsNone(compute_radar_ground_range_bounds(0.5, 0.5, 230000, 3048, 0))

    def test_compute_radar_beam_height_sanity_bound_at_max_range(self):
        """
        Sanity check (self-derived from the stated beam-height equation,
        not sourced from an external published chart) that NEXRAD's lowest
        (0.5 degree) tilt is several kilometers above ground by the time it
        reaches the conventional 230 km maximum range -- consistent with
        the well-known qualitative behavior (commonly cited as ~15,000-
        18,000 ft AGL) that NEXRAD's lowest tilt significantly overshoots
        low-level targets at long range.
        """
        height = compute_radar_beam_height(230000, 0.5)
        self.assertGreater(height, 3500)
        self.assertLess(height, 6500)

    def test_compute_terrain_elevation_angle_zero_at_station_height(self):
        """
        Test that a terrain point at exactly the station's own elevation
        and nonzero distance reads slightly below zero, due to the
        curvature-drop correction (a flat horizon appears to dip below
        local level at any positive distance).
        """
        angle = compute_terrain_elevation_angle(10000, 0, 0)
        self.assertLess(angle, 0)
        self.assertGreater(angle, -1)

    def test_compute_terrain_elevation_angle_zero_distance(self):
        """
        Test that a terrain point at zero distance but above the station
        elevation reads a 90 degree (straight up) elevation angle.
        """
        self.assertAlmostEqual(
            compute_terrain_elevation_angle(0, 100, 0), 90.0, delta=1e-6
        )

    def test_compute_terrain_elevation_angle_increases_with_height(self):
        """
        Test that, for a fixed ground distance, the elevation angle
        increases monotonically with terrain height.
        """
        angles = [
            compute_terrain_elevation_angle(50000, height, 0)
            for height in (0, 500, 1000, 2000, 5000)
        ]
        self.assertEqual(angles, sorted(angles))

    def test_compute_terrain_elevation_angle_consistent_with_beam_height(self):
        """
        Test that compute_terrain_elevation_angle agrees (to within the
        small-angle approximation it uses) with the exact beam-height
        model: for a beam at a shallow elevation angle, the height it
        reaches at a given range should be recovered as approximately the
        same elevation angle when fed back through the terrain-angle
        formula.
        """
        for elevation_angle in (0.1, 0.5, 1.0, 2.0):
            slant_range = 50000
            height = compute_radar_beam_height(slant_range, elevation_angle, 0)
            recovered = compute_terrain_elevation_angle(slant_range, height, 0)
            self.assertAlmostEqual(recovered, elevation_angle, delta=0.01)
