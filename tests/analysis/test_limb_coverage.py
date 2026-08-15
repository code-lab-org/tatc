"""
Unit tests for the limb sounding coverage analysis functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from shapely.geometry import MultiPoint
from shapely.geometry import Point as ShapelyPoint
from skyfield.positionlib import Geocentric
from skyfield.units import Distance, Velocity

from tatc.analysis import ScanDirection, collect_limb_observations
from tatc.analysis.limb_coverage import (
    _constant_rate_scan_fractions,
    _default_scan_direction,
    _interpolate_limb_point,
    _limb_tangent_point,
    _sample_limb_scan,
)
from tatc.constants import EARTH_MEAN_RADIUS, timescale

from .common import IssConstellationTestCase


def _make_geocentric(position_m, velocity_m_per_s, t):
    """
    Builds a synthetic Geocentric position/velocity for direct,
    deterministic control over the geometry helper functions, bypassing
    real orbit propagation.
    """
    position_m = np.array(position_m, dtype=float)
    velocity_m_per_s = np.array(velocity_m_per_s, dtype=float)
    return Geocentric(
        Distance(m=position_m).au,
        Velocity(km_per_s=velocity_m_per_s / 1000.0).au_per_d,
        t,
    )


class TestLimbTangentPoint(unittest.TestCase):
    """
    Unit tests for `_limb_tangent_point`.
    """

    def setUp(self):
        self.t = timescale.from_datetime(datetime(2024, 1, 1, tzinfo=timezone.utc))
        # circular-orbit-like configuration (position along +x, velocity
        # along +y): V along +y, N (orbit normal) along +z, B along +x,
        # coinciding with the position direction exactly in this case.
        self.sat_altitude = 700e3
        self.r_sat = EARTH_MEAN_RADIUS + self.sat_altitude
        self.sat = _make_geocentric(
            [[self.r_sat], [0], [0]], [[0], [7.5e3], [0]], self.t
        )

    def test_zero_elevation_azimuth_looks_along_velocity(self):
        """
        Test that scanning for a target elevation equal to the satellite's
        own altitude (a purely horizontal view, el=0) with scan_azimuth=0
        (looking along velocity) yields a tangent point at the satellite's
        own position -- the degenerate case where the ray is tangent to
        the sphere of the satellite's own radius exactly at the satellite.
        """
        tp_p, in_domain = _limb_tangent_point(self.sat, 0.0, [self.sat_altitude])
        self.assertTrue(in_domain[0])
        np.testing.assert_allclose(
            tp_p.ravel(), [self.r_sat, 0, 0], atol=1e-6
        )

    def test_tangent_point_radius_matches_analytic_formula(self):
        """
        Cross-checks the tangent point's geocentric radius against the
        analytic spherical-Earth relationship r_tp = R_sat * cos(el) for
        several target elevations, independent of scan_azimuth.
        """
        for azimuth in [0.0, 45.0, 90.0, 180.0, 270.0]:
            for target_elevation in [0.0, 100e3, 300e3]:
                with self.subTest(azimuth=azimuth, target_elevation=target_elevation):
                    tp_p, in_domain = _limb_tangent_point(
                        self.sat, azimuth, [target_elevation]
                    )
                    self.assertTrue(in_domain[0])
                    expected_radius = (EARTH_MEAN_RADIUS + target_elevation)
                    np.testing.assert_allclose(
                        np.linalg.norm(tp_p.ravel()), expected_radius, rtol=1e-9
                    )

    def test_target_above_satellite_altitude_is_out_of_domain(self):
        """
        Test that a requested elevation at or above the satellite's own
        altitude has no valid viewing angle and is flagged out of domain.
        """
        _, in_domain = _limb_tangent_point(
            self.sat, 90.0, [self.sat_altitude + 1e3]
        )
        self.assertFalse(in_domain[0])

    def test_cross_track_azimuth_stays_perpendicular_to_velocity(self):
        """
        Test that a scan_azimuth=90 (cross-track) tangent point lies in
        the plane orthogonal to velocity (here the x-z plane, since
        velocity is along +y), i.e. has no y-component.
        """
        tp_p, in_domain = _limb_tangent_point(self.sat, 90.0, [0.0])
        self.assertTrue(in_domain[0])
        self.assertAlmostEqual(tp_p.ravel()[1], 0.0, places=6)

    def test_vectorizes_across_multiple_samples(self):
        """
        Test that passing multiple satellite samples and target elevations
        (one column of `sat_pv` per target elevation) produces the same
        results as calling the function once per sample.
        """
        sat_multi = _make_geocentric(
            [[self.r_sat, self.r_sat], [0, 0], [0, 0]],
            [[0, 0], [7.5e3, 7.5e3], [0, 0]],
            timescale.from_datetimes([self.t.utc_datetime(), self.t.utc_datetime()]),
        )
        tp_p, in_domain = _limb_tangent_point(sat_multi, 90.0, [0.0, 100e3])
        tp_p_0, in_domain_0 = _limb_tangent_point(self.sat, 90.0, [0.0])
        tp_p_1, in_domain_1 = _limb_tangent_point(self.sat, 90.0, [100e3])
        np.testing.assert_allclose(tp_p[:, 0], tp_p_0.ravel(), atol=1e-6)
        np.testing.assert_allclose(tp_p[:, 1], tp_p_1.ravel(), atol=1e-6)
        self.assertEqual(list(in_domain), [in_domain_0[0], in_domain_1[0]])


class TestConstantRateScanFractions(unittest.TestCase):
    """
    Unit tests for `_constant_rate_scan_fractions`.
    """

    def test_uniform_values_give_uniform_fractions(self):
        """
        Test that evenly-spaced values are reached at evenly-spaced
        fractions of elapsed time.
        """
        fractions = _constant_rate_scan_fractions(np.array([0.0, 1.0, 2.0, 3.0]))
        np.testing.assert_allclose(fractions, [0.0, 1 / 3, 2 / 3, 1.0])

    def test_nonuniform_values_give_fractions_proportional_to_change(self):
        """
        Test that unevenly-spaced values are reached at fractions
        proportional to their cumulative absolute change, not their
        position in the sequence -- a value reached after a large jump
        arrives at a proportionally later fraction than one reached after
        a small jump.
        """
        # jumps of 1, 4, 5 -> total 10 -> cumulative fractions 0, .1, .5, 1
        fractions = _constant_rate_scan_fractions(np.array([0.0, 1.0, 5.0, 10.0]))
        np.testing.assert_allclose(fractions, [0.0, 0.1, 0.5, 1.0])

    def test_uses_path_length_not_net_displacement(self):
        """
        Test that a non-monotonic sequence is timed by its total absolute
        path length (as a constant-speed process would take), not by net
        displacement -- e.g. a value that doubles back over already-
        covered ground still takes time to do so.
        """
        # up 5, down 5, up 5 -> total path length 15
        fractions = _constant_rate_scan_fractions(np.array([0.0, 5.0, 0.0, 5.0]))
        np.testing.assert_allclose(fractions, [0.0, 1 / 3, 2 / 3, 1.0])

    def test_all_equal_values_falls_back_to_even_spacing(self):
        """
        Test that a sequence with no change at all (a degenerate case with
        no meaningful rate to infer) falls back to even time spacing
        rather than dividing by zero.
        """
        fractions = _constant_rate_scan_fractions(np.array([5.0, 5.0, 5.0]))
        np.testing.assert_allclose(fractions, [0.0, 0.5, 1.0])

    def test_single_value_is_a_single_zero_fraction(self):
        """
        Test that a single-value sequence is trivially reached at
        fraction zero.
        """
        fractions = _constant_rate_scan_fractions(np.array([5.0]))
        np.testing.assert_allclose(fractions, [0.0])


class TestDefaultScanDirection(unittest.TestCase):
    """
    Unit tests for `_default_scan_direction`.
    """

    def test_forward_azimuth_defaults_upward(self):
        """
        Test that azimuths closer to 0 deg (forward-looking) than 180 deg
        default to an upward scan.
        """
        for azimuth in [0.0, 1.0, 45.0, 89.0, -45.0, 271.0]:
            with self.subTest(azimuth=azimuth):
                self.assertEqual(
                    _default_scan_direction(azimuth), ScanDirection.UPWARD
                )

    def test_rearward_azimuth_defaults_downward(self):
        """
        Test that azimuths closer to 180 deg (rearward-looking) than 0 deg
        default to a downward scan.
        """
        for azimuth in [180.0, 179.0, 91.0, 135.0, 225.0, -135.0]:
            with self.subTest(azimuth=azimuth):
                self.assertEqual(
                    _default_scan_direction(azimuth), ScanDirection.DOWNWARD
                )

    def test_exactly_sideways_defaults_upward(self):
        """
        Test that azimuths exactly equidistant from 0 and 180 deg (90 and
        270 deg) default to upward, per the documented tie-break.
        """
        self.assertEqual(_default_scan_direction(90.0), ScanDirection.UPWARD)
        self.assertEqual(_default_scan_direction(270.0), ScanDirection.UPWARD)
        self.assertEqual(_default_scan_direction(-90.0), ScanDirection.UPWARD)


class TestInterpolateLimbPoint(unittest.TestCase):
    """
    Unit tests for `_interpolate_limb_point`.
    """

    def test_interpolates_at_exact_crossing(self):
        """
        Test that the interpolated point falls exactly halfway between two
        samples when the sample elevation is exactly halfway between their
        elevations.
        """
        points = [
            {
                "time": datetime(2024, 1, 1, tzinfo=timezone.utc),
                "longitude": 0.0,
                "latitude": 0.0,
                "elevation": 0.0,
            },
            {
                "time": datetime(2024, 1, 1, 0, 0, 10, tzinfo=timezone.utc),
                "longitude": 10.0,
                "latitude": 10.0,
                "elevation": 100.0,
            },
        ]
        result = _interpolate_limb_point(points, 50.0)
        self.assertAlmostEqual(result["longitude"], 5.0)
        self.assertAlmostEqual(result["latitude"], 5.0)
        self.assertAlmostEqual(result["elevation"], 50.0)
        self.assertEqual(
            result["time"], datetime(2024, 1, 1, 0, 0, 5, tzinfo=timezone.utc)
        )

    def test_clamps_to_nearest_endpoint_when_never_crossed(self):
        """
        Test that when the sample elevation is never crossed by the
        profile, the result clamps to the nearest endpoint rather than
        extrapolating.
        """
        points = [
            {
                "time": datetime(2024, 1, 1, tzinfo=timezone.utc),
                "longitude": 0.0,
                "latitude": 0.0,
                "elevation": 100.0,
            },
            {
                "time": datetime(2024, 1, 1, 0, 0, 10, tzinfo=timezone.utc),
                "longitude": 10.0,
                "latitude": 10.0,
                "elevation": 80.0,
            },
        ]
        result = _interpolate_limb_point(points, 50.0)
        self.assertEqual(result["longitude"], 10.0)
        self.assertEqual(result["elevation"], 80.0)

    def test_longitude_interpolates_along_shortest_path_across_antimeridian(self):
        """
        Test that longitude interpolation wraps across the antimeridian
        along the shortest path, rather than the long way around.
        """
        points = [
            {
                "time": datetime(2024, 1, 1, tzinfo=timezone.utc),
                "longitude": 170.0,
                "latitude": 0.0,
                "elevation": 100.0,
            },
            {
                "time": datetime(2024, 1, 1, 0, 0, 10, tzinfo=timezone.utc),
                "longitude": -170.0,
                "latitude": 0.0,
                "elevation": 0.0,
            },
        ]
        result = _interpolate_limb_point(points, 50.0)
        self.assertAlmostEqual(abs(result["longitude"]), 180.0)


class TestSampleLimbScanTiming(IssConstellationTestCase):
    """
    Integration tests verifying `_sample_limb_scan`'s constant-angular-rate
    timing model.
    """

    def test_larger_angular_gap_gets_more_elapsed_time(self):
        """
        Test that, for elevations with a much larger jump in viewing angle
        between the first two samples than the last two, the elapsed time
        for that first jump is correspondingly much larger -- unlike an
        "equal time per requested elevation" model, which would have split
        the scan duration evenly regardless of angular spacing.
        """
        start = datetime(2022, 6, 20, tzinfo=timezone.utc)
        scan_duration = timedelta(seconds=30)
        # a large angular gap from 0 to 300 km, followed by a much smaller
        # one from 300 to 350 km (for this ISS-like orbit: roughly 9.3 deg
        # vs 2.5 deg)
        points = _sample_limb_scan(
            self.satellite, start, 90.0, [0.0, 300e3, 350e3], scan_duration
        )
        self.assertEqual(len(points), 3)
        first_gap = (points[1]["time"] - points[0]["time"]).total_seconds()
        second_gap = (points[2]["time"] - points[1]["time"]).total_seconds()
        self.assertGreater(first_gap, 3.5 * second_gap)
        self.assertAlmostEqual(
            first_gap + second_gap, scan_duration.total_seconds(), places=6
        )


class TestCollectLimbObservations(IssConstellationTestCase):
    """
    Integration tests for `collect_limb_observations`.
    """

    def setUp(self):
        super().setUp()
        self.times = [
            datetime(2022, 6, 20, hour, tzinfo=timezone.utc) for hour in range(3)
        ]
        self.scan_elevations = [0.0, 20e3, 40e3, 60e3, 80e3, 100e3]
        self.scan_duration = timedelta(seconds=5)

    def test_collect_limb_observations_returns_results(self):
        """
        Test that limb observations are found for a satellite in an
        ISS-like orbit.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        self.assertFalse(results.empty)
        self.assertTrue((results.satellite == "Test").all())
        self.assertEqual(len(results), len(self.times))

    def test_collect_limb_observations_geometry_types(self):
        """
        Test that each scan's `geometry` is a MultiPoint (the full swept
        track) and `position` is a single Point (the interpolated
        representative sample), both 3D (lon/lat/elevation).
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        self.assertGreater(len(results), 0)
        for geometry in results.geometry:
            self.assertIsInstance(geometry, MultiPoint)
        for position in results.position:
            self.assertIsInstance(position, ShapelyPoint)
            self.assertTrue(position.has_z)

    def test_collect_limb_observations_point_count_matches_scan_elevations(self):
        """
        Test that each scan's track has one point per requested elevation,
        since every requested elevation is well below an ISS-altitude
        satellite's own altitude (all in domain).
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        for geometry in results.geometry:
            self.assertEqual(len(geometry.geoms), len(self.scan_elevations))

    def test_collect_limb_observations_start_end_bracket_time(self):
        """
        Test that the representative `time` always falls within the
        scan's `start`/`end` bounds.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        self.assertTrue((results.start <= results.time).all())
        self.assertTrue((results.time <= results.end).all())

    def test_collect_limb_observations_default_sample_elevation_is_midpoint(self):
        """
        Test that omitting `sample_elevation` selects a representative
        point near the midpoint of the requested `scan_elevations` range.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        midpoint = (min(self.scan_elevations) + max(self.scan_elevations)) / 2
        for position in results.position:
            self.assertAlmostEqual(position.z, midpoint, delta=1e3)

    def test_collect_limb_observations_all_elevations_above_altitude_is_empty(self):
        """
        Test that a scan requesting only elevations above the satellite's
        own altitude produces no rows at all.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=[10_000e3, 20_000e3],
            scan_duration=self.scan_duration,
        )
        self.assertTrue(results.empty)

    def test_collect_limb_observations_empty_for_no_times(self):
        """
        Test that requesting no scan start times produces an empty result.
        """
        results = collect_limb_observations(
            self.satellite,
            [],
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        self.assertTrue(results.empty)

    def test_collect_limb_observations_sorted_by_time(self):
        """
        Test that results are sorted by the representative `time`.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        self.assertTrue((results.time.values == np.sort(results.time.values)).all())

    def test_collect_limb_observations_bottom_to_top_ascends(self):
        """
        Test that a UPWARD scan visits tangent points
        in ascending elevation order.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
            scan_direction=ScanDirection.UPWARD,
        )
        for geometry in results.geometry:
            elevations = [point.z for point in geometry.geoms]
            self.assertEqual(elevations, sorted(elevations))

    def test_collect_limb_observations_top_to_bottom_descends(self):
        """
        Test that a DOWNWARD scan visits tangent points in descending
        elevation order -- the reverse of UPWARD.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
            scan_direction=ScanDirection.DOWNWARD,
        )
        for geometry in results.geometry:
            elevations = [point.z for point in geometry.geoms]
            self.assertEqual(elevations, sorted(elevations, reverse=True))

    def test_collect_limb_observations_scan_direction_ignores_input_order(self):
        """
        Test that scan_direction determines sweep order regardless of the
        order scan_elevations was originally given in -- an unsorted input
        list still produces a properly ordered scan.
        """
        shuffled = [60e3, 0.0, 100e3, 20e3, 80e3, 40e3]
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=shuffled,
            scan_duration=self.scan_duration,
            scan_direction=ScanDirection.UPWARD,
        )
        for geometry in results.geometry:
            elevations = [point.z for point in geometry.geoms]
            self.assertEqual(elevations, sorted(elevations))

    def test_collect_limb_observations_start_is_earliest_regardless_of_direction(self):
        """
        Test that `start`/`end` always reflect the chronologically first
        and last sample -- for DOWNWARD this is the highest elevation
        first, not the lowest.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=90.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
            scan_direction=ScanDirection.DOWNWARD,
        )
        for _, row in results.iterrows():
            points = list(row.geometry.geoms)
            self.assertEqual(points[0].z, max(p.z for p in points))
            self.assertEqual(points[-1].z, min(p.z for p in points))

    def test_collect_limb_observations_default_direction_from_forward_azimuth(self):
        """
        Test that omitting scan_direction for a forward-looking azimuth
        (0 deg, e.g. MLS) defaults to an upward (ascending) scan.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=0.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        for geometry in results.geometry:
            elevations = [point.z for point in geometry.geoms]
            self.assertEqual(elevations, sorted(elevations))

    def test_collect_limb_observations_default_direction_from_rearward_azimuth(self):
        """
        Test that omitting scan_direction for a rearward-looking azimuth
        (180 deg) defaults to a downward (descending) scan.
        """
        results = collect_limb_observations(
            self.satellite,
            self.times,
            scan_azimuth=180.0,
            scan_elevations=self.scan_elevations,
            scan_duration=self.scan_duration,
        )
        for geometry in results.geometry:
            elevations = [point.z for point in geometry.geoms]
            self.assertEqual(elevations, sorted(elevations, reverse=True))


if __name__ == "__main__":
    unittest.main()
