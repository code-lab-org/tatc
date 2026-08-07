"""
Unit tests for the radio occultation (RO) coverage analysis functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from shapely.geometry import MultiPoint
from shapely.geometry import Point as ShapelyPoint
from skyfield.positionlib import Geocentric
from skyfield.units import Distance, Velocity

from tatc.analysis import collect_ro_observations
from tatc.analysis.ro_coverage import (
    _interpolate_ro_point,
    _receiver_frame_vectors,
    _sample_ro_arc,
    _tangent_point_geometry,
)
from tatc.constants import timescale
from tatc.schemas import GeneralPerturbationsOrbit, Instrument, Satellite


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


class TestReceiverFrameVectors(unittest.TestCase):
    """
    Unit tests for `_receiver_frame_vectors`.
    """

    def setUp(self):
        self.t = timescale.from_datetime(datetime(2024, 1, 1, tzinfo=timezone.utc))

    def test_circular_orbit_matches_hand_computed_frame(self):
        """
        Test the VNB frame for a simple circular-orbit-like configuration
        (position along +x, velocity along +y): V should point along the
        velocity (+y), N (orbit normal, r x v) along +z, and B (V x N,
        completing the right-handed frame) along +x -- coinciding with the
        position direction exactly in this case, since velocity is purely
        tangential.
        """
        rx = _make_geocentric([[7000e3], [0], [0]], [[0], [7.5e3], [0]], self.t)
        v_u, n_u, b_u = _receiver_frame_vectors(rx)
        np.testing.assert_allclose(v_u.ravel(), [0, 1, 0], atol=1e-12)
        np.testing.assert_allclose(n_u.ravel(), [0, 0, 1], atol=1e-12)
        np.testing.assert_allclose(b_u.ravel(), [1, 0, 0], atol=1e-12)

    def test_frame_is_exactly_orthonormal_for_eccentric_orbit(self):
        """
        Regression test: V, N, B must be an exact orthonormal frame, even
        when velocity has a radial component (i.e. is not purely
        tangential) -- unlike a previous implementation, which used the
        position unit vector (r-hat) in place of B, only exactly
        orthogonal to V for a circular orbit (for an eccentric orbit the
        two can differ by tens of degrees).
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 25544U 98067A   22171.11255782  .00008307  00000+0  15444-3 0  9992",
                "2 25544  51.6448 322.0970 0003980 282.3738 231.6559 15.49798078345636",
            ]
        )
        t = timescale.from_datetimes(
            [datetime(2022, 6, 20, i, tzinfo=timezone.utc) for i in range(5)]
        )
        track = orbit.to_gp_orbit().get_orbit_track_at_time(t)
        v_u, n_u, b_u = _receiver_frame_vectors(track)
        np.testing.assert_allclose(np.einsum("ij,ij->j", v_u, n_u), 0, atol=1e-9)
        np.testing.assert_allclose(np.einsum("ij,ij->j", v_u, b_u), 0, atol=1e-9)
        np.testing.assert_allclose(np.einsum("ij,ij->j", n_u, b_u), 0, atol=1e-9)
        np.testing.assert_allclose(np.linalg.norm(b_u, axis=0), 1, atol=1e-9)


class TestTangentPointGeometry(unittest.TestCase):
    """
    Unit tests for `_tangent_point_geometry`.
    """

    def setUp(self):
        self.t = timescale.from_datetime(datetime(2024, 1, 1, tzinfo=timezone.utc))

    def test_collinear_through_earth_center_is_intersecting(self):
        """
        Test that a receiver and transmitter on exactly opposite sides of
        the Earth (collinear through the center) produce a tangent point
        at the origin, correctly flagged as "intersecting" (sign < 0).
        """
        rx = _make_geocentric([[7000e3], [0], [0]], [[0], [7.5e3], [0]], self.t)
        tx = _make_geocentric([[-7000e3], [0], [0]], [[0], [-3.0e3], [0]], self.t)
        v_u, n_u, b_u = _receiver_frame_vectors(rx)
        tp_p, tp_v, tp_sign, _, _ = _tangent_point_geometry(tx, rx, v_u, n_u, b_u)
        np.testing.assert_allclose(tp_p.ravel(), [0, 0, 0], atol=1e-6)
        self.assertEqual(tp_sign[0], -1)
        self.assertIsNone(tp_v)

    def test_same_side_is_not_intersecting(self):
        """
        Test that a receiver and transmitter on the same side of the Earth
        (tangent point not between them) are correctly flagged as
        "parallel"/non-intersecting (sign > 0).
        """
        rx = _make_geocentric([[7000e3], [0], [0]], [[0], [7.5e3], [0]], self.t)
        tx = _make_geocentric([[8000e3], [2000e3], [0]], [[0], [3.0e3], [0]], self.t)
        v_u, n_u, b_u = _receiver_frame_vectors(rx)
        _, _, tp_sign, _, _ = _tangent_point_geometry(tx, rx, v_u, n_u, b_u)
        self.assertEqual(tp_sign[0], 1)

    def test_tangent_point_matches_independent_closest_point_formula(self):
        """
        Cross-checks the tangent point position against an independently
        (re-)implemented closest-point-on-a-line-to-the-origin formula:
        for the line through the transmitter with direction d = tx - rx,
        the point closest to the origin is tx - d*(tx.d)/(d.d).
        """
        rx = _make_geocentric([[7000e3], [500e3], [0]], [[0], [7.5e3], [1e3]], self.t)
        tx = _make_geocentric(
            [[-6000e3], [-1000e3], [3000e3]], [[1e3], [-3.0e3], [0]], self.t
        )
        v_u, n_u, b_u = _receiver_frame_vectors(rx)
        tp_p, _, _, _, _ = _tangent_point_geometry(tx, rx, v_u, n_u, b_u)
        tx_p = np.array(tx.position.m).ravel()
        rx_p = np.array(rx.position.m).ravel()
        d = tx_p - rx_p
        expected_tp = tx_p - d * np.dot(tx_p, d) / np.dot(d, d)
        np.testing.assert_allclose(tp_p.ravel(), expected_tp, rtol=1e-9)

    def test_tangent_point_velocity_matches_finite_difference(self):
        """
        Cross-checks the analytic tangent point velocity
        (`compute_velocity=True`, otherwise unused by any caller in the
        codebase) against a central finite difference of the tangent point
        position, using synthetic constant-velocity (straight-line) motion
        so the finite difference is essentially exact for a small enough
        time step.
        """
        dt = 0.001  # seconds
        rx_p = np.array([7000e3, 500e3, 0.0])
        rx_v = np.array([0.0, 7.5e3, 1e3])
        tx_p = np.array([-6000e3, -1000e3, 3000e3])
        tx_v = np.array([1e3, -3.0e3, 0.0])

        def tangent_point_position(offset):
            rx = _make_geocentric(
                (rx_p + rx_v * offset).reshape(3, 1), rx_v.reshape(3, 1), self.t
            )
            tx = _make_geocentric(
                (tx_p + tx_v * offset).reshape(3, 1), tx_v.reshape(3, 1), self.t
            )
            v_u, n_u, b_u = _receiver_frame_vectors(rx)
            return _tangent_point_geometry(tx, rx, v_u, n_u, b_u)[0]

        finite_diff_velocity = (
            tangent_point_position(dt) - tangent_point_position(-dt)
        ) / (2 * dt)

        rx0 = _make_geocentric(rx_p.reshape(3, 1), rx_v.reshape(3, 1), self.t)
        tx0 = _make_geocentric(tx_p.reshape(3, 1), tx_v.reshape(3, 1), self.t)
        v_u, n_u, b_u = _receiver_frame_vectors(rx0)
        _, tp_v, _, _, _ = _tangent_point_geometry(
            tx0, rx0, v_u, n_u, b_u, compute_velocity=True
        )
        np.testing.assert_allclose(
            tp_v.ravel(), finite_diff_velocity.ravel(), rtol=1e-6
        )


class TestInterpolateRoPoint(unittest.TestCase):
    """
    Unit tests for `_interpolate_ro_point`.
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
                "elevation": 100.0,
                "rx_tx_pitch": 0.0,
                "rx_tx_yaw": 0.0,
                "tp_tx_azimuth": 0.0,
            },
            {
                "time": datetime(2024, 1, 1, 0, 0, 10, tzinfo=timezone.utc),
                "longitude": 10.0,
                "latitude": 10.0,
                "elevation": 0.0,
                "rx_tx_pitch": 10.0,
                "rx_tx_yaw": 10.0,
                "tp_tx_azimuth": 10.0,
            },
        ]
        result = _interpolate_ro_point(points, 50.0)
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
                "rx_tx_pitch": 0.0,
                "rx_tx_yaw": 0.0,
                "tp_tx_azimuth": 0.0,
            },
            {
                "time": datetime(2024, 1, 1, 0, 0, 10, tzinfo=timezone.utc),
                "longitude": 10.0,
                "latitude": 10.0,
                "elevation": 80.0,
                "rx_tx_pitch": 10.0,
                "rx_tx_yaw": 10.0,
                "tp_tx_azimuth": 10.0,
            },
        ]
        # both elevations (100, 80) are above sample_elevation=50, so the
        # sign of (elevation - sample_elevation) never changes -- a
        # genuine non-crossing case (unlike sample_elevation=80, which
        # would exactly touch the second point's elevation and register as
        # a crossing via the sign-based diff detection, not exercising
        # this fallback at all). 50 is closer to the second point's diff
        # (30) than the first's (50), so it should clamp to the second.
        result = _interpolate_ro_point(points, 50.0)
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
                "rx_tx_pitch": 0.0,
                "rx_tx_yaw": 0.0,
                "tp_tx_azimuth": 0.0,
            },
            {
                "time": datetime(2024, 1, 1, 0, 0, 10, tzinfo=timezone.utc),
                "longitude": -170.0,
                "latitude": 0.0,
                "elevation": 0.0,
                "rx_tx_pitch": 0.0,
                "rx_tx_yaw": 0.0,
                "tp_tx_azimuth": 0.0,
            },
        ]
        result = _interpolate_ro_point(points, 50.0)
        # shortest path from 170 to -170 crosses the antimeridian (180),
        # so the midpoint should be +/-180, not 0
        self.assertAlmostEqual(abs(result["longitude"]), 180.0)


class TestCollectRoObservations(unittest.TestCase):
    """
    Integration tests for `collect_ro_observations`.
    """

    def setUp(self):
        self.receiver = Satellite(
            name="RX",
            orbit=GeneralPerturbationsOrbit.from_tle(
                [
                    "1 25544U 98067A   22171.11255782  .00008307  00000+0  15444-3 0  9992",
                    "2 25544  51.6448 322.0970 0003980 282.3738 231.6559 15.49798078345636",
                ]
            ),
            instruments=[Instrument(name="RO Receiver", field_of_regard=180.0)],
        )
        self.transmitter = Satellite(
            name="TX",
            orbit=GeneralPerturbationsOrbit.from_tle(
                [
                    "1 24876U 97035A   22171.50000000  .00000015  00000-0  00000-0 0  9990",
                    "2 24876  55.0000 100.0000 0100000  90.0000 270.0000  2.00561000123456",
                ]
            ),
            instruments=[Instrument(name="RO Transmitter", field_of_regard=180.0)],
        )
        self.start = datetime(2022, 6, 20, tzinfo=timezone.utc)
        self.end = self.start + timedelta(hours=6)

    def test_collect_ro_observations_returns_results(self):
        """
        Test that RO observations are found for a several-hour analysis
        period between a LEO receiver and a MEO transmitter.
        """
        results = collect_ro_observations(
            self.receiver, self.transmitter, self.start, self.end
        )
        self.assertFalse(results.empty)
        self.assertTrue((results.receiver == "RX").all())
        self.assertTrue((results.transmitter == "TX").all())

    def test_collect_ro_observations_empty_for_short_period(self):
        """
        Test that no observations are found for a period too short to
        contain any RO profile.
        """
        results = collect_ro_observations(
            self.receiver,
            self.transmitter,
            self.start,
            self.start + timedelta(seconds=1),
        )
        self.assertTrue(results.empty)

    def test_collect_ro_observations_geometry_types(self):
        """
        Test that each observation's `geometry` is a MultiPoint (the full
        sampled profile) and `position` is a single Point (the
        interpolated representative sample), both 3D (lon/lat/elevation).
        """
        results = collect_ro_observations(
            self.receiver, self.transmitter, self.start, self.end
        )
        self.assertGreater(len(results), 0)
        for geometry in results.geometry:
            self.assertIsInstance(geometry, MultiPoint)
        for position in results.position:
            self.assertIsInstance(position, ShapelyPoint)
            self.assertTrue(position.has_z)

    def test_collect_ro_observations_smaller_max_yaw_is_more_restrictive(self):
        """
        Test that decreasing `max_yaw` never increases the number of
        observations found (a tighter yaw bound can only exclude
        profiles, not add new ones).
        """
        permissive = collect_ro_observations(
            self.receiver, self.transmitter, self.start, self.end, max_yaw=65
        )
        restrictive = collect_ro_observations(
            self.receiver, self.transmitter, self.start, self.end, max_yaw=5
        )
        self.assertLessEqual(len(restrictive), len(permissive))

    def test_collect_ro_observations_respects_range_elevation(self):
        """
        Test that every returned profile's sampled tangent point track
        stays within the requested elevation range.
        """
        range_elevation = (-100e3, 40e3)
        results = collect_ro_observations(
            self.receiver,
            self.transmitter,
            self.start,
            self.end,
            range_elevation=range_elevation,
        )
        self.assertGreater(len(results), 0)
        for geometry in results.geometry:
            for point in geometry.geoms:
                self.assertGreaterEqual(point.z, range_elevation[0])
                self.assertLessEqual(point.z, range_elevation[1])

    def test_collect_ro_observations_multiple_transmitters(self):
        """
        Test that observations from multiple transmitters are combined
        and sorted by time.
        """
        transmitter_2 = Satellite(
            name="TX2",
            orbit=GeneralPerturbationsOrbit.from_tle(
                [
                    "1 24876U 97035A   22171.50000000  .00000015  00000-0  00000-0 0  9990",
                    "2 24876  55.0000  50.0000 0100000 200.0000  90.0000  2.00561000123456",
                ]
            ),
            instruments=[Instrument(name="RO Transmitter", field_of_regard=180.0)],
        )
        results = collect_ro_observations(
            self.receiver, [self.transmitter, transmitter_2], self.start, self.end
        )
        self.assertTrue(set(results.transmitter.unique()) <= {"TX", "TX2"})
        self.assertTrue((results.time.values == np.sort(results.time.values)).all())

    def test_collect_ro_observations_is_rising_is_boolean(self):
        """
        Test that `is_rising` is a boolean flag on every observation.
        """
        results = collect_ro_observations(
            self.receiver, self.transmitter, self.start, self.end
        )
        self.assertGreater(len(results), 0)
        self.assertTrue(results.is_rising.isin([True, False]).all())

    def test_sample_ro_arc_closes_at_arc_boundary_when_still_in_range(self):
        """
        Test that an observation correctly closes at the sampled arc's own
        boundary, not just via exiting the elevation range -- using an
        effectively unbounded elevation range to guarantee the tangent
        point never leaves it, so the only way the sampled window's single
        observation can end is by reaching the last sample.
        """
        observations = _sample_ro_arc(
            self.transmitter,
            self.receiver,
            self.start,
            self.start + timedelta(minutes=5),
            timedelta(seconds=30),
            (-1e9, 1e9),
        )
        self.assertEqual(len(observations), 1)
        self.assertEqual(len(observations[0]["points"]), 11)


if __name__ == "__main__":
    unittest.main()
