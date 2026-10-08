"""
Unit tests for the space-based (satellite-to-satellite) coverage analysis
functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
from shapely.geometry import MultiLineString
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric
from skyfield.units import Distance

from tatc.analysis import collect_space_observations
from tatc.analysis.space_sampling import (
    _compute_min_altitude,
    _get_space_residual,
)
from tatc.constants import EARTH_EQUATORIAL_RADIUS, EARTH_POLAR_RADIUS, timescale
from tatc.schemas import CircularOrbit, Satellite, WalkerConstellation
from tatc.utils import hash_geometry
from tatc.utils.ellipsoid import geodetic_to_rectangular


def _column(*positions):
    """
    Builds an array of positions (meters, shape (3, N)).
    """
    return np.array(positions, dtype=float).T


class TestSpaceResidual(unittest.TestCase):
    """
    Unit tests for `_get_space_residual` and `_compute_min_altitude`.
    """

    def test_opposite_sides_are_occluded(self):
        r = EARTH_EQUATORIAL_RADIUS + 500e3
        residual = _get_space_residual(_column([r, 0, 0]), _column([-r, 0, 0]))
        self.assertGreater(residual[0], 0)

    def test_nearby_are_visible(self):
        r = EARTH_EQUATORIAL_RADIUS + 500e3
        residual = _get_space_residual(_column([r, 0, 0]), _column([r, 100e3, 0]))
        self.assertLess(residual[0], 0)

    def test_equatorial_chord_altitude(self):
        # an equatorial chord between points 500 km above the equator,
        # separated by 30 degrees of longitude, has its tangent point midway
        # at a geocentric distance of r cos(15 deg)
        r = EARTH_EQUATORIAL_RADIUS + 500e3
        from_p = _column([r, 0, 0])
        to_p = _column([r * np.cos(np.radians(30)), r * np.sin(np.radians(30)), 0])
        expected = r * np.cos(np.radians(15)) - EARTH_EQUATORIAL_RADIUS
        self.assertAlmostEqual(
            _compute_min_altitude(from_p, to_p)[0], expected, delta=1
        )
        # near the minimum altitude, the constraint is violated by the difference
        for min_grazing_altitude in [expected - 5e3, expected + 5e3]:
            self.assertAlmostEqual(
                _get_space_residual(
                    from_p, to_p, min_grazing_altitude=min_grazing_altitude
                )[0],
                min_grazing_altitude - expected,
                delta=1,
            )
        # and elsewhere, by a bound of it with the correct sign
        self.assertLess(_get_space_residual(from_p, to_p)[0], 0)
        self.assertGreater(
            _get_space_residual(from_p, to_p, min_grazing_altitude=300e3)[0], 0
        )

    def test_tangent_point_beyond_segment_uses_endpoint(self):
        # radially aligned positions: the line's tangent point (the geocenter)
        # is beyond the segment, whose minimum altitude is the lower position
        from_p = geodetic_to_rectangular([10], [45], [400e3])
        to_p = geodetic_to_rectangular([10], [45], [2000e3])
        self.assertAlmostEqual(
            _compute_min_altitude(from_p, to_p)[0], 400e3, delta=1e-3
        )
        self.assertAlmostEqual(
            _get_space_residual(from_p, to_p, min_grazing_altitude=405e3)[0],
            5e3,
            delta=1e-3,
        )
        self.assertGreater(
            _get_space_residual(from_p, to_p, min_grazing_altitude=500e3)[0], 0
        )
        self.assertLess(
            _get_space_residual(from_p, to_p, min_grazing_altitude=300e3)[0], 0
        )

    def test_occlusion_disabled(self):
        r = EARTH_EQUATORIAL_RADIUS + 500e3
        residual = _get_space_residual(
            _column([r, 0, 0]), _column([-r, 0, 0]), min_grazing_altitude=None
        )
        np.testing.assert_array_equal(residual, [-1])

    def test_range_constraints(self):
        r = EARTH_EQUATORIAL_RADIUS + 500e3
        from_p = _column([r, 0, 0], [r, 0, 0])
        to_p = _column([r, 1000e3, 0], [r, 3000e3, 0])
        np.testing.assert_allclose(
            _get_space_residual(
                from_p, to_p, max_range=2000e3, min_grazing_altitude=None
            ),
            [-1000e3, 1000e3],
        )
        np.testing.assert_allclose(
            _get_space_residual(
                from_p, to_p, min_range=2000e3, min_grazing_altitude=None
            ),
            [1000e3, -1000e3],
        )
        # the largest violation of all constraints
        np.testing.assert_allclose(
            _get_space_residual(
                from_p,
                to_p,
                min_range=1500e3,
                max_range=2000e3,
                min_grazing_altitude=None,
            ),
            [500e3, 1000e3],
        )

    def test_bounds_decide_sign(self):
        # the residual has the sign of the exact minimum altitude constraint,
        # whether or not decided by the bounds of the closest approach
        rng = np.random.default_rng(0)
        size = 2000
        from_p = geodetic_to_rectangular(
            rng.uniform(-180, 180, size),
            rng.uniform(-90, 90, size),
            rng.uniform(300e3, 2000e3, size),
        )
        to_p = geodetic_to_rectangular(
            rng.uniform(-180, 180, size),
            rng.uniform(-90, 90, size),
            rng.uniform(300e3, 2000e3, size),
        )
        for min_grazing_altitude in [0, 100e3]:
            exact = min_grazing_altitude - _compute_min_altitude(from_p, to_p)
            residual = _get_space_residual(
                from_p, to_p, min_grazing_altitude=min_grazing_altitude
            )
            np.testing.assert_array_equal(residual <= 0, exact <= 0)
            # and is exact where not decided by the bounds
            band = np.abs(exact) < EARTH_EQUATORIAL_RADIUS - EARTH_POLAR_RADIUS
            self.assertTrue(np.any(band))
            self.assertTrue(np.all(np.abs(residual - exact) < 22e3))


def _make_satellite(name, altitude, inclination, raan=0, true_anomaly=0):
    return Satellite(
        name=name,
        orbit=CircularOrbit(
            altitude=altitude,
            inclination=inclination,
            right_ascension_ascending_node=raan,
            true_anomaly=true_anomaly,
        ),
    )


def _brute_force_visible(from_satellite, to_satellite, times, min_grazing_altitude):
    """
    Evaluates whether the line of sight between two satellites stays above a
    geodetic altitude at each time, independently of TAT-C's ellipsoid
    utilities: by sampling 2001 points along it with Skyfield.
    """
    t = timescale.from_datetimes(times)
    from_p = from_satellite.orbit.to_gp_orbit().get_orbit_track(times).position.m
    to_p = to_satellite.orbit.to_gp_orbit().get_orbit_track(times).position.m
    s = np.linspace(0, 1, 2001)
    points = from_p[:, :, np.newaxis] + (to_p - from_p)[:, :, np.newaxis] * s
    altitude = np.array(
        [
            wgs84.geographic_position_of(
                Geocentric(Distance(m=points[:, :, k]).au, None, t)
            ).elevation.m
            for k in range(len(s))
        ]
    )
    return np.min(altitude, axis=0) > min_grazing_altitude


class TestCollectSpaceObservations(unittest.TestCase):
    """
    Unit tests for `collect_space_observations`.
    """

    def setUp(self):
        self.start = datetime(2025, 1, 1, tzinfo=timezone.utc)
        self.end = self.start + timedelta(hours=3)
        self.leader = _make_satellite("leader", 500e3, 50)
        self.near = _make_satellite("near", 500e3, 50, true_anomaly=30)
        self.far = _make_satellite("far", 500e3, 50, true_anomaly=60)
        self.polar = _make_satellite("polar", 700e3, 97, raan=60)

    def test_same_orbit_phased_visible(self):
        observations = collect_space_observations(
            self.leader, self.start, self.end, self.near
        )
        self.assertEqual(len(observations), 1)
        self.assertEqual(observations.start.iloc[0], pd.Timestamp(self.start))
        self.assertEqual(observations.end.iloc[0], pd.Timestamp(self.end))
        self.assertEqual(observations.from_satellite.iloc[0], "leader")
        self.assertEqual(observations.to_satellite.iloc[0], "near")

    def test_same_orbit_phased_occluded(self):
        # occluded by the Earth
        self.assertEqual(
            len(
                collect_space_observations(self.leader, self.start, self.end, self.far)
            ),
            0,
        )
        # occluded by the atmosphere (the chord grazes at about 260-280 km)
        self.assertEqual(
            len(
                collect_space_observations(
                    self.leader,
                    self.start,
                    self.end,
                    self.near,
                    min_grazing_altitude=300e3,
                )
            ),
            0,
        )
        # beyond the maximum range (about 3.56 Mm)
        self.assertEqual(
            len(
                collect_space_observations(
                    self.leader, self.start, self.end, self.near, max_range=3000e3
                )
            ),
            0,
        )

    def test_matches_brute_force(self):
        min_grazing_altitude = 100e3
        observations = collect_space_observations(
            self.leader,
            self.start,
            self.end,
            self.polar,
            min_grazing_altitude=min_grazing_altitude,
        )
        self.assertGreater(len(observations), 1)
        times = [
            self.start + timedelta(seconds=float(x))
            for x in np.arange(0, (self.end - self.start).total_seconds(), 10)
        ]
        visible = _brute_force_visible(
            self.leader, self.polar, times, min_grazing_altitude
        )
        timestamps = pd.DatetimeIndex(times)
        inside = np.zeros(len(times), dtype=bool)
        near_boundary = np.zeros(len(times), dtype=bool)
        tolerance = pd.Timedelta(seconds=1)
        for _, row in observations.iterrows():
            inside |= (timestamps >= row.start) & (timestamps <= row.end)
            for boundary in (row.start, row.end):
                near_boundary |= np.abs(timestamps - boundary) < tolerance
        np.testing.assert_array_equal(inside[~near_boundary], visible[~near_boundary])

    def test_default_to_satellites(self):
        satellites = [self.leader, self.near, self.polar]
        observations = collect_space_observations(satellites, self.start, self.end)
        pairs = set(zip(observations.from_satellite, observations.to_satellite))
        # no satellite observes itself
        self.assertTrue(all(a != b for a, b in pairs))
        # every pair is observed in both directions, at the same times
        self.assertEqual(len(pairs), 6)
        for a, b in pairs:
            forward = observations[
                (observations.from_satellite == a) & (observations.to_satellite == b)
            ]
            backward = observations[
                (observations.from_satellite == b) & (observations.to_satellite == a)
            ]
            np.testing.assert_array_equal(forward.start, backward.start)
            np.testing.assert_array_equal(forward.end, backward.end)

    def test_equal_satellites_not_paired(self):
        copy = self.leader.model_copy(deep=True)
        observations = collect_space_observations(
            self.leader, self.start, self.end, [copy, self.near]
        )
        self.assertEqual(set(observations.to_satellite), {"near"})

    def test_output_format(self):
        observations = collect_space_observations(
            [self.leader, self.polar], self.start, self.end
        )
        self.assertGreater(len(observations), 0)
        self.assertEqual(
            list(observations.columns),
            [
                "target_hash",
                "geometry",
                "from_satellite",
                "to_satellite",
                "start",
                "end",
            ],
        )
        self.assertEqual(observations.crs, "EPSG:4326")
        self.assertEqual(str(observations.start.dtype), "datetime64[ns, UTC]")
        self.assertTrue(observations.start.is_monotonic_increasing)
        self.assertTrue(np.all(observations.end > observations.start))
        for _, row in observations.iterrows():
            self.assertIsInstance(row.geometry, MultiLineString)
            self.assertEqual(len(row.geometry.geoms), 2)
            self.assertTrue(row.geometry.has_z)
            self.assertEqual(row.target_hash, hash_geometry(row.geometry))
        # the lines of sight join the satellites' positions at the start and end
        row = observations.iloc[0]
        satellites = {"leader": self.leader, "polar": self.polar}
        for line, time in zip(row.geometry.geoms, [row.start, row.end]):
            for point, name in zip(line.coords, [row.from_satellite, row.to_satellite]):
                track = (
                    satellites[name]
                    .orbit.to_gp_orbit()
                    .get_orbit_track([time.to_pydatetime(warn=False)])
                )
                position = wgs84.geographic_position_of(track)
                np.testing.assert_allclose(
                    point,
                    [
                        position.longitude.degrees[0],
                        position.latitude.degrees[0],
                        position.elevation.m[0],
                    ],
                    atol=1e-3,
                )

    def test_empty(self):
        for observations in [
            collect_space_observations(self.leader, self.start, self.end, self.far),
            collect_space_observations(self.leader, self.start, self.start),
            collect_space_observations(self.leader, self.start, self.end),
        ]:
            self.assertEqual(len(observations), 0)
            self.assertEqual(
                list(observations.columns),
                [
                    "target_hash",
                    "geometry",
                    "from_satellite",
                    "to_satellite",
                    "start",
                    "end",
                ],
            )

    def test_input_errors(self):
        constellation = WalkerConstellation(
            name="walker",
            orbit=self.leader.orbit,
            number_satellites=4,
            number_planes=2,
        )
        with self.assertRaises(TypeError):
            collect_space_observations(constellation, self.start, self.end)
        with self.assertRaises(TypeError):
            collect_space_observations(self.leader, self.start, self.end, constellation)
        with self.assertRaises(ValueError):
            collect_space_observations(
                self.leader, self.start.replace(tzinfo=None), self.end
            )
        with self.assertRaises(ValueError):
            collect_space_observations(
                self.leader,
                self.start,
                self.end,
                self.near,
                min_range=2000e3,
                max_range=1000e3,
            )
        with self.assertRaises(ValueError):
            collect_space_observations(
                self.leader, self.start, self.end, self.near, min_duration=timedelta(0)
            )
