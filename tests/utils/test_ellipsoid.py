"""
Unit tests for the WGS 84 ellipsoid geometry functions in tatc.utils.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import numpy as np
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc.utils.ellipsoid import (
    _ellipsoidal_tangent_point,
    _itrs_rotation,
    compute_ellipsoid_intersection,
    compute_tangent_point,
    geodetic_to_rectangular,
    rectangular_to_geodetic,
)
from tatc.constants import timescale


class TestGeodeticAltitude(unittest.TestCase):
    """
    Unit tests for the geodetic altitude from `rectangular_to_geodetic`.
    """

    def test_matches_skyfield_wgs84(self):
        """
        Test that geodetic altitudes match Skyfield's WGS 84 conversion
        for points spanning all latitudes (including the poles) and
        altitudes from well below the surface (as for the straight-line
        tangent points of radio occultations) to low Earth orbit.
        """
        t = timescale.from_datetime(datetime(2024, 1, 1, tzinfo=timezone.utc))
        rng = np.random.default_rng(0)
        latitudes = np.concatenate([rng.uniform(-90, 90, 200), [-90, 90, 0]])
        longitudes = rng.uniform(-180, 180, latitudes.size)
        elevations = rng.uniform(-250e3, 1000e3, latitudes.size)
        position_itrs = np.array(
            [
                wgs84.latlon(lat, lon, elevation_m=h).itrs_xyz.m
                for lat, lon, h in zip(latitudes, longitudes, elevations)
            ]
        ).T
        np.testing.assert_allclose(
            rectangular_to_geodetic(position_itrs)[2], elevations, atol=1e-6
        )
        # a single point converted via Skyfield's GCRS path also agrees
        point = wgs84.latlon(45.0, 10.0, elevation_m=80e3).at(t)
        np.testing.assert_allclose(
            rectangular_to_geodetic(itrs.rotation_at(t) @ point.position.m)[2],
            80e3,
            atol=1e-6,
        )


class TestItrsRotation(unittest.TestCase):
    """
    Unit tests for `_itrs_rotation`.
    """

    def test_scalar_time_broadcasts_like_array_time(self):
        """
        Test that a scalar time gets a trailing axis of length one, giving
        the same rotated positions as an equivalent one-element time array.
        """
        dt = datetime(2024, 1, 1, tzinfo=timezone.utc)
        p = np.array([[7000e3], [1000e3], [-2000e3]])
        scalar = _itrs_rotation(timescale.from_datetime(dt))
        array = _itrs_rotation(timescale.from_datetimes([dt]))
        self.assertEqual(scalar.shape, (3, 3, 1))
        np.testing.assert_allclose(
            np.einsum("ij...,j...->i...", scalar, p),
            np.einsum("ij...,j...->i...", array, p),
        )


class TestEllipsoidalTangentPoint(unittest.TestCase):
    """
    Unit tests for `_ellipsoidal_tangent_point`.
    """

    def test_matches_brute_force_minimum_geodetic_altitude(self):
        """
        Test that, for random LEO-GNSS-like lines with tangent heights from
        -250 to +100 km (the range of radio occultation straight-line
        heights) at random locations, the tangent point is the point of
        minimum geodetic altitude found by densely sampling each line.
        """
        t = timescale.from_datetimes([datetime(2024, 1, 1, tzinfo=timezone.utc)] * 200)
        rotation = _itrs_rotation(t)
        rng = np.random.default_rng(1)
        up = rng.normal(size=(3, 200))
        up /= np.linalg.norm(up, axis=0)
        side = np.cross(up, rng.normal(size=(3, 200)), axis=0)
        side /= np.linalg.norm(side, axis=0)
        guess = up * (6371e3 + rng.uniform(-250e3, 100e3, 200))
        rx, tx = guess - side * 2800e3, guess + side * 25500e3
        tp = _ellipsoidal_tangent_point(rx, tx - rx, t)

        d = (tx - rx) / np.linalg.norm(tx - rx, axis=0)
        offsets = np.linspace(-30e3, 30e3, 6001)
        line = tp[:, :, np.newaxis] + d[:, :, np.newaxis] * offsets
        line_itrs = np.einsum("ijn,jnk->ink", rotation, line)
        altitudes = rectangular_to_geodetic(line_itrs.reshape(3, -1))[2].reshape(
            200, -1
        )
        tp_altitude = rectangular_to_geodetic(
            np.einsum("ij...,j...->i...", rotation, tp)
        )[2]
        # within one 10 m sample of the minimum, and no higher than it
        self.assertLess(
            np.max(np.abs(offsets[np.argmin(altitudes, axis=1)])), 10.0 + 1e-6
        )
        self.assertLess(np.max(tp_altitude - altitudes.min(axis=1)), 1e-3)
        # the tangent point lies on the line
        residual = (tp - rx) - np.einsum("ij,ij->j", tp - rx, d) * d
        np.testing.assert_allclose(np.linalg.norm(residual, axis=0), 0, atol=1e-6)


class TestGeodeticConversions(unittest.TestCase):
    """
    Unit tests for `geodetic_to_rectangular` and `rectangular_to_geodetic`.
    """

    def setUp(self):
        rng = np.random.default_rng(0)
        self.longitude = rng.uniform(-180, 180, 500)
        self.latitude = np.concatenate([rng.uniform(-90, 90, 496), [-90, 90, 0, 45]])
        self.elevation = rng.uniform(-1000, 900e3, 500)

    def test_geodetic_to_rectangular_matches_skyfield(self):
        """
        Test that rectangular coordinates match Skyfield's WGS 84 positions.
        """
        position = geodetic_to_rectangular(
            self.longitude, self.latitude, self.elevation
        )
        expected = wgs84.latlon(
            self.latitude, self.longitude, self.elevation
        ).itrs_xyz.m
        np.testing.assert_allclose(position, expected, atol=1e-6)

    def test_round_trip(self):
        """
        Test that converting to rectangular and back to geodetic coordinates
        recovers the coordinates, including at the poles.
        """
        longitude, latitude, elevation = rectangular_to_geodetic(
            geodetic_to_rectangular(self.longitude, self.latitude, self.elevation)
        )
        polar = np.abs(self.latitude) == 90
        np.testing.assert_allclose(longitude[~polar], self.longitude[~polar], atol=1e-9)
        np.testing.assert_allclose(latitude, self.latitude, atol=1e-9)
        np.testing.assert_allclose(elevation, self.elevation, atol=1e-6)

    def test_scalar_position(self):
        """
        Test that a single position (shape (3,)) converts to scalars.
        """
        longitude, latitude, elevation = rectangular_to_geodetic(
            geodetic_to_rectangular(-74.0, 40.7, 10.0)
        )
        self.assertAlmostEqual(float(longitude), -74.0)
        self.assertAlmostEqual(float(latitude), 40.7)
        self.assertAlmostEqual(float(elevation), 10.0, places=6)


class TestComputeEllipsoidIntersection(unittest.TestCase):
    """
    Unit tests for `compute_ellipsoid_intersection`.
    """

    def test_ray_toward_surface_point(self):
        """
        Test that rays toward points on the surface (at an elevation)
        intersect it at those points, and rays away from it miss.
        """
        target = geodetic_to_rectangular([0, 100, -60], [0, 45, -80], 2000)
        origin = geodetic_to_rectangular([5, 95, -50], [10, 40, -70], 700e3)
        points, found = compute_ellipsoid_intersection(origin, target - origin, 2000)
        self.assertTrue(np.all(found))
        _, _, elevation = rectangular_to_geodetic(points)
        # the ellipsoid of semi-axes extended by an elevation approximates
        # the surface at that elevation to within centimeters
        np.testing.assert_allclose(elevation, 2000, atol=0.1)
        np.testing.assert_allclose(points, target, atol=1.0)
        points, found = compute_ellipsoid_intersection(origin, origin - target, 2000)
        self.assertFalse(np.any(found))
        np.testing.assert_allclose(points, origin)


class TestComputeTangentPoint(unittest.TestCase):
    """
    Unit tests for `compute_tangent_point`.
    """

    def test_minimum_geodetic_altitude_on_line(self):
        """
        Test that the tangent point (Earth-fixed) lies on its line, at the
        line's minimum geodetic altitude.
        """
        rng = np.random.default_rng(1)
        guess = geodetic_to_rectangular(
            rng.uniform(-180, 180, 50), rng.uniform(-85, 85, 50), -50e3
        )
        side = np.cross(guess, rng.normal(size=(3, 50)), axis=0)
        side /= np.linalg.norm(side, axis=0)
        position, direction = guess - side * 2800e3, side * 28300e3
        tp = compute_tangent_point(position, direction)
        d = direction / np.linalg.norm(direction, axis=0)
        offsets = np.linspace(-30e3, 30e3, 6001)
        line = tp[:, :, np.newaxis] + d[:, :, np.newaxis] * offsets
        altitudes = rectangular_to_geodetic(line.reshape(3, -1))[2].reshape(50, -1)
        self.assertLess(
            np.max(np.abs(offsets[np.argmin(altitudes, axis=1)])), 10.0 + 1e-6
        )
        self.assertLess(
            np.max(rectangular_to_geodetic(tp)[2] - altitudes.min(axis=1)), 1e-3
        )
        residual = (tp - position) - np.einsum("ij,ij->j", tp - position, d) * d
        np.testing.assert_allclose(np.linalg.norm(residual, axis=0), 0, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
