"""
Unit tests for the tangent point geometry functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import numpy as np
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc.utils.tangent_point import (
    _ellipsoidal_tangent_point,
    _geodetic_altitude,
    _itrs_rotation,
)
from tatc.constants import timescale


class TestGeodeticAltitude(unittest.TestCase):
    """
    Unit tests for `_geodetic_altitude`.
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
            _geodetic_altitude(position_itrs), elevations, atol=1e-6
        )
        # a single point converted via Skyfield's GCRS path also agrees
        point = wgs84.latlon(45.0, 10.0, elevation_m=80e3).at(t)
        np.testing.assert_allclose(
            _geodetic_altitude(itrs.rotation_at(t) @ point.position.m), 80e3, atol=1e-6
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
        altitudes = _geodetic_altitude(line_itrs.reshape(3, -1)).reshape(200, -1)
        tp_altitude = _geodetic_altitude(np.einsum("ij...,j...->i...", rotation, tp))
        # within one 10 m sample of the minimum, and no higher than it
        self.assertLess(
            np.max(np.abs(offsets[np.argmin(altitudes, axis=1)])), 10.0 + 1e-6
        )
        self.assertLess(np.max(tp_altitude - altitudes.min(axis=1)), 1e-3)
        # the tangent point lies on the line
        residual = (tp - rx) - np.einsum("ij,ij->j", tp - rx, d) * d
        np.testing.assert_allclose(np.linalg.norm(residual, axis=0), 0, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
