"""
Unit tests for the sun glint utilities in tatc.utils.glint.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from skyfield.api import wgs84

from tatc.constants import de421
from tatc.schemas import CircularOrbit
from tatc.utils import compute_glint_angle, compute_specular_point


class TestGlint(unittest.TestCase):
    """
    Unit tests for compute_specular_point and compute_glint_angle.
    """

    def setUp(self):
        self.start = datetime(2026, 3, 20, tzinfo=timezone.utc)
        self.orbit = CircularOrbit(
            altitude=705e3, inclination=98.2, epoch=self.start
        ).to_gp_orbit()
        self.times = [self.start + timedelta(minutes=m) for m in range(0, 100, 2)]
        self.track = self.orbit.get_orbit_track(self.times)

    def test_specular_point_has_zero_glint_angle(self):
        """
        Test that the glint angle vanishes at the specular point, which exists
        while the satellite is over the day side, and not otherwise.
        """
        point = compute_specular_point(self.track)
        found = np.isfinite(point.latitude.degrees)
        self.assertGreater(np.sum(found), 0)
        self.assertGreater(np.sum(~found), 0)
        angle = compute_glint_angle(self.track, point)
        np.testing.assert_array_less(angle[found], 1e-4)
        subpoint = wgs84.subpoint_of(self.track)
        solar_elevation = (
            (de421["earth"] + subpoint)
            .at(self.track.t)
            .observe(de421["sun"])
            .apparent()
            .altaz()[0]
            .degrees
        )
        self.assertTrue(np.all(found[solar_elevation > 5]))
        self.assertFalse(np.any(found[solar_elevation < -30]))

    def test_specular_point_minimizes_glint_angle(self):
        """
        Test that no nearby point has a smaller glint angle than the
        specular point.
        """
        point = compute_specular_point(self.track)
        k = int(np.flatnonzero(np.isfinite(point.latitude.degrees))[3])
        lat, lon = point.latitude.degrees[k], point.longitude.degrees[k]
        offsets = np.linspace(-1, 1, 21)
        grid_lat, grid_lon = np.meshgrid(lat + offsets, lon + offsets)
        track = self.orbit.get_orbit_track([self.times[k]] * grid_lat.size)
        angles = compute_glint_angle(
            track, wgs84.latlon(grid_lat.ravel(), grid_lon.ravel())
        )
        # the smallest glint angle of the grid is at its center (the specular point)
        self.assertEqual(np.argmin(angles), grid_lat.size // 2)

    def test_equal_elevations(self):
        """
        Test that the Sun and the satellite are at the same elevation angle
        seen from the specular point, in opposite azimuths.
        """
        point = compute_specular_point(self.track)
        found = np.flatnonzero(np.isfinite(point.latitude.degrees))
        track = self.track[found]
        target = wgs84.latlon(
            point.latitude.degrees[found], point.longitude.degrees[found]
        )
        sat_alt, sat_az, _ = (track - target.at(track.t)).altaz()
        sun_alt, sun_az, _ = (
            (de421["earth"] + target)
            .at(track.t)
            .observe(de421["sun"])
            .apparent()
            .altaz()
        )
        np.testing.assert_allclose(sat_alt.degrees, sun_alt.degrees, atol=0.01)
        # azimuths (ill-conditioned near the zenith) below 80 degrees elevation
        low = sun_alt.degrees < 80
        azimuth_difference = (sat_az.degrees[low] - sun_az.degrees[low]) % 360
        np.testing.assert_allclose(azimuth_difference, 180, atol=0.05)

    def test_scalar_inputs(self):
        """
        Test that a single time and target give scalar results.
        """
        track = self.track[0]
        point = compute_specular_point(track)
        self.assertEqual(np.ndim(point.latitude.degrees), 0)
        angle = compute_glint_angle(track, wgs84.subpoint_of(track))
        self.assertEqual(np.ndim(angle), 0)
        self.assertTrue(0 <= angle <= 180)


if __name__ == "__main__":
    unittest.main()
