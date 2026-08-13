"""
Unit tests for the tatc.utils.observation module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest

from tatc.utils import (
    compute_field_of_regard,
    compute_max_access_time,
    compute_max_transit_time,
    compute_min_along_track_distance,
    compute_min_elevation_angle,
    field_of_regard_to_swath_width,
    swath_width_to_field_of_regard,
    swath_width_to_field_of_view,
)


class TestObservation(unittest.TestCase):  # pylint: disable=too-many-public-methods
    """
    Unit tests for the tatc.utils.observation module.
    """

    def test_swath_width_to_field_of_regard(self):
        """
        Test that the swath width can be converted to field of regard.
        """
        self.assertAlmostEqual(
            swath_width_to_field_of_regard(705000, 185815), 15.0, delta=0.001
        )

    def test_swath_width_to_field_of_regard_zero_swath(self):
        """
        Test that a zero swath width (nadir-only observation) requires no
        field of regard.
        """
        self.assertEqual(swath_width_to_field_of_regard(705000, 0), 0.0)

    def test_swath_width_to_field_of_regard_elevation(self):
        """
        Test that a positive elevation (observing a point above the mean
        Earth radius) changes the field of regard required for the same
        nominal swath width and altitude.
        """
        self.assertAlmostEqual(
            swath_width_to_field_of_regard(705000, 185815, 5000),
            15.105822,
            delta=1e-5,
        )

    def test_swath_width_to_field_of_regard_inverts_field_of_regard_to_swath_width(
        self,
    ):
        """
        Test that swath_width_to_field_of_regard is the inverse of
        field_of_regard_to_swath_width across a range of altitudes and
        field of regard values (as fractions of the horizon-limited maximum
        field of regard at each altitude, so every combination stays within
        the physically achievable range).
        """
        for altitude in (500000, 705000, 800000, 20000000):
            max_field_of_regard = compute_field_of_regard(altitude, 0)
            for fraction in (0.05, 0.2, 0.5, 0.9):
                field_of_regard = max_field_of_regard * fraction
                swath_width = field_of_regard_to_swath_width(altitude, field_of_regard)
                self.assertAlmostEqual(
                    swath_width_to_field_of_regard(altitude, swath_width),
                    field_of_regard,
                    delta=1e-6,
                )

    def test_swath_width_to_field_of_view_matches_field_of_regard_at_nadir(self):
        """
        Test that, at a look angle of zero (swath centered on nadir), the
        field of view matches the symmetric field of regard exactly.
        """
        self.assertAlmostEqual(
            swath_width_to_field_of_view(705000, 185815, 0),
            swath_width_to_field_of_regard(705000, 185815),
            delta=1e-9,
        )

    def test_swath_width_to_field_of_view_zero_swath(self):
        """
        Test that a zero swath width requires no field of view, regardless
        of look angle.
        """
        self.assertEqual(swath_width_to_field_of_view(705000, 0, 20), 0.0)

    def test_swath_width_to_field_of_view_decreases_with_look_angle(self):
        """
        Test that, for a fixed altitude and swath width, the field of view
        decreases as the off-nadir look angle increases (a fixed ground
        swath subtends a smaller angle as viewing becomes more oblique).
        """
        field_of_views = [
            swath_width_to_field_of_view(705000, 185815, look_angle)
            for look_angle in (0, 5, 10, 20, 30)
        ]
        self.assertEqual(field_of_views, sorted(field_of_views, reverse=True))

    def test_swath_width_to_field_of_view_off_nadir(self):
        """
        Test that the field of view can be computed for a specified
        off-nadir look angle.
        """
        self.assertAlmostEqual(
            swath_width_to_field_of_view(705000, 185815, 20),
            12.992488,
            delta=1e-5,
        )

    def test_swath_width_to_field_of_view_saturates_beyond_horizon(self):
        """
        Test that a look angle at or beyond the horizon-limited maximum
        saturates to the same field of view rather than raising an error.
        """
        max_look_angle = compute_field_of_regard(705000, 0) / 2
        self.assertAlmostEqual(
            swath_width_to_field_of_view(705000, 185815, max_look_angle),
            swath_width_to_field_of_view(705000, 185815, max_look_angle + 10),
            delta=1e-12,
        )

    def test_field_of_regard_to_swath_width(self):
        """
        Test that the field of regard can be converted to swath width.
        """
        self.assertAlmostEqual(
            field_of_regard_to_swath_width(705000, 15.0), 185815, delta=1.0
        )

    def test_field_of_regard_to_swath_width_modis(self):
        """
        Test against the published MODIS instrument specification (Terra/
        Aqua, 705 km altitude, +/-55 degree scan angle for a 110 degree
        field of regard, 2330 km swath width).
        """
        self.assertAlmostEqual(
            field_of_regard_to_swath_width(705000, 110.0),
            2330000,
            delta=100,
        )

    def test_field_of_regard_to_swath_width_zero(self):
        """
        Test that a zero field of regard (nadir-only observation) yields a
        zero-width swath.
        """
        self.assertEqual(field_of_regard_to_swath_width(705000, 0), 0.0)

    def test_field_of_regard_to_swath_width_elevation(self):
        """
        Test that a positive elevation (observing a point above the mean
        Earth radius) changes the swath width for the same field of regard
        and altitude.
        """
        self.assertAlmostEqual(
            field_of_regard_to_swath_width(705000, 15.0, 5000),
            184495.638,
            delta=1e-3,
        )

    def test_field_of_regard_to_swath_width_increases_with_field_of_regard(self):
        """
        Test that, for a fixed altitude, the swath width increases
        monotonically with the field of regard.
        """
        swath_widths = [
            field_of_regard_to_swath_width(705000, field_of_regard)
            for field_of_regard in (1, 5, 10, 15, 20, 30, 50)
        ]
        self.assertEqual(swath_widths, sorted(swath_widths))

    def test_field_of_regard_to_swath_width_saturates_beyond_horizon(self):
        """
        Test that a field of regard at or beyond the horizon-limited
        maximum saturates to the maximum observable swath width rather
        than raising an error or growing without bound.
        """
        max_field_of_regard = compute_field_of_regard(705000, 0)
        self.assertEqual(
            field_of_regard_to_swath_width(705000, max_field_of_regard),
            field_of_regard_to_swath_width(705000, max_field_of_regard + 50),
        )
        self.assertEqual(
            field_of_regard_to_swath_width(705000, max_field_of_regard),
            field_of_regard_to_swath_width(705000, 180),
        )

    def test_compute_field_of_regard(self):
        """
        Test that the field of regard can be computed for a given altitude
        and minimum elevation angle.
        """
        self.assertAlmostEqual(
            compute_field_of_regard(705000, 81.66446), 15.0, delta=0.001
        )

    def test_compute_min_elevation_angle(self):
        """
        Test that the minimum elevation angle can be computed for a given
        altitude and field of regard.
        """
        self.assertAlmostEqual(
            compute_min_elevation_angle(705000, 15.0), 81.66446, delta=0.001
        )

    def test_compute_min_elevation_angle_modis(self):
        """
        Test against the published MODIS instrument specification
        (Terra/Aqua, 705 km altitude, 110 degree field of regard). The
        elevation angle at the swath edge should be roughly consistent
        with MODIS's documented maximum view zenith angle of about 65
        degrees (elevation = 90 - view zenith angle = ~25 degrees).
        """
        self.assertAlmostEqual(
            compute_min_elevation_angle(705000, 110.0), 25.0, delta=1.0
        )

    def test_compute_min_elevation_angle_saturated(self):
        """
        Test that the minimum elevation angle is saturated at 0 degrees for
        a given altitude and field of regard.
        """
        self.assertEqual(compute_min_elevation_angle(30000000, 180.0), 0.0)

    def test_compute_max_access_time(self):
        """
        Test that the maximum access time can be computed for a given
        altitude and minimum elevation angle.
        """
        self.assertAlmostEqual(compute_max_access_time(705000, 81.66446), 28, delta=1)

    def test_compute_max_access_time_iss(self):
        """
        Test against the commonly cited maximum ISS visibility duration:
        at its typical ~408 km altitude, a directly overhead pass (0
        degree minimum elevation, horizon-to-horizon) lasts up to about
        10 minutes (per NASA's Spot The Station and similar references).
        """
        self.assertAlmostEqual(compute_max_access_time(408000, 0), 600, delta=30)

    def test_compute_max_transit_time(self):
        """
        Test that the maximum transit time can be computed for a given
        altitude, inclination, and along track distance.
        """
        self.assertAlmostEqual(
            compute_max_transit_time(705000, 51.6, 100000),
            15.458200,
            delta=1e-5,
        )

    def test_compute_max_transit_time_zero_along_track(self):
        """
        Test that a zero along track distance requires zero transit time.
        """
        self.assertEqual(compute_max_transit_time(705000, 51.6, 0), 0.0)

    def test_compute_max_transit_time_linear_in_along_track(self):
        """
        Test that the maximum transit time scales linearly with the along
        track distance, for a fixed altitude and inclination.
        """
        self.assertAlmostEqual(
            compute_max_transit_time(705000, 51.6, 200000),
            2 * compute_max_transit_time(705000, 51.6, 100000),
            delta=1e-9,
        )

    def test_compute_max_transit_time_equatorial_slower_than_polar(self):
        """
        Test that, for the same along track distance and altitude, an
        equatorial orbit requires more transit time than a polar orbit.
        Each orbit's slowest (worst-case) ground velocity is evaluated at
        its own extreme latitude: for the equatorial orbit this is the
        equator itself, where Earth's eastward rotation partially cancels
        the ground-relative velocity (a subtraction); for the polar orbit
        this is the pole, where rotation contributes nothing and the
        ground-relative velocity equals the full inertial ground speed.
        The equatorial case remains slower overall.
        """
        self.assertGreater(
            compute_max_transit_time(705000, 0, 100000),
            compute_max_transit_time(705000, 90, 100000),
        )

    def test_compute_min_along_track_distance(self):
        """
        Test that the minimum along track distance can be computed for a
        given altitude, inclination, and access time.
        """
        self.assertAlmostEqual(
            compute_min_along_track_distance(705000, 51.6, 20),
            129381.170,
            delta=1e-3,
        )

    def test_compute_min_along_track_distance_zero_access_time(self):
        """
        Test that a zero access time yields zero along track distance.
        """
        self.assertEqual(compute_min_along_track_distance(705000, 51.6, 0), 0.0)

    def test_compute_min_along_track_distance_linear_in_access_time(self):
        """
        Test that the minimum along track distance scales linearly with
        the access time, for a fixed altitude and inclination.
        """
        self.assertAlmostEqual(
            compute_min_along_track_distance(705000, 51.6, 40),
            2 * compute_min_along_track_distance(705000, 51.6, 20),
            delta=1e-9,
        )

    def test_compute_min_along_track_distance_inverts_compute_max_transit_time(self):
        """
        Test that compute_min_along_track_distance is the exact inverse of
        compute_max_transit_time across a range of altitudes, inclinations,
        and along track distances.
        """
        for altitude in (500000, 705000, 800000):
            for inclination in (0, 30, 51.6, 90, 98):
                for along_track in (1000, 50000, 200000):
                    transit_time = compute_max_transit_time(
                        altitude, inclination, along_track
                    )
                    self.assertAlmostEqual(
                        compute_min_along_track_distance(
                            altitude, inclination, transit_time
                        ),
                        along_track,
                        delta=1e-6,
                    )
