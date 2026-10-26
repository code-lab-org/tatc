"""
Unit tests for the tatc.utils.observation module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from pyproj import Geod

from tatc.analysis import collect_orbit_track
from tatc.constants import EARTH_MEAN_RADIUS
from tatc.schemas import Satellite, SunSynchronousOrbit
from tatc.utils import (
    compute_along_track_field_of_view,
    compute_ground_inertial_velocity,
    compute_ground_surface_velocity,
    compute_field_of_regard,
    compute_local_incidence_angle,
    compute_max_access_time,
    compute_max_transit_time,
    compute_min_along_track_distance,
    compute_min_elevation_angle,
    field_of_regard_to_swath_width,
    geodetic_to_rectangular,
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

    def test_compute_along_track_field_of_view(self):
        """
        Test that the along-track field of view subtends the distance the
        ground track advances in one time step over the ellipsoid at the
        given latitude, from the height above the ellipsoid, and by default
        at the fastest ground track velocity over the orbit, which exceeds
        the inertial ground velocity for a retrograde orbit.
        """
        altitude, time_step = 700e3, 10
        for inclination in (51.6, 98.2):
            for latitude in (0, 30, -45):
                velocity = compute_ground_surface_velocity(
                    altitude, 0, inclination, latitude
                )
                surface_radius = np.linalg.norm(geodetic_to_rectangular(0, latitude))
                advance = velocity * time_step * surface_radius / EARTH_MEAN_RADIUS
                height = EARTH_MEAN_RADIUS + altitude - surface_radius
                self.assertAlmostEqual(
                    compute_along_track_field_of_view(
                        altitude, time_step, inclination, latitude
                    ),
                    np.degrees(2 * np.arctan(advance / 2 / height)),
                    delta=1e-9,
                )
                self.assertGreaterEqual(
                    compute_along_track_field_of_view(altitude, time_step, inclination),
                    compute_along_track_field_of_view(
                        altitude, time_step, inclination, latitude
                    ),
                )
        # retrograde: faster than the ground track of a non-rotating Earth
        inertial = compute_ground_inertial_velocity(altitude)
        self.assertGreater(
            compute_along_track_field_of_view(altitude, time_step, 98.2),
            np.degrees(2 * np.arctan(inertial * time_step / 2 / altitude)),
        )
        # an eccentric orbit is evaluated at perigee (lower and faster)
        self.assertGreater(
            compute_along_track_field_of_view(altitude, time_step, 98.2, None, 0.01),
            compute_along_track_field_of_view(altitude, time_step, 98.2),
        )
        # proportional to the time step for small angles
        self.assertAlmostEqual(
            compute_along_track_field_of_view(altitude, 2, 98.2)
            / compute_along_track_field_of_view(altitude, 1, 98.2),
            2,
            delta=1e-4,
        )

    def test_along_track_field_of_view_tiles_propagated_orbit(self):
        """
        Test that nadir footprints with the tailored along-track field of
        view, at the satellite's height above the ellipsoid, are at least as
        long as the ground track advances in a time step around a
        propagated orbit, including near the equator, where the satellite is
        lowest above the ellipsoid (about 7 km below its mean altitude).
        """
        start = datetime(2022, 7, 14, tzinfo=timezone.utc)
        orbit = SunSynchronousOrbit(
            altitude=834e3, equator_crossing_time="13:30", epoch=start
        )
        satellite = Satellite(name="Test", orbit=orbit)
        time_step = 1.7864
        field_of_view = compute_along_track_field_of_view(
            orbit.get_mean_altitude(), time_step, orbit.get_inclination()
        )
        times = [start + timedelta(seconds=10 * k) for k in range(610)]
        track = collect_orbit_track(satellite, times)
        later = collect_orbit_track(
            satellite, [t + timedelta(seconds=time_step) for t in times]
        )
        advance = Geod(ellps="WGS84").inv(
            track.geometry.x, track.geometry.y, later.geometry.x, later.geometry.y
        )[2]
        height = np.array([point.z for point in track.geometry])
        length = 2 * height * np.tan(np.radians(field_of_view / 2))
        self.assertTrue(np.all(length >= advance))


class TestLocalIncidenceAngle(unittest.TestCase):
    """
    Unit tests for compute_local_incidence_angle.
    """

    def test_flat_terrain(self):
        """
        Test that the local incidence angle on flat terrain is 90 degrees
        less the elevation angle, at any azimuth and aspect.
        """
        for elevation_angle in [10, 45, 90]:
            for azimuth in [0, 123, 270]:
                self.assertAlmostEqual(
                    compute_local_incidence_angle(elevation_angle, azimuth, 0, 77),
                    90 - elevation_angle,
                )

    def test_slope_facing_toward_and_away(self):
        """
        Test that a slope facing the satellite reduces the incidence angle
        by the slope, and one facing away increases it by the slope.
        """
        self.assertAlmostEqual(compute_local_incidence_angle(40, 270, 20, 270), 30)
        self.assertAlmostEqual(compute_local_incidence_angle(40, 270, 20, 90), 70)

    def test_slope_facing_across(self):
        """
        Test the angle for a slope facing perpendicular to the line of
        sight, where cos(angle) = sin(elevation) cos(slope).
        """
        self.assertAlmostEqual(
            compute_local_incidence_angle(40, 0, 20, 90),
            np.degrees(np.arccos(np.sin(np.radians(40)) * np.cos(np.radians(20)))),
        )

    def test_hidden_slope(self):
        """
        Test that a slope facing away more steeply than the satellite's
        elevation angle gives an angle above 90 degrees (hidden).
        """
        self.assertGreater(compute_local_incidence_angle(30, 0, 45, 180), 90)
        self.assertAlmostEqual(compute_local_incidence_angle(30, 0, 45, 180), 105)

    def test_arrays_broadcast(self):
        """
        Test that arrays broadcast elementwise and return an array, while
        scalars return a float.
        """
        angles = compute_local_incidence_angle(
            np.array([40, 40, 90]), 270, np.array([20, 0, 30]), 270
        )
        self.assertIsInstance(angles, np.ndarray)
        np.testing.assert_allclose(angles, [30, 50, 30])
        self.assertIsInstance(compute_local_incidence_angle(40, 270), float)
