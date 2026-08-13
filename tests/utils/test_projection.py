"""
Unit tests for the tatc.utils.projection module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import numpy as np
from pyproj import Transformer
from shapely.geometry import Point
from skyfield.api import EarthSatellite

from tatc import config, constants
from tatc.constants import timescale
from tatc.schemas import CircularOrbit
from tatc.utils.observation import (
    compute_field_of_regard,
    field_of_regard_to_swath_width,
)
from tatc.utils.orbital import compute_ground_surface_velocity
from tatc.utils.projection import (
    buffer_footprint,
    buffer_target,
    compute_footprint,
    compute_limb,
    compute_projected_ray_position,
)


def _great_circle_distance(lat1, lon1, lat2, lon2):
    earth_radius = 6371000
    phi1, phi2 = np.radians(lat1), np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2) ** 2
    return 2 * earth_radius * np.arcsin(np.sqrt(a))


class TestProjection(unittest.TestCase):  # pylint: disable=too-many-public-methods
    """
    Unit tests for the tatc.utils.projection module.
    """

    def setUp(self):
        noon_utc = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        self.orbit = CircularOrbit(
            mean_altitude=705000,
            true_anomaly=0,
            epoch=noon_utc,
            inclination=51.6,
            right_ascension_ascending_node=0.0,
        )
        self.satellite = EarthSatellite.from_satrec(
            self.orbit.to_gp_orbit().elements[0].to_satrec(), timescale
        )
        self.orbit_track = self.satellite.at(timescale.from_datetime(noon_utc))
        self.subpoint = self.orbit_track.subpoint()

    def test_compute_projected_ray_position_nadir_matches_subpoint(self):
        """
        Test that a pure nadir ray (zero field of view, roll, and pitch)
        lands at essentially the same location as Skyfield's own
        independently-computed sub-satellite point. The nadir ray must
        follow the geodetic vertical (the WGS 84 ellipsoid surface normal),
        not the geocentric direction (straight toward the Earth's center):
        those two only coincide at the equator and poles, and this test
        point (true_anomaly=0, near the ascending node) sits close to the
        equator, where the two would be hard to tell apart. See
        `test_compute_projected_ray_position_nadir_matches_subpoint_off_equator`
        for a test point where the distinction actually matters.
        """
        position = compute_projected_ray_position(
            self.orbit_track, 0, 0, 0, 0, False, 0, 0
        )
        distance = _great_circle_distance(
            self.subpoint.latitude.degrees,
            self.subpoint.longitude.degrees,
            position.latitude.degrees,
            position.longitude.degrees,
        )
        self.assertLess(distance, 1)

    def test_compute_projected_ray_position_nadir_matches_subpoint_off_equator(self):
        """
        Test that a pure nadir ray still matches Skyfield's sub-satellite
        point away from the equator (near maximum latitude for this
        orbit's inclination), where the geocentric and geodetic nadir
        directions diverge by kilometers if conflated. This is the
        regression test for a bug where the nadir ray pointed toward the
        Earth's center (geocentric) instead of along the local WGS 84
        ellipsoid normal (geodetic), which agreed with Skyfield's subpoint
        only near the equator/poles and was off by ~2 km at 45 degrees
        latitude for a 700 km altitude orbit.
        """
        noon_utc = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        orbit = CircularOrbit(
            mean_altitude=705000,
            true_anomaly=90,
            epoch=noon_utc,
            inclination=51.6,
            right_ascension_ascending_node=0.0,
        )
        satellite = EarthSatellite.from_satrec(
            orbit.to_gp_orbit().elements[0].to_satrec(), timescale
        )
        orbit_track = satellite.at(timescale.from_datetime(noon_utc))
        subpoint = orbit_track.subpoint()
        position = compute_projected_ray_position(orbit_track, 0, 0, 0, 0, False, 0, 0)
        # confirm this test point is actually away from the equator
        self.assertGreater(abs(subpoint.latitude.degrees), 45)
        distance = _great_circle_distance(
            subpoint.latitude.degrees,
            subpoint.longitude.degrees,
            position.latitude.degrees,
            position.longitude.degrees,
        )
        self.assertLess(distance, 1)

    def test_compute_projected_ray_position_elevation(self):
        """
        Test that the projected position's elevation matches the
        requested elevation.
        """
        position = compute_projected_ray_position(
            self.orbit_track, 0, 0, 0, 0, False, 0, 1000
        )
        self.assertAlmostEqual(position.elevation.m, 1000, delta=1e-3)

    def test_compute_projected_ray_position_roll_moves_away_from_subpoint(self):
        """
        Test that increasing the magnitude of the roll angle moves the
        projected position monotonically farther from the sub-satellite
        point, in both directions.
        """
        distances = [
            _great_circle_distance(
                self.subpoint.latitude.degrees,
                self.subpoint.longitude.degrees,
                *self._latlon(roll_angle=roll),
            )
            for roll in (0, 5, 10, 20)
        ]
        self.assertEqual(distances, sorted(distances))
        distances_negative = [
            _great_circle_distance(
                self.subpoint.latitude.degrees,
                self.subpoint.longitude.degrees,
                *self._latlon(roll_angle=roll),
            )
            for roll in (0, -5, -10, -20)
        ]
        self.assertEqual(distances_negative, sorted(distances_negative))

    def test_compute_projected_ray_position_pitch_moves_away_from_subpoint(self):
        """
        Test that increasing the magnitude of the pitch angle moves the
        projected position monotonically farther from the sub-satellite
        point, in both directions.
        """
        distances = [
            _great_circle_distance(
                self.subpoint.latitude.degrees,
                self.subpoint.longitude.degrees,
                *self._latlon(pitch_angle=pitch),
            )
            for pitch in (0, 5, 10, 20)
        ]
        self.assertEqual(distances, sorted(distances))

    def test_compute_projected_ray_position_rectangular_matches_elliptical_at_principal_angles(
        self,
    ):
        """
        Test that rectangular and elliptical field of view shapes produce
        identical rays at the four axis-aligned angles (0, 90, 180, 270
        degrees), where an ellipse touches its bounding rectangle.
        """
        for angle in (0, 90, 180, 270):
            rectangular = compute_projected_ray_position(
                self.orbit_track, 10, 4, 0, 0, True, angle, 0
            )
            elliptical = compute_projected_ray_position(
                self.orbit_track, 10, 4, 0, 0, False, angle, 0
            )
            self.assertAlmostEqual(
                rectangular.latitude.degrees, elliptical.latitude.degrees, delta=1e-9
            )
            self.assertAlmostEqual(
                rectangular.longitude.degrees, elliptical.longitude.degrees, delta=1e-9
            )

    def test_compute_projected_ray_position_rectangular_differs_off_axis(self):
        """
        Test that rectangular and elliptical field of view shapes produce
        different rays at an off-axis angle, confirming the rectangular
        shape is not silently ignored.
        """
        rectangular = compute_projected_ray_position(
            self.orbit_track, 10, 4, 0, 0, True, 45, 0
        )
        elliptical = compute_projected_ray_position(
            self.orbit_track, 10, 4, 0, 0, False, 45, 0
        )
        self.assertNotAlmostEqual(
            rectangular.latitude.degrees, elliptical.latitude.degrees, delta=1e-6
        )

    def test_compute_projected_ray_position_rectangular_all_edge_segments(self):
        """
        Test that sampling one angle from each of the rectangle's 8
        boundary segments (split at the 4 corners and, within the top and
        bottom edges, at the axes) produces 8 distinct, well-defined
        positions, none of which collapse to the sub-satellite point.
        """
        theta = np.degrees(np.arctan(4 / 10))
        sample_angles = [
            theta / 2,  # right edge, upper half
            45,  # top edge, right half
            (theta + (180 - theta)) / 2,  # top edge, left half
            179,  # left edge, upper half
            180 + theta / 2,  # left edge, lower half
            225,  # bottom edge, right half
            304,  # bottom edge, left half
            350,  # right edge, lower half
        ]
        positions = [
            compute_projected_ray_position(self.orbit_track, 10, 4, 0, 0, True, a, 0)
            for a in sample_angles
        ]
        latlons = [(p.latitude.degrees, p.longitude.degrees) for p in positions]
        self.assertEqual(len(set(latlons)), len(latlons))
        for lat, lon in latlons:
            self.assertLess(
                _great_circle_distance(
                    self.subpoint.latitude.degrees,
                    self.subpoint.longitude.degrees,
                    lat,
                    lon,
                ),
                200000,
            )

    def test_compute_projected_ray_position_saturates_beyond_horizon(self):
        """
        Test that a roll angle beyond the horizon-limited maximum (where
        the ray misses the WGS 84 geoid entirely) falls back to a stable
        limb-edge position rather than raising an error, and that this
        fallback position no longer changes with further increases in
        roll angle.
        """
        max_look_angle = compute_field_of_regard(705000, 0) / 2
        beyond_horizon = compute_projected_ray_position(
            self.orbit_track, 0, 0, max_look_angle + 5, 0, False, 0, 0
        )
        further_beyond = compute_projected_ray_position(
            self.orbit_track, 0, 0, max_look_angle + 20, 0, False, 0, 0
        )
        self.assertAlmostEqual(
            beyond_horizon.latitude.degrees, further_beyond.latitude.degrees, delta=1e-9
        )
        self.assertAlmostEqual(
            beyond_horizon.longitude.degrees,
            further_beyond.longitude.degrees,
            delta=1e-9,
        )

    def test_compute_projected_ray_position_vectorized_matches_scalar(self):
        """
        Test that passing a vector of times produces the same results as
        calling the function separately for each individual time.
        """
        times = timescale.utc(2020, 3, 20, 12, [0, 1, 2])
        orbit_track = self.satellite.at(times)
        vectorized = compute_projected_ray_position(
            orbit_track, 0, 0, 5, 0, False, 0, 0
        )
        for i in range(3):
            scalar = compute_projected_ray_position(
                self.satellite.at(times[i]), 0, 0, 5, 0, False, 0, 0
            )
            self.assertAlmostEqual(
                vectorized.latitude.degrees[i], scalar.latitude.degrees, delta=1e-9
            )
            self.assertAlmostEqual(
                vectorized.longitude.degrees[i], scalar.longitude.degrees, delta=1e-9
            )

    def test_compute_footprint_scalar_orbit_track(self):
        """
        Test that a scalar (single-time) orbit_track produces a single
        valid footprint, rather than raising an error.
        """
        footprints = compute_footprint(self.orbit_track, 10, 10, 0, 0, False)
        self.assertEqual(len(footprints), 1)
        self.assertTrue(footprints[0].is_valid)

    def test_compute_footprint_scalar_matches_single_element_vector(self):
        """
        Test that a scalar orbit_track produces the same footprint as a
        vectorized orbit_track containing that single time.
        """
        vector_track = self.satellite.at(timescale.utc(2020, 3, 20, 12, [0]))
        scalar_footprint = compute_footprint(self.orbit_track, 10, 10, 0, 0, False)
        vector_footprint = compute_footprint(vector_track, 10, 10, 0, 0, False)
        self.assertTrue(scalar_footprint[0].equals(vector_footprint[0]))

    def test_compute_footprint_vectorized_multiple_times(self):
        """
        Test that a vectorized orbit_track produces one valid footprint per time.
        """
        times = timescale.utc(2020, 3, 20, 12, [0, 1, 2])
        orbit_track = self.satellite.at(times)
        footprints = compute_footprint(orbit_track, 10, 10, 0, 0, False)
        self.assertEqual(len(footprints), 3)
        for footprint in footprints:
            self.assertTrue(footprint.is_valid)

    def test_compute_footprint_default_number_points_elliptical(self):
        """
        Test that the default number of points for an elliptical footprint
        matches the runtime configuration.
        """
        footprints = compute_footprint(self.orbit_track, 10, 10, 0, 0, False)
        self.assertEqual(
            len(footprints[0].exterior.coords) - 1,
            config.rc.footprint_points_elliptical,
        )

    def test_compute_footprint_default_number_points_rectangular(self):
        """
        Test that the default number of points for a rectangular footprint
        is 4 times the per-side runtime configuration (one segment per side).
        """
        footprints = compute_footprint(self.orbit_track, 10, 4, 0, 0, True)
        self.assertEqual(
            len(footprints[0].exterior.coords) - 1,
            4 * config.rc.footprint_points_rectangular_side,
        )

    def test_compute_footprint_explicit_number_points(self):
        """
        Test that an explicit number_points overrides the runtime
        configuration default.
        """
        footprints = compute_footprint(
            self.orbit_track, 10, 10, 0, 0, False, number_points=10
        )
        self.assertEqual(len(footprints[0].exterior.coords) - 1, 10)

    def test_compute_footprint_contains_subpoint(self):
        """
        Test that both elliptical and rectangular nadir footprints contain
        the sub-satellite point.
        """
        subpoint_geom = Point(
            self.subpoint.longitude.degrees, self.subpoint.latitude.degrees
        )
        elliptical = compute_footprint(self.orbit_track, 10, 10, 0, 0, False)
        rectangular = compute_footprint(self.orbit_track, 10, 4, 0, 0, True)
        self.assertTrue(elliptical[0].contains(subpoint_geom))
        self.assertTrue(rectangular[0].contains(subpoint_geom))

    def test_compute_limb_scalar_orbit_track(self):
        """
        Test that a scalar (single-time) orbit_track produces a single
        valid limb, rather than raising an error.
        """
        limbs = compute_limb(self.orbit_track)
        self.assertEqual(len(limbs), 1)
        self.assertTrue(limbs[0].is_valid)

    def test_compute_limb_scalar_matches_single_element_vector(self):
        """
        Test that a scalar orbit_track produces the same limb as a
        vectorized orbit_track containing that single time.
        """
        vector_track = self.satellite.at(timescale.utc(2020, 3, 20, 12, [0]))
        scalar_limb = compute_limb(self.orbit_track)
        vector_limb = compute_limb(vector_track)
        self.assertTrue(scalar_limb[0].equals(vector_limb[0]))

    def test_compute_limb_vectorized_multiple_times(self):
        """
        Test that a vectorized orbit_track produces one valid limb per time.
        """
        times = timescale.utc(2020, 3, 20, 12, [0, 1, 2])
        orbit_track = self.satellite.at(times)
        limbs = compute_limb(orbit_track)
        self.assertEqual(len(limbs), 3)
        for limb in limbs:
            self.assertTrue(limb.is_valid)

    def test_compute_limb_matches_expected_horizon_angle(self):
        """
        Test that every point on the limb boundary is at approximately the
        expected Earth central angle (subpoint to horizon, for a zero
        minimum elevation angle) from the sub-satellite point.
        """
        altitude = self.orbit.mean_altitude
        expected_angle = np.degrees(
            np.arccos(
                constants.EARTH_MEAN_RADIUS / (constants.EARTH_MEAN_RADIUS + altitude)
            )
        )
        limb = compute_limb(self.orbit_track)[0]
        for x, y, _ in limb.exterior.coords:
            distance = np.degrees(
                _great_circle_distance(
                    self.subpoint.latitude.degrees,
                    self.subpoint.longitude.degrees,
                    y,
                    x,
                )
                / constants.EARTH_MEAN_RADIUS
            )
            self.assertAlmostEqual(distance, expected_angle, delta=0.5)

    def test_compute_limb_contains_subpoint(self):
        """
        Test that the limb contains the sub-satellite point (the closest,
        and definitely visible, point on the Earth's surface).
        """
        subpoint_geom = Point(
            self.subpoint.longitude.degrees, self.subpoint.latitude.degrees
        )
        limb = compute_limb(self.orbit_track)[0]
        self.assertTrue(limb.contains(subpoint_geom))

    def test_compute_limb_number_points(self):
        """
        Test that the number_points parameter controls the number of
        vertices in the limb polygon.
        """
        limb = compute_limb(self.orbit_track, number_points=8)[0]
        self.assertEqual(len(limb.exterior.coords) - 1, 8)

    def test_buffer_footprint_contains_center(self):
        """
        Test that buffering a point produces a polygon containing that
        original point.
        """
        point = Point(0, 0)
        to_crs = Transformer.from_crs("EPSG:4326", "EPSG:4087", always_xy=True)
        from_crs = Transformer.from_crs("EPSG:4087", "EPSG:4326", always_xy=True)
        result = buffer_footprint(point, to_crs, from_crs, 100000, 0)
        self.assertTrue(result.contains(point))

    def test_buffer_footprint_excludes_far_point(self):
        """
        Test that a point well outside the swath width is not contained
        in the buffered footprint.
        """
        point = Point(0, 0)
        to_crs = Transformer.from_crs("EPSG:4326", "EPSG:4087", always_xy=True)
        from_crs = Transformer.from_crs("EPSG:4087", "EPSG:4326", always_xy=True)
        result = buffer_footprint(point, to_crs, from_crs, 100000, 0)
        self.assertFalse(result.contains(Point(5, 5)))

    def test_buffer_footprint_matches_expected_radius(self):
        """
        Test that the buffered polygon's boundary is approximately
        swath_width / 2 (great-circle distance) from the center point.
        """
        point = Point(0, 0)
        swath_width = 200000
        to_crs = Transformer.from_crs("EPSG:4326", "EPSG:4087", always_xy=True)
        from_crs = Transformer.from_crs("EPSG:4087", "EPSG:4326", always_xy=True)
        result = buffer_footprint(point, to_crs, from_crs, swath_width, 0)
        for x, y, _ in result.exterior.coords:
            distance = _great_circle_distance(0, 0, y, x)
            self.assertAlmostEqual(distance, swath_width / 2, delta=1000)

    def test_buffer_footprint_increases_with_swath_width(self):
        """
        Test that the buffered footprint's area increases monotonically
        with the swath width.
        """
        point = Point(0, 0)
        to_crs = Transformer.from_crs("EPSG:4326", "EPSG:4087", always_xy=True)
        from_crs = Transformer.from_crs("EPSG:4087", "EPSG:4326", always_xy=True)
        smaller = buffer_footprint(point, to_crs, from_crs, 100000, 0)
        larger = buffer_footprint(point, to_crs, from_crs, 200000, 0)
        self.assertGreater(larger.area, smaller.area)

    def test_buffer_footprint_elevation(self):
        """
        Test that the buffered footprint is projected to the requested elevation.
        """
        point = Point(0, 0)
        to_crs = Transformer.from_crs("EPSG:4326", "EPSG:4087", always_xy=True)
        from_crs = Transformer.from_crs("EPSG:4087", "EPSG:4326", always_xy=True)
        result = buffer_footprint(point, to_crs, from_crs, 100000, 500)
        self.assertAlmostEqual(next(iter(result.exterior.coords))[2], 500, delta=1e-6)

    def test_buffer_footprint_zero_swath_width_is_empty(self):
        """
        Test that a zero swath width produces an empty geometry rather
        than raising an error (buffering a point by zero distance).
        """
        point = Point(0, 0)
        to_crs = Transformer.from_crs("EPSG:4326", "EPSG:4087", always_xy=True)
        from_crs = Transformer.from_crs("EPSG:4087", "EPSG:4326", always_xy=True)
        result = buffer_footprint(point, to_crs, from_crs, 0, 0)
        self.assertTrue(result.is_empty)

    def test_buffer_target_contains_original_geometry(self):
        """
        Test that the buffered target contains the original geometry.
        """
        point = Point(0, 0)
        result = buffer_target(point, 705000, 51.6, 20, 60)
        self.assertTrue(result.contains(point))

    def test_buffer_target_matches_expected_distance(self):
        """
        Test that the buffer distance matches the independently-computed
        expected distance: half the swath width (from the field of
        regard) plus the ground distance traveled in one time step (using
        the ground velocity at the orbit's extreme latitude, the fastest,
        most conservative point).
        """
        point = Point(0, 0)
        altitude, inclination, field_of_regard, time_step = 705000, 51.6, 20, 60
        swath_width = field_of_regard_to_swath_width(altitude, field_of_regard)
        extreme_latitude = min(inclination, 180 - inclination)
        ground_velocity = compute_ground_surface_velocity(
            altitude, 0, inclination, extreme_latitude
        )
        expected_distance = ground_velocity * time_step + swath_width / 2
        result = buffer_target(point, altitude, inclination, field_of_regard, time_step)
        for x, y, *_ in result.exterior.coords:
            distance = _great_circle_distance(0, 0, y, x)
            self.assertAlmostEqual(distance, expected_distance, delta=1000)

    def test_buffer_target_is_conservative_at_high_latitude(self):
        """
        Test that the default `distance_crs` (an equidistant cylindrical
        projection whose standard parallel tracks the target's own
        latitude) keeps the buffer close to the intended distance in every
        direction, even far from the equator -- unlike a fixed low-latitude
        standard parallel (e.g. `EPSG:4087`), whose east-west ground
        distance for a fixed buffer shrinks by `cos(latitude)` away from
        the equator (see
        `test_buffer_target_fixed_low_latitude_crs_under_buffers_at_high_latitude`).
        """
        altitude, inclination, field_of_regard, time_step = 705000, 51.6, 20, 60
        swath_width = field_of_regard_to_swath_width(altitude, field_of_regard)
        extreme_latitude = min(inclination, 180 - inclination)
        ground_velocity = compute_ground_surface_velocity(
            altitude, 0, inclination, extreme_latitude
        )
        expected_distance = ground_velocity * time_step + swath_width / 2
        point = Point(0, 51.6)
        result = buffer_target(point, altitude, inclination, field_of_regard, time_step)
        distances = [
            _great_circle_distance(51.6, 0, y, x) for x, y, *_ in result.exterior.coords
        ]
        # every direction (including the worst-case, most-compressed one)
        # should still reach at least the intended distance
        self.assertGreater(min(distances), expected_distance * 0.95)

    def test_buffer_target_fixed_low_latitude_crs_under_buffers_at_high_latitude(self):
        """
        Regression test documenting why `distance_crs` now defaults to a
        latitude-tracking projection instead of a fixed `EPSG:4087` (true
        scale only at the equator): forcing `EPSG:4087` at a high latitude
        under-buffers in the east-west direction by roughly `cos(latitude)`,
        which could silently exclude a target that is still within reach.
        """
        altitude, inclination, field_of_regard, time_step = 705000, 51.6, 20, 60
        swath_width = field_of_regard_to_swath_width(altitude, field_of_regard)
        extreme_latitude = min(inclination, 180 - inclination)
        ground_velocity = compute_ground_surface_velocity(
            altitude, 0, inclination, extreme_latitude
        )
        expected_distance = ground_velocity * time_step + swath_width / 2
        point = Point(0, 51.6)
        result = buffer_target(
            point,
            altitude,
            inclination,
            field_of_regard,
            time_step,
            distance_crs="EPSG:4087",
        )
        distances = [
            _great_circle_distance(51.6, 0, y, x) for x, y, *_ in result.exterior.coords
        ]
        # the east-west (worst-case) direction falls well short of the
        # intended distance -- roughly cos(51.6 deg) = 0.62x
        self.assertLess(min(distances), expected_distance * 0.7)

    def test_buffer_target_increases_with_time_step(self):
        """
        Test that the buffered target's area increases monotonically with
        the time step (more distance traveled).
        """
        point = Point(0, 0)
        smaller = buffer_target(point, 705000, 51.6, 20, 10)
        larger = buffer_target(point, 705000, 51.6, 20, 200)
        self.assertGreater(larger.area, smaller.area)

    def test_buffer_target_increases_with_field_of_regard(self):
        """
        Test that the buffered target's area increases monotonically with
        the field of regard (wider swath).
        """
        point = Point(0, 0)
        smaller = buffer_target(point, 705000, 51.6, 5, 60)
        larger = buffer_target(point, 705000, 51.6, 60, 60)
        self.assertGreater(larger.area, smaller.area)

    def test_buffer_target_distance_scaling(self):
        """
        Test that distance_scaling scales the buffer distance: doubling
        it should roughly double how far the boundary extends beyond the
        original point.
        """
        point = Point(0, 0)
        result_1x = buffer_target(point, 705000, 51.6, 20, 60, distance_scaling=1.0)
        result_2x = buffer_target(point, 705000, 51.6, 20, 60, distance_scaling=2.0)
        coords_1x = list(result_1x.exterior.coords)
        coords_2x = list(result_2x.exterior.coords)
        max_dist_1x = max(_great_circle_distance(0, 0, y, x) for x, y in coords_1x)
        max_dist_2x = max(_great_circle_distance(0, 0, y, x) for x, y in coords_2x)
        self.assertAlmostEqual(max_dist_2x / max_dist_1x, 2.0, delta=0.05)

    def _latlon(self, roll_angle=0, pitch_angle=0):
        position = compute_projected_ray_position(
            self.orbit_track, 0, 0, roll_angle, pitch_angle, False, 0, 0
        )
        return position.latitude.degrees, position.longitude.degrees
