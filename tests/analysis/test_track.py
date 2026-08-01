"""
Unit tests for the track analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

import numpy as np
from pyproj import Transformer
from shapely.geometry import MultiPolygon, Point as ShapelyPoint, Polygon
from skyfield.api import wgs84

from tatc.analysis import (
    OrbitCoordinate,
    OrbitOutput,
    collect_ground_pixels,
    collect_ground_track,
    collect_orbit_track,
    compute_ground_track,
)
from tatc.schemas import GroundStation, Instrument, Point, PointedInstrument, Satellite
from tatc.utils.geometry import geodesic_distance
from tatc.utils.observation import field_of_regard_to_swath_width

from .common import IssConstellationTestCase


class TestGroundTrackAnalysis(IssConstellationTestCase):
    def setUp(self):
        super().setUp()
        self.point = Point(id=0, latitude=0, longitude=0, min_elevation_angle=10)
        self.station = GroundStation(
            name="Station 1", latitude=0, longitude=180, min_elevation_angle=10
        )

    def test_collect_orbit_track(self):
        """
        Test that orbit track collection works for a single satellite and a list of times.
        """
        collect_orbit_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
        )

    def test_collect_orbit_track_empty(self):
        """
        Test that orbit track collection returns an empty DataFrame when no times are provided.
        """
        collect_orbit_track(
            self.satellite,
            [],
        )

    def test_collect_orbit_track_with_mask(self):
        """
        Test that orbit track collection works for a single satellite, a list of times, and a mask.
        """
        mask = Polygon([[-90, 45], [-90, 45], [90, 45], [90, -45], [-90, -45]])
        collect_orbit_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
                for i in range(10)
            ],
            mask=mask,
        )

    def test_collect_orbit_track_wgs84_altitude_matches_published_iss_range(self):
        """
        Test that the WGS84 orbit track altitude falls within the ISS's
        publicly reported operating altitude range (~330-460 km), confirming
        the reported height is the satellite's own altitude rather than a
        zero-elevation ground point or the `elevation` parameter.
        """
        results = collect_orbit_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
        )
        for point in results.geometry:
            self.assertGreater(point.z, 330e3)
            self.assertLess(point.z, 460e3)

    def test_collect_orbit_track_elevation_only_affects_swath_width(self):
        """
        Test that the `elevation` parameter changes the computed swath width
        but does not change the reported WGS84 position, which always
        reflects the satellite's true altitude.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        results_low = collect_orbit_track(self.satellite, times, elevation=0)
        results_high = collect_orbit_track(self.satellite, times, elevation=100e3)
        for p_low, p_high in zip(results_low.geometry, results_high.geometry):
            self.assertAlmostEqual(p_low.x, p_high.x)
            self.assertAlmostEqual(p_low.y, p_high.y)
            self.assertAlmostEqual(p_low.z, p_high.z)
        self.assertTrue(
            (results_low.swath_width.values != results_high.swath_width.values).all()
        )

    def test_collect_orbit_track_ecef_roundtrips_to_wgs84(self):
        """
        Test that the ECEF coordinate output, converted back to geodetic
        longitude/latitude/height using an independent pyproj transform,
        matches the directly-computed WGS84 output for the same times.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        wgs84_results = collect_orbit_track(
            self.satellite, times, coordinates=OrbitCoordinate.WGS84
        )
        ecef_results = collect_orbit_track(
            self.satellite, times, coordinates=OrbitCoordinate.ECEF
        )
        to_wgs84 = Transformer.from_crs("EPSG:4978", "EPSG:4326", always_xy=True)
        for wgs84_point, ecef_point in zip(wgs84_results.geometry, ecef_results.geometry):
            lon, lat, height = to_wgs84.transform(
                ecef_point.x, ecef_point.y, ecef_point.z
            )
            self.assertAlmostEqual(lon, wgs84_point.x, places=6)
            self.assertAlmostEqual(lat, wgs84_point.y, places=6)
            self.assertAlmostEqual(height, wgs84_point.z, places=2)

    def test_collect_orbit_track_crs_by_coordinates(self):
        """
        Test that the returned GeoDataFrame's CRS matches the requested
        coordinate system: geographic degrees for WGS84, geocentric meters
        for ECEF, and unset for ECI (an inertial, time-varying frame with
        no fixed EPSG code).
        """
        times = [datetime(2022, 6, 1, tzinfo=timezone.utc)]
        self.assertEqual(
            collect_orbit_track(
                self.satellite, times, coordinates=OrbitCoordinate.WGS84
            ).crs.to_string(),
            "EPSG:4326",
        )
        self.assertEqual(
            collect_orbit_track(
                self.satellite, times, coordinates=OrbitCoordinate.ECEF
            ).crs.to_string(),
            "EPSG:4978",
        )
        self.assertIsNone(
            collect_orbit_track(
                self.satellite, times, coordinates=OrbitCoordinate.ECI
            ).crs
        )

    def test_collect_orbit_track_mask_filters_consistently_across_coordinates(self):
        """
        Test that a mask (always interpreted in WGS84 lon/lat) filters the
        same set of times regardless of the requested output `coordinates`.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
            for i in range(20)
        ]
        mask = Polygon([[-90, 45], [90, 45], [90, -45], [-90, -45]])
        wgs84_times = list(
            collect_orbit_track(
                self.satellite, times, mask=mask, coordinates=OrbitCoordinate.WGS84
            ).time
        )
        ecef_times = list(
            collect_orbit_track(
                self.satellite, times, mask=mask, coordinates=OrbitCoordinate.ECEF
            ).time
        )
        eci_times = list(
            collect_orbit_track(
                self.satellite, times, mask=mask, coordinates=OrbitCoordinate.ECI
            ).time
        )
        self.assertTrue(0 < len(wgs84_times) < len(times))
        self.assertEqual(wgs84_times, ecef_times)
        self.assertEqual(wgs84_times, eci_times)

    def test_collect_orbit_track_velocity_eci_matches_published_iss_speed(self):
        """
        Test that the ECI velocity magnitude matches the ISS's publicly
        reported orbital speed of approximately 7.66 km/s.
        """
        results = collect_orbit_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
            coordinates=OrbitCoordinate.ECI,
            orbit_output=OrbitOutput.POSITION_VELOCITY,
        )
        for velocity in results.velocity:
            speed = np.linalg.norm([velocity.x, velocity.y, velocity.z])
            self.assertGreater(speed, 7.5e3)
            self.assertLess(speed, 7.8e3)

    def test_collect_orbit_track_velocity_eci_faster_than_ecef(self):
        """
        Test that the inertial (ECI) speed exceeds the Earth-fixed (ECEF)
        speed, as expected for this prograde (51.6 deg inclination) orbit
        where Earth's rotation partially cancels the satellite's ground-
        relative motion.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        eci_results = collect_orbit_track(
            self.satellite,
            times,
            coordinates=OrbitCoordinate.ECI,
            orbit_output=OrbitOutput.POSITION_VELOCITY,
        )
        ecef_results = collect_orbit_track(
            self.satellite,
            times,
            coordinates=OrbitCoordinate.ECEF,
            orbit_output=OrbitOutput.POSITION_VELOCITY,
        )
        for eci_velocity, ecef_velocity in zip(eci_results.velocity, ecef_results.velocity):
            eci_speed = np.linalg.norm([eci_velocity.x, eci_velocity.y, eci_velocity.z])
            ecef_speed = np.linalg.norm(
                [ecef_velocity.x, ecef_velocity.y, ecef_velocity.z]
            )
            self.assertGreater(eci_speed, ecef_speed)

    def test_collect_orbit_track_velocity_wgs84_enu_matches_ecef_magnitude(self):
        """
        Test that the WGS84 East/North/Up velocity has the same magnitude as
        the ECEF velocity it is rotated from, since a rotation preserves
        vector length.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        wgs84_results = collect_orbit_track(
            self.satellite,
            times,
            coordinates=OrbitCoordinate.WGS84,
            orbit_output=OrbitOutput.POSITION_VELOCITY,
        )
        ecef_results = collect_orbit_track(
            self.satellite,
            times,
            coordinates=OrbitCoordinate.ECEF,
            orbit_output=OrbitOutput.POSITION_VELOCITY,
        )
        for enu_velocity, ecef_velocity in zip(
            wgs84_results.velocity, ecef_results.velocity
        ):
            enu_speed = np.linalg.norm([enu_velocity.x, enu_velocity.y, enu_velocity.z])
            ecef_speed = np.linalg.norm(
                [ecef_velocity.x, ecef_velocity.y, ecef_velocity.z]
            )
            self.assertAlmostEqual(enu_speed, ecef_speed, places=3)

    def test_collect_orbit_track_swath_width_independent_of_coordinates(self):
        """
        Test that `swath_width` is the same regardless of the requested
        output `coordinates`, since it depends only on the satellite's
        altitude and the `elevation` parameter, not the output frame.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        wgs84_widths = collect_orbit_track(
            self.satellite, times, coordinates=OrbitCoordinate.WGS84
        ).swath_width.values
        ecef_widths = collect_orbit_track(
            self.satellite, times, coordinates=OrbitCoordinate.ECEF
        ).swath_width.values
        eci_widths = collect_orbit_track(
            self.satellite, times, coordinates=OrbitCoordinate.ECI
        ).swath_width.values
        np.testing.assert_allclose(wgs84_widths, ecef_widths)
        np.testing.assert_allclose(wgs84_widths, eci_widths)

    def test_collect_ground_track(self):
        """
        Test that ground track collection works for a single satellite and a list of times.
        """
        collect_ground_track(
            self.satellite,
            times=[
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
        )

    def test_collect_ground_track_empty(self):
        """
        Test that ground track collection returns an empty DataFrame when
        no times are provided.
        """
        collect_ground_track(
            self.satellite,
            [],
        )

    def test_collect_ground_track_with_mask(self):
        """
        Test that ground track collection works for a single satellite, a list
        of times, and a mask.
        """
        mask = Polygon([[-90, 45], [-90, 45], [90, 45], [90, -45], [-90, -45]])
        collect_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
                for i in range(10)
            ],
            mask=mask,
        )

    def test_collect_ground_track_returns_one_row_per_time(self):
        """
        Test that ground track collection returns one polygon per requested
        time, each a valid, non-empty Polygon or MultiPolygon.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        results = collect_ground_track(self.satellite, times)
        self.assertEqual(len(results.index), len(times))
        for geometry in results.geometry:
            self.assertIsInstance(geometry, (Polygon, MultiPolygon))
            self.assertFalse(geometry.is_empty)

    def test_collect_ground_track_footprint_contains_subsatellite_point(self):
        """
        Test that the (nadir-pointing) instrument's footprint contains the
        satellite's true sub-satellite point.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        results = collect_ground_track(self.satellite, times)
        orbit_track = self.orbit.to_gp_orbit().get_orbit_track(times)
        for i, geometry in enumerate(results.geometry):
            subpoint = wgs84.subpoint_of(orbit_track[i])
            self.assertTrue(
                geometry.contains(
                    ShapelyPoint(subpoint.longitude.degrees, subpoint.latitude.degrees)
                )
            )

    def test_collect_ground_track_footprint_radius_matches_swath_width(self):
        """
        Test that the SPICE-computed footprint radius (geodesic distance
        from the sub-satellite point to the farthest footprint vertex)
        matches the analytic swath width formula from
        `field_of_regard_to_swath_width`, for a narrow field of regard
        (which keeps the footprint a single, near-circular polygon).
        """
        instrument = Instrument(name="Narrow", field_of_regard=10.0)
        satellite = Satellite(
            name="Narrow", orbit=self.orbit, instruments=[instrument]
        )
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(5)
        ]
        results = collect_ground_track(satellite, times)
        orbit_track = self.orbit.to_gp_orbit().get_orbit_track(times)
        gp_orbit = self.orbit.to_gp_orbit()
        expected_radius = (
            field_of_regard_to_swath_width(
                gp_orbit.get_mean_altitude(), instrument.field_of_regard
            )
            / 2
        )
        for i, geometry in enumerate(results.geometry):
            self.assertIsInstance(geometry, Polygon)
            subpoint = wgs84.subpoint_of(orbit_track[i])
            max_distance = max(
                geodesic_distance(
                    subpoint.longitude.degrees, subpoint.latitude.degrees, x, y
                )
                for x, y, *_ in geometry.exterior.coords
            )
            self.assertAlmostEqual(
                max_distance, expected_radius, delta=expected_radius * 0.05
            )

    def test_collect_ground_track_sat_altaz_is_zenith_at_footprint_center(self):
        """
        Test that the satellite appears at zenith (90 deg altitude) as seen
        from its own footprint center, since this is a nadir-pointing
        instrument.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        results = collect_ground_track(self.satellite, times, sat_altaz=True)
        for sat_alt in results.sat_alt:
            self.assertAlmostEqual(sat_alt, 90.0, places=1)

    def test_collect_ground_track_solar_altaz_valid_range(self):
        """
        Test that the solar altitude/azimuth angles fall within their valid
        ranges.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(10)
        ]
        results = collect_ground_track(self.satellite, times, solar_altaz=True)
        for solar_alt, solar_az in zip(results.solar_alt, results.solar_az):
            self.assertGreaterEqual(solar_alt, -90.0)
            self.assertLessEqual(solar_alt, 90.0)
            self.assertGreaterEqual(solar_az, 0.0)
            self.assertLess(solar_az, 360.0)

    def test_collect_ground_track_mask_limits_footprints_to_region(self):
        """
        Test that every returned footprint intersects the provided mask,
        and that the mask actually excludes some of the requested times.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
            for i in range(10)
        ]
        mask = Polygon([[-90, 45], [90, 45], [90, -45], [-90, -45]])
        results = collect_ground_track(self.satellite, times, mask=mask)
        self.assertTrue(0 < len(results.index) < len(times))
        for geometry in results.geometry:
            self.assertTrue(mask.intersects(geometry))

    def test_compute_ground_track(self):
        """
        Test that ground track computation works for a single satellite and a list
        of times.
        """
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), Polygon)

    def test_compute_ground_track_no_instr_index(self):
        """
        Test that ground track computation works for a single satellite and a list
        of times with no instrument index.
        """
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), Polygon)

    def test_compute_ground_track_multipolygon(self):
        """
        Test that ground track computation works for a single satellite and a list
        of times with a multipolygon result.
        """
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, 40, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), MultiPolygon)

    def test_collect_ground_pixels_requires_rectangular_pointed_instrument(self):
        """
        Test that ground pixel collection raises a ValueError for an
        instrument that is not a rectangular PointedInstrument.
        """
        with self.assertRaises(ValueError):
            collect_ground_pixels(
                self.satellite,
                [datetime(2022, 6, 1, tzinfo=timezone.utc)],
            )

    def test_collect_ground_pixels_returns_grid_per_time(self):
        """
        Test that ground pixel collection returns one point per pixel per
        requested time, for a rectangular pixel array.
        """
        instrument = PointedInstrument(
            name="Pixels",
            cross_track_field_of_view=10.0,
            along_track_field_of_view=10.0,
            is_rectangular=True,
            cross_track_pixels=3,
            along_track_pixels=3,
        )
        satellite = Satellite(name="Pixels", orbit=self.orbit, instruments=[instrument])
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(3)
        ]
        results = collect_ground_pixels(satellite, times)
        self.assertEqual(len(results.index), len(times) * 9)
        self.assertTrue((results.groupby("time").size() == 9).all())

    def test_collect_ground_pixels_center_pixel_matches_subpoint(self):
        """
        Test that, for a nadir-pointing instrument with an odd-sized pixel
        grid, one pixel coincides with the satellite's true sub-satellite
        point (the exact center of the grid has zero cone-angle offset from
        boresight).
        """
        instrument = PointedInstrument(
            name="Pixels",
            cross_track_field_of_view=10.0,
            along_track_field_of_view=10.0,
            is_rectangular=True,
            cross_track_pixels=3,
            along_track_pixels=3,
        )
        satellite = Satellite(name="Pixels", orbit=self.orbit, instruments=[instrument])
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(3)
        ]
        results = collect_ground_pixels(satellite, times)
        orbit_track = self.orbit.to_gp_orbit().get_orbit_track(times)
        for i, time in enumerate(results.time.unique()):
            subpoint = wgs84.subpoint_of(orbit_track[i])
            rows = results[results.time == time]
            min_distance = min(
                geodesic_distance(
                    subpoint.longitude.degrees, subpoint.latitude.degrees, p.x, p.y
                )
                for p in rows.geometry
            )
            self.assertAlmostEqual(min_distance, 0.0, delta=1.0)

    def test_collect_ground_pixels_sat_altaz_center_pixel_is_zenith(self):
        """
        Test that the center pixel's satellite altitude is at zenith
        (90 deg), since it coincides with the nadir sub-satellite point.
        """
        instrument = PointedInstrument(
            name="Pixels",
            cross_track_field_of_view=10.0,
            along_track_field_of_view=10.0,
            is_rectangular=True,
            cross_track_pixels=3,
            along_track_pixels=3,
        )
        satellite = Satellite(name="Pixels", orbit=self.orbit, instruments=[instrument])
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(3)
        ]
        results = collect_ground_pixels(satellite, times, sat_altaz=True)
        for _, rows in results.groupby("time"):
            self.assertAlmostEqual(rows.sat_alt.max(), 90.0, places=1)

    def test_collect_ground_pixels_solar_altaz_valid_range(self):
        """
        Test that the solar altitude/azimuth angles fall within their valid
        ranges.
        """
        instrument = PointedInstrument(
            name="Pixels",
            cross_track_field_of_view=10.0,
            along_track_field_of_view=10.0,
            is_rectangular=True,
            cross_track_pixels=2,
            along_track_pixels=2,
        )
        satellite = Satellite(name="Pixels", orbit=self.orbit, instruments=[instrument])
        times = [datetime(2022, 6, 1, tzinfo=timezone.utc)]
        results = collect_ground_pixels(satellite, times, solar_altaz=True)
        for solar_alt, solar_az in zip(results.solar_alt, results.solar_az):
            self.assertGreaterEqual(solar_alt, -90.0)
            self.assertLessEqual(solar_alt, 90.0)
            self.assertGreaterEqual(solar_az, 0.0)
            self.assertLess(solar_az, 360.0)

    def test_collect_ground_pixels_mask_limits_pixels_to_region(self):
        """
        Test that every returned pixel falls within the provided mask, and
        that the mask actually excludes some of the requested times.
        """
        instrument = PointedInstrument(
            name="Pixels",
            cross_track_field_of_view=10.0,
            along_track_field_of_view=10.0,
            is_rectangular=True,
            cross_track_pixels=3,
            along_track_pixels=3,
        )
        satellite = Satellite(name="Pixels", orbit=self.orbit, instruments=[instrument])
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
            for i in range(10)
        ]
        mask = Polygon([[-90, 45], [90, 45], [90, -45], [-90, -45]])
        results = collect_ground_pixels(satellite, times, mask=mask)
        self.assertTrue(0 < results.time.nunique() < len(times))
        for geometry in results.geometry:
            self.assertTrue(mask.contains(geometry))

    def test_collect_ground_pixels_empty(self):
        """
        Test that ground pixel collection returns an empty DataFrame when
        no times are provided.
        """
        instrument = PointedInstrument(
            name="Pixels",
            cross_track_field_of_view=10.0,
            along_track_field_of_view=10.0,
            is_rectangular=True,
        )
        satellite = Satellite(name="Pixels", orbit=self.orbit, instruments=[instrument])
        results = collect_ground_pixels(satellite, [])
        self.assertTrue(results.empty)
