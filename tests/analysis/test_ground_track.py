"""
Unit tests for the ground track analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import MultiPolygon, Point as ShapelyPoint, Polygon, box
from skyfield.api import wgs84

from tatc.analysis import (
    collect_ground_pixels,
    collect_ground_track,
    collect_orbit_track,
    compute_ground_track,
)
from tatc.schemas import (
    GeneralPerturbationsOrbit,
    GroundStation,
    Instrument,
    Point,
    PointedInstrument,
    Satellite,
)
from tatc.utils.geometry import geodesic_distance, split_polygon
from tatc.utils.observation import field_of_regard_to_swath_width

from .common import IssConstellationTestCase


class TestGroundTrackAnalysis(IssConstellationTestCase):
    """
    Unit tests for the ground track analysis functions.
    """

    def setUp(self):
        super().setUp()
        self.point = Point(id=0, latitude=0, longitude=0, min_elevation_angle=10)
        self.station = GroundStation(
            name="Station 1", latitude=0, longitude=180, min_elevation_angle=10
        )

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
        satellite = Satellite(name="Narrow", orbit=self.orbit, instruments=[instrument])
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

    def test_collect_ground_track_mask_as_geodataframe(self):
        """
        Test that a mask passed as a GeoDataFrame (rather than a bare
        Polygon) works the same way.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
            for i in range(10)
        ]
        polygon = Polygon([[-90, 45], [90, 45], [90, -45], [-90, -45]])
        mask = gpd.GeoDataFrame(geometry=[polygon], crs="EPSG:4326")
        results = collect_ground_track(self.satellite, times, mask=mask)
        self.assertTrue(0 < len(results.index) < len(times))
        for geometry in results.geometry:
            self.assertTrue(polygon.intersects(geometry))

    def test_collect_ground_track_mask_as_geoseries(self):
        """
        Test that a mask passed as a GeoSeries (rather than a bare
        Polygon) works the same way.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
            for i in range(10)
        ]
        polygon = Polygon([[-90, 45], [90, 45], [90, -45], [-90, -45]])
        mask = gpd.GeoSeries([polygon], crs="EPSG:4326")
        results = collect_ground_track(self.satellite, times, mask=mask)
        self.assertTrue(0 < len(results.index) < len(times))
        for geometry in results.geometry:
            self.assertTrue(polygon.intersects(geometry))

    def test_collect_ground_track_mask_excludes_everything_at_coarse_stage(self):
        """
        Test that a mask far enough away to exclude every point even at
        the (conservative) observable-period culling stage returns an empty
        result.
        """
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(seconds=i)
            for i in range(10)
        ]
        mask = Polygon([[10, 89], [11, 89], [11, 89.5], [10, 89.5]])
        results = collect_ground_track(self.satellite, times, mask=mask)
        self.assertTrue(results.empty)

    def test_collect_ground_track_mask_excludes_everything_at_footprint_stage(self):
        """
        Test that a mask close enough to pass the (conservative)
        observable-period culling stage, but too far for any actual instrument footprint to
        intersect, still returns an empty result (rather than incorrectly
        returning the period-culled, unfiltered points).
        """
        instrument = Instrument(name="Narrow", field_of_regard=10.0)
        satellite = Satellite(name="Narrow", orbit=self.orbit, instruments=[instrument])
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(seconds=i)
            for i in range(10)
        ]
        subpoint = collect_orbit_track(satellite, [times[0]]).geometry.iloc[0]
        # 0.5 deg (~55 km) away: within the observable-period stage's
        # conservative field of regard (including its margins), but
        # well beyond the true footprint
        mask = ShapelyPoint(subpoint.x + 0.5, subpoint.y).buffer(0.01)
        results = collect_ground_track(satellite, times, mask=mask)
        self.assertTrue(results.empty)

    def test_collect_ground_track_mask_matches_unmasked_intersections(self):
        """
        Test that a mask culls exactly the times whose footprint does not
        intersect it, for masks across the anti-meridian (given past 180
        degrees longitude) and across a pole (spanning all longitudes),
        which a culling stage must not wrongly exclude.
        """
        instrument = Instrument(name="Narrow", field_of_regard=60.0)
        satellite = Satellite(name="Narrow", orbit=self.orbit, instruments=[instrument])
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
            for i in range(1440)
        ]
        unmasked = collect_ground_track(satellite, times)
        for mask in [
            Polygon([(170, -10), (190, -10), (190, 10), (170, 10)]),
            box(-180, -90, 180, -50),
        ]:
            with self.subTest(mask=mask.wkt):
                masked = collect_ground_track(satellite, times, mask=mask)
                expected = unmasked[
                    [split_polygon(mask).intersects(g) for g in unmasked.geometry]
                ]
                self.assertGreater(len(expected.index), 0)
                # rows are sorted by time
                self.assertEqual(list(masked.time), list(expected.time))

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

    def test_collect_ground_pixels_mask_sorts_rows_by_time(self):
        """
        Test that masked ground pixels are sorted by time, with the pixels at
        each time in the same order as without the mask.
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
        masked = collect_ground_pixels(satellite, times, mask=mask)
        unmasked = collect_ground_pixels(satellite, times)
        expected = unmasked[[mask.intersects(g) for g in unmasked.geometry]]
        self.assertTrue(masked.time.is_monotonic_increasing)
        self.assertEqual(list(masked.time), list(expected.time))
        # the same pixels, to within rounding (the masked orbit track is
        # indexed from a longer one)
        np.testing.assert_allclose(
            [(g.x, g.y) for g in masked.geometry],
            [(g.x, g.y) for g in expected.geometry],
            rtol=0,
            atol=1e-9,
        )

    def test_collect_ground_pixels_mask_as_geodataframe(self):
        """
        Test that a mask passed as a GeoDataFrame (rather than a bare
        Polygon) works the same way.
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
        polygon = Polygon([[-90, 45], [90, 45], [90, -45], [-90, -45]])
        mask = gpd.GeoDataFrame(geometry=[polygon], crs="EPSG:4326")
        results = collect_ground_pixels(satellite, times, mask=mask)
        self.assertTrue(0 < results.time.nunique() < len(times))
        for geometry in results.geometry:
            self.assertTrue(polygon.contains(geometry))

    def test_collect_ground_pixels_mask_as_geoseries(self):
        """
        Test that a mask passed as a GeoSeries (rather than a bare
        Polygon) works the same way.
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
        polygon = Polygon([[-90, 45], [90, 45], [90, -45], [-90, -45]])
        mask = gpd.GeoSeries([polygon], crs="EPSG:4326")
        results = collect_ground_pixels(satellite, times, mask=mask)
        self.assertTrue(0 < results.time.nunique() < len(times))
        for geometry in results.geometry:
            self.assertTrue(polygon.contains(geometry))

    def test_collect_ground_pixels_mask_excludes_everything_at_coarse_stage(self):
        """
        Test that a mask far enough away to exclude every point even at
        the (conservative) observable-period culling stage returns an empty
        result.
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
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(seconds=i)
            for i in range(10)
        ]
        mask = Polygon([[10, 89], [11, 89], [11, 89.5], [10, 89.5]])
        results = collect_ground_pixels(satellite, times, mask=mask)
        self.assertTrue(results.empty)

    def test_collect_ground_pixels_mask_excludes_everything_at_footprint_stage(self):
        """
        Test that a mask close enough to pass the (conservative)
        observable-period culling stage, but too far for any actual pixel footprint to
        intersect, still returns an empty result.
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
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(seconds=i)
            for i in range(10)
        ]
        subpoint = collect_orbit_track(satellite, [times[0]]).geometry.iloc[0]
        mask = ShapelyPoint(subpoint.x + 0.5, subpoint.y).buffer(0.01)
        results = collect_ground_pixels(satellite, times, mask=mask)
        self.assertTrue(results.empty)

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

    def test_collect_ground_track_repeat_cycle_far_from_epoch(self):
        """
        Test that, far from the epoch of an orbit with a repeat cycle, the
        ground track is propagated with the orbit's repeat track (as are
        observations), so that a mask culls the same footprints as clipping
        the unmasked ground track, and footprints repeat with the cycle.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(
            [
                "1 39084U 13008A   26213.27824675  .00000294  00000+0  75333-4 0  9990",
                "2 39084  98.2277 282.8718 0001275  92.4910 267.6434 14.57104473704466",
            ],
            remove_drag=True,
            repeat_cycle="auto",
        )
        repeat_cycle = orbit.get_repeat_cycle()
        satellite = Satellite(
            name="Landsat 8",
            orbit=orbit,
            instruments=[Instrument(name="Imager", field_of_regard=15)],
        )
        mask = box(-110, 30, -90, 50)
        start = orbit.get_epoch() + 6 * repeat_cycle + timedelta(days=4)
        times = [start + timedelta(seconds=10 * i) for i in range(int(86400 / 10))]
        masked = collect_ground_track(satellite, times, mask=mask)
        clipped = gpd.clip(collect_ground_track(satellite, times), mask)
        self.assertGreater(len(masked), 0)
        self.assertEqual(sorted(masked.time), sorted(clipped.time))
        first = collect_ground_track(
            satellite, [t - 6 * repeat_cycle for t in times], mask=mask
        )
        self.assertEqual(
            sorted(masked.time), sorted(t + 6 * repeat_cycle for t in first.time)
        )


class TestGroundTrackOfSatellites(IssConstellationTestCase):
    """
    Unit tests for the ground tracks and pixels of several satellites.
    """

    def test_satellites_equal_each_satellite(self):
        """
        Test that the ground track and ground pixels of several satellites
        equal those of each satellite, concatenated and sorted by time, with
        and without a mask.
        """
        instrument = PointedInstrument(
            name="Pixels",
            cross_track_field_of_view=10.0,
            along_track_field_of_view=10.0,
            is_rectangular=True,
            cross_track_pixels=3,
            along_track_pixels=3,
        )
        members = [
            member.model_copy(update={"instruments": [instrument]})
            for member in self.constellation.generate_members()
        ]
        times = [
            datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=2 * i)
            for i in range(30)
        ]
        mask = box(-90, -45, 90, 45)
        for collect in (collect_ground_track, collect_ground_pixels):
            for region in (None, mask):
                with self.subTest(collect=collect.__name__, mask=region is not None):
                    expected = (
                        pd.concat([collect(m, times, mask=region) for m in members])
                        .sort_values("time", kind="stable")
                        .reset_index(drop=True)
                    )
                    pd.testing.assert_frame_equal(
                        collect(members, times, mask=region), expected
                    )
