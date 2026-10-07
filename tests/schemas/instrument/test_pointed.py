"""
Unit tests for the PointedInstrument schema.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from pydantic import ValidationError
from skyfield.api import EarthSatellite, wgs84

from tatc.constants import timescale
from tatc.schemas import CircularOrbit, PointedInstrument
from tatc.utils import field_of_regard_to_swath_width, geodesic_distance
from tatc.utils.projection import (
    VelocityFrame,
    ViewGeometry,
    compute_projected_ray_position,
    compute_view_angles,
)


class TestPointedInstrument(unittest.TestCase):
    """
    Unit tests for the PointedInstrument schema.
    """

    def setUp(self):
        noon_utc = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        self.test_time = timescale.from_datetime(noon_utc)
        self.test_sat = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=500000,
                true_anomaly=0,
                epoch=noon_utc,
                inclination=0.0,
                right_ascension_ascending_node=0.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        self.orbit_track = self.test_sat.at(self.test_time)  # type: ignore

    def test_good_data(self):
        """
        Test that a PointedInstrument can be created with valid data.
        """
        good_data = {
            "name": "Test Instrument",
            "cross_track_field_of_view": 20.0,
            "along_track_field_of_view": 10.0,
            "roll_angle": 5.0,
            "pitch_angle": -5.0,
            "is_rectangular": True,
            "cross_track_pixels": 4,
            "along_track_pixels": 2,
            "cross_track_oversampling": 0.1,
            "along_track_oversampling": 0.2,
        }
        o = PointedInstrument(**good_data)
        self.assertEqual(
            o.cross_track_field_of_view, good_data["cross_track_field_of_view"]
        )
        self.assertEqual(
            o.along_track_field_of_view, good_data["along_track_field_of_view"]
        )
        self.assertEqual(o.roll_angle, good_data["roll_angle"])
        self.assertEqual(o.pitch_angle, good_data["pitch_angle"])
        self.assertEqual(o.is_rectangular, good_data["is_rectangular"])
        self.assertEqual(o.cross_track_pixels, good_data["cross_track_pixels"])
        self.assertEqual(o.along_track_pixels, good_data["along_track_pixels"])
        self.assertEqual(
            o.cross_track_oversampling, good_data["cross_track_oversampling"]
        )
        self.assertEqual(
            o.along_track_oversampling, good_data["along_track_oversampling"]
        )

    def test_defaults(self):
        """
        Test that optional fields default to a single non-rectangular,
        boresight-pointed pixel with no oversampling.
        """
        o = PointedInstrument(
            name="Test Instrument",
            cross_track_field_of_view=20.0,
            along_track_field_of_view=10.0,
        )
        self.assertEqual(o.roll_angle, 0)
        self.assertEqual(o.pitch_angle, 0)
        self.assertFalse(o.is_rectangular)
        self.assertEqual(o.cross_track_pixels, 1)
        self.assertEqual(o.along_track_pixels, 1)
        self.assertEqual(o.cross_track_oversampling, 0)
        self.assertEqual(o.along_track_oversampling, 0)
        self.assertEqual(o.velocity_frame, VelocityFrame.EARTH_FIXED)

    def test_field_of_view_bounds(self):
        """
        Test that cross/along track field of view must be in (0, 180].
        """
        PointedInstrument(
            name="t", cross_track_field_of_view=180.0, along_track_field_of_view=180.0
        )
        with self.assertRaises(ValidationError):
            PointedInstrument(
                name="t", cross_track_field_of_view=0.0, along_track_field_of_view=10.0
            )
        with self.assertRaises(ValidationError):
            PointedInstrument(
                name="t",
                cross_track_field_of_view=10.0,
                along_track_field_of_view=180.1,
            )

    def test_angle_bounds(self):
        """
        Test that roll/pitch angle must be in [-180, 180].
        """
        PointedInstrument(
            name="t",
            cross_track_field_of_view=10.0,
            along_track_field_of_view=10.0,
            roll_angle=180.0,
            pitch_angle=-180.0,
        )
        with self.assertRaises(ValidationError):
            PointedInstrument(
                name="t",
                cross_track_field_of_view=10.0,
                along_track_field_of_view=10.0,
                roll_angle=180.1,
            )
        with self.assertRaises(ValidationError):
            PointedInstrument(
                name="t",
                cross_track_field_of_view=10.0,
                along_track_field_of_view=10.0,
                pitch_angle=-180.1,
            )

    def test_pixel_count_bounds(self):
        """
        Test that cross/along track pixel counts must be at least 1.
        """
        with self.assertRaises(ValidationError):
            PointedInstrument(
                name="t",
                cross_track_field_of_view=10.0,
                along_track_field_of_view=10.0,
                cross_track_pixels=0,
            )
        with self.assertRaises(ValidationError):
            PointedInstrument(
                name="t",
                cross_track_field_of_view=10.0,
                along_track_field_of_view=10.0,
                along_track_pixels=0,
            )

    def test_oversampling_bounds(self):
        """
        Test that cross/along track oversampling must be in [0, 1).
        """
        with self.assertRaises(ValidationError):
            PointedInstrument(
                name="t",
                cross_track_field_of_view=10.0,
                along_track_field_of_view=10.0,
                cross_track_oversampling=-0.1,
            )
        with self.assertRaises(ValidationError):
            PointedInstrument(
                name="t",
                cross_track_field_of_view=10.0,
                along_track_field_of_view=10.0,
                along_track_oversampling=1.0,
            )

    def test_get_cross_track_instantaneous_field_of_view_single_pixel(self):
        """
        Test that a single pixel's instantaneous field of view equals the
        full cross track field of view when there is no oversampling.
        """
        o = PointedInstrument(
            name="t", cross_track_field_of_view=20.0, along_track_field_of_view=10.0
        )
        self.assertAlmostEqual(o.get_cross_track_instantaneous_field_of_view(), 20.0)

    def test_get_cross_track_instantaneous_field_of_view_with_oversampling(self):
        """
        Test that oversampling (fractional pixel overlap) inflates each
        pixel's instantaneous field of view beyond the non-overlapping
        pixel pitch (field of view / pixel count), consistent with the
        pixel footprints overlapping by the requested fraction.
        """
        o = PointedInstrument(
            name="t",
            cross_track_field_of_view=20.0,
            along_track_field_of_view=10.0,
            cross_track_pixels=4,
            cross_track_oversampling=0.5,
        )
        # pitch = 20 / 4 = 5 degrees; ifov = pitch / (1 - 0.5) = 10 degrees
        self.assertAlmostEqual(o.get_cross_track_instantaneous_field_of_view(), 10.0)

    def test_get_along_track_instantaneous_field_of_view_with_oversampling(self):
        """
        Test that oversampling inflates the along track instantaneous
        field of view analogously to the cross track case.
        """
        o = PointedInstrument(
            name="t",
            cross_track_field_of_view=20.0,
            along_track_field_of_view=10.0,
            along_track_pixels=2,
            along_track_oversampling=0.5,
        )
        # pitch = 10 / 2 = 5 degrees; ifov = pitch / (1 - 0.5) = 10 degrees
        self.assertAlmostEqual(o.get_along_track_instantaneous_field_of_view(), 10.0)

    def test_get_pixel_cone_and_clock_angle_single_pixel(self):
        """
        Test that a single pixel (spanning the full field of view) is
        centered on boresight: zero cone angle.
        """
        o = PointedInstrument(
            name="t", cross_track_field_of_view=20.0, along_track_field_of_view=10.0
        )
        cone, clock = o.get_pixel_cone_and_clock_angle(0, 0)
        self.assertAlmostEqual(cone, 0.0)

    def test_get_pixel_cone_and_clock_angle_cross_track(self):
        """
        Test the cone and clock angles for a 2-pixel cross-track-only
        array: each pixel is offset a quarter of the full field of view
        from boresight, in opposite cross-track (clock 0/180) directions.
        """
        o = PointedInstrument(
            name="t",
            cross_track_field_of_view=20.0,
            along_track_field_of_view=10.0,
            cross_track_pixels=2,
        )
        cone_0, clock_0 = o.get_pixel_cone_and_clock_angle(0, 0)
        cone_1, clock_1 = o.get_pixel_cone_and_clock_angle(1, 0)
        self.assertAlmostEqual(cone_0, 5.0)
        self.assertAlmostEqual(cone_1, 5.0)
        self.assertAlmostEqual(clock_0, 180.0)
        self.assertAlmostEqual(clock_1, 0.0)

    def test_get_pixel_cone_and_clock_angle_along_track(self):
        """
        Test the cone and clock angles for a 2-pixel along-track-only
        array: each pixel is offset a quarter of the full field of view
        from boresight, in opposite along-track (clock 90/-90) directions.
        """
        o = PointedInstrument(
            name="t",
            cross_track_field_of_view=20.0,
            along_track_field_of_view=10.0,
            along_track_pixels=2,
        )
        cone_0, clock_0 = o.get_pixel_cone_and_clock_angle(0, 0)
        cone_1, clock_1 = o.get_pixel_cone_and_clock_angle(0, 1)
        self.assertAlmostEqual(cone_0, 2.5)
        self.assertAlmostEqual(cone_1, 2.5)
        self.assertAlmostEqual(clock_0, 90.0)
        self.assertAlmostEqual(clock_1, -90.0)

    def test_compute_footprint_center_matches_subpoint_when_unpointed(self):
        """
        Test that with zero roll/pitch, the footprint center coincides
        with Skyfield's own WGS 84 sub-satellite point (same as a plain
        nadir Instrument).
        """
        o = PointedInstrument(
            name="t", cross_track_field_of_view=20.0, along_track_field_of_view=10.0
        )
        center = o.compute_footprint_center(self.orbit_track)
        subpoint = wgs84.subpoint_of(self.orbit_track)
        self.assertAlmostEqual(
            geodesic_distance(
                center.longitude.degrees,
                center.latitude.degrees,
                subpoint.longitude.degrees,
                subpoint.latitude.degrees,
            ),
            0,
            delta=1e-3,
        )

    def test_compute_footprint_center_moves_with_roll(self):
        """
        Test that increasing the magnitude of the roll angle moves the
        footprint center monotonically farther from the sub-satellite
        point, in both directions.
        """
        subpoint = wgs84.subpoint_of(self.orbit_track)

        def distance_for_roll(roll):
            o = PointedInstrument(
                name="t",
                cross_track_field_of_view=1.0,
                along_track_field_of_view=1.0,
                roll_angle=roll,
            )
            center = o.compute_footprint_center(self.orbit_track)
            return geodesic_distance(
                subpoint.longitude.degrees,
                subpoint.latitude.degrees,
                center.longitude.degrees,
                center.latitude.degrees,
            )

        distances = [distance_for_roll(r) for r in (0, 5, 10, 20)]
        self.assertEqual(distances, sorted(distances))
        distances_negative = [distance_for_roll(r) for r in (0, -5, -10, -20)]
        self.assertEqual(distances_negative, sorted(distances_negative))

    def test_compute_footprint_center_velocity_frame(self):
        """
        Test that the velocity frame passes through to projections: in an
        inclined orbit, a pitched footprint center moves with an inertial
        velocity frame, while an unpointed center does not.
        """
        orbit_track = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=500000,
                true_anomaly=0,
                epoch=self.test_time.utc_datetime(),
                inclination=51.6,
                right_ascension_ascending_node=0.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        ).at(self.test_time)

        def center(velocity_frame, pitch):
            o = PointedInstrument(
                name="t",
                cross_track_field_of_view=1.0,
                along_track_field_of_view=1.0,
                pitch_angle=pitch,
                velocity_frame=velocity_frame,
            )
            return o.compute_footprint_center(orbit_track)

        for pitch, moves in [(0, False), (20, True)]:
            earth_fixed = center(VelocityFrame.EARTH_FIXED, pitch)
            inertial = center("inertial", pitch)
            distance = geodesic_distance(
                earth_fixed.longitude.degrees,
                earth_fixed.latitude.degrees,
                inertial.longitude.degrees,
                inertial.latitude.degrees,
            )
            if moves:
                self.assertGreater(distance, 1e3)
            else:
                self.assertAlmostEqual(distance, 0, delta=1e-3)

    def test_is_in_field_of_view_rolled_edges(self):
        """
        Test that a rolled view is rotated rigidly: a 40 deg view rolled by 10
        deg spans scan angles from -10 to 30 deg.
        """
        o = PointedInstrument(
            name="t",
            cross_track_field_of_view=40.0,
            along_track_field_of_view=10.0,
            roll_angle=10,
            is_rectangular=True,
        )
        for roll, expected in [
            (29.9, True),
            (30.1, False),
            (-9.9, True),
            (-10.1, False),
        ]:
            target = compute_projected_ray_position(self.orbit_track, 0, 0, roll)
            self.assertEqual(
                o.is_in_field_of_view(self.orbit_track, target).tolist(),
                [expected],
                f"roll {roll}",
            )

    def test_is_in_field_of_view(self):
        """
        Test that targets just inside (outside) the edges of rectangular and
        elliptical rolled and pitched views are (are not) in the field of
        view, with targets placed relative to the view center.
        """
        for is_rectangular in (True, False):
            o = PointedInstrument(
                name="t",
                cross_track_field_of_view=40.0,
                along_track_field_of_view=10.0,
                roll_angle=10,
                pitch_angle=-5,
                is_rectangular=is_rectangular,
            )

            def target(cross, along, rectangular=False, angle=0):
                # a ray offset from the view center by the given half widths
                return compute_projected_ray_position(
                    self.orbit_track, 2 * cross, 2 * along, 10, -5, rectangular, angle
                )

            for position, expected in [
                (target(0, 0), True),
                (target(19.9, 0), True),
                (target(20.1, 0), False),
                (target(19.9, 0, angle=180), True),
                (target(20.1, 0, angle=180), False),
                (target(0, 4.9, angle=90), True),
                (target(0, 5.1, angle=90), False),
                (target(0, 4.9, angle=270), True),
                (target(0, 5.1, angle=270), False),
                # near a corner: inside the rectangle, outside the ellipse
                (target(19.5, 4.8, True, 15), is_rectangular),
            ]:
                self.assertEqual(
                    o.is_in_field_of_view(self.orbit_track, position).tolist(),
                    [expected],
                    f"rectangular {is_rectangular}",
                )

    def test_compute_projected_pixel_position_matches_direct_cone_offset(self):
        """
        Test that a pixel's projected position matches an independent
        computation using its cone angle as a direct roll offset (the
        angular displacement convention used elsewhere, e.g. `roll_angle`,
        which is not halved). This is a regression test for a bug where
        `compute_projected_pixel_position` passed the pixel's cone angle
        directly as a field of view, which `compute_projected_ray_position`
        halves internally, landing pixels at roughly half their intended
        offset from boresight.
        """
        o = PointedInstrument(
            name="t",
            cross_track_field_of_view=20.0,
            along_track_field_of_view=20.0,
            cross_track_pixels=2,
        )
        cone, _ = o.get_pixel_cone_and_clock_angle(1, 0)
        pixel_position = o.compute_projected_pixel_position(self.orbit_track, 1, 0)
        # clock = 0 for this pixel, so its offset is purely in the roll direction
        direct_position = compute_projected_ray_position(
            self.orbit_track,
            cross_track_field_of_view=0,
            along_track_field_of_view=0,
            roll_angle=cone,
            pitch_angle=0,
            is_rectangular=False,
            angle=0,
            elevation=0,
        )
        self.assertAlmostEqual(
            geodesic_distance(
                pixel_position.longitude.degrees,
                pixel_position.latitude.degrees,
                direct_position.longitude.degrees,
                direct_position.latitude.degrees,
            ),
            0,
            delta=50.0,
        )

    def test_compute_footprint_pixel_array_shape(self):
        """
        Test that the pixel array contains one point per cross/along
        track pixel, for each orbit track time.
        """
        o = PointedInstrument(
            name="t",
            cross_track_field_of_view=20.0,
            along_track_field_of_view=10.0,
            cross_track_pixels=3,
            along_track_pixels=2,
        )
        times = timescale.utc(2020, 3, 20, 12, 0, [0, 1, 2])  # type: ignore
        orbit_track = self.test_sat.at(times)  # type: ignore
        pixel_arrays = o.compute_footprint_pixel_array(orbit_track)
        self.assertEqual(len(pixel_arrays), 3)
        for pixel_array in pixel_arrays:
            self.assertEqual(len(pixel_array.geoms), 6)

    def test_compute_footprint_cross_track_extent_matches_swath_width(self):
        """
        Test that the cross-track and along-track extents of a computed
        footprint match `field_of_regard_to_swath_width` evaluated at the
        satellite's actual height above the WGS 84 ellipsoid, cross-checking
        the closed-form swath width formula against the WGS 84 ellipsoid
        footprint geometry (same validation as for the nadir Instrument,
        extended to distinct cross/along track fields of view).
        """
        o = PointedInstrument(
            name="t", cross_track_field_of_view=30.0, along_track_field_of_view=10.0
        )
        height = wgs84.geographic_position_of(self.orbit_track).elevation.m
        expected_cross_track = field_of_regard_to_swath_width(
            height, o.cross_track_field_of_view
        )
        expected_along_track = field_of_regard_to_swath_width(
            height, o.along_track_field_of_view
        )
        # angle = [0, 90, 180, 270, 360] degrees, sampling both axes exactly
        footprint = o.compute_footprint(self.orbit_track, number_points=5)
        coords = list(footprint[0].exterior.coords)
        cross_track_extent = geodesic_distance(
            coords[0][0], coords[0][1], coords[2][0], coords[2][1]
        )
        along_track_extent = geodesic_distance(
            coords[1][0], coords[1][1], coords[3][0], coords[3][1]
        )
        self.assertAlmostEqual(cross_track_extent, expected_cross_track, delta=50.0)
        self.assertAlmostEqual(along_track_extent, expected_along_track, delta=50.0)


if __name__ == "__main__":
    unittest.main()


class TestRollAngleProfile(unittest.TestCase):
    """
    Unit tests for a PointedInstrument roll angle that varies around the
    orbit (roll_angle_profile).
    """

    def setUp(self):
        epoch = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        self.satellite = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=700000,
                true_anomaly=0,
                epoch=epoch,
                inclination=98.0,
                right_ascension_ascending_node=0.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        period = 2 * np.pi * np.sqrt((6378137.0 + 700000) ** 3 / 3.986004418e14)
        # a quarter orbit apart: arguments of latitude near 0, 90, 180, 270 deg
        self.track = self.satellite.at(
            timescale.from_datetimes(
                [epoch + timedelta(seconds=float(k * period / 4)) for k in range(4)]
            )
        )
        self.base = dict(
            name="radar",
            field_of_regard=90,
            cross_track_field_of_view=10,
            along_track_field_of_view=1,
            roll_angle=-30,
            is_rectangular=True,
        )

    def test_default_uses_roll_angle(self):
        """
        Test that, without a profile, the roll angle is constant.
        """
        instrument = PointedInstrument(**self.base)
        self.assertEqual(instrument.get_roll_angle(self.track), -30)

    def test_profile_sorted_and_interpolated_periodically(self):
        """
        Test that the profile is sorted by argument of latitude and
        interpolated linearly, wrapping around 360 degrees.
        """
        instrument = PointedInstrument(
            **self.base, roll_angle_profile=[(270, -32), (90, -28)]
        )
        self.assertEqual(instrument.roll_angle_profile, [(90, -28), (270, -32)])
        roll = instrument.get_roll_angle(self.track)
        # 0 and 180 deg are midway between the profile points (wrapping at 360)
        np.testing.assert_allclose(roll, [-30, -28, -30, -32], atol=0.05)

    def test_profile_validation(self):
        """
        Test that arguments of latitude outside [0, 360) and roll angles
        outside [-180, 180] are rejected.
        """
        for profile in ([(360, -30)], [(-1, -30)], [(10, 181)], []):
            with self.assertRaises(ValidationError):
                PointedInstrument(**self.base, roll_angle_profile=profile)

    def test_constant_profile_matches_fixed_roll(self):
        """
        Test that a constant profile gives the same footprints and fields
        of view as the fixed roll angle.
        """
        fixed = PointedInstrument(**self.base)
        steered = PointedInstrument(**self.base, roll_angle_profile=[(0, -30)])
        for a, b in zip(
            fixed.compute_footprint(self.track), steered.compute_footprint(self.track)
        ):
            self.assertTrue(a.equals(b))
        target = fixed.compute_footprint_center(self.track[1])
        self.assertTrue(steered.is_in_field_of_view(self.track[1], target)[0])

    def test_profile_moves_view(self):
        """
        Test that the view center follows the profile: rolling 5 degrees
        farther to the right in the northern part of the orbit moves the
        footprint center farther from the ground track there only.
        """
        fixed = PointedInstrument(**self.base)
        steered = PointedInstrument(
            **self.base,
            roll_angle_profile=[(0, -30), (90, -35), (180, -30), (270, -30)],
        )
        sub = wgs84.subpoint_of(self.track)
        distance = {}
        for name, instrument in [("fixed", fixed), ("steered", steered)]:
            center = instrument.compute_footprint_center(self.track)
            distance[name] = [
                geodesic_distance(
                    sub.longitude.degrees[k],
                    sub.latitude.degrees[k],
                    center.longitude.degrees[k],
                    center.latitude.degrees[k],
                )
                for k in range(4)
            ]
        self.assertGreater(distance["steered"][1], distance["fixed"][1] + 50e3)
        for k in (0, 2, 3):
            self.assertAlmostEqual(
                distance["steered"][k], distance["fixed"][k], delta=1e3
            )


class TestTiltedScanViewGeometry(unittest.TestCase):
    """
    Unit tests for a pitched PointedInstrument in scan geometry: the pitch
    angle tilts the scan plane, as for a scanner that tilts fore or aft
    (such as PACE OCI).
    """

    def setUp(self):
        epoch = datetime(2026, 10, 4, tzinfo=timezone.utc)
        satellite = EarthSatellite.from_satrec(
            CircularOrbit(mean_altitude=676e3, inclination=98.1, epoch=epoch)
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        self.track = satellite.at(timescale.from_datetime(epoch))
        # an OCI-like scanner tilted 20 deg forward
        self.base = dict(
            name="tilted scanner",
            field_of_regard=120,
            cross_track_field_of_view=112.9,
            along_track_field_of_view=0.1,
            pitch_angle=20,
            is_rectangular=True,
            cross_track_pixels=11,
        )

    def test_scan_line_matches_tilted_frame(self):
        """
        Test that the pixels of a tilted scan lie on the same line as those
        of a frame view pitched by the same angle (the plane containing the
        cross-track axis and the tilted boresight), rather than on a cone of
        constant along-track angle, which would be up to hundreds of
        kilometers farther forward at the swath edges.
        """
        scan = PointedInstrument(**self.base, view_geometry="scan")
        for index in range(11):
            pixel = scan.compute_projected_pixel_position(self.track, index, 0)
            cross_offset, _ = scan._get_pixel_offsets(index, 0)
            plane = compute_projected_ray_position(
                self.track,
                2 * abs(cross_offset),
                2 * abs(cross_offset),
                0,
                20,
                False,
                0 if cross_offset >= 0 else 180,
            )
            self.assertLess(
                geodesic_distance(
                    pixel.longitude.degrees,
                    pixel.latitude.degrees,
                    plane.longitude.degrees,
                    plane.latitude.degrees,
                ),
                1,
            )

    def test_field_of_view_follows_tilted_scan_plane(self):
        """
        Test that targets on the tilted scan line are in the field of view
        at all scan angles, and that targets displaced along track by more
        than half the along-track field of view are not.
        """
        scan = PointedInstrument(**self.base, view_geometry="scan")
        for roll in (-55, -30, 0, 30, 55):
            for pitch, inside in [(0, True), (0.04, True), (0.06, False)]:
                ray = compute_projected_ray_position(
                    self.track, 0, 0, roll, pitch, tilt_angle=20
                )
                target = wgs84.latlon(ray.latitude.degrees, ray.longitude.degrees)
                self.assertEqual(
                    bool(scan.is_in_field_of_view(self.track, target)[0]), inside
                )


class TestPitchAngleProfile(unittest.TestCase):
    """
    Unit tests for a PointedInstrument pitch angle that varies around the
    orbit (pitch_angle_profile).
    """

    def setUp(self):
        epoch = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        self.satellite = EarthSatellite.from_satrec(
            CircularOrbit(
                mean_altitude=700000,
                true_anomaly=0,
                epoch=epoch,
                inclination=98.0,
                right_ascension_ascending_node=0.0,
            )
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        period = 2 * np.pi * np.sqrt((6378137.0 + 700000) ** 3 / 3.986004418e14)
        # a quarter orbit apart: arguments of latitude near 0, 90, 180, 270 deg
        self.track = self.satellite.at(
            timescale.from_datetimes(
                [epoch + timedelta(seconds=float(k * period / 4)) for k in range(4)]
            )
        )
        self.base = dict(
            name="imager",
            field_of_regard=120,
            cross_track_field_of_view=100,
            along_track_field_of_view=1,
            pitch_angle=20,
            is_rectangular=True,
        )

    def test_default_uses_pitch_angle(self):
        """
        Test that, without a profile, the pitch angle is constant.
        """
        instrument = PointedInstrument(**self.base)
        self.assertEqual(instrument.get_pitch_angle(self.track), 20)
        self.assertEqual(instrument.get_pitch_angle(), 20)

    def test_profile_sorted_and_interpolated_periodically(self):
        """
        Test that the profile is sorted by argument of latitude and
        interpolated linearly, wrapping around 360 degrees.
        """
        instrument = PointedInstrument(
            **self.base, pitch_angle_profile=[(270, -20), (90, 20)]
        )
        self.assertEqual(instrument.pitch_angle_profile, [(90, 20), (270, -20)])
        pitch = instrument.get_pitch_angle(self.track)
        # 0 and 180 deg are midway between the profile points (wrapping at 360)
        np.testing.assert_allclose(pitch, [0, 20, 0, -20], atol=0.2)
        # the roll angle is unaffected
        self.assertEqual(instrument.get_roll_angle(self.track), 0)

    def test_profile_validation(self):
        """
        Test that arguments of latitude outside [0, 360) and pitch angles
        outside [-180, 180] are rejected.
        """
        for profile in ([(360, 20)], [(-1, 20)], [(10, -181)], []):
            with self.assertRaises(ValidationError) as context:
                PointedInstrument(**self.base, pitch_angle_profile=profile)
            if profile == [(10, -181)]:
                self.assertIn("Pitch angles", str(context.exception))

    def test_constant_profile_matches_fixed_pitch(self):
        """
        Test that a constant profile gives the same footprints, footprint
        centers, pixel positions, and fields of view as the fixed pitch angle,
        in both view geometries.
        """
        for view_geometry in ("frame", "scan"):
            fixed = PointedInstrument(**self.base, view_geometry=view_geometry)
            steered = PointedInstrument(
                **{**self.base, "pitch_angle": 0},
                view_geometry=view_geometry,
                pitch_angle_profile=[(0, 20)],
            )
            for a, b in zip(
                fixed.compute_footprint(self.track),
                steered.compute_footprint(self.track),
            ):
                self.assertTrue(a.equals(b))
            target = fixed.compute_footprint_center(self.track[1])
            self.assertAlmostEqual(
                steered.compute_footprint_center(self.track[1]).latitude.degrees,
                target.latitude.degrees,
            )
            self.assertTrue(steered.is_in_field_of_view(self.track[1], target)[0])
            pixel = fixed.compute_projected_pixel_position(self.track, 0, 0)
            np.testing.assert_allclose(
                steered.compute_projected_pixel_position(
                    self.track, 0, 0
                ).latitude.degrees,
                pixel.latitude.degrees,
            )

    def test_profile_moves_view(self):
        """
        Test that the view center follows the profile: pitching forward in
        the northern part of the orbit and aft in the southern part moves the
        footprint center ahead of the satellite in the north and behind it in
        the south, by about 700 km * tan(20 deg) = 255 km.
        """
        instrument = PointedInstrument(
            **{**self.base, "pitch_angle": 0},
            pitch_angle_profile=[(0, 0), (90, 20), (180, 0), (270, -20)],
        )
        sub = wgs84.subpoint_of(self.track)
        center = instrument.compute_footprint_center(self.track)
        offsets = [
            geodesic_distance(
                sub.longitude.degrees[k],
                sub.latitude.degrees[k],
                center.longitude.degrees[k],
                center.latitude.degrees[k],
            )
            for k in range(4)
        ]
        np.testing.assert_allclose(offsets, [0, 255e3, 0, 255e3], atol=15e3)
        # the center is near the subsatellite point 36 s (255 km) later in
        # the north (forward) and 36 s earlier in the south (aft)
        times = self.track.t.utc_datetime()
        for k, sign in [(1, 1), (3, -1)]:
            ahead, behind = (
                wgs84.subpoint_of(
                    self.satellite.at(
                        timescale.from_datetime(
                            times[k] + timedelta(seconds=float(direction * 36))
                        )
                    )
                )
                for direction in (sign, -sign)
            )
            to_ahead, to_behind = (
                geodesic_distance(
                    position.longitude.degrees,
                    position.latitude.degrees,
                    center.longitude.degrees[k],
                    center.latitude.degrees[k],
                )
                for position in (ahead, behind)
            )
            self.assertLess(to_ahead, 20e3)
            self.assertGreater(to_behind, 400e3)


class TestScanViewGeometry(unittest.TestCase):
    """
    Unit tests for PointedInstrument views defined in scan (angular)
    geometry, as for cross-track scanners.
    """

    def setUp(self):
        epoch = datetime(2026, 10, 1, tzinfo=timezone.utc)
        satellite = EarthSatellite.from_satrec(
            CircularOrbit(mean_altitude=834e3, inclination=98.7, epoch=epoch)
            .to_gp_orbit()
            .elements[0]
            .to_satrec(),
            timescale,
        )
        self.track = satellite.at(timescale.from_datetime(epoch))
        # a VIIRS-like scanner: 112.1 deg across track, 0.814 deg along track
        self.base = dict(
            name="scanner",
            field_of_regard=115,
            cross_track_field_of_view=112.1,
            along_track_field_of_view=0.814,
            is_rectangular=True,
        )

    def along_track_extent(self, instrument, cross_angle):
        """
        Distance (km) between the footprint's fore and aft edges at a
        cross-track angle (degrees): between the polygon vertices nearest
        that cross-track angle ahead of and behind the cross-track plane.
        """
        footprint = instrument.compute_footprint(self.track, number_points=200)[0]
        best = {}
        for lon, lat in np.array(footprint.exterior.coords)[:, :2]:
            roll, pitch = compute_view_angles(self.track, wgs84.latlon(lat, lon))
            side = pitch > 0
            if side not in best or abs(roll - cross_angle) < best[side][0]:
                best[side] = (abs(roll - cross_angle), lon, lat)
        (_, lon_0, lat_0), (_, lon_1, lat_1) = best[True], best[False]
        return geodesic_distance(lon_0, lat_0, lon_1, lat_1) / 1e3

    def test_default_frame_geometry(self):
        """
        Test that views are defined in frame geometry by default.
        """
        self.assertEqual(
            PointedInstrument(**self.base).view_geometry, ViewGeometry.FRAME
        )

    def test_scan_footprint_widens_away_from_nadir(self):
        """
        Test that, in scan geometry, the footprint's along-track extent is
        that of the rays at the edges of the along-track field of view: at
        nadir, about 11.8 km in both geometries; 50 deg across track, wider
        in scan geometry (where the along-track angular extent is constant)
        than in frame geometry (where it narrows with the cosine of the
        cross-track angle) by about 1 / cos(50 deg).
        """
        scan = PointedInstrument(**self.base, view_geometry=ViewGeometry.SCAN)
        frame = PointedInstrument(**self.base)
        half = self.base["along_track_field_of_view"] / 2
        fore = compute_projected_ray_position(self.track, 0, 0, 50, half)
        aft = compute_projected_ray_position(self.track, 0, 0, 50, -half)
        expected = (
            geodesic_distance(
                fore.longitude.degrees,
                fore.latitude.degrees,
                aft.longitude.degrees,
                aft.latitude.degrees,
            )
            / 1e3
        )
        self.assertAlmostEqual(self.along_track_extent(scan, 0), 11.8, delta=0.5)
        self.assertAlmostEqual(self.along_track_extent(frame, 0), 11.8, delta=0.5)
        self.assertAlmostEqual(self.along_track_extent(scan, 50), expected, delta=0.5)
        self.assertAlmostEqual(
            self.along_track_extent(scan, 50) / self.along_track_extent(frame, 50),
            1 / np.cos(np.radians(50)),
            delta=0.1,
        )

    def test_scan_field_of_view_bounds(self):
        """
        Test that, in scan geometry, a target just inside the angular
        bounds at the swath edge is in the field of view, and one just
        outside is not (in frame geometry, the along-track bound narrows
        there, so the inside target is outside).
        """
        scan = PointedInstrument(**self.base, view_geometry=ViewGeometry.SCAN)
        frame = PointedInstrument(**self.base)
        for pitch, inside in [(0.39, True), (0.42, False)]:
            ray = compute_projected_ray_position(self.track, 0, 0, 55, pitch)
            target = wgs84.latlon(ray.latitude.degrees, ray.longitude.degrees)
            self.assertEqual(scan.is_in_field_of_view(self.track, target)[0], inside)
            self.assertFalse(frame.is_in_field_of_view(self.track, target)[0])

    def test_scan_pixel_positions(self):
        """
        Test that, in scan geometry, pixel centers are evenly spaced in
        cross-track angle (the rays at the pixels' angular offsets).
        """
        scan = PointedInstrument(
            **{**self.base, "cross_track_field_of_view": 100},
            cross_track_pixels=5,
            view_geometry=ViewGeometry.SCAN,
        )
        for index, angle in enumerate([-40, -20, 0, 20, 40]):
            pixel = scan.compute_projected_pixel_position(self.track, index, 0)
            ray = compute_projected_ray_position(self.track, 0, 0, angle, 0)
            self.assertAlmostEqual(
                pixel.latitude.degrees, ray.latitude.degrees, places=8
            )
            self.assertAlmostEqual(
                pixel.longitude.degrees, ray.longitude.degrees, places=8
            )
