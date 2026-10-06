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
from tatc.utils.projection import VelocityFrame, compute_projected_ray_position


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
