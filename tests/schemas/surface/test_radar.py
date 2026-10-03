"""
Unit tests for the RadarStation schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import math
import unittest

from pydantic import ValidationError
from shapely.geometry import MultiPolygon, Polygon

from tatc.schemas import RadarBand, RadarStation, TerrainMask
from tatc.utils import compute_radar_ground_range_bounds


class TestRadarStation(unittest.TestCase):
    """
    Unit tests for the RadarStation schema.
    """

    def test_good_data(self):
        """
        Test that the RadarStation schema correctly initializes with valid data.
        """
        good_data = {
            "name": "KOUN",
            "latitude": 35.236,
            "longitude": -97.463,
            "elevation": 370,
            "max_range": 230000,
            "min_elevation_angle": 0.5,
            "max_elevation_angle": 19.5,
        }
        o = RadarStation(**good_data)
        self.assertEqual(o.name, good_data.get("name"))
        self.assertEqual(o.latitude, good_data.get("latitude"))
        self.assertEqual(o.longitude, good_data.get("longitude"))
        self.assertEqual(o.max_range, good_data.get("max_range"))
        self.assertEqual(o.min_elevation_angle, good_data.get("min_elevation_angle"))
        self.assertEqual(o.max_elevation_angle, good_data.get("max_elevation_angle"))

    def test_defaults(self):
        """
        Test that max_range, min_elevation_angle, max_elevation_angle, and
        beam_width default to conventional NEXRAD WSR-88D values when omitted.
        """
        o = RadarStation(name="test", latitude=35.236, longitude=-97.463)
        self.assertEqual(o.max_range, 230000)
        self.assertEqual(o.min_elevation_angle, 0.5)
        self.assertEqual(o.max_elevation_angle, 19.5)
        self.assertEqual(o.beam_width, 0.95)

    def test_bad_beam_width_negative(self):
        """
        Test that a negative beam_width is rejected.
        """
        with self.assertRaises(ValidationError):
            RadarStation(name="test", latitude=0, longitude=0, beam_width=-0.1)

    def test_zero_beam_width_bounds_beam_centers(self):
        """
        Test that a zero beam_width bounds coverage by the lowest and
        highest beam centers.
        """
        o = RadarStation(
            name="test", latitude=0, longitude=0, elevation=400, beam_width=0
        )
        self.assertEqual(
            o.compute_ground_ranges(3048),
            compute_radar_ground_range_bounds(0.5, 19.5, 230000, 3048, 400),
        )

    def test_beam_width_extends_ground_ranges(self):
        """
        Test that half the beam_width extends the observed elevation angles
        beyond the lowest and highest beam centers, widening the annulus
        at both edges.
        """
        o = RadarStation(name="test", latitude=0, longitude=0, elevation=400)
        centers = RadarStation(
            name="test", latitude=0, longitude=0, elevation=400, beam_width=0
        )
        self.assertEqual(
            o.compute_ground_ranges(3048),
            compute_radar_ground_range_bounds(0.025, 19.975, 230000, 3048, 400),
        )
        inner, outer = o.compute_ground_ranges(3048)
        self.assertLess(inner, centers.compute_ground_ranges(3048)[0])
        self.assertGreater(outer, centers.compute_ground_ranges(3048)[1])

    def test_get_effective_max_elevation_angle(self):
        """
        Test that the effective maximum elevation angle adds half the beam
        width and is capped at 90 degrees.
        """
        o = RadarStation(name="test", latitude=0, longitude=0)
        self.assertAlmostEqual(
            o.get_effective_max_elevation_angle(), 19.975, delta=1e-9
        )
        o = RadarStation(name="test", latitude=0, longitude=0, max_elevation_angle=90)
        self.assertEqual(o.get_effective_max_elevation_angle(), 90)

    def test_bad_name_missing(self):
        """
        Test that the RadarStation schema raises a ValidationError when the
        required name field is missing.
        """
        bad_data = {"latitude": 35.236, "longitude": -97.463}
        with self.assertRaises(ValidationError):
            RadarStation(**bad_data)

    def test_bad_max_range_not_positive(self):
        """
        Test that a non-positive max_range is rejected.
        """
        with self.assertRaises(ValidationError):
            RadarStation(name="test", latitude=0, longitude=0, max_range=0)

    def test_bad_elevation_angle_out_of_range(self):
        """
        Test that elevation angles outside [-90, 90] (minimum) or [0, 90]
        (maximum) are rejected.
        """
        with self.assertRaises(ValidationError):
            RadarStation(
                name="test", latitude=0, longitude=0, min_elevation_angle=-90.1
            )
        with self.assertRaises(ValidationError):
            RadarStation(name="test", latitude=0, longitude=0, max_elevation_angle=-0.1)
        with self.assertRaises(ValidationError):
            RadarStation(name="test", latitude=0, longitude=0, max_elevation_angle=90.1)

    def test_negative_min_elevation_angle(self):
        """
        Test that a negative minimum elevation angle (scanning below local
        horizontal from an elevated site) is accepted and extends coverage
        at a target height above the station.
        """
        o = RadarStation(
            name="test",
            latitude=0,
            longitude=0,
            elevation=2290,
            min_elevation_angle=-0.2,
        )
        self.assertEqual(o.min_elevation_angle, -0.2)
        nominal = RadarStation(name="test", latitude=0, longitude=0, elevation=2290)
        self.assertGreater(
            o.compute_ground_ranges(3048)[1], nominal.compute_ground_ranges(3048)[1]
        )

    def test_elevation_angle_boundary_values(self):
        """
        Test that elevation angle values exactly at bounds (0 and 90) are accepted.
        """
        o = RadarStation(
            name="test",
            latitude=0,
            longitude=0,
            min_elevation_angle=0,
            max_elevation_angle=90,
        )
        self.assertEqual(o.min_elevation_angle, 0)
        self.assertEqual(o.max_elevation_angle, 90)

    def test_bad_max_elevation_angle_less_than_min(self):
        """
        Test that a max_elevation_angle below min_elevation_angle is rejected.
        """
        with self.assertRaises(ValidationError):
            RadarStation(
                name="test",
                latitude=0,
                longitude=0,
                min_elevation_angle=10,
                max_elevation_angle=5,
            )

    def test_inherits_point_latitude_validation(self):
        """
        Test that RadarStation inherits Point's latitude validation.
        """
        with self.assertRaises(ValidationError):
            RadarStation(name="test", latitude=100, longitude=0)

    def test_compute_footprint_empty_at_and_below_station_elevation(self):
        """
        Test that a target at or below the station's own elevation is not
        observable when every beam departs above local horizontal.
        """
        station = RadarStation(name="test", latitude=0, longitude=0, elevation=2000)
        self.assertTrue(station.compute_footprint(elevation=2000).is_empty)
        self.assertTrue(station.compute_footprint(elevation=1500).is_empty)

    def test_compute_footprint_below_station_negative_tilt(self):
        """
        Test that a beam departing below local horizontal observes a target
        below the station in an annulus between its descending and climbing
        crossings, and nothing below the beam's lowest point.
        """
        station = RadarStation(
            name="test",
            latitude=0,
            longitude=0,
            elevation=2290,
            min_elevation_angle=-0.2,
        )
        footprint = station.compute_footprint(elevation=2000)
        self.assertIsInstance(footprint, Polygon)
        self.assertEqual(len(footprint.interiors), 1)
        self.assertTrue(station.compute_footprint(elevation=1500).is_empty)

    def test_compute_footprint_annulus_above_station(self):
        """
        Test that a target well above the station's elevation produces an
        annulus (one interior "cone of silence" ring).
        """
        station = RadarStation(name="test", latitude=0, longitude=0, elevation=0)
        footprint = station.compute_footprint(elevation=3048)
        self.assertIsInstance(footprint, Polygon)
        self.assertFalse(footprint.is_empty)
        self.assertEqual(len(footprint.interiors), 1)

    def test_compute_footprint_empty_when_unreachable(self):
        """
        Test that a target too high to be reached within max_range yields
        an empty footprint.
        """
        station = RadarStation(
            name="test", latitude=0, longitude=0, elevation=0, max_range=1000
        )
        footprint = station.compute_footprint(elevation=100000)
        self.assertTrue(footprint.is_empty)

    def test_compute_footprint_annulus_area_matches_ground_ranges(self):
        """
        Test that the computed footprint's area (in a local equal-distance
        projection, approximated here via the planar area near the
        equator) is consistent with the analytic annulus area
        pi * (outer^2 - inner^2), within a tolerance that accounts for
        oblate-Earth/projection distortion.
        """
        station = RadarStation(name="test", latitude=0, longitude=0, elevation=0)
        inner, outer = station.compute_ground_ranges(elevation=3048)
        footprint = station.compute_footprint(elevation=3048)
        # reproject to the same local equidistant CRS used internally to
        # compare areas on a comparable (meters) basis
        import pyproj
        from shapely.ops import transform

        to_crs = pyproj.Transformer.from_crs(
            "EPSG:4326", "+proj=eqc +lat_ts=0 +datum=WGS84 +units=m", always_xy=True
        )
        projected = transform(to_crs.transform, footprint)
        expected_area = math.pi * (outer**2 - inner**2)
        self.assertAlmostEqual(projected.area / expected_area, 1.0, delta=0.05)


class TestTerrainMask(unittest.TestCase):
    """
    Unit tests for the TerrainMask schema.
    """

    def test_good_data(self):
        """
        Test that the TerrainMask schema correctly initializes with valid data.
        """
        mask = TerrainMask(
            azimuth=[0, 90, 180, 270], min_elevation_angle=[0.5, 5, 0.5, 2]
        )
        self.assertEqual(mask.azimuth, [0, 90, 180, 270])
        self.assertEqual(mask.min_elevation_angle, [0.5, 5, 0.5, 2])

    def test_bad_mismatched_lengths(self):
        """
        Test that mismatched azimuth/min_elevation_angle lengths are rejected.
        """
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[0, 90, 180], min_elevation_angle=[0.5, 5])

    def test_bad_azimuth_out_of_range(self):
        """
        Test that azimuth values outside [0, 360) are rejected.
        """
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[-1, 90], min_elevation_angle=[0.5, 5])
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[0, 360], min_elevation_angle=[0.5, 5])

    def test_bad_azimuth_not_increasing(self):
        """
        Test that non-strictly-increasing azimuth values are rejected.
        """
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[90, 0], min_elevation_angle=[0.5, 5])
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[0, 0], min_elevation_angle=[0.5, 5])

    def test_bad_too_few_samples(self):
        """
        Test that a single-sample mask is rejected (at least two required).
        """
        with self.assertRaises(ValidationError):
            TerrainMask(azimuth=[0], min_elevation_angle=[0.5])

    def test_get_min_elevation_angle_at_samples(self):
        """
        Test that the interpolated value exactly at a sample azimuth
        matches that sample.
        """
        mask = TerrainMask(
            azimuth=[0, 90, 180, 270], min_elevation_angle=[0.5, 5, 0.5, 2]
        )
        for azimuth, expected in zip([0, 90, 180, 270], [0.5, 5, 0.5, 2]):
            self.assertAlmostEqual(
                mask.get_min_elevation_angle(azimuth), expected, delta=1e-9
            )

    def test_get_min_elevation_angle_interpolates_between_samples(self):
        """
        Test that the interpolated value midway between two samples is the
        linear average.
        """
        mask = TerrainMask(azimuth=[0, 90], min_elevation_angle=[0, 10])
        self.assertAlmostEqual(mask.get_min_elevation_angle(45), 5, delta=1e-9)

    def test_get_min_elevation_angle_wraps_around_zero(self):
        """
        Test that interpolation wraps around the 0/360 degree boundary
        (e.g. between the last sample and the first).
        """
        mask = TerrainMask(azimuth=[0, 270], min_elevation_angle=[0, 10])
        # halfway from 270 to 360 (0), wrapping, should be ~5
        self.assertAlmostEqual(mask.get_min_elevation_angle(315), 5, delta=1e-9)

    def test_get_min_elevation_angle_accepts_azimuth_outside_0_360(self):
        """
        Test that an azimuth outside [0, 360) is wrapped before lookup.
        """
        mask = TerrainMask(
            azimuth=[0, 90, 180, 270], min_elevation_angle=[0.5, 5, 0.5, 2]
        )
        self.assertAlmostEqual(
            mask.get_min_elevation_angle(90),
            mask.get_min_elevation_angle(450),
            delta=1e-9,
        )
        self.assertAlmostEqual(
            mask.get_min_elevation_angle(0),
            mask.get_min_elevation_angle(-360),
            delta=1e-9,
        )


class TestRadarStationTerrainMask(unittest.TestCase):
    """
    Unit tests for RadarStation's terrain_mask-aware behavior.
    """

    def setUp(self):
        self.mask = TerrainMask(
            azimuth=[0, 45, 90, 180, 270, 315],
            min_elevation_angle=[0.5, 5, 0.5, 0.5, 0.5, 5],
        )
        self.station = RadarStation(
            name="test", latitude=0, longitude=0, elevation=0, terrain_mask=self.mask
        )
        self.unmasked = RadarStation(name="test", latitude=0, longitude=0, elevation=0)

    def test_get_effective_min_elevation_angle_uses_terrain(self):
        """
        Test that the effective minimum elevation angle is raised at a
        blocked azimuth and matches the nominal value at a clear azimuth.
        """
        self.assertAlmostEqual(
            self.station.get_effective_min_elevation_angle(45), 5, delta=1e-6
        )
        self.assertAlmostEqual(
            self.station.get_effective_min_elevation_angle(180),
            self.station.min_elevation_angle,
            delta=1e-6,
        )

    def test_get_effective_min_elevation_angle_capped_at_max(self):
        """
        Test that the effective minimum elevation angle never exceeds the
        effective maximum elevation angle, even if terrain blockage would
        otherwise be higher.
        """
        mask = TerrainMask(azimuth=[0, 180], min_elevation_angle=[0.5, 50])
        station = RadarStation(
            name="test",
            latitude=0,
            longitude=0,
            max_elevation_angle=19.5,
            terrain_mask=mask,
        )
        self.assertAlmostEqual(
            station.get_effective_min_elevation_angle(180), 19.975, delta=1e-9
        )

    def test_get_effective_min_elevation_angle_no_mask(self):
        """
        Test that, without a terrain_mask, the effective minimum elevation
        angle is the lower half-power edge of the lowest beam.
        """
        self.assertAlmostEqual(
            self.unmasked.get_effective_min_elevation_angle(45),
            self.unmasked.min_elevation_angle - self.unmasked.beam_width / 2,
            delta=1e-9,
        )

    def test_get_effective_min_elevation_angle_terrain_within_beam(self):
        """
        Test that terrain below the lowest beam center but above its lower
        half-power edge raises the effective minimum elevation angle to the
        terrain angle.
        """
        mask = TerrainMask(azimuth=[0, 180], min_elevation_angle=[0.3, 0.3])
        station = RadarStation(name="test", latitude=0, longitude=0, terrain_mask=mask)
        self.assertAlmostEqual(
            station.get_effective_min_elevation_angle(90), 0.3, delta=1e-9
        )

    def test_compute_ground_range_profile_shape(self):
        """
        Test that the ground range profile returns the requested number of
        azimuth samples, evenly spaced.
        """
        profile = self.station.compute_ground_range_profile(
            elevation=3048, number_points=36
        )
        self.assertEqual(len(profile), 36)
        azimuths = [azimuth for azimuth, _, _ in profile]
        self.assertEqual(azimuths, sorted(azimuths))
        self.assertAlmostEqual(azimuths[1] - azimuths[0], 10, delta=1e-9)

    def test_compute_ground_range_profile_shrinks_outer_at_blocked_azimuth(self):
        """
        Test that, for a target elevation well above the station, the outer
        ground range is smaller at a terrain-blocked azimuth than at a
        clear one.
        """
        profile = {
            azimuth: (inner, outer)
            for azimuth, inner, outer in self.station.compute_ground_range_profile(
                elevation=3048, number_points=360
            )
        }
        self.assertLess(profile[45][1], profile[180][1])

    def test_compute_ground_range_profile_inner_unaffected_by_terrain(self):
        """
        Test that the inner (cone of silence) ground range is the same at
        every non-blocked azimuth, since it depends only on
        max_elevation_angle, not terrain.
        """
        profile = self.station.compute_ground_range_profile(
            elevation=3048, number_points=360
        )
        inner_values = {round(inner, 3) for _, inner, outer in profile if outer > 0}
        self.assertEqual(len(inner_values), 1)

    def test_compute_footprint_with_terrain_mask_smaller_than_symmetric(self):
        """
        Test that a terrain-masked footprint at an elevated target is
        smaller in area than the azimuthally symmetric (unmasked) footprint.
        """
        masked = self.station.compute_footprint(elevation=3048)
        unmasked = self.unmasked.compute_footprint(elevation=3048)
        self.assertLess(masked.area, unmasked.area)

    def test_compute_footprint_with_terrain_mask_is_valid(self):
        """
        Test that a terrain-masked footprint is a valid, non-empty geometry.
        """
        footprint = self.station.compute_footprint(elevation=3048, number_points=72)
        self.assertIsInstance(footprint, (Polygon, MultiPolygon))
        self.assertTrue(footprint.is_valid)
        self.assertFalse(footprint.is_empty)

    def test_compute_footprint_with_terrain_mask_below_station(self):
        """
        Test that a terrain-masked footprint of a target below the station,
        whose inner bound varies with the terrain-raised lowest angle, is
        valid, non-empty, and smaller than the unmasked footprint.
        """
        mask = TerrainMask(
            azimuth=[0, 90, 180, 270], min_elevation_angle=[-1, -0.3, 0.5, -1]
        )
        masked = RadarStation(
            name="test",
            latitude=0,
            longitude=0,
            elevation=2290,
            min_elevation_angle=-0.2,
            terrain_mask=mask,
        )
        unmasked = RadarStation(
            name="test",
            latitude=0,
            longitude=0,
            elevation=2290,
            min_elevation_angle=-0.2,
        )
        footprint = masked.compute_footprint(elevation=2000, number_points=72)
        self.assertTrue(footprint.is_valid)
        self.assertFalse(footprint.is_empty)
        self.assertLess(footprint.area, unmasked.compute_footprint(elevation=2000).area)

    def test_compute_footprint_fully_blocked_is_empty(self):
        """
        Test that a station with every azimuth blocked beyond what the
        target elevation requires returns an empty footprint.
        """
        mask = TerrainMask(azimuth=[0, 180], min_elevation_angle=[19.5, 19.5])
        station = RadarStation(
            name="test", latitude=0, longitude=0, max_range=1000, terrain_mask=mask
        )
        footprint = station.compute_footprint(elevation=100000, number_points=36)
        self.assertTrue(footprint.is_empty)


class TestRadarBand(unittest.TestCase):
    """
    Unit tests for the RadarBand enum and RadarStation.from_band.
    """

    def test_band_defaults_to_none(self):
        """
        Test that a plain RadarStation has no band tagged by default.
        """
        station = RadarStation(name="test", latitude=0, longitude=0)
        self.assertIsNone(station.band)

    def test_from_band_tags_band(self):
        """
        Test that from_band tags the constructed station with the
        requested band.
        """
        for band in (RadarBand.S, RadarBand.C, RadarBand.X):
            station = RadarStation.from_band(band, name="test", latitude=0, longitude=0)
            self.assertEqual(station.band, band)

    def test_from_band_max_range_decreases_s_to_x(self):
        """
        Test that the nominal max_range decreases from S-band to C-band to
        X-band, reflecting greater susceptibility to attenuation at
        shorter wavelengths.
        """
        s = RadarStation.from_band(RadarBand.S, name="s", latitude=0, longitude=0)
        c = RadarStation.from_band(RadarBand.C, name="c", latitude=0, longitude=0)
        x = RadarStation.from_band(RadarBand.X, name="x", latitude=0, longitude=0)
        self.assertGreater(s.max_range, c.max_range)
        self.assertGreater(c.max_range, x.max_range)

    def test_from_band_s_matches_default_max_range(self):
        """
        Test that the S-band nominal max_range matches RadarStation's own
        default (both represent the same NEXRAD WSR-88D convention).
        """
        s = RadarStation.from_band(RadarBand.S, name="s", latitude=0, longitude=0)
        default = RadarStation(name="default", latitude=0, longitude=0)
        self.assertEqual(s.max_range, default.max_range)

    def test_from_band_max_range_override(self):
        """
        Test that an explicit max_range keyword argument overrides the
        band's nominal default.
        """
        station = RadarStation.from_band(
            RadarBand.X, name="test", latitude=0, longitude=0, max_range=12345
        )
        self.assertEqual(station.max_range, 12345)

    def test_from_band_passes_through_other_fields(self):
        """
        Test that other RadarStation fields (e.g. location) pass through
        from_band's keyword arguments unchanged.
        """
        station = RadarStation.from_band(
            RadarBand.C, name="test", latitude=12.5, longitude=-34.5, elevation=100
        )
        self.assertEqual(station.latitude, 12.5)
        self.assertEqual(station.longitude, -34.5)
        self.assertEqual(station.elevation, 100)

    def test_band_does_not_affect_footprint_geometry(self):
        """
        Test that tagging a band (beyond its effect on max_range) does not
        change footprint computation: a plain RadarStation and one built
        via from_band with an identical, explicit max_range produce the
        same footprint.
        """
        plain = RadarStation(name="plain", latitude=0, longitude=0, max_range=60000)
        tagged = RadarStation.from_band(
            RadarBand.X, name="tagged", latitude=0, longitude=0, max_range=60000
        )
        self.assertTrue(
            plain.compute_footprint(elevation=3048).equals(
                tagged.compute_footprint(elevation=3048)
            )
        )
