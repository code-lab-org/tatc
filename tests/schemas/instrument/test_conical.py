"""
Unit tests for the ConicalInstrument schema.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import numpy as np
from pydantic import ValidationError
from shapely.geometry import Point as ShapelyPoint
from skyfield.api import EarthSatellite, wgs84

from tatc.constants import timescale
from tatc.schemas import (
    CircularOrbit,
    ConicalInstrument,
    Instrument,
    PointedInstrument,
    Satellite,
)
from tatc.utils import (
    compute_cone_and_azimuth,
    compute_projected_ray_position,
    field_of_regard_to_swath_width,
)


class TestConicalInstrument(unittest.TestCase):
    """
    Unit tests for the ConicalInstrument schema.
    """

    def setUp(self):
        noon_utc = datetime(2020, 3, 20, 12, tzinfo=timezone.utc)
        self.orbit = CircularOrbit(
            mean_altitude=700000,
            true_anomaly=0,
            epoch=noon_utc,
            inclination=51.6,
            right_ascension_ascending_node=0.0,
        )
        self.satellite = EarthSatellite.from_satrec(
            self.orbit.to_gp_orbit().elements[0].to_satrec(), timescale
        )
        self.orbit_track = self.satellite.at(timescale.from_datetime(noon_utc))

    def ray(self, cone_angle, azimuth, seconds=0):
        """Ground position of a ray at a cone angle and scan azimuth, seconds after the test time."""
        c, a = np.radians(cone_angle), np.radians(azimuth)
        roll = np.degrees(np.arctan2(np.sin(c) * np.sin(a), np.cos(c)))
        pitch = np.degrees(np.arcsin(np.sin(c) * np.cos(a)))
        orbit_track = self.satellite.at(timescale.utc(2020, 3, 20, 12, 0, seconds))
        return compute_projected_ray_position(orbit_track, 0, 0, roll, pitch)

    def test_defaults(self):
        """
        Test that a conical instrument defaults to a full forward-centered
        rotation and a field of regard containing the scanned arc.
        """
        o = ConicalInstrument(name="t", cone_angle=45, along_track_field_of_view=2)
        self.assertEqual(o.scan_center_azimuth, 0)
        self.assertEqual(o.scan_half_width, 180)
        self.assertEqual(o.field_of_regard, 2 * (45 + 1 + 1))
        o = ConicalInstrument(
            name="t", cone_angle=45, along_track_field_of_view=2, field_of_regard=100
        )
        self.assertEqual(o.field_of_regard, 100)

    def test_bounds(self):
        """
        Test that the cone angle, field of view, and scan sector are bounded.
        """
        for kwargs in [
            {"cone_angle": 0, "along_track_field_of_view": 1},
            {"cone_angle": 90, "along_track_field_of_view": 1},
            {"cone_angle": 45, "along_track_field_of_view": 0},
            {"cone_angle": 45, "along_track_field_of_view": 1, "scan_half_width": 0},
            {"cone_angle": 45, "along_track_field_of_view": 1, "scan_half_width": 181},
            {
                "cone_angle": 45,
                "along_track_field_of_view": 1,
                "scan_center_azimuth": 181,
            },
        ]:
            with self.assertRaises(ValidationError):
                ConicalInstrument(name="t", **kwargs)

    def test_satellite_parses_conical_instrument(self):
        """
        Test that a satellite parses a conical instrument specification.
        """
        satellite = Satellite(
            name="t",
            orbit=self.orbit,
            instruments=[
                {"cone_angle": 45, "along_track_field_of_view": 1},
                {"cross_track_field_of_view": 10, "along_track_field_of_view": 1},
                {"field_of_regard": 50},
            ],
        )
        self.assertEqual(
            [type(i) for i in satellite.instruments],
            [ConicalInstrument, PointedInstrument, Instrument],
        )

    def test_get_swath_width(self):
        """
        Test that the swath width spans the full cone for a sector including
        the cross-track directions, and the sector's ends otherwise.
        """
        full = ConicalInstrument(name="t", cone_angle=45, along_track_field_of_view=1)
        self.assertAlmostEqual(
            full.get_swath_width(700e3),
            field_of_regard_to_swath_width(700e3, 90),
            delta=1,
        )
        forward = full.model_copy(update={"scan_half_width": 60})
        self.assertLess(forward.get_swath_width(700e3), full.get_swath_width(700e3))
        # a sector entirely to the left of the track
        left = full.model_copy(
            update={"scan_center_azimuth": 90, "scan_half_width": 30}
        )
        self.assertLess(left.get_swath_width(700e3), full.get_swath_width(700e3) / 2)

    def test_is_in_field_of_view(self):
        """
        Test that points on the arc within (outside) the scan sector are (are
        not) in the field of view, and that points on the arc a few seconds
        earlier or later are in the field of view if, and only if, the
        satellite's motion in that time is within the along-track distance
        swept by the footprint.
        """
        # 2 deg subtends 2 * 700 km * tan(1 deg) = 24 km at nadir, about 3.6 s
        o = ConicalInstrument(
            name="t", cone_angle=45, along_track_field_of_view=2, scan_half_width=60
        )
        for azimuth, seconds, expected in [
            (0, 0, True),
            (59, 0, True),
            (-59, 0, True),
            (61, 0, False),
            (180, 0, False),
            (0, 1.5, True),
            (0, -1.5, True),
            (0, 2.5, False),
            (0, -2.5, False),
            (50, 1.5, True),
            (50, -2.5, False),
        ]:
            self.assertEqual(
                o.is_in_field_of_view(
                    self.orbit_track, self.ray(45, azimuth, seconds)
                ).tolist(),
                [expected],
                f"azimuth {azimuth}, seconds {seconds}",
            )
        aft = o.model_copy(update={"scan_center_azimuth": 180})
        self.assertTrue(
            aft.is_in_field_of_view(self.orbit_track, self.ray(45, 170)).all()
        )
        self.assertFalse(
            aft.is_in_field_of_view(self.orbit_track, self.ray(45, 0)).all()
        )

    def test_compute_footprint(self):
        """
        Test that the footprint of a forward sector contains points on its arc
        and not nadir, and that of a full rotation contains points around the
        cone (except near the cross-track directions, where the swept band
        narrows to the arc itself).
        """
        subpoint = wgs84.subpoint_of(self.orbit_track)
        nadir = ShapelyPoint(subpoint.longitude.degrees, subpoint.latitude.degrees)
        point = lambda p: ShapelyPoint(p.longitude.degrees, p.latitude.degrees)
        forward = ConicalInstrument(
            name="t", cone_angle=45, along_track_field_of_view=2, scan_half_width=60
        )
        footprint = forward.compute_footprint(self.orbit_track)[0]
        self.assertTrue(
            footprint.contains(
                point(forward.compute_footprint_center(self.orbit_track))
            )
        )
        self.assertTrue(footprint.contains(point(self.ray(45, 50, 1.5))))
        self.assertFalse(footprint.contains(point(self.ray(45, 0, 2.5))))
        self.assertFalse(footprint.contains(point(self.ray(45, 70))))
        self.assertFalse(footprint.contains(nadir))
        full = forward.model_copy(update={"scan_half_width": 180})
        footprint = full.compute_footprint(self.orbit_track)[0]
        for azimuth in (0, 60, 120, 180, -120, -60):
            self.assertTrue(footprint.contains(point(self.ray(45, azimuth))), azimuth)
        self.assertFalse(footprint.contains(nadir))

    def test_wide_field_of_view_does_not_widen_swath(self):
        """
        Test that a wide along-track field of view extends the footprint
        along track everywhere on the arc, but not across track at its ends.
        """
        point = lambda p: ShapelyPoint(p.longitude.degrees, p.latitude.degrees)
        # 10 deg subtends about 18 s of motion
        o = ConicalInstrument(
            name="t", cone_angle=45, along_track_field_of_view=10, scan_half_width=75
        )
        footprint = o.compute_footprint(self.orbit_track)[0]
        for cone, azimuth, seconds, expected in [
            (45, 0, 8, True),
            (45, 74, 8, True),
            (45, 74, -8, True),
            (48, 74, 0, False),
        ]:
            target = self.ray(cone, azimuth, seconds)
            self.assertEqual(
                footprint.contains(point(target)), expected, (cone, azimuth)
            )
            self.assertEqual(
                o.is_in_field_of_view(self.orbit_track, target).tolist(), [expected]
            )

    def test_compute_footprint_vectorized(self):
        """
        Test that a vectorized orbit track produces one footprint per time.
        """
        o = ConicalInstrument(name="t", cone_angle=45, along_track_field_of_view=2)
        footprints = o.compute_footprint(
            self.satellite.at(timescale.utc(2020, 3, 20, 12, [0, 10, 20]))
        )
        self.assertEqual(len(footprints), 3)

    def test_compute_cone_and_azimuth(self):
        """
        Test that the cone angle and scan azimuth of a ray are recovered.
        """
        for cone, azimuth in [(30, 0), (45, 75), (45, -120), (20, 180)]:
            c, a = compute_cone_and_azimuth(self.orbit_track, self.ray(cone, azimuth))
            self.assertAlmostEqual(float(c), cone, delta=1e-6)
            self.assertAlmostEqual(
                ((float(a) - azimuth + 180) % 360) - 180, 0, delta=1e-6
            )
