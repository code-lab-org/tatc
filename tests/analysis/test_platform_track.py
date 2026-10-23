"""
Unit tests for surface platform analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import numpy as np
from pyproj import Geod
from shapely.geometry import box
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc.analysis import (
    collect_ground_track,
    collect_platform_track,
    compute_platform_track,
)
from tatc.schemas import (
    Instrument,
    PointedInstrument,
    Satellite,
    SurfacePlatform,
)
from tatc.utils import ViewGeometry

from .common import IssConstellationTestCase

GEOD = Geod(ellps="WGS84")


def _distances(platform, geometry):
    """
    Distances (meters) from a platform to the exterior points of a footprint.
    """
    lon, lat, *_ = zip(*geometry.exterior.coords)
    return GEOD.inv(
        np.full(len(lon), platform.longitude),
        np.full(len(lat), platform.latitude),
        lon,
        lat,
    )[2]


class TestPlatformAnalysis(IssConstellationTestCase):
    """
    Unit tests for surface platform analysis functions.
    """

    def setUp(self):
        super().setUp()
        self.camera = Instrument(name="Camera", field_of_regard=60)
        self.oblique = PointedInstrument(
            name="Oblique",
            cross_track_field_of_view=40,
            along_track_field_of_view=30,
            pitch_angle=60,
            is_rectangular=True,
        )
        self.platform = SurfacePlatform(
            name="Balloon",
            latitude=33.45,
            longitude=-112.07,
            elevation=300,
            height=500,
            heading=45,
            instruments=[self.camera, self.oblique],
        )
        self.platforms = [
            self.platform,
            SurfacePlatform(
                name="Tower",
                latitude=33.46,
                longitude=-112.06,
                elevation=310,
                height=30,
                instruments=[self.camera],
            ),
        ]

    def test_collect_platform_track_empty(self):
        """
        Test that an empty platform list produces an empty data frame.
        """
        track = collect_platform_track([])
        self.assertTrue(track.empty)
        self.assertListEqual(
            list(track.columns), ["platform", "instrument", "elevation", "geometry"]
        )

    def test_collect_platform_track_singular_matches_list(self):
        """
        Test that passing a single platform produces the same result as
        passing a list containing that one platform.
        """
        singular = collect_platform_track(self.platform)
        as_list = collect_platform_track([self.platform])
        self.assertEqual(len(singular), 1)
        self.assertTrue(singular.geometry.iloc[0].equals(as_list.geometry.iloc[0]))

    def test_collect_platform_track_instruments(self):
        """
        Test that collecting every instrument of several platforms produces
        one row per platform and instrument, in order.
        """
        track = collect_platform_track(self.platforms, instrument_index=None)
        self.assertListEqual(list(track.platform), ["Balloon", "Balloon", "Tower"])
        self.assertListEqual(list(track.instrument), ["Camera", "Oblique", "Camera"])
        self.assertListEqual(list(track.elevation), [300, 300, 310])
        track = collect_platform_track(self.platform, instrument_index=1)
        self.assertListEqual(list(track.instrument), ["Oblique"])

    def test_collect_platform_track_instrument_index_out_of_range(self):
        """
        Test that an instrument index out of range names the platform.
        """
        with self.assertRaisesRegex(IndexError, "platform 'Tower'"):
            collect_platform_track(self.platforms, instrument_index=1)

    def test_collect_platform_track_nadir_footprint(self):
        """
        Test that a nadir cone projects to a circle of radius height times the
        tangent of the half angle at the ground or a target elevation.
        """
        for elevation, height in [(None, 500), (600, 200)]:
            track = collect_platform_track(self.platform, elevation=elevation)
            distances = _distances(self.platform, track.geometry.iloc[0])
            np.testing.assert_allclose(
                distances, height * np.tan(np.radians(30)), atol=0.1
            )
            elevations = [z for *_, z in track.geometry.iloc[0].exterior.coords]
            np.testing.assert_allclose(elevations, track.elevation.iloc[0])

    def test_collect_platform_track_horizon(self):
        """
        Test that a hemispheric view projects to the horizon.
        """
        platform = self.platform.model_copy(
            update={"instruments": [Instrument(field_of_regard=180)]}
        )
        track = collect_platform_track(platform)
        distances = _distances(platform, track.geometry.iloc[0])
        np.testing.assert_allclose(distances, np.sqrt(2 * 6371e3 * 500), rtol=0.005)

    def test_collect_platform_track_oblique_footprint(self):
        """
        Test that a forward-pitched view projects ahead of the platform,
        toward its heading, between the near and far edges of the view.
        """
        footprint = collect_platform_track(
            self.platform, instrument_index=1
        ).geometry.iloc[0]
        azimuth, _, distance = GEOD.inv(
            self.platform.longitude,
            self.platform.latitude,
            footprint.centroid.x,
            footprint.centroid.y,
        )
        self.assertAlmostEqual(azimuth, 45, delta=0.1)
        self.assertGreater(distance, 500 * np.tan(np.radians(45)))
        self.assertLess(distance, 500 * np.tan(np.radians(75)))
        self.assertAlmostEqual(min(_distances(self.platform, footprint)), 500, delta=1)

    def test_collect_platform_track_scan_geometry(self):
        """
        Test that a scanning view rolled to the left of the heading projects
        to the left.
        """
        platform = self.platform.model_copy(
            update={
                "instruments": [
                    PointedInstrument(
                        cross_track_field_of_view=20,
                        along_track_field_of_view=2,
                        roll_angle=45,
                        view_geometry=ViewGeometry.SCAN,
                        is_rectangular=True,
                    )
                ]
            }
        )
        footprint = collect_platform_track(platform).geometry.iloc[0]
        azimuth, _, _ = GEOD.inv(
            platform.longitude,
            platform.latitude,
            footprint.centroid.x,
            footprint.centroid.y,
        )
        self.assertAlmostEqual(azimuth % 360, 315, delta=0.1)

    def test_collect_platform_track_elevation_above_instruments(self):
        """
        Test that a target elevation at or above a platform's instruments
        raises a ValueError.
        """
        with self.assertRaisesRegex(ValueError, "Tower"):
            collect_platform_track(self.platforms, elevation=340)
        with self.assertRaises(ValueError):
            collect_platform_track(self.platform, elevation=800)

    def test_collect_platform_track_mask(self):
        """
        Test that a mask clips results to the intersecting platform(s) only,
        in order.
        """
        far = self.platform.model_copy(update={"name": "Far", "longitude": -100})
        mask = box(-113, 33, -112, 34)
        track = collect_platform_track([far, *self.platforms], None, mask=mask)
        self.assertListEqual(list(track.platform), ["Balloon", "Balloon", "Tower"])
        self.assertListEqual(list(track.index), [0, 1, 2])

    def test_matches_satellite_ground_track(self):
        """
        Test that a platform at a satellite's position, heading along its
        ground track, projects the same footprint as the satellite.
        """
        time = datetime(2022, 6, 20, tzinfo=timezone.utc)
        instrument = PointedInstrument(
            cross_track_field_of_view=20,
            along_track_field_of_view=10,
            roll_angle=15,
            pitch_angle=-5,
            is_rectangular=True,
        )
        satellite = Satellite(name="ISS", orbit=self.orbit, instruments=[instrument])
        expected = collect_ground_track(satellite, [time]).geometry.iloc[0]
        orbit_track = self.orbit.to_gp_orbit().get_orbit_track([time])
        _, velocity = orbit_track.frame_xyz_and_velocity(itrs)
        subpoint = wgs84.geographic_position_of(orbit_track)
        lat = np.radians(subpoint.latitude.degrees[0])
        lon = np.radians(subpoint.longitude.degrees[0])
        east = [-np.sin(lon), np.cos(lon), 0]
        north = [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)]
        v = np.array(velocity.m_per_s)[:, 0]
        platform = SurfacePlatform(
            name="ISS",
            latitude=subpoint.latitude.degrees[0],
            longitude=subpoint.longitude.degrees[0],
            height=subpoint.elevation.m[0],
            heading=np.degrees(np.arctan2(np.dot(v, east), np.dot(v, north))) % 360,
            instruments=[instrument],
        )
        footprint = collect_platform_track(platform).geometry.iloc[0]
        np.testing.assert_allclose(
            np.array(footprint.exterior.coords),
            np.array(expected.exterior.coords),
            atol=1e-6,
        )

    def test_compute_platform_track_dissolves(self):
        """
        Test that compute_platform_track dissolves footprints into a single
        geometry.
        """
        merged = compute_platform_track(self.platforms, instrument_index=None)
        self.assertEqual(len(merged), 1)
        self.assertListEqual(list(merged.columns), ["geometry"])
        individual = collect_platform_track(self.platforms, instrument_index=None)
        self.assertTrue(merged.geometry.iloc[0].covers(individual.union_all()))

    def test_compute_platform_track_empty(self):
        """
        Test that an empty platform list produces an empty data frame.
        """
        self.assertTrue(compute_platform_track([]).empty)


if __name__ == "__main__":
    unittest.main()
