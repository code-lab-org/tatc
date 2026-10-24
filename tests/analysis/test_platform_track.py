"""
Unit tests for surface platform analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import numpy as np
from pyproj import Geod
from shapely.geometry import Polygon, box
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
    TerrainMask,
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

    def test_collect_platform_track_zenith_footprint(self):
        """
        Test that a zenith cone projects to a circle of radius the height of
        a target above the instruments times the tangent of the half angle.
        """
        platform = self.platforms[1].model_copy(
            update={
                "instruments": [
                    PointedInstrument(
                        cross_track_field_of_view=60,
                        along_track_field_of_view=60,
                        roll_angle=180,
                    )
                ]
            }
        )
        track = collect_platform_track(platform, elevation=1340)
        distances = _distances(platform, track.geometry.iloc[0])
        # geodesic distances on the ellipsoid, below the target elevation, to
        # the points on the circle and between them (on great circles)
        radius = 1000 * np.tan(np.radians(30))
        self.assertAlmostEqual(max(distances), radius, delta=1e-3 * radius)
        self.assertGreater(min(distances), 0.99 * radius)
        elevations = [z for *_, z in track.geometry.iloc[0].exterior.coords]
        np.testing.assert_allclose(elevations, 1340)

    def test_collect_platform_track_ground_horizon(self):
        """
        Test that a view toward a target above the instruments that is partly
        blocked by the ground is clipped to the horizon, spanning the view's
        azimuths along the lines of sight tangent to the ground.
        """
        platform = self.platforms[1].model_copy(
            update={
                "instruments": [
                    PointedInstrument(
                        cross_track_field_of_view=40,
                        along_track_field_of_view=40,
                        roll_angle=90,
                        is_rectangular=True,
                    )
                ]
            }
        )
        footprint = collect_platform_track(platform, elevation=1310).geometry.iloc[0]
        self.assertEqual(footprint.geom_type, "Polygon")
        self.assertTrue(footprint.is_valid)
        distances = _distances(platform, footprint)
        # the upper edge of the view, 20 degrees above horizontal
        self.assertAlmostEqual(
            min(distances), 970 / np.tan(np.radians(20)), delta=0.01 * 970
        )
        # the line of sight tangent to the ground, extended to 1 km above it
        horizon = np.sqrt(2 * 6371e3 * 30) + np.sqrt(2 * 6371e3 * 1000)
        np.testing.assert_allclose(max(distances), horizon, rtol=0.005)
        # the horizon spans the view's azimuths, to the left of the heading
        lon, lat, *_ = zip(*footprint.exterior.coords)
        azimuths = GEOD.inv(
            np.full(len(lon), platform.longitude),
            np.full(len(lat), platform.latitude),
            lon,
            lat,
        )[0][distances > 0.995 * horizon]
        self.assertAlmostEqual(min(azimuths), -110, delta=0.5)
        self.assertAlmostEqual(max(azimuths), -70, delta=0.5)

    def test_collect_platform_track_ground_horizon_tilted(self):
        """
        Test that tilted views toward a target above the instruments that are
        partly blocked by the ground project to valid polygons.
        """
        for roll_angle, pitch_angle in [(95, -60), (100, 30), (90, 45), (120, -15)]:
            platform = self.platforms[1].model_copy(
                update={
                    "instruments": [
                        PointedInstrument(
                            cross_track_field_of_view=30,
                            along_track_field_of_view=20,
                            roll_angle=roll_angle,
                            pitch_angle=pitch_angle,
                            is_rectangular=True,
                        )
                    ]
                }
            )
            footprint = collect_platform_track(platform, elevation=1310).geometry
            self.assertEqual(footprint.iloc[0].geom_type, "Polygon")
            self.assertTrue(footprint.iloc[0].is_valid)

    def test_collect_platform_track_ground_blocked(self):
        """
        Test that a view about nadir toward a target above the instruments
        is blocked by the ground: entirely, with an empty footprint, or
        except for an annulus between the horizontal and the horizon.
        """
        platform = self.platforms[1].model_copy(
            update={
                "instruments": [
                    Instrument(field_of_regard=60),
                    Instrument(field_of_regard=180),
                ]
            }
        )
        footprints = collect_platform_track(
            platform, instrument_index=None, elevation=1310
        ).geometry
        self.assertTrue(footprints.iloc[0].is_empty)
        annulus = footprints.iloc[1]
        self.assertTrue(annulus.is_valid)
        self.assertEqual(len(annulus.interiors), 1)
        # horizontal lines of sight reach 970 m above the instruments
        hole = Polygon(annulus.interiors[0])
        np.testing.assert_allclose(
            _distances(platform, hole), np.sqrt(2 * 6371e3 * 970), rtol=0.005
        )
        np.testing.assert_allclose(
            _distances(platform, annulus),
            np.sqrt(2 * 6371e3 * 30) + np.sqrt(2 * 6371e3 * 1000),
            rtol=0.005,
        )

    def test_collect_platform_track_terrain_mask(self):
        """
        Test that a terrain mask clips footprints toward a target above the
        instruments to the region visible over it, but not footprints toward
        a target below the instruments.
        """
        mask = TerrainMask(
            azimuth=[0, 240, 250, 290, 300], min_elevation_angle=[-1, -1, 2, 2, -1]
        )
        platform = self.platforms[1].model_copy(
            update={
                "terrain_mask": mask,
                "instruments": [
                    PointedInstrument(
                        cross_track_field_of_view=40,
                        along_track_field_of_view=40,
                        roll_angle=90,
                        is_rectangular=True,
                    )
                ],
            }
        )
        unmasked = platform.model_copy(update={"terrain_mask": None})
        footprint = collect_platform_track(platform, elevation=1310).geometry.iloc[0]
        self.assertTrue(footprint.is_valid)
        self.assertTrue(footprint.has_z)
        self.assertTrue(
            footprint.equals(
                collect_platform_track(unmasked, elevation=1310)
                .geometry.iloc[0]
                .intersection(platform.compute_visible_region(1310))
            )
        )
        # the ridge (2 degrees above horizontal) limits the view to the west
        self.assertLess(max(_distances(platform, footprint)), 30000)
        # but not toward a target below the instruments
        self.assertTrue(
            collect_platform_track(platform)
            .geometry.iloc[0]
            .equals(collect_platform_track(unmasked).geometry.iloc[0])
        )

    def test_collect_platform_track_elevation_at_instruments(self):
        """
        Test that a target elevation at a platform's instruments raises a
        ValueError.
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
