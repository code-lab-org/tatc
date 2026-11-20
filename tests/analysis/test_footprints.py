"""
Unit tests for the analysis of the footprints of any platforms.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import box

from tatc.analysis import (
    aggregate_observations,
    collect_airborne_ground_track,
    collect_footprint_observations,
    collect_footprints,
    collect_ground_pixels,
    collect_ground_track,
    collect_observations,
    collect_radar_track,
    collect_surface_track,
    dissolve_footprints,
    find_coincident_ground_track,
    find_coincident_observations,
)
from tatc.schemas import (
    Aircraft,
    CircularOrbit,
    Instrument,
    PointedInstrument,
    RadarStation,
    Satellite,
    SurfacePlatform,
    Trajectory,
)

from ..schemas.air.common import START, straight_trajectory
from .common import IssConstellationTestCase

T0 = datetime(2022, 6, 20, tzinfo=timezone.utc)
"""A time near the epoch of the ISS orbit's elements."""


def _camera(name="camera", **kwargs):
    """
    Gets a rectangular camera with a grid of 3 by 2 pixels.
    """
    return PointedInstrument(
        name=name,
        cross_track_field_of_view=40,
        along_track_field_of_view=30,
        is_rectangular=True,
        cross_track_pixels=3,
        along_track_pixels=2,
        **kwargs,
    )


class TestFootprints(IssConstellationTestCase):
    """
    Unit tests for the footprints of any platforms.
    """

    def setUp(self):
        super().setUp()
        self.data = straight_trajectory(n=120, latitude=35.0, longitude=-110.0)
        self.aircraft = Aircraft(
            name="P-3",
            trajectory=Trajectory.from_dataframe(self.data),
            instruments=[_camera()],
        )
        self.tower = SurfacePlatform(
            name="Tower",
            latitude=35.02,
            longitude=-109.99,
            height=100,
            instruments=[_camera(), Instrument(name="wide", field_of_regard=120)],
        )
        self.radar = RadarStation(name="KXXX", latitude=35.0, longitude=-110.0)
        self.times = list(pd.date_range(START, periods=120, freq="1s").to_pydatetime())
        # a satellite with an orbit near the times
        self.leo = Satellite(
            name="LEO",
            orbit=CircularOrbit(mean_altitude=500e3, inclination=50, epoch=START),
            instruments=[Instrument(name="wide", field_of_regard=120.0)],
        )

    def test_collect_footprints_matches_each_kind(self):
        """
        Test that the footprints of a single platform are those of its kind's
        analysis.
        """
        self.assertTrue(
            collect_footprints(self.aircraft, self.times[:5]).equals(
                collect_airborne_ground_track(self.aircraft, self.times[:5])
            )
        )
        self.assertTrue(
            collect_footprints(self.leo, self.times[:5]).equals(
                collect_ground_track(self.leo, self.times[:5])
            )
        )
        self.assertTrue(
            collect_footprints(self.tower, instrument_index=None).equals(
                collect_surface_track(self.tower, instrument_index=None)
            )
        )
        self.assertTrue(
            collect_footprints(self.radar, elevation=3000).equals(
                collect_radar_track(self.radar, 3000)
            )
        )

    def test_collect_footprints_mixed(self):
        """
        Test that the footprints of platforms of every kind are combined,
        sorted by time, with static footprints last.
        """
        footprints = collect_footprints(
            [self.radar, self.tower, self.aircraft, self.leo],
            self.times[:10],
            elevation=0,
            platform_altaz=True,
        )
        counts = {
            column: int(footprints[column].notna().sum())
            for column in ("satellite", "aircraft", "surface_platform", "radar_station")
        }
        self.assertDictEqual(
            counts,
            {
                "satellite": 10,
                "aircraft": 10,
                "surface_platform": 10,
                "radar_station": 1,
            },
        )
        for column in ("sat_alt", "aircraft_alt", "surface_platform_alt"):
            self.assertIn(column, footprints.columns)
        # sorted by time, with the radar (without a time) last
        timed = footprints["time"].iloc[:-1]
        self.assertTrue(timed.is_monotonic_increasing)
        self.assertTrue(pd.isna(footprints["time"].iloc[-1]))
        self.assertEqual(footprints["radar_station"].iloc[-1], "KXXX")
        # no beam of the radar reaches the ground at its own elevation
        self.assertFalse(footprints["valid_obs"].iloc[-1])
        # the aircraft passes over the tower's view
        pairs = find_coincident_ground_track(footprints, timedelta(seconds=5))
        self.assertTrue(
            (
                (pairs["aircraft_1"] == "P-3")
                & (pairs["surface_platform_2"] == "Tower")
            ).any()
        )

    def test_collect_footprints_radar_overpass(self):
        """
        Test that the overpass of a satellite over a radar's coverage at an
        elevation is coincident with it.
        """
        track = collect_ground_track(self.leo, self.times, elevation=3000)
        center = track.geometry.iloc[60].centroid
        radar = self.radar.model_copy(
            update={"latitude": center.y, "longitude": center.x}
        )
        footprints = collect_footprints([self.leo, radar], self.times, elevation=3000)
        pairs = find_coincident_ground_track(footprints)
        self.assertIn(self.times[60], set(pairs["time_1"]))
        self.assertTrue((pairs["radar_station_2"] == "KXXX").all())

    def test_collect_footprints_arguments(self):
        """
        Test that arguments that do not apply to a kind of platform raise
        errors.
        """
        with self.assertRaisesRegex(ValueError, "times are required"):
            collect_footprints(self.satellite)
        with self.assertRaisesRegex(ValueError, "elevation is required"):
            collect_footprints(self.radar)
        with self.assertRaisesRegex(ValueError, "no pixels"):
            collect_footprints(self.radar, elevation=0, pixels=True)
        with self.assertRaisesRegex(ValueError, "solar_altaz requires times"):
            collect_footprints(self.tower, solar_altaz=True)
        with self.assertRaisesRegex(TypeError, r"platforms\[1\] must be"):
            collect_footprints([self.tower, "tower"])
        with self.assertRaisesRegex(TypeError, "generate_members"):
            collect_footprints(self.constellation, self.times)
        self.assertTrue(collect_footprints([]).empty)

    def test_collect_footprints_pixels(self):
        """
        Test that the pixels of an aircraft's and a fixed platform's cameras
        are their grids of pixel centers, with angles seen from each pixel.
        """
        pixels = collect_footprints(
            [self.aircraft, self.tower],
            self.times[:2],
            pixels=True,
            platform_altaz=True,
            solar_altaz=True,
        )
        self.assertTrue((pixels.geom_type == "Point").all())
        self.assertEqual(int(pixels["aircraft"].notna().sum()), 2 * 6)
        self.assertEqual(int(pixels["surface_platform"].notna().sum()), 2 * 6)
        self.assertTrue(pixels["valid_obs"].all())
        self.assertTrue(pixels["solar_alt"].notna().all())
        # the pixels of a camera looking down lie within its footprint
        footprint = collect_surface_track(self.tower).geometry.iloc[0]
        tower_pixels = pixels[pixels["surface_platform"] == "Tower"]
        self.assertTrue(tower_pixels.within(footprint.buffer(1e-9)).all())
        # each pixel's elevation angle is below the zenith
        self.assertTrue((tower_pixels["surface_platform_alt"] < 90).all())
        # without times, the pixels of a fixed platform are static and valid
        static = collect_footprints(self.tower, pixels=True)
        self.assertNotIn("time", static.columns)
        self.assertEqual(len(static), 6)
        self.assertTrue(static["valid_obs"].all())

    def test_collect_ground_pixels_vectorized(self):
        """
        Test that the angles of a satellite's pixels are seen from each pixel.
        """
        satellite = self.leo.model_copy(update={"instruments": [_camera("pixels")]})
        pixels = collect_ground_pixels(
            satellite, self.times[:2], sat_altaz=True, solar_altaz=True
        )
        self.assertEqual(len(pixels), 2 * 6)
        # the angles are seen from each pixel, not shared at each time
        self.assertTrue((pixels.groupby("time")["sat_alt"].nunique() > 1).all())
        self.assertTrue(pixels["glint_angle"].notna().all())

    def test_dissolve_footprints(self):
        """
        Test that the valid footprints dissolve into one geometry, or one for
        each group.
        """
        footprints = gpd.GeoDataFrame(
            {
                "aircraft": ["A", "A", "B", None],
                "surface_platform": [None, None, None, "T"],
                "valid_obs": [True, True, True, False],
            },
            geometry=[
                box(0, 0, 1, 1),
                box(1, 0, 2, 1),
                box(5, 5, 6, 6),
                box(9, 9, 10, 10),
            ],
            crs="EPSG:4326",
        )
        dissolved = dissolve_footprints(footprints)
        self.assertListEqual(list(dissolved.columns), ["geometry"])
        self.assertAlmostEqual(dissolved.geometry.iloc[0].area, 3)
        by_aircraft = dissolve_footprints(footprints, by="aircraft")
        self.assertListEqual(list(by_aircraft.columns), ["aircraft", "geometry"])
        self.assertListEqual(list(by_aircraft["aircraft"]), ["A", "B"])
        self.assertAlmostEqual(by_aircraft.geometry.iloc[0].area, 2)
        empty = dissolve_footprints(footprints[~footprints["valid_obs"]], by="aircraft")
        self.assertTrue(empty.empty)
        self.assertListEqual(list(empty.columns), ["aircraft", "geometry"])


class TestFootprintObservations(IssConstellationTestCase):
    """
    Unit tests for observations of points by footprints.
    """

    @staticmethod
    def _footprints(rows, column="aircraft"):
        """
        Gets footprints from (name, instrument, seconds, valid, geometry).
        """
        return gpd.GeoDataFrame(
            {
                "time": [START + timedelta(seconds=row[2]) for row in rows],
                column: [row[0] for row in rows],
                "instrument": [row[1] for row in rows],
                "valid_obs": [row[3] for row in rows],
            },
            geometry=[row[4] for row in rows],
            crs="EPSG:4326",
        )

    def test_runs(self):
        """
        Test that each observation spans a run of consecutive valid
        footprints of an instrument that contain the point.
        """
        inside, outside = box(0, 0, 2, 2), box(5, 5, 6, 6)
        footprints = self._footprints(
            [
                ("A", "camera", 0, True, inside),
                ("A", "camera", 1, True, inside),
                ("A", "camera", 2, True, outside),
                ("A", "camera", 3, True, inside),
                ("A", "camera", 4, False, inside),
                ("A", "camera", 5, True, inside),
                ("A", "lidar", 1, True, inside),
            ]
        )
        observations = collect_footprint_observations(ShapelyPoint(1, 1), footprints)
        self.assertListEqual(
            list(observations.columns),
            [
                "target_hash",
                "geometry",
                "aircraft",
                "instrument",
                "start",
                "epoch",
                "end",
            ],
        )
        seconds = lambda column: list(  # noqa: E731
            (observations[column] - START).dt.total_seconds()
        )
        self.assertListEqual(
            list(observations["instrument"]), ["camera", "lidar", "camera", "camera"]
        )
        self.assertListEqual(seconds("start"), [0, 1, 3, 5])
        self.assertListEqual(seconds("end"), [1, 1, 3, 5])
        self.assertListEqual(seconds("epoch"), [0.5, 1, 3, 5])

    def test_empty_and_static(self):
        """
        Test that no points or footprints give no observations, and that
        static footprints (without times) raise an error.
        """
        footprints = self._footprints([("A", "camera", 0, True, box(0, 0, 2, 2))])
        self.assertTrue(collect_footprint_observations([], footprints).empty)
        empty = collect_footprint_observations(ShapelyPoint(9, 9), footprints)
        self.assertTrue(empty.empty)
        self.assertIn("aircraft", empty.columns)
        with self.assertRaisesRegex(ValueError, "must have a time"):
            collect_footprint_observations(
                ShapelyPoint(1, 1), footprints.drop(columns="time")
            )
        with self.assertRaises(TypeError):
            collect_footprint_observations(box(0, 0, 1, 1), footprints)

    def test_coverage_metrics(self):
        """
        Test that the observations of different kinds of platforms aggregate,
        and coincide, with the platforms named in the columns of their kinds.
        """
        target = box(0, 0, 2, 2)
        footprints = pd.concat(
            [
                self._footprints([("A", "camera", 0, True, target)]),
                self._footprints(
                    [("T", "camera", s, True, target) for s in (0, 1, 2)],
                    column="surface_platform",
                ),
            ],
            ignore_index=True,
        )
        observations = collect_footprint_observations(ShapelyPoint(1, 1), footprints)
        aggregated = aggregate_observations(observations)
        self.assertEqual(len(aggregated), 1)
        self.assertEqual(aggregated["aircraft"].iloc[0], "A")
        self.assertEqual(aggregated["surface_platform"].iloc[0], "T")
        self.assertEqual(aggregated["access"].iloc[0], timedelta(seconds=2))
        coincident = find_coincident_observations(observations)
        self.assertEqual(len(coincident), 1)
        self.assertEqual(coincident["aircraft"].iloc[0], "A")
        self.assertEqual(coincident["surface_platform"].iloc[0], "T")
        self.assertEqual(coincident["instrument"].iloc[0], "camera, camera")

    def test_matches_collect_observations_targets(self):
        """
        Test that a point observed by a satellite's footprints has the target
        hash of its observations by `collect_observations`, so that both
        combine.
        """
        satellite = Satellite(
            name="Test",
            orbit=self.orbit,
            instruments=[Instrument(name="Test", field_of_regard=120.0)],
        )
        times = [T0 + timedelta(minutes=i) for i in range(0, 180)]
        footprints = collect_ground_track(satellite, times)
        point = ShapelyPoint(footprints.geometry.iloc[30].centroid.coords[0])
        sampled = collect_footprint_observations(point, footprints)
        self.assertFalse(sampled.empty)
        observed = collect_observations(point, satellite, times[0], times[-1])
        self.assertSetEqual(set(sampled["target_hash"]), set(observed["target_hash"]))
        # each sampled period lies within an observed period (to a time step)
        for row in sampled.itertuples():
            self.assertTrue(
                (
                    (observed["start"] <= row.start + timedelta(minutes=1))
                    & (observed["end"] >= row.end - timedelta(minutes=1))
                ).any()
            )
        combined = aggregate_observations(pd.concat([sampled, observed]))
        self.assertTrue(np.all(combined["satellite"].str.startswith("Test")))
