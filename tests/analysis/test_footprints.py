"""
Unit tests for the analysis of the footprints of any platforms.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import warnings
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import LineString, box

from tatc.analysis import (
    aggregate_observations,
    collect_airborne_ground_track,
    collect_footprint_observations,
    collect_footprints,
    collect_ground_pixels,
    collect_ground_track,
    collect_observations,
    collect_radar_track,
    collect_region_observations,
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
from tatc.utils.geometry import hash_geometry

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

    def test_default_elevations(self):
        """
        Test that, without an elevation, the footprints of satellites and
        aircraft are on the ellipsoid and those of a surface platform at its
        ground elevation, by each kind's analysis and together, and that a
        radar station requires one.
        """
        times = self.times[:3]
        tower = self.tower.model_copy(update={"elevation": 250})

        def heights(frame):
            return set(
                np.round(shapely.get_coordinates(frame.geometry, include_z=True)[:, 2])
            )

        self.assertEqual(
            heights(collect_ground_track(self.leo, times, elevation=None)), {0}
        )
        self.assertEqual(
            heights(
                collect_airborne_ground_track(self.aircraft, times, elevation=None)
            ),
            {0},
        )
        self.assertEqual(heights(collect_surface_track(tower)), {250})
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            footprints = collect_footprints([self.aircraft, tower], times)
        self.assertEqual(heights(footprints[footprints.aircraft.notna()]), {0})
        self.assertEqual(
            heights(footprints[footprints.surface_platform.notna()]), {250}
        )
        with self.assertRaisesRegex(ValueError, "radar stations"):
            collect_radar_track(self.radar, None)  # type: ignore

    def test_collect_footprints_mixed(self):
        """
        Test that the footprints of platforms of every kind are combined,
        sorted by time, with the footprints of fixed platforms at each time,
        and static without times.
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
                "radar_station": 10,
            },
        )
        # one set of angle columns for every kind of platform, with no angles
        # for radar stations
        self.assertListEqual(
            [column for column in footprints.columns if column.endswith("_alt")],
            ["platform_alt"],
        )
        angled = footprints["radar_station"].isna() & footprints["valid_obs"]
        self.assertTrue(footprints.loc[angled, "platform_alt"].notna().all())
        self.assertTrue(
            footprints.loc[footprints["radar_station"].notna(), "platform_alt"]
            .isna()
            .all()
        )
        # sorted by time, with a time for every footprint
        self.assertFalse(footprints["time"].isna().any())
        self.assertTrue(footprints["time"].is_monotonic_increasing)
        # no beam of the radar reaches the ground at its own elevation
        radar = footprints[footprints["radar_station"].notna()]
        self.assertFalse(radar["valid_obs"].any())
        # without times, the footprints of fixed platforms are static, after
        # those of moving platforms
        static = collect_footprints(
            [self.radar, self.aircraft, self.tower], elevation=0
        )
        self.assertTrue(static["time"].iloc[-2:].isna().all())
        self.assertListEqual(
            list(static["radar_station"].iloc[-2:].fillna(static["surface_platform"])),
            ["Tower", "KXXX"],
        )
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
        # with times, the radar's coverage samples the points that the
        # satellite observes, at the same times (outside the radar's cone of
        # silence)
        point = ShapelyPoint(center.x + 0.5, center.y, 3000)
        observations = collect_footprint_observations(point, footprints)
        self.assertSetEqual(set(observations["radar_station"].dropna()), {"KXXX"})
        radar_rows = observations[observations["radar_station"].notna()]
        self.assertTrue(radar_rows["instrument"].isna().all())
        coincident = find_coincident_observations(observations)
        self.assertEqual(len(coincident), 1)
        self.assertEqual(coincident["satellite"].iloc[0], "LEO")
        self.assertEqual(coincident["radar_station"].iloc[0], "KXXX")
        self.assertEqual(coincident["instrument"].iloc[0], "wide")
        self.assertLessEqual(coincident["start"].iloc[0], self.times[60])
        self.assertGreaterEqual(coincident["end"].iloc[0], self.times[60])

    def test_collect_footprints_elevation_by_kind(self):
        """
        Test that an elevation may be given for each kind of platform, and
        that footprints of different kinds at different elevations warn.
        """
        with self.assertWarnsRegex(UserWarning, "different elevations"):
            footprints = collect_footprints(
                [self.leo, self.radar], self.times, elevation={"radar_station": 3000}
            )
        heights = footprints.geometry.apply(lambda g: g.exterior.coords[0][2])
        self.assertTrue((heights[footprints.satellite.notna()] == 0).all())
        self.assertTrue((heights[footprints.radar_station.notna()] == 3000).all())
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            same = collect_footprints(
                [self.leo, self.radar], self.times, elevation=3000
            )
        self.assertEqual(len(same), len(footprints))
        with self.assertRaisesRegex(ValueError, "platform columns"):
            collect_footprints(self.radar, elevation={"radar": 3000})
        with self.assertRaisesRegex(ValueError, "elevation is required"):
            collect_footprints(
                [self.leo, self.radar], self.times, elevation={"satellite": 0}
            )

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
        # platforms of a kind must have distinct names (of different kinds,
        # they need not)
        twin = self.tower.model_copy(update={"latitude": 35.1})
        with self.assertRaisesRegex(ValueError, "more than one surface platform"):
            collect_footprints([self.tower, twin])
        with self.assertRaisesRegex(ValueError, "more than one aircraft"):
            collect_airborne_ground_track([self.aircraft, self.aircraft])
        renamed = self.radar.model_copy(update={"name": "Tower"})
        self.assertEqual(
            len(
                collect_footprints(
                    [self.tower, renamed], elevation=0, instrument_index=None
                )
            ),
            3,
        )
        # a tuple of platforms, of any kind, is a list of them
        self.assertEqual(
            len(collect_ground_track((self.leo,), self.times[:2])),
            len(collect_ground_track([self.leo], self.times[:2])),
        )

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
        self.assertTrue((tower_pixels["platform_alt"] < 90).all())
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
            satellite, self.times[:2], platform_altaz=True, solar_altaz=True
        )
        self.assertEqual(len(pixels), 2 * 6)
        # the angles are seen from each pixel, not shared at each time
        self.assertTrue((pixels.groupby("time")["platform_alt"].nunique() > 1).all())
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
    Unit tests for observations of points and regions by footprints.
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

    def test_runs_break_at_gaps(self):
        """
        Test that a run does not span the times without footprints of its
        track (as those a mask culled), unless they are within `max_gap`.
        """
        inside = box(0, 0, 2, 2)
        footprints = self._footprints(
            [("A", "camera", s, True, inside) for s in (0, 1, 2, 100, 101)]
            + [("B", "camera", s, True, inside) for s in (0, 10, 20)]
        )
        observations = collect_footprint_observations(ShapelyPoint(1, 1), footprints)
        seconds = lambda column: list(  # noqa: E731
            (observations[column] - START).dt.total_seconds()
        )
        # the time step of A is 1 second, and of B 10 seconds
        self.assertListEqual(list(observations["aircraft"]), ["A", "B", "A"])
        self.assertListEqual(seconds("start"), [0, 0, 100])
        self.assertListEqual(seconds("end"), [2, 20, 101])
        merged = collect_footprint_observations(
            ShapelyPoint(1, 1), footprints, max_gap=timedelta(seconds=100)
        )
        self.assertEqual(len(merged), 2)
        with self.assertRaisesRegex(ValueError, "non-negative"):
            collect_footprint_observations(
                ShapelyPoint(1, 1), footprints, max_gap=timedelta(seconds=-1)
            )

    def test_elevation_warning(self):
        """
        Test that a point within footprints at another elevation warns.
        """
        footprints = self._footprints(
            [("A", "camera", 0, True, shapely.force_3d(box(0, 0, 2, 2), 500))]
        )
        with self.assertWarnsRegex(UserWarning, "other elevations"):
            collect_footprint_observations(ShapelyPoint(1, 1), footprints)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            observations = collect_footprint_observations(
                ShapelyPoint(1, 1, 500), footprints
            )
        self.assertEqual(len(observations), 1)

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
        with self.assertRaisesRegex(TypeError, "Polygon"):
            collect_footprint_observations(LineString([(0, 0), (1, 1)]), footprints)

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

    def test_regions(self):
        """
        Test that each observation of a region spans a run of consecutive
        valid footprints of an instrument that overlap it in a positive area,
        with the part of the region they cover as its swath.
        """
        region = box(0.5, 0.5, 1.5, 1.5)
        footprints = self._footprints(
            [
                ("A", "camera", 0, True, box(0, 0, 1, 2)),
                ("A", "camera", 1, True, box(1, 0, 2, 2)),
                ("A", "camera", 2, True, box(5, 5, 6, 6)),
                ("A", "camera", 3, True, box(0, 0, 1, 2)),
                # touches the region along an edge: not an observation
                ("A", "camera", 4, True, box(1.5, 0, 3, 2)),
                ("A", "camera", 5, True, box(0, 0, 2, 2)),
            ]
        )
        observations = collect_footprint_observations(region, footprints)
        seconds = lambda column: list(  # noqa: E731
            (observations[column] - START).dt.total_seconds()
        )
        self.assertListEqual(seconds("start"), [0, 3, 5])
        self.assertListEqual(seconds("end"), [1, 3, 5])
        self.assertListEqual(seconds("epoch"), [0.5, 3, 5])
        self.assertTrue(observations.geometry.iloc[0].equals(region))
        self.assertTrue(observations.geometry.iloc[1].equals(box(0.5, 0.5, 1, 1.5)))
        self.assertSetEqual(set(observations["target_hash"]), {hash_geometry(region)})
        # a region at an elevation has its swaths at its elevation, and one at
        # another elevation than the footprints warns
        raised = shapely.force_3d(region, 500)
        footprints["geometry"] = shapely.force_3d(footprints.geometry.to_numpy(), 500)
        with self.assertWarnsRegex(UserWarning, "other elevations"):
            collect_footprint_observations(region, footprints)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            observations = collect_footprint_observations(raised, footprints)
        self.assertEqual(len(observations), 3)
        elevations = shapely.get_coordinates(
            observations.geometry.iloc[0], include_z=True
        )[:, 2]
        self.assertTrue(np.allclose(elevations, 500))
        self.assertSetEqual(set(observations["target_hash"]), {hash_geometry(raised)})

    def test_points_and_regions(self):
        """
        Test that points and regions are observed together, as each alone.
        """
        point, region = ShapelyPoint(1, 1), box(0.5, 0.5, 1.5, 1.5)
        footprints = self._footprints(
            [("A", "camera", s, True, box(0, 0, 2, 2)) for s in (0, 1, 2)]
        )
        observations = collect_footprint_observations(
            gpd.GeoSeries([point, region]), footprints
        )
        self.assertEqual(len(observations), 2)
        for target in (point, region):
            alone = collect_footprint_observations(target, footprints)
            self.assertIn(
                alone["target_hash"].iloc[0], set(observations["target_hash"])
            )
        aggregated = aggregate_observations(observations)
        self.assertEqual(len(aggregated), 2)
        self.assertTrue((aggregated["access"] == timedelta(seconds=2)).all())

    def test_point_angles(self):
        """
        Test that, with the platforms, the observations of points have the
        angles of their platforms seen from the points at their epochs, and
        those of regions have none.
        """
        aircraft = Aircraft(
            name="P-3",
            trajectory=Trajectory.from_dataframe(straight_trajectory(n=120)),
            instruments=[_camera()],
        )
        tower = SurfacePlatform(
            name="Tower",
            latitude=35,
            longitude=-110,
            height=100,
            instruments=[_camera()],
        )
        times = list(pd.date_range(START, periods=120, freq="1s").to_pydatetime())
        footprints = collect_footprints(aircraft, times)
        center = footprints.geometry.iloc[60].centroid
        point = ShapelyPoint(center.x, center.y)
        region = center.buffer(0.001)
        observations = collect_footprint_observations(
            [point, region], footprints, platforms=[aircraft, tower]
        )
        self.assertListEqual(
            list(observations.columns[-2:]), ["platform_alt", "platform_az"]
        )
        points = observations[observations.geometry.geom_type == "Point"]
        self.assertEqual(len(points), 1)
        # the nadir camera's run over the point centers on the time of its
        # footprint, when the aircraft is overhead
        self.assertLess(abs((points["epoch"].iloc[0] - times[60]).total_seconds()), 1)
        self.assertGreater(points["platform_alt"].iloc[0], 85)
        regions = observations[observations.geometry.geom_type != "Point"]
        self.assertEqual(len(regions), 1)
        self.assertTrue(regions[["platform_alt", "platform_az"]].isna().all(axis=None))
        # without the platforms, the observations have no angles
        self.assertNotIn(
            "platform_alt", collect_footprint_observations(point, footprints).columns
        )
        empty = collect_footprint_observations(
            ShapelyPoint(0, 0), footprints, platforms=aircraft
        )
        self.assertTrue(empty.empty)
        self.assertIn("platform_alt", empty.columns)
        with self.assertRaisesRegex(ValueError, "P-3"):
            collect_footprint_observations(point, footprints, platforms=tower)

    def test_matches_collect_region_observations_targets(self):
        """
        Test that a region observed by a satellite's footprints has the target
        hash of its observations by `collect_region_observations`, so that
        both combine.
        """
        satellite = Satellite(
            name="Test",
            orbit=self.orbit,
            instruments=[Instrument(name="Test", field_of_regard=120.0)],
        )
        times = [T0 + timedelta(minutes=i) for i in range(0, 180)]
        footprints = collect_ground_track(satellite, times)
        center = footprints.geometry.iloc[30].centroid
        region = box(center.x - 1, center.y - 1, center.x + 1, center.y + 1)
        sampled = collect_footprint_observations(region, footprints)
        self.assertFalse(sampled.empty)
        observed = collect_region_observations(region, satellite, times[0], times[-1])
        self.assertSetEqual(set(sampled["target_hash"]), set(observed["target_hash"]))
        # each sampled swath lies within the region
        self.assertTrue(
            all(region.buffer(1e-9).contains(swath) for swath in sampled.geometry)
        )
        combined = aggregate_observations(pd.concat([sampled, observed]))
        self.assertEqual(len(combined), 1)
