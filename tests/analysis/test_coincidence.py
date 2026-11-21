"""
Unit tests for the coincident observation analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

import geopandas as gpd
import pandas as pd
from shapely.geometry import Point as ShapelyPoint, box

from tatc.analysis import (
    aggregate_observations,
    collect_ground_track,
    collect_observations,
    collect_surface_track,
    find_coincident_ground_track,
    find_coincident_observations,
)
from tatc.schemas import Instrument, Satellite, SurfacePlatform

from .common import IssConstellationTestCase

T0 = datetime(2022, 6, 1, tzinfo=timezone.utc)


def _minutes(value):
    return T0 + timedelta(minutes=value)


def _make_observation(satellite, start, end, instrument="Test", target_hash="0"):
    """
    Build a synthetic observation record matching the schema
    `collect_observations` produces, with start and end in minutes after T0.
    """
    return {
        "target_hash": target_hash,
        "geometry": ShapelyPoint(0, 0),
        "satellite": satellite,
        "instrument": instrument,
        "start": pd.Timestamp(_minutes(start)),
        "epoch": pd.Timestamp(_minutes((start + end) / 2)),
        "end": pd.Timestamp(_minutes(end)),
    }


def _make_footprint(satellite, time, geometry, instrument="Test"):
    """
    Build a synthetic footprint record matching the schema
    `collect_ground_track` produces, with time in minutes after T0.
    """
    return {
        "time": pd.Timestamp(_minutes(time)),
        "satellite": satellite,
        "instrument": instrument,
        "valid_obs": True,
        "geometry": geometry,
    }


class TestFindCoincidentObservations(IssConstellationTestCase):
    """
    Unit tests for `find_coincident_observations`.
    """

    def test_overlap(self):
        """
        Test that the coincident period of two satellites is the overlap of
        their observation periods.
        """
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 6, 15)],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite, "A, B")
        self.assertEqual(results.iloc[0].instrument, "Test, Test")
        self.assertEqual(results.iloc[0].start, _minutes(6))
        self.assertEqual(results.iloc[0].end, _minutes(10))
        self.assertEqual(results.iloc[0].epoch, _minutes(8))

    def test_no_overlap(self):
        """
        Test that separate observation periods are not coincident, giving an
        empty data frame with the result columns.
        """
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 12, 15)],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations)
        self.assertTrue(results.empty)
        self.assertListEqual(
            list(results.columns),
            [
                "target_hash",
                "geometry",
                "satellite",
                "instrument",
                "start",
                "epoch",
                "end",
            ],
        )

    def test_max_time_difference(self):
        """
        Test that observations coincide if the gap between their periods is
        at most the maximum time difference, widening each period by half.
        """
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 12, 15)],
            crs="EPSG:4326",
        )
        self.assertTrue(
            find_coincident_observations(
                observations, max_time_difference=timedelta(minutes=1)
            ).empty
        )
        results = find_coincident_observations(
            observations, max_time_difference=timedelta(minutes=4)
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].start, _minutes(10))
        self.assertEqual(results.iloc[0].end, _minutes(12))

    def test_instantaneous(self):
        """
        Test that instantaneous observations coincide if their epochs differ
        by at most the maximum time difference, regardless of their periods,
        and that each coincident observation spans the contributing epochs.
        """
        # periods overlap (minutes 8 to 10), epochs (minutes 5 and 19) do not
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 8, 30)],
            crs="EPSG:4326",
        )
        self.assertEqual(len(find_coincident_observations(observations).index), 1)
        self.assertTrue(
            find_coincident_observations(
                observations,
                max_time_difference=timedelta(minutes=10),
                instantaneous=True,
            ).empty
        )
        results = find_coincident_observations(
            observations, max_time_difference=timedelta(minutes=15), instantaneous=True
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite, "A, B")
        self.assertEqual(results.iloc[0].start, _minutes(5))
        self.assertEqual(results.iloc[0].end, _minutes(19))
        self.assertEqual(results.iloc[0].epoch, _minutes(12))

    def test_same_satellite(self):
        """
        Test that observations by different instruments on the same
        satellite are not coincident with each other.
        """
        observations = gpd.GeoDataFrame(
            [
                _make_observation("A", 0, 10, instrument="X"),
                _make_observation("A", 5, 15, instrument="Y"),
            ],
            crs="EPSG:4326",
        )
        self.assertTrue(find_coincident_observations(observations).empty)

    def test_maximal_period(self):
        """
        Test that a coincident period spans the time that at least two
        satellites observe, even as the satellites observing change, and
        records each satellite and instrument once.
        """
        observations = gpd.GeoDataFrame(
            [
                _make_observation("A", 0, 10),
                _make_observation("B", 2, 4, instrument="X"),
                _make_observation("B", 2, 4, instrument="X"),
                _make_observation("C", 3, 12),
                _make_observation("D", 20, 25),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite, "A, B, C")
        self.assertEqual(results.iloc[0].instrument, "Test, X, Test")
        self.assertEqual(results.iloc[0].start, _minutes(2))
        self.assertEqual(results.iloc[0].end, _minutes(10))

    def test_min_platforms(self):
        """
        Test that a minimum of three satellites requires all three to
        observe at once.
        """
        observations = gpd.GeoDataFrame(
            [
                _make_observation("A", 0, 10),
                _make_observation("B", 2, 8),
                _make_observation("C", 6, 12),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations, min_platforms=3)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite, "A, B, C")
        self.assertEqual(results.iloc[0].start, _minutes(6))
        self.assertEqual(results.iloc[0].end, _minutes(8))

    def test_isolates_targets(self):
        """
        Test that observations of different targets are not coincident.
        """
        observations = gpd.GeoDataFrame(
            [
                _make_observation("A", 0, 10, target_hash="0"),
                _make_observation("B", 5, 15, target_hash="1"),
                _make_observation("C", 5, 15, target_hash="1"),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_observations(observations)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].target_hash, "1")
        self.assertEqual(results.iloc[0].satellite, "B, C")

    def test_invalid_arguments(self):
        """
        Test that invalid arguments and region observations raise errors.
        """
        observations = gpd.GeoDataFrame(
            [_make_observation("A", 0, 10), _make_observation("B", 6, 15)],
            crs="EPSG:4326",
        )
        with self.assertRaises(ValueError):
            find_coincident_observations(observations, min_platforms=1)
        with self.assertRaises(ValueError):
            find_coincident_observations(
                observations, max_time_difference=timedelta(minutes=-1)
            )
        observations.geometry = [box(0, 0, 1, 1), box(0, 0, 1, 1)]
        with self.assertRaises(ValueError):
            find_coincident_observations(observations)

    def test_collect_observations(self):
        """
        Test that the coincident observations of a constellation lie within
        observations by at least two of its satellites, and aggregate.
        """
        satellites = self.constellation.generate_members()
        observations = collect_observations(
            ShapelyPoint(0, 45), satellites, T0, T0 + timedelta(days=1)
        )
        # successive passes of the constellation are about 40 minutes apart
        max_time_difference = timedelta(minutes=45)
        results = find_coincident_observations(
            observations, max_time_difference=max_time_difference
        )
        self.assertFalse(results.empty)
        half = max_time_difference / 2
        for result in results.itertuples():
            contributors = observations[
                (observations.start - half <= result.end)
                & (observations.end + half >= result.start)
            ]
            self.assertGreaterEqual(contributors.satellite.nunique(), 2)
            self.assertSetEqual(
                set(result.satellite.split(", ")), set(contributors.satellite)
            )
        self.assertEqual(len(aggregate_observations(results).index), len(results.index))


class TestFindCoincidentGroundTrack(IssConstellationTestCase):
    """
    Unit tests for `find_coincident_ground_track`.
    """

    def test_overlap(self):
        """
        Test that the overlapping footprints of two satellites at the same
        time are coincident, with their overlap as the geometry.
        """
        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("A", 0, box(0, 0, 2, 2)),
                _make_footprint("B", 0, box(1, 1, 3, 3)),
                _make_footprint("B", 1, box(1, 1, 3, 3)),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_ground_track(ground_track)
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite_1, "A")
        self.assertEqual(results.iloc[0].satellite_2, "B")
        self.assertEqual(results.iloc[0].time_1, _minutes(0))
        self.assertEqual(results.iloc[0].time_2, _minutes(0))
        self.assertEqual(results.iloc[0].time_difference, timedelta(0))
        self.assertTrue(results.iloc[0].valid_obs_1 and results.iloc[0].valid_obs_2)
        self.assertTrue(results.iloc[0].geometry.equals(box(1, 1, 2, 2)))
        self.assertEqual(results.crs, ground_track.crs)

    def test_max_time_difference(self):
        """
        Test that footprints coincide if their times differ by at most the
        maximum time difference, ordered by time.
        """
        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("B", 0, box(1, 1, 3, 3)),
                _make_footprint("A", 5, box(0, 0, 2, 2)),
            ],
            crs="EPSG:4326",
        )
        self.assertTrue(find_coincident_ground_track(ground_track).empty)
        results = find_coincident_ground_track(
            ground_track, max_time_difference=timedelta(minutes=5)
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite_1, "B")
        self.assertEqual(results.iloc[0].satellite_2, "A")
        self.assertEqual(results.iloc[0].time_difference, timedelta(minutes=5))

    def test_same_satellite_and_touching(self):
        """
        Test that footprints of the same satellite, or that only touch,
        are not coincident.
        """
        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("A", 0, box(0, 0, 2, 2), instrument="X"),
                _make_footprint("A", 0, box(1, 1, 3, 3), instrument="Y"),
                _make_footprint("B", 0, box(3, 3, 4, 4)),
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_ground_track(ground_track)
        self.assertTrue(results.empty)
        self.assertIn("time_difference", results.columns)
        self.assertIn("satellite_2", results.columns)

    def test_aircraft(self):
        """
        Test that footprints of a satellite and an aircraft (an underflight)
        are coincident, even if the aircraft has the same name, and that a
        footprint must name its platform.
        """
        aircraft = _make_footprint(None, 1, box(1, 1, 3, 3))
        aircraft["aircraft"] = aircraft.pop("satellite")
        aircraft["aircraft"] = "A"
        ground_track = pd.concat(
            [
                gpd.GeoDataFrame(
                    [_make_footprint("A", 0, box(0, 0, 2, 2))], crs="EPSG:4326"
                ),
                gpd.GeoDataFrame([aircraft], crs="EPSG:4326"),
            ],
            ignore_index=True,
        )
        results = find_coincident_ground_track(
            ground_track, max_time_difference=timedelta(minutes=1)
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(results.iloc[0].satellite_1, "A")
        self.assertTrue(pd.isna(results.iloc[0].aircraft_1))
        self.assertEqual(results.iloc[0].aircraft_2, "A")
        self.assertTrue(pd.isna(results.iloc[0].satellite_2))
        with self.assertRaisesRegex(
            ValueError, "satellite, aircraft, surface_platform, or radar_station"
        ):
            find_coincident_ground_track(ground_track.drop(columns="aircraft"))

    def test_static_footprints(self):
        """
        Test that the footprints of fixed platforms (without a time) are
        coincident with overlapping footprints at any time, and with each
        other, after the footprints with times.
        """
        tower = {
            "surface_platform": "Tower",
            "instrument": "Camera",
            "valid_obs": True,
            "geometry": box(1, 1, 3, 3),
        }
        balloon = {**tower, "surface_platform": "Balloon", "geometry": box(2, 2, 4, 4)}
        ground_track = pd.concat(
            [
                gpd.GeoDataFrame(
                    [
                        _make_footprint("A", 0, box(0, 0, 2, 2)),
                        _make_footprint("B", 60, box(0, 0, 2, 2)),
                        _make_footprint("C", 0, box(5, 5, 6, 6)),
                    ],
                    crs="EPSG:4326",
                ),
                gpd.GeoDataFrame([tower, balloon], crs="EPSG:4326"),
            ],
            ignore_index=True,
        )
        results = find_coincident_ground_track(ground_track)
        # A and B, an hour apart, are each coincident with the tower, but not
        # with each other
        self.assertListEqual(
            results["satellite_1"].fillna(results["surface_platform_1"]).tolist(),
            ["A", "B", "Tower"],
        )
        self.assertListEqual(
            results["satellite_2"].fillna(results["surface_platform_2"]).tolist(),
            ["Tower", "Tower", "Balloon"],
        )
        self.assertTrue(results["time_difference"].isna().all())
        self.assertTrue(results["time_2"].isna().all())
        self.assertTrue(results.geometry.iloc[0].equals(box(1, 1, 2, 2)))
        self.assertTrue(results.geometry.iloc[-1].equals(box(2, 2, 3, 3)))
        # only static footprints, without a time column
        static = find_coincident_ground_track(
            gpd.GeoDataFrame([tower, balloon], crs="EPSG:4326")
        )
        self.assertEqual(len(static.index), 1)
        self.assertTrue(static["time_difference"].isna().all())

    def test_collect_surface_track(self):
        """
        Test that the overpasses of a satellite over the view of a fixed
        platform are coincident with it.
        """
        satellite = Satellite(
            name="Satellite",
            orbit=self.orbit,
            instruments=[Instrument(name="Test", field_of_regard=60.0)],
        )
        times = [T0 + timedelta(minutes=i) for i in range(0, 100)]
        ground_track = collect_ground_track(satellite, times)
        center = ground_track.geometry.iloc[50].centroid
        tower = SurfacePlatform(
            name="Tower",
            latitude=center.y,
            longitude=center.x,
            height=1000,
            instruments=[Instrument(name="Camera", field_of_regard=120.0)],
        )
        results = find_coincident_ground_track(
            pd.concat([ground_track, collect_surface_track(tower)], ignore_index=True)
        )
        self.assertFalse(results.empty)
        self.assertIn(times[50], set(results["time_1"]))
        self.assertTrue((results["surface_platform_2"] == "Tower").all())
        self.assertTrue(results["time_difference"].isna().all())

    def test_empty(self):
        """
        Test that an empty ground track has no coincident footprints.
        """
        ground_track = collect_ground_track([], [])
        results = find_coincident_ground_track(ground_track)
        self.assertTrue(results.empty)
        with self.assertRaises(ValueError):
            find_coincident_ground_track(
                ground_track, max_time_difference=timedelta(minutes=-1)
            )

    def test_collect_ground_track(self):
        """
        Test that the coincident footprints of a constellation lie within
        the footprints of both satellites.
        """
        instrument = Instrument(name="Test", field_of_regard=60.0)
        satellites = [
            Satellite(name=member.name, orbit=member.orbit, instruments=[instrument])
            for member in self.constellation.generate_members()
        ]
        times = [T0 + timedelta(minutes=i) for i in range(0, 24 * 60, 2)]
        ground_track = collect_ground_track(satellites, times)
        results = find_coincident_ground_track(
            ground_track, max_time_difference=timedelta(minutes=10)
        )
        self.assertFalse(results.empty)
        self.assertTrue((results.satellite_1 != results.satellite_2).all())
        self.assertTrue((results.time_difference <= timedelta(minutes=10)).all())
        self.assertTrue((results.time_1 <= results.time_2).all())
        footprints = ground_track.set_index(["time", "satellite"]).geometry
        for result in results.itertuples():
            for time, satellite in (
                (result.time_1, result.satellite_1),
                (result.time_2, result.satellite_2),
            ):
                footprint = footprints.loc[(time, satellite)]
                self.assertTrue(footprint.buffer(1e-9).contains(result.geometry))

    def test_between(self):
        """
        Test that `between` keeps only the pairs of two kinds of platforms,
        with the first kind first, whatever its time.
        """
        aircraft = _make_footprint(None, 0, box(1, 1, 3, 3))
        aircraft["aircraft"] = "P-3"
        ground_track = gpd.GeoDataFrame(
            [
                aircraft,
                _make_footprint("A", 1, box(0, 0, 2, 2)),
                _make_footprint("B", 1, box(0, 0, 2, 2)),
            ],
            crs="EPSG:4326",
        )
        difference = timedelta(minutes=1)
        self.assertEqual(len(find_coincident_ground_track(ground_track, difference)), 3)
        results = find_coincident_ground_track(
            ground_track, difference, between=("satellite", "aircraft")
        )
        self.assertListEqual(list(results.satellite_1), ["A", "B"])
        self.assertListEqual(list(results.aircraft_2), ["P-3", "P-3"])
        self.assertTrue((results.time_difference == -difference).all())
        reverse = find_coincident_ground_track(
            ground_track, difference, between=("aircraft", "satellite")
        )
        self.assertListEqual(list(reverse.satellite_2), ["A", "B"])
        self.assertTrue((reverse.time_difference == difference).all())
        same = find_coincident_ground_track(
            ground_track, difference, between=("satellite", "satellite")
        )
        self.assertListEqual(
            list(zip(same.satellite_1, same.satellite_2)), [("A", "B")]
        )
        self.assertTrue(
            find_coincident_ground_track(
                ground_track, difference, between=("satellite", "radar_station")
            ).empty
        )
        for between in ("satellite", ("satellite",), ("satellite", "station")):
            with self.subTest(between=between):
                with self.assertRaisesRegex(ValueError, "between"):
                    find_coincident_ground_track(ground_track, between=between)

    def test_many_platforms(self):
        """
        Test that the footprints of many platforms pair as those of each pair
        of platforms do.
        """
        footprints = [
            _make_footprint(name, time, box(x, 0, x + 2, 2))
            for name, x in (("A", 0), ("B", 1), ("C", 1.5), ("D", 10))
            for time in (0, 1)
        ]
        ground_track = gpd.GeoDataFrame(footprints, crs="EPSG:4326")
        results = find_coincident_ground_track(ground_track)
        self.assertListEqual(
            sorted(zip(results.satellite_1, results.satellite_2)),
            sorted([("A", "B"), ("A", "C"), ("B", "C")] * 2),
        )
        self.assertTrue((results.time_difference == timedelta(0)).all())

    def test_fixed_platform_copies(self):
        """
        Test that a footprint pairs only with the copy nearest in time of
        each instrument of a fixed platform collected at times (its footprint
        is the same at every time), the earlier of two as near, and each
        footprint of one fixed platform with the nearest copy of another.
        """

        def radar(time, name="KIWA", geometry=box(1, 1, 3, 3)):
            footprint = _make_footprint(None, time, geometry)
            del footprint["satellite"], footprint["instrument"]
            return {**footprint, "radar_station": name}

        def tower(time, instrument):
            footprint = _make_footprint(None, time, box(0, 0, 2, 2), instrument)
            del footprint["satellite"]
            return {**footprint, "surface_platform": "Tower"}

        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("A", 5, box(0, 0, 2, 2)),
                _make_footprint("A", 6.5, box(0, 0, 2, 2)),
                *[radar(time) for time in range(10)],
                *[tower(time, instrument) for time in (0, 11) for instrument in "XY"],
            ],
            crs="EPSG:4326",
        )
        difference = timedelta(minutes=10)
        results = find_coincident_ground_track(
            ground_track, difference, between=("satellite", "radar_station")
        )
        # A's footprint at 6.5 minutes is as near to the copies at 6 and 7
        # minutes: the earlier is kept
        self.assertListEqual(list(results.time_1), [_minutes(5), _minutes(6.5)])
        self.assertListEqual(list(results.time_2), [_minutes(5), _minutes(6)])
        towers = find_coincident_ground_track(
            ground_track, difference, between=("satellite", "surface_platform")
        )
        self.assertListEqual(
            list(zip(towers.time_1, towers.time_2, towers.instrument_2)),
            [
                (_minutes(5), _minutes(0), "X"),
                (_minutes(5), _minutes(0), "Y"),
                (_minutes(6.5), _minutes(11), "X"),
                (_minutes(6.5), _minutes(11), "Y"),
            ],
        )
        fixed = find_coincident_ground_track(
            ground_track, difference, between=("radar_station", "surface_platform")
        )
        # each radar footprint with the nearest copy of each tower instrument
        self.assertEqual(len(fixed), 20)
        self.assertTrue((abs(fixed.time_difference) <= timedelta(minutes=5)).all())

    def test_fixed_platform_copies_of_other_geometries(self):
        """
        Test that a footprint pairs with the nearest copy of a fixed
        platform's instrument among copies of different geometries (as of
        footprints collected at other elevations), and with its static
        footprints at any time.
        """

        def tower(time, geometry, name="Tower"):
            footprint = _make_footprint(None, 0 if time is None else time, geometry)
            del footprint["satellite"]
            if time is None:
                footprint["time"] = pd.NaT
            return {**footprint, "surface_platform": name}

        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("A", 5, box(0, 0, 2, 2)),
                *[tower(time, box(0, 0, 2, 2)) for time in (0, 1, 2)],
                *[tower(time, box(1, 1, 3, 3)) for time in (4, 8)],
                tower(None, box(1, 0, 3, 2), name="Static"),
            ],
            crs="EPSG:4326",
        )
        ground_track["time"] = pd.to_datetime(ground_track["time"], utc=True)
        results = find_coincident_ground_track(
            ground_track,
            timedelta(minutes=10),
            between=("satellite", "surface_platform"),
        )
        self.assertListEqual(list(results.surface_platform_2), ["Tower", "Static"])
        # the copy at 4 minutes is nearer than those at 0 to 2 minutes
        self.assertEqual(results.time_2.iloc[0], pd.Timestamp(_minutes(4)))
        self.assertTrue(pd.isna(results.time_2.iloc[1]))
        # beyond the maximum time difference of every copy, none
        none = find_coincident_ground_track(
            ground_track.iloc[[0, 1, 2, 3]], timedelta(minutes=2)
        )
        self.assertTrue(none.empty)

    def test_valid_pairs(self):
        """
        Test that a pair is a valid observation if both footprints are, so
        that dissolving the pairs gives the area validly observed in common.
        """
        from tatc.analysis import dissolve_footprints

        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("A", 0, box(0, 0, 2, 2)),
                _make_footprint("B", 0, box(1, 1, 3, 3)),
                {**_make_footprint("C", 0, box(0, 0, 1.5, 1.5)), "valid_obs": False},
            ],
            crs="EPSG:4326",
        )
        results = find_coincident_ground_track(ground_track)
        self.assertListEqual(
            list(zip(results.satellite_1, results.satellite_2, results.valid_obs)),
            [("A", "B", True), ("A", "C", False), ("B", "C", False)],
        )
        common = dissolve_footprints(results)
        self.assertTrue(common.geometry.iloc[0].equals(box(1, 1, 2, 2)))
        self.assertNotIn(
            "valid_obs",
            find_coincident_ground_track(ground_track.drop(columns="valid_obs")),
        )

    def test_bounds_at_other_times(self):
        """
        Test that footprints whose bounds overlap but whose times differ by
        more than the maximum time difference are not coincident, and that
        footprints whose bounds overlap without overlapping are not either.
        """
        ground_track = gpd.GeoDataFrame(
            [
                _make_footprint("A", 0, box(0, 0, 2, 2)),
                _make_footprint("B", 30, box(0, 0, 2, 2)),
                _make_footprint(
                    "B", 0, box(0, 0, 3, 3).difference(box(0, 0, 2.5, 2.5))
                ),
            ],
            crs="EPSG:4326",
        )
        self.assertTrue(find_coincident_ground_track(ground_track).empty)
        self.assertEqual(
            len(find_coincident_ground_track(ground_track, timedelta(minutes=30))), 1
        )
