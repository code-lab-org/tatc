"""
Unit tests for the region coverage analysis functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import Polygon, box

from tatc.analysis import (
    aggregate_observations,
    reduce_observations,
    collect_multi_region_observations,
    collect_observations,
    collect_region_observations,
    compute_region_access_periods,
)
from tatc.schemas import ConicalInstrument, Instrument, PointedInstrument, Satellite
from tatc.utils import hash_geometry, split_polygon

from .common import IssConstellationTestCase


class TestVisiblePolygonIntervalSeries(IssConstellationTestCase):
    """
    Unit tests for `compute_region_access_periods`.
    """

    def setUp(self):
        super().setUp()
        self.narrow = Instrument(name="Narrow", field_of_regard=60)
        self.narrow_satellite = Satellite(
            name="Narrow", orbit=self.orbit, instruments=[self.narrow]
        )
        self.start = datetime(2022, 6, 1, tzinfo=timezone.utc)

    def test_periods_contain_footprint_intersections(self):
        """
        Test that the periods contain every sampled time at which the
        instrument's footprint intersects the region, including a region
        across the anti-meridian and one across a pole.
        """
        times = [self.start + timedelta(seconds=30 * i) for i in range(2880)]
        footprints = self.narrow.compute_footprint(self.orbit.get_orbit_track(times))
        for region in [
            box(-114.8, 31.3, -109.0, 37.0),
            Polygon([(170, -10), (190, -10), (190, 10), (170, 10)]),
            box(-180, -90, 180, -50),
        ]:
            with self.subTest(region=region.wkt):
                periods = compute_region_access_periods(
                    region,
                    self.narrow_satellite,
                    times[0],
                    times[-1],
                    self.narrow.field_of_regard,
                )
                observed = [
                    time
                    for time, footprint in zip(times, footprints)
                    if split_polygon(region).intersects(footprint)
                ]
                self.assertGreater(len(observed), 0)
                for time in observed:
                    self.assertTrue(
                        any(pd.Timestamp(time) in period for period in periods)
                    )

    def test_periods_are_empty_for_unobservable_region(self):
        """
        Test that there are no periods for a region beyond the reach of the
        instrument's field of regard (near the pole, for the ISS orbit).
        """
        periods = compute_region_access_periods(
            box(-180, 80, 180, 90),
            self.narrow_satellite,
            self.start,
            self.start + timedelta(days=1),
            self.narrow.field_of_regard,
        )
        self.assertTrue(periods.empty)

    def test_periods_span_analysis_for_global_region(self):
        """
        Test that a region covering the globe is observable over a single
        period spanning the whole analysis period.
        """
        end = self.start + timedelta(hours=6)
        periods = compute_region_access_periods(
            box(-180, -90, 180, 90),
            self.narrow_satellite,
            self.start,
            end,
            self.narrow.field_of_regard,
        )
        self.assertEqual(len(periods), 1)
        self.assertEqual(periods.iloc[0].left, pd.Timestamp(self.start))
        self.assertEqual(periods.iloc[0].right, pd.Timestamp(end))

    def test_periods_use_region_elevation(self):
        """
        Test that the periods of a region with z coordinates use its
        elevation by default, as when the elevation is given.
        """
        flat = box(-114.8, 31.3, -109.0, 37.0)
        raised = Polygon([(x, y, 3000) for x, y in flat.exterior.coords])
        end = self.start + timedelta(days=1)
        by_z = compute_region_access_periods(
            raised, self.narrow_satellite, self.start, end, self.narrow.field_of_regard
        )
        given = compute_region_access_periods(
            flat,
            self.narrow_satellite,
            self.start,
            end,
            self.narrow.field_of_regard,
            3000,
        )
        self.assertGreater(len(by_z), 0)
        self.assertEqual(list(by_z), list(given))

    def test_point_is_not_a_region(self):
        """
        Test that a point raises a TypeError.
        """
        with self.assertRaises(TypeError):
            compute_region_access_periods(
                ShapelyPoint(0, 0),
                self.satellite,
                self.start,
                self.start + timedelta(hours=1),
            )


class TestCollectRegionObservations(IssConstellationTestCase):
    """
    Unit tests for collecting observations of regions.
    """

    def setUp(self):
        super().setUp()
        self.narrow = Instrument(name="Narrow", field_of_regard=60)
        self.narrow_satellite = Satellite(
            name="Narrow", orbit=self.orbit, instruments=[self.narrow]
        )
        self.start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        self.end = self.start + timedelta(days=1)

    def test_small_region_matches_point(self):
        """
        Test that a region about 100 m across is observed over nearly the
        same periods as the point at its center.
        """
        point = collect_observations(
            ShapelyPoint(-74.03, 40.74), self.narrow_satellite, self.start, self.end
        )
        region = collect_region_observations(
            box(-74.0305, 40.7395, -74.0295, 40.7405),
            self.narrow_satellite,
            self.start,
            self.end,
        )
        self.assertGreater(len(point.index), 0)
        self.assertEqual(len(point.index), len(region.index))
        for column in ["start", "end"]:
            difference = (region[column] - point[column]).dt.total_seconds().abs()
            self.assertLess(difference.max(), 1)

    def test_region_periods_contain_footprint_intersections(self):
        """
        Test that the observation periods of a region contain every sampled
        time at which the instrument's footprint intersects it.
        """
        region = box(-114.8, 31.3, -109.0, 37.0)
        observations = collect_region_observations(
            region, self.narrow_satellite, self.start, self.end
        )
        times = [self.start + timedelta(seconds=30 * i) for i in range(2880)]
        footprints = self.narrow.compute_footprint(self.orbit.get_orbit_track(times))
        observed = [t for t, f in zip(times, footprints) if region.intersects(f)]
        self.assertGreater(len(observed), 0)
        for time in observed:
            self.assertTrue(
                ((observations.start <= time) & (time <= observations.end)).any()
            )

    def test_region_swath_is_split(self):
        """
        Test that observations of a region across the anti-meridian record
        swaths within the region split along the anti-meridian, and the hash
        of the region.
        """
        region = Polygon([(170, -10), (190, -10), (190, 10), (170, 10)])
        observations = collect_region_observations(
            region, self.narrow_satellite, self.start, self.end
        )
        self.assertGreater(len(observations.index), 0)
        self.assertTrue((observations.target_hash == hash_geometry(region)).all())
        bounds = observations.total_bounds
        self.assertGreaterEqual(bounds[0], -180)
        self.assertLessEqual(bounds[2], 180)
        split = split_polygon(region).buffer(1e-9)
        for swath in observations.geometry:
            self.assertGreater(swath.area, 0)
            self.assertTrue(split.contains(swath))

    def test_multi_observations_of_region(self):
        """
        Test that multiple satellite observations of a region report the
        region's hash, the swath, satellite and instrument names, and start and
        end, sorted by start, for each member of a constellation.
        """
        members = self.constellation.generate_members()
        region = box(-114.8, 31.3, -109.0, 37.0)
        observations = collect_multi_region_observations(
            region, members, self.start, self.end
        )
        self.assertGreater(len(observations.index), 0)
        self.assertEqual(
            list(observations.columns),
            ["target_hash", "geometry", "satellite", "instrument", "start", "end"],
        )
        self.assertTrue((observations.target_hash == hash_geometry(region)).all())
        self.assertTrue(observations.start.is_monotonic_increasing)
        self.assertTrue(set(observations.satellite) <= {m.name for m in members})

    def test_observations_of_pointed_and_conical_instruments(self):
        """
        Test that the observation periods of pointed and conical instruments
        contain every sampled time at which their footprint intersects the
        region, and start and end within a sample of the first and last, and
        that their swaths contain the parts of the region within the sampled
        footprints.
        """
        region = box(-114.8, 31.3, -109.0, 37.0)
        times = pd.date_range(self.start, self.end, freq="2s", inclusive="left")
        orbit_track = self.orbit.get_orbit_track(list(times))
        for instrument in [
            PointedInstrument(
                name="Pushbroom",
                field_of_regard=40,
                cross_track_field_of_view=20,
                along_track_field_of_view=0.5,
                roll_angle=8,
                is_rectangular=True,
            ),
            ConicalInstrument(
                name="Conical",
                cone_angle=45,
                along_track_field_of_view=1,
                scan_half_width=65,
            ),
        ]:
            with self.subTest(instrument=instrument.name):
                satellite = Satellite(
                    name="Test", orbit=self.orbit, instruments=[instrument]
                )
                observations = collect_region_observations(
                    region, satellite, self.start, self.end
                )
                footprints = np.array(
                    instrument.compute_footprint(orbit_track), dtype=object
                )
                observed = shapely.intersects(footprints, region)
                self.assertGreater(observed.sum(), 0)
                inside = np.zeros(len(times), dtype=bool)
                for start, end, swath in zip(
                    observations.start, observations.end, observations.geometry
                ):
                    inside |= (times >= start) & (times <= end)
                    sampled = observed & (times >= start) & (times <= end)
                    self.assertGreater(sampled.sum(), 0)
                    self.assertLess((times[sampled][0] - start).total_seconds(), 2)
                    self.assertLess((end - times[sampled][-1]).total_seconds(), 2)
                    viewed = shapely.intersection(
                        shapely.union_all(footprints[sampled]), region
                    )
                    self.assertLess(
                        shapely.difference(viewed, swath).area, 1e-3 * viewed.area
                    )
                self.assertFalse(np.any(observed & ~inside))

    def test_aggregate_and_reduce_observations(self):
        """
        Test that region observations (without a point identifier) are
        aggregated and reduced for each region, by its hash, to the union of
        its swaths.
        """
        regions = [box(-10, 35, -5, 40), box(10, 40, 15, 45)]
        observations = pd.concat(
            [
                collect_region_observations(
                    region, self.narrow_satellite, self.start, self.end
                )
                for region in regions
            ],
            ignore_index=True,
        )
        reduced = reduce_observations(aggregate_observations(observations))
        self.assertEqual(len(reduced.index), 2)
        self.assertNotIn("point_id", reduced.columns)
        self.assertEqual(
            set(reduced.target_hash), {hash_geometry(region) for region in regions}
        )
        for row in reduced.itertuples():
            swaths = observations[observations.target_hash == row.target_hash]
            self.assertEqual(row.samples, len(swaths.index))
            self.assertAlmostEqual(
                row.geometry.area, shapely.union_all(swaths.geometry).area
            )

    def test_region_elevation_from_z(self):
        """
        Test that a region's z coordinates set its elevation: its recorded
        swath has them, and its observed points are at that elevation,
        as for a point at the same elevation.
        """
        flat = box(-74.0305, 40.7395, -74.0295, 40.7405)
        region = Polygon([(x, y, 2000) for x, y in flat.exterior.coords])
        observations = collect_region_observations(
            region, self.narrow_satellite, self.start, self.end
        )
        point = collect_observations(
            ShapelyPoint(-74.03, 40.74, 2000),
            self.narrow_satellite,
            self.start,
            self.end,
        )
        self.assertGreater(len(observations.index), 0)
        self.assertTrue(observations.geometry.iloc[0].has_z)
        self.assertEqual(len(point.index), len(observations.index))
        for column in ["start", "end"]:
            difference = (observations[column] - point[column]).dt.total_seconds().abs()
            self.assertLess(difference.max(), 1)

    def test_point_is_not_a_region(self):
        """
        Test that a point raises a TypeError (see collect_observations).
        """
        with self.assertRaises(TypeError):
            collect_region_observations(
                ShapelyPoint(0, 0), self.satellite, self.start, self.end
            )
