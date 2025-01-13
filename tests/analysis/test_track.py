import unittest

from datetime import datetime, timezone, timedelta
from shapely.geometry import Polygon, MultiPolygon

from tatc.analysis import (
    collect_orbit_track,
    collect_ground_track,
    compute_ground_track,
)
from tatc.schemas import (
    Point,
    GroundStation,
    Satellite,
    Instrument,
    TwoLineElements,
    WalkerConstellation,
)


class TestGroundTrackAnalysis(unittest.TestCase):
    def setUp(self):
        self.point = Point(id=0, latitude=0, longitude=0, min_elevation_angle=10)
        self.station = GroundStation(
            name="Station 1", latitude=0, longitude=180, min_elevation_angle=10
        )
        self.instrument = Instrument(name="Test", field_of_regard=180.0)
        self.orbit = TwoLineElements(
            tle=[
                "1 25544U 98067A   22171.11255782  .00008307  00000+0  15444-3 0  9992",
                "2 25544  51.6448 322.0970 0003980 282.3738 231.6559 15.49798078345636",
            ]
        )
        self.satellite = Satellite(
            name="Test", orbit=self.orbit, instruments=[self.instrument]
        )
        self.constellation = WalkerConstellation(
            name="Test",
            orbit=self.orbit,
            instruments=[self.instrument],
            number_satellites=4,
            number_planes=2,
        )

    def test_collect_orbit_track(self):
        results = collect_orbit_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
        )

    def test_collect_orbit_track_empty(self):
        results = collect_orbit_track(
            self.satellite,
            [],
        )

    def test_collect_orbit_track_with_mask(self):
        mask = Polygon([[-90, 45], [-90, 45], [90, 45], [90, -45], [-90, -45]])
        results = collect_orbit_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
                for i in range(10)
            ],
            mask=mask,
        )

    def test_collect_ground_track(self):
        results = collect_ground_track(
            self.satellite,
            times=[
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
        )

    def test_collect_ground_track_empty(self):
        results = collect_ground_track(
            self.satellite,
            [],
        )

    def test_collect_ground_track_utm(self):
        results = collect_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
            crs="utm",
        )

    def test_collect_ground_track_with_mask(self):
        mask = Polygon([[-90, 45], [-90, 45], [90, 45], [90, -45], [-90, -45]])
        results = collect_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
                for i in range(10)
            ],
            mask=mask,
        )

    def test_compute_ground_track_point(self):
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
            method="point",
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), Polygon)

    def test_compute_ground_track_point_no_instr_index(self):
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
            method="point",
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), Polygon)

    def test_compute_ground_track_point_multipolygon(self):
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, 40, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
            method="point",
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), MultiPolygon)

    def test_compute_ground_track_line_short(self):
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
            method="line",
            crs="EPSG:4087",
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), Polygon)

    def test_compute_ground_track_line_long(self):
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, tzinfo=timezone.utc) + timedelta(minutes=5 * i)
                for i in range(12)
            ],
            method="line",
            crs="EPSG:4087",
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), MultiPolygon)

    def test_compute_ground_track_line_multipolygon(self):
        results = compute_ground_track(
            self.satellite,
            [
                datetime(2022, 6, 1, 1, 40, tzinfo=timezone.utc) + timedelta(minutes=i)
                for i in range(10)
            ],
            method="line",
            crs="EPSG:4087",
        )
        self.assertEqual(len(results.index), 1)
        self.assertEqual(type(results.iloc[0].geometry), MultiPolygon)
