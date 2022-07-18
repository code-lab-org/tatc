import unittest

from datetime import datetime, timezone

from tatc.analysis import (
    collect_observations,
    collect_downlinks,
    compute_latencies,
    reduce_latencies,
)
from tatc.schemas import (
    Point,
    GroundStation,
    Satellite,
    Instrument,
    TwoLineElements,
    WalkerConstellation,
)


class TestLatencyAnalysis(unittest.TestCase):
    def setUp(self):
        self.point = Point(id=0, latitude=0, longitude=0, min_elevation_angle=10)
        self.station = GroundStation(
            name="Station 1", latitude=0, longitude=180, min_elevation_angle=10
        )
        self.stations = [
            self.station,
            GroundStation(
                name="Station 2", latitude=50, longitude=0, min_elevation_angle=10
            ),
            GroundStation(
                name="Station 3", latitude=50, longitude=90, min_elevation_angle=10
            ),
            GroundStation(
                name="Station 4", latitude=50, longitude=-90, min_elevation_angle=10
            ),
        ]
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

    def test_collect_downlinks(self):
        results = collect_downlinks(
            self.station,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 10, tzinfo=timezone.utc),
        )

    def test_collect_downlinks_empty(self):
        results = collect_downlinks(
            self.station,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 1, tzinfo=timezone.utc),
        )
        self.assertTrue(results.empty)

    def test_collect_multi_downlinks(self):
        results = collect_downlinks(
            self.stations,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 10, tzinfo=timezone.utc),
        )

    def test_collect_multi_downlinks_empty(self):
        results = collect_downlinks(
            self.stations,
            self.satellite,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 1, tzinfo=timezone.utc),
        )
        self.assertTrue(results.empty)

    def test_compute_latency(self):
        results = compute_latencies(
            collect_observations(
                self.point,
                self.satellite,
                self.instrument,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
            collect_downlinks(
                self.station,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
        )

    def test_compute_latency_empty(self):
        results = compute_latencies(
            collect_observations(
                self.point,
                self.satellite,
                self.instrument,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
            ),
            collect_downlinks(
                self.station,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
        )
        self.assertTrue(results.empty)

    def test_compute_latency_no_downlinks(self):
        results = compute_latencies(
            collect_observations(
                self.point,
                self.satellite,
                self.instrument,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
            collect_downlinks(
                self.station,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 1, 1, tzinfo=timezone.utc),
            ),
        )
        self.assertTrue(results.dropna().empty)

    def test_compute_latency_multi_station(self):
        results = compute_latencies(
            collect_observations(
                self.point,
                self.satellite,
                self.instrument,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
            collect_downlinks(
                self.stations,
                self.satellite,
                datetime(2022, 6, 1, tzinfo=timezone.utc),
                datetime(2022, 6, 10, tzinfo=timezone.utc),
            ),
        )

    def test_reduce_latency(self):
        results = reduce_latencies(
            compute_latencies(
                collect_observations(
                    self.point,
                    self.satellite,
                    self.instrument,
                    datetime(2022, 6, 1, tzinfo=timezone.utc),
                    datetime(2022, 6, 10, tzinfo=timezone.utc),
                ),
                collect_downlinks(
                    self.stations,
                    self.satellite,
                    datetime(2022, 6, 1, tzinfo=timezone.utc),
                    datetime(2022, 6, 10, tzinfo=timezone.utc),
                ),
            )
        )

    def test_reduce_latency_empty(self):
        results = reduce_latencies(
            compute_latencies(
                collect_observations(
                    self.point,
                    self.satellite,
                    self.instrument,
                    datetime(2022, 6, 1, tzinfo=timezone.utc),
                    datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
                ),
                collect_downlinks(
                    self.stations,
                    self.satellite,
                    datetime(2022, 6, 1, tzinfo=timezone.utc),
                    datetime(2022, 6, 10, tzinfo=timezone.utc),
                ),
            )
        )
        self.assertTrue(results.empty)
