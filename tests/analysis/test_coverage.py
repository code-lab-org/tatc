import unittest

from datetime import datetime, timezone

from tatc.analysis import (
    collect_observations,
    collect_multi_observations,
    aggregate_observations,
    reduce_observations,
)
from tatc.schemas import (
    Point,
    Satellite,
    Instrument,
    TwoLineElements,
    WalkerConstellation,
)


class TestCoverageAnalysis(unittest.TestCase):
    def setUp(self):
        self.point = Point(id=0, latitude=0, longitude=0)
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

    def test_collect_observations(self):
        results = collect_observations(
            self.point,
            self.satellite,
            self.instrument,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
            omit_solar=True,
        )

    def test_collect_observations_with_solar(self):
        results = collect_observations(
            self.point,
            self.satellite,
            self.instrument,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
            omit_solar=False,
        )

    def test_collect_observations_all_culminate(self):
        start = datetime(2022, 6, 1, 0, 43, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 45, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            self.instrument,
            start,
            end,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_miss_first_rise(self):
        start = datetime(2022, 6, 1, 0, 43, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 1, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            self.instrument,
            start,
            end,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].start, start)

    def test_collect_observations_miss_last_set(self):
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 45, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            self.instrument,
            start,
            end,
        )
        self.assertEqual(len(results), 1)
        self.assertEqual(results.iloc[0].end, end)

    def test_collect_observations_null(self):
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point,
            self.satellite,
            self.instrument,
            start,
            end,
        )
        self.assertTrue(results.empty)

    def test_collect_observations_null_with_solar(self):
        start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        end = datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc)
        results = collect_observations(
            self.point, self.satellite, self.instrument, start, end, omit_solar=False
        )
        self.assertTrue(results.empty)

    def test_collect_multi_observations(self):
        results = collect_multi_observations(
            self.point,
            [self.constellation],
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )

    def test_aggregate_observations(self):
        results = collect_multi_observations(
            self.point,
            [self.constellation],
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )
        results = aggregate_observations(results)
        for i in range(len(results.index)):
            self.assertEqual(
                results.iloc[i].end - results.iloc[i].start, results.iloc[i].access
            )
        for i in range(1, len(results.index)):
            self.assertEqual(
                results.iloc[i].start - results.iloc[i-1].end, results.iloc[i].revisit
            )

    def test_aggregate_observations_null(self):
        results = collect_multi_observations(
            self.point,
            [self.constellation],
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
        )
        results = aggregate_observations(results)
        self.assertTrue(results.empty)

    def test_reduce_observations(self):
        results = collect_multi_observations(
            self.point,
            [self.constellation],
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
        )
        results = aggregate_observations(results)
        results = reduce_observations(results)

    def test_reduce_observations_null(self):
        results = collect_multi_observations(
            self.point,
            [self.constellation],
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 1, 0, 30, tzinfo=timezone.utc),
        )
        results = aggregate_observations(results)
        results = reduce_observations(results)
        self.assertTrue(results.empty)
