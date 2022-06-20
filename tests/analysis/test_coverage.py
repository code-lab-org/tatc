import unittest

from datetime import datetime, timezone

from tatc.analysis import (
    collect_observations,
    collect_multi_observations,
    aggregate_observations,
    reduce_observations,
)
from tatc.schemas import Point, Satellite, Instrument, TwoLineElements


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

    def test_collect_observations_omit_solar(self):
        collect_observations(
            self.point,
            self.satellite,
            self.instrument,
            datetime(2022, 6, 1, tzinfo=timezone.utc),
            datetime(2022, 6, 2, tzinfo=timezone.utc),
            omit_solar=True,
        )

    def test_collect_observations_no_omit_solar(self):
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
