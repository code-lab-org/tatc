import unittest

from datetime import datetime, timezone

import numpy as np
from tatc.analysis import compute_dop, DopMethod
from tatc.schemas import (
    Point,
    CircularOrbit,
    WalkerConstellation,
)


class TestDopAnalysis(unittest.TestCase):
    def setUp(self):
        self.null_island = Point(id=0, latitude=0, longitude=0)
        self.orbit = CircularOrbit(
            altitude=20180e3,
            inclination=55,
            epoch=datetime(2024, 1, 1, tzinfo=timezone.utc),
        )
        self.times = [
            datetime(2024, 1, 1, tzinfo=timezone.utc),
            datetime(2024, 1, 2, tzinfo=timezone.utc),
        ]

        self.gps_constellation = WalkerConstellation(
            name="GPS",
            orbit=self.orbit,
            instruments=[],
            number_satellites=24,
            number_planes=6,
        )

    def test_compute_gdop(self):
        results = compute_dop(
            self.times,
            self.null_island,
            self.gps_constellation.generate_members(),
            10,
            DopMethod.GDOP,
        )
        self.assertEqual(len(results), len(self.times))
        self.assertNotIn(np.nan, results.dop.values)

    def test_compute_pdop(self):
        results = compute_dop(
            self.times,
            self.null_island,
            self.gps_constellation.generate_members(),
            10,
            DopMethod.PDOP,
        )
        self.assertEqual(len(results), len(self.times))
        self.assertNotIn(np.nan, results.dop.values)

    def test_compute_hdop(self):
        results = compute_dop(
            self.times,
            self.null_island,
            self.gps_constellation.generate_members(),
            10,
            DopMethod.HDOP,
        )
        self.assertEqual(len(results), len(self.times))
        self.assertNotIn(np.nan, results.dop.values)

    def test_compute_vdop(self):
        results = compute_dop(
            self.times,
            self.null_island,
            self.gps_constellation.generate_members(),
            10,
            DopMethod.VDOP,
        )
        self.assertEqual(len(results), len(self.times))
        self.assertNotIn(np.nan, results.dop.values)

    def test_compute_tdop(self):
        results = compute_dop(
            self.times,
            self.null_island,
            self.gps_constellation.generate_members(),
            10,
            DopMethod.TDOP,
        )
        self.assertEqual(len(results), len(self.times))
        self.assertNotIn(np.nan, results.dop.values)

    def test_compute_gdop_nan(self):
        results = compute_dop(
            self.times,
            self.null_island,
            self.gps_constellation.generate_members(),
            80,
            DopMethod.GDOP,
        )
        self.assertEqual(len(results), len(self.times))
        # minimum elevation angle of 80 is too high to yield enough visible satellites
        self.assertTrue(np.all(np.isnan(results.dop.values)))
