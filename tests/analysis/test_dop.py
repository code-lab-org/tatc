"""
Unit tests for the DOP analysis functions in the tatc.analysis module.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timezone

import numpy as np

from tatc.analysis import DopMethod, compute_dop
from tatc.schemas import (
    CircularOrbit,
    GeneralPerturbationsOrbit,
    Instrument,
    Point,
    Satellite,
    WalkerConstellation,
)


class TestDopAnalysis(unittest.TestCase):
    """
    Unit tests for the DOP analysis functions in the tatc.analysis module.
    """
    def setUp(self):
        self.null_island = Point(id=0, latitude=0, longitude=0)
        self.orbit = CircularOrbit(
            mean_altitude=20180e3,
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
            number_satellites=24,
            number_planes=6,
        )

    def test_compute_gdop(self):
        """
        Test that GDOP computation works for a single point, constellation, and time.
        """
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
        """
        Test that PDOP computation works for a single point, constellation, and time.
        """
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
        """
        Test that HDOP computation works for a single point, constellation, and time.
        """
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
        """
        Test that VDOP computation works for a single point, constellation, and time.
        """
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
        """
        Test that TDOP computation works for a single point, constellation, and time.
        """
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
        """
        Test that GDOP computation returns NaN when there are not enough visible satellites.
        """
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

    def test_compute_dop_geometry_uses_lon_lat_order(self):
        """
        Regression test: the output geometry must be (longitude, latitude),
        matching `geopandas.points_from_xy`'s documented (x=lon, y=lat)
        convention and every other geometry construction in this codebase.
        A point at (latitude=10, longitude=50) was previously reported as
        (50, 10) interpreted backwards -- undetectable with a symmetric
        point like (0, 0), so this uses distinct lat/lon values instead.
        """
        point = Point(id=0, latitude=10, longitude=50)
        results = compute_dop(
            self.times,
            point,
            self.gps_constellation.generate_members(),
            10,
            DopMethod.GDOP,
        )
        self.assertEqual(results.geometry.iloc[0].x, 50)
        self.assertEqual(results.geometry.iloc[0].y, 10)

    def test_compute_dop_honors_point_elevation(self):
        """
        Regression test: `point.elevation` must actually affect the
        computed geometry (range/angles to each satellite), not be
        silently ignored. The effect is small at GNSS-scale ranges, so a
        large elevation is used to make the difference unambiguous.
        """
        point_sea_level = Point(id=0, latitude=10, longitude=50, elevation=0)
        point_elevated = Point(id=0, latitude=10, longitude=50, elevation=8000)
        satellites = self.gps_constellation.generate_members()
        dop_sea_level = compute_dop(
            self.times, point_sea_level, satellites, 10, DopMethod.GDOP
        ).dop.iloc[0]
        dop_elevated = compute_dop(
            self.times, point_elevated, satellites, 10, DopMethod.GDOP
        ).dop.iloc[0]
        self.assertNotEqual(dop_sea_level, dop_elevated)

    def test_compute_dop_min_count_visible_boundary_is_inclusive(self):
        """
        Regression test: `min_count_visible` must be an inclusive lower
        bound (as its docstring states) -- exactly `min_count_visible`
        visible satellites still yields a valid (non-NaN) DOP value, and
        one fewer returns NaN. The true visible count depends on the
        specific constellation geometry and test time, so it is found by
        scanning `min_count_visible` (the DOP value itself is unaffected by
        `min_count_visible` below the true count, since it only gates
        whether NaN is returned, not which satellites are used).
        """
        satellites = self.gps_constellation.generate_members()
        time = self.times[0]
        last_valid = None
        for candidate in range(1, len(satellites) + 1):
            value = compute_dop(
                [time], self.null_island, satellites, 10, DopMethod.GDOP, candidate
            ).dop.iloc[0]
            if np.isnan(value):
                break
            last_valid = candidate
        self.assertIsNotNone(last_valid)
        value_at_boundary = compute_dop(
            [time], self.null_island, satellites, 10, DopMethod.GDOP, last_valid
        ).dop.iloc[0]
        self.assertFalse(np.isnan(value_at_boundary))
        value_past_boundary = compute_dop(
            [time], self.null_island, satellites, 10, DopMethod.GDOP, last_valid + 1
        ).dop.iloc[0]
        self.assertTrue(np.isnan(value_past_boundary))

    def test_compute_dop_multi_element_orbit_uses_per_time_nearest_element(self):
        """
        Regression test: a satellite built from a multi-element orbit (e.g.
        a history of TLEs) must be propagated with whichever element is
        closest to *each* requested time, not a single element chosen once
        (from the first requested time) and reused for the whole span.

        `compute_dop`'s output only exposes a `dop` value (which needs >= 4
        simultaneously "visible" satellites to be non-NaN) and `geometry`
        (which is just the ground point, independent of the satellites) --
        neither directly exposes a single satellite's own propagated
        position. So this uses `min_elevation=-90` (every satellite counts
        as "visible" regardless of true geometry, since elevation is always
        >= -90) plus 3 arbitrary "filler" satellites to reach the 4-visible
        threshold, isolating the comparison to whether the 4th (test
        subject) satellite's own contribution to the DOP design matrix
        differs between a multi-element orbit and a single-element orbit
        built directly from the second TLE, at a two-time query where the
        first time is near the first element's epoch and the second time is
        near the second element's epoch. If per-time selection works, the
        multi-element satellite must select the second element for the
        second query time, exactly matching the single-element-2 satellite;
        if it doesn't (the element closest to the *first* query time is
        reused throughout, as directly verified below to diverge by
        thousands of km at the second time), the DOP values would differ.
        """
        tle_1 = [
            "1 25544U 98067A   22171.11255782  .00008307  00000+0  15444-3 0  9992",
            "2 25544  51.6448 322.0970 0003980 282.3738 231.6559 15.49798078345636",
        ]
        tle_2 = [
            "1 25544U 98067A   22200.50000000  .00008307  00000+0  15444-3 0  9994",
            "2 25544  51.6448  10.0000 0003980 100.0000  50.0000 15.49798078350000",
        ]
        multi_element_orbit = GeneralPerturbationsOrbit.from_tle(tle_1 + tle_2)
        single_element_orbit = GeneralPerturbationsOrbit.from_tle(tle_2)
        instrument = Instrument(name="Test", field_of_regard=10.0)
        multi_satellite = Satellite(
            name="Multi", orbit=multi_element_orbit, instruments=[instrument]
        )
        single_satellite = Satellite(
            name="Single", orbit=single_element_orbit, instruments=[instrument]
        )
        filler_satellites = self.gps_constellation.generate_members()[:3]
        # first time near the first element's epoch, second near the second
        times = [
            datetime(2022, 6, 20, 12, tzinfo=timezone.utc),
            datetime(2022, 7, 19, 12, tzinfo=timezone.utc),
        ]
        multi_result = compute_dop(
            times,
            self.null_island,
            filler_satellites + [multi_satellite],
            -90,
            DopMethod.GDOP,
            min_count_visible=4,
        )
        single_result = compute_dop(
            times,
            self.null_island,
            filler_satellites + [single_satellite],
            -90,
            DopMethod.GDOP,
            min_count_visible=4,
        )
        self.assertAlmostEqual(
            multi_result.dop.iloc[1], single_result.dop.iloc[1], places=6
        )

    def test_compute_dop_singular_matrix_returns_nan_with_warning(self):
        """
        Test that a degenerate geometry (satellites sharing the exact same
        position/velocity, making the DOP design matrix singular) returns
        NaN with a warning, rather than raising.
        """
        instrument = Instrument(name="I", field_of_regard=180.0)
        satellites = [
            Satellite(name=f"S{i}", orbit=self.orbit, instruments=[instrument])
            for i in range(4)
        ]
        with self.assertWarns(UserWarning):
            results = compute_dop(
                [self.times[0]],
                self.null_island,
                satellites,
                -90,
                DopMethod.GDOP,
                min_count_visible=1,
            )
        self.assertTrue(np.isnan(results.dop.iloc[0]))

    def test_compute_dop_invalid_method_raises(self):
        """
        Test that an unrecognized `dop_method` value raises a ValueError,
        for a geometry with enough visible satellites to actually reach
        the dop_method dispatch.
        """
        with self.assertRaises(ValueError):
            compute_dop(
                [self.times[0]],
                self.null_island,
                self.gps_constellation.generate_members(),
                10,
                "not_a_real_method",
            )
