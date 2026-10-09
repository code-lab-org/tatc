"""
Unit tests for the input validation shared by tatc.analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import warnings
from datetime import datetime, time, timedelta, timezone

from tatc.analysis import (
    DopMethod,
    collect_downlinks,
    collect_ground_pixels,
    collect_ground_track,
    collect_limb_observations,
    collect_multi_observations,
    collect_observations,
    collect_orbit_track,
    collect_ro_observations,
    compute_dop,
    collect_region_observations,
    compute_access_periods,
    compute_ground_track,
)
from shapely.geometry import Polygon

from tatc.analysis.check import (
    _check_instrument_index,
    _check_satellite,
    _check_satellites,
    _check_time_window,
)
from tatc.schemas import (
    GroundStation,
    Instrument,
    Point,
    Satellite,
    SunSynchronousOrbit,
)

from .common import IssConstellationTestCase


class TestValidation(IssConstellationTestCase):
    """
    Unit tests for the input validation shared by tatc.analysis functions.
    """

    def setUp(self):
        super().setUp()
        self.point = Point(id=0, latitude=0, longitude=0)
        self.start = datetime(2022, 6, 1, tzinfo=timezone.utc)
        self.end = self.start + timedelta(hours=1)
        self.times = [self.start, self.end]

    def test_check_satellite(self):
        """
        Test that a satellite passes and a constellation raises a TypeError
        directing to its generate_members method.
        """
        self.assertIs(_check_satellite(self.satellite), self.satellite)
        with self.assertRaisesRegex(TypeError, "generate_members"):
            _check_satellite(self.constellation)
        with self.assertRaises(TypeError):
            _check_satellite(self.point)

    def test_check_satellites(self):
        """
        Test that a satellite or list of satellites is returned as a list and
        that constellations (alone or as list members) raise a TypeError.
        """
        self.assertEqual(_check_satellites(self.satellite), [self.satellite])
        members = self.constellation.generate_members()
        self.assertEqual(_check_satellites(members), members)
        with self.assertRaisesRegex(TypeError, "generate_members"):
            _check_satellites(self.constellation)
        with self.assertRaisesRegex(TypeError, r"satellites\[1\]"):
            _check_satellites([self.satellite, self.constellation])
        with self.assertRaisesRegex(TypeError, "list of Satellites"):
            _check_satellites(self.satellite, allow_single=False)
        with self.assertRaisesRegex(TypeError, "generate_members"):
            _check_satellites(self.constellation, allow_single=False)

    def test_analysis_functions_reject_constellations(self):
        """
        Test that analysis functions limited to satellites reject a
        constellation rather than silently analyzing its lead member.
        """
        calls = {
            "collect_observations": lambda c: collect_observations(
                self.point, c, self.start, self.end
            ),
            "collect_multi_observations": lambda c: collect_multi_observations(
                self.point, c, self.start, self.end
            ),
            "collect_downlinks": lambda c: collect_downlinks(
                GroundStation(name="Test", latitude=0, longitude=0),
                c,
                self.start,
                self.end,
            ),
            "collect_limb_observations": lambda c: collect_limb_observations(
                c, self.times, 0, [10e3, 50e3], timedelta(seconds=10)
            ),
            "collect_ro_observations (receiver)": lambda c: collect_ro_observations(
                c, self.satellite, self.start, self.end
            ),
            "collect_ro_observations (transmitters)": lambda c: collect_ro_observations(
                self.satellite, c, self.start, self.end
            ),
            "collect_orbit_track": lambda c: collect_orbit_track(c, self.times),
            "collect_ground_track": lambda c: collect_ground_track(c, self.times),
            "compute_ground_track": lambda c: compute_ground_track(c, self.times),
            "collect_ground_pixels": lambda c: collect_ground_pixels(c, self.times),
            "compute_dop": lambda c: compute_dop(
                self.times, self.point, c, 0, DopMethod.PDOP
            ),
        }
        for name, call in calls.items():
            with self.subTest(name):
                with self.assertRaisesRegex(TypeError, "generate_members"):
                    call(self.constellation)

    def test_check_time_window(self):
        """
        Test that a time window requires timezone-aware datetimes and an end
        no earlier than its start (an empty window is allowed).
        """
        _check_time_window(self.start, self.end)
        _check_time_window(self.start, self.start)
        naive = self.start.replace(tzinfo=None)
        with self.assertRaisesRegex(ValueError, "start must be a timezone-aware"):
            _check_time_window(naive, self.end)
        with self.assertRaisesRegex(ValueError, "end must be a timezone-aware"):
            _check_time_window(self.start, self.end.replace(tzinfo=None))
        with self.assertRaisesRegex(ValueError, "is before start"):
            _check_time_window(self.end, self.start)

    def test_analysis_functions_check_time_windows(self):
        """
        Test that analysis functions over a time window reject a reversed
        window rather than silently returning an empty result.
        """
        region = Polygon([(-10, -10), (10, -10), (10, 10), (-10, 10)])
        station = GroundStation(name="Test", latitude=0, longitude=0)
        calls = {
            "collect_observations": lambda s, e: collect_observations(
                self.point, self.satellite, s, e
            ),
            "compute_access_periods": lambda s, e: compute_access_periods(
                self.point, self.satellite, s, e
            ),
            "collect_region_observations": lambda s, e: collect_region_observations(
                region, self.satellite, s, e
            ),
            "collect_downlinks": lambda s, e: collect_downlinks(
                station, self.satellite, s, e
            ),
            "collect_ro_observations": lambda s, e: collect_ro_observations(
                self.satellite, self.satellite, s, e
            ),
        }
        for name, call in calls.items():
            with self.subTest(name):
                with self.assertRaisesRegex(ValueError, "is before start"):
                    call(self.end, self.start)
                with self.assertRaisesRegex(ValueError, "timezone-aware"):
                    call(self.start.replace(tzinfo=None), self.end)

    def test_check_instrument_index(self):
        """
        Test that an instrument index out of range raises an IndexError
        naming the satellite (negative indices count from the end).
        """
        self.assertEqual(_check_instrument_index(self.satellite, 0), 0)
        self.assertEqual(_check_instrument_index(self.satellite, -1), -1)
        for index in (1, -2):
            with self.subTest(index=index):
                with self.assertRaisesRegex(
                    IndexError, "'Test', which has 1 instrument$"
                ):
                    _check_instrument_index(self.satellite, index)
        calls = {
            "collect_observations": lambda: collect_observations(
                self.point, self.satellite, self.start, self.end, instrument_index=1
            ),
            "collect_orbit_track": lambda: collect_orbit_track(
                self.satellite, self.times, instrument_index=1
            ),
            "collect_ground_track": lambda: collect_ground_track(
                self.satellite, self.times, instrument_index=1
            ),
        }
        for name, call in calls.items():
            with self.subTest(name):
                with self.assertRaisesRegex(IndexError, "out of range for satellite"):
                    call()

    def test_warn_low_perigees(self):
        """
        Test that an analysis warns once of satellites with a perigee below
        100 km, as from an altitude mistakenly specified in kilometers, but
        not of satellites in valid orbits.
        """
        low = Satellite(
            name="Low",
            orbit=SunSynchronousOrbit(
                altitude=700, equator_crossing_time=time(10, 30), epoch=self.start
            ),
            instruments=[Instrument(name="Test")],
        )
        with self.assertWarnsRegex(UserWarning, r"\('Low'\).*0\.7 km.*in meters"):
            collect_orbit_track([low, self.satellite], self.times)
        with self.assertWarnsRegex(UserWarning, "in meters"):
            compute_access_periods(self.point, low, self.start, self.end)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            _check_satellites([self.satellite])
