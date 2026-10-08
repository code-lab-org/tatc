"""
Unit tests for the input validation shared by tatc.analysis functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta, timezone

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
    compute_ground_track,
)
from tatc.analysis.check import _check_satellite, _check_satellites
from tatc.schemas import GroundStation, Point

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
