"""
Unit tests for the LimbSounder and OccultationReceiver schemas.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

from pydantic import TypeAdapter, ValidationError

from tatc.schemas import (
    Aircraft,
    LimbSounder,
    OccultationReceiver,
    Satellite,
    ScanDirection,
    ScanStart,
    SunSynchronousOrbit,
    SurfacePlatform,
)
from tatc.schemas.air import Trajectory
from tatc.schemas.instrument import AllInstruments

START = datetime(2026, 1, 1, tzinfo=timezone.utc)
ORBIT = SunSynchronousOrbit(altitude=700e3, equator_crossing_time="13:30")


def _limb(**kwargs) -> LimbSounder:
    fields = {
        "name": "Limb",
        "scan_azimuth": 0,
        "scan_elevations": [10e3, 50e3, 90e3],
        "scan_duration": timedelta(seconds=20),
    }
    return LimbSounder(**{**fields, **kwargs})


class TestLimbSounder(unittest.TestCase):
    """
    Unit tests for the LimbSounder schema.
    """

    def test_defaults(self):
        """
        Test the default scans, cycle period, and sample elevation.
        """
        forward = _limb()
        self.assertEqual(forward.type, "limb")
        self.assertEqual(
            forward.get_scans(), [ScanStart(direction=ScanDirection.UPWARD)]
        )
        self.assertEqual(forward.get_cycle_period(), timedelta(seconds=20))
        self.assertEqual(forward.get_sample_elevation(), 50e3)
        rearward = _limb(scan_azimuth=180)
        self.assertEqual(
            rearward.get_scans(), [ScanStart(direction=ScanDirection.DOWNWARD)]
        )

    def test_parsed_by_type_or_fields(self):
        """
        Test that data parses to a limb sounder by its type or its fields, and
        round-trips through JSON.
        """
        adapter = TypeAdapter(AllInstruments)
        sounder = _limb(
            scans=[
                ScanStart(direction="upward"),
                ScanStart(offset=timedelta(seconds=30), direction="downward"),
            ],
            cycle_period=timedelta(seconds=60),
        )
        data = sounder.model_dump(mode="json")
        self.assertEqual(adapter.validate_python(data), sounder)
        del data["type"]
        self.assertEqual(adapter.validate_python(data), sounder)

    def test_bad_data(self):
        """
        Test that invalid scans, sample elevations, and unused fields raise.
        """
        with self.assertRaisesRegex(ValidationError, "increasing offsets"):
            _limb(
                scans=[
                    ScanStart(offset=timedelta(seconds=10), direction="upward"),
                    ScanStart(offset=timedelta(seconds=5), direction="downward"),
                ],
                cycle_period=timedelta(seconds=60),
            )
        with self.assertRaisesRegex(ValidationError, "increasing offsets"):
            _limb(scans=[ScanStart(offset=timedelta(seconds=20), direction="upward")])
        with self.assertRaisesRegex(ValidationError, "within the scan"):
            _limb(sample_elevation=100e3)
        with self.assertRaises(ValidationError):
            _limb(scan_elevations=[10e3])
        with self.assertRaises(ValidationError):
            _limb(scan_duration=timedelta(0))
        for field, value in [
            ("field_of_regard", 50),
            ("req_target_sunlit", True),
            ("min_target_glint_angle", 10),
        ]:
            with self.assertRaisesRegex(ValidationError, "does not use"):
                _limb(**{field: value})

    def test_no_view(self):
        """
        Test that a limb sounder has no swath or footprint.
        """
        sounder = _limb()
        with self.assertRaisesRegex(TypeError, "collect_limb_soundings"):
            sounder.get_swath_width(700e3)
        with self.assertRaisesRegex(TypeError, "collect_limb_soundings"):
            sounder.compute_footprint(None)

    def test_scan_times(self):
        """
        Test the start times of a cycle of alternating scans.
        """
        sounder = _limb(
            scans=[
                ScanStart(direction="upward"),
                ScanStart(offset=timedelta(seconds=30), direction="downward"),
            ],
            cycle_period=timedelta(seconds=60),
        )
        times = sounder.get_scan_times(START, START + timedelta(seconds=150))
        self.assertEqual(
            times[ScanDirection.UPWARD],
            [START + timedelta(seconds=s) for s in (0, 60, 120)],
        )
        self.assertEqual(
            times[ScanDirection.DOWNWARD],
            [START + timedelta(seconds=s) for s in (30, 90)],
        )
        self.assertEqual(
            sounder.get_scan_times(START, START),
            {ScanDirection.UPWARD: [], ScanDirection.DOWNWARD: []},
        )


class TestOccultationReceiver(unittest.TestCase):
    """
    Unit tests for the OccultationReceiver schema.
    """

    def test_defaults_and_parsing(self):
        """
        Test the defaults, which are those of collect_ro_observations, and
        parsing by type or fields.
        """
        receiver = OccultationReceiver(name="RO", constellations=["gps"])
        self.assertEqual(receiver.type, "occultation")
        self.assertEqual(receiver.max_yaw, 65)
        self.assertEqual(receiver.range_elevation, (-200e3, 60e3))
        self.assertEqual(receiver.sample_elevation, -80e3)
        self.assertEqual(receiver.min_profile_duration, timedelta(seconds=30))
        adapter = TypeAdapter(AllInstruments)
        data = receiver.model_dump(mode="json")
        self.assertEqual(adapter.validate_python(data), receiver)
        del data["type"]
        self.assertEqual(adapter.validate_python(data), receiver)

    def test_bad_data(self):
        """
        Test that invalid yaws, elevations, and unused fields raise.
        """
        with self.assertRaises(ValidationError):
            OccultationReceiver(max_yaw=0)
        with self.assertRaises(ValidationError):
            OccultationReceiver(max_yaw=95)
        with self.assertRaisesRegex(ValidationError, "range of elevations"):
            OccultationReceiver(sample_elevation=100e3)
        with self.assertRaisesRegex(ValidationError, "does not use"):
            OccultationReceiver(req_self_sunlit=True)
        with self.assertRaisesRegex(TypeError, "collect_occultations"):
            OccultationReceiver().get_min_elevation_angle(700e3)


class TestSounderPlatforms(unittest.TestCase):
    """
    Unit tests for the platforms that carry sounders.
    """

    def test_satellite_carries_sounders(self):
        """
        Test that a satellite carries sounders beside other instruments and
        round-trips through JSON.
        """
        satellite = Satellite(
            name="S",
            orbit=ORBIT,
            instruments=[_limb(), OccultationReceiver(name="RO")],
        )
        self.assertEqual(
            Satellite.model_validate_json(satellite.model_dump_json()), satellite
        )

    def test_others_reject_sounders(self):
        """
        Test that aircraft and surface platforms reject sounders.
        """
        trajectory = Trajectory.from_waypoints(
            [(-100, 35), (-99, 35)], start=START, speed=200, altitude=8000
        )
        for sounder in (_limb(), OccultationReceiver(name="RO")):
            with self.assertRaisesRegex(ValidationError, "only supported on a sat"):
                Aircraft(name="A", trajectory=trajectory, instruments=[sounder])
            with self.assertRaisesRegex(ValidationError, "only supported on a sat"):
                SurfacePlatform(
                    name="P",
                    latitude=35,
                    longitude=-100,
                    elevation=800,
                    height=1.5,
                    instruments=[sounder],
                )


if __name__ == "__main__":
    unittest.main()
