"""
Unit tests for loading catalog entries and building their instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import copy
import unittest
import warnings
from datetime import datetime, timedelta, timezone

import yaml
from pydantic import ValidationError

from tatc.analysis import collect_observations
from tatc.catalog import (
    STANDARD_VERSION,
    InstrumentEntry,
    Tier,
    get_entry,
    get_entry_ids,
    list_entries,
)
from tatc.catalog import loader
from tatc.schemas import (
    Attitude,
    AttitudeSchedule,
    AttitudeSegment,
    ConicalInstrument,
    GeneralPerturbationsOrbit,
    Instrument,
    Point,
    Satellite,
)
from tatc.utils import VelocityFrame

GPM_TLE = [
    "1 39574U 14009C   26277.54227627  .00004528  00000+0  10004-3 0  9997",
    "2 39574  64.9701 301.6724 0010793 268.6572  91.3332 15.45217995714704",
]


class TestLoading(unittest.TestCase):
    """
    Unit tests for loading entries.
    """

    def test_entries(self):
        """
        Test that the cataloged entries load, each with a current record.
        """
        self.assertIn("abi", get_entry_ids())
        self.assertIn("gmi", get_entry_ids())
        for entry_id in get_entry_ids():
            self.assertEqual(get_entry(entry_id).id, entry_id)

    def test_records_current(self):
        """
        Test that each record tested its configuration as cataloged, so a
        configuration cannot change without being validated again.
        """
        for entry in list_entries():
            self.assertIsNotNone(entry.record, entry.id)
            self.assertEqual(entry.record.standard_version, STANDARD_VERSION, entry.id)
            for key, record in entry.record.configurations.items():
                self.assertEqual(
                    record.digest, entry.configurations[key].digest, (entry.id, key)
                )

    def test_alias_lookup(self):
        """
        Test lookup by alias and name, ignoring case.
        """
        self.assertEqual(get_entry("GPM Microwave Imager").id, "gmi")
        self.assertEqual(get_entry("GMI").id, "gmi")

    def test_unknown_entry(self):
        """
        Test that an unknown entry is an error that lists those available.
        """
        with self.assertRaisesRegex(KeyError, "abi"):
            get_entry("unknown")

    def test_list_filters(self):
        """
        Test filtering by class and tier.
        """
        self.assertEqual(
            [e.id for e in list_entries(instrument_class="conical")], ["gmi"]
        )
        self.assertEqual(
            [e.id for e in list_entries(instrument_class="geostationary")], ["abi"]
        )
        self.assertEqual([e.id for e in list_entries("conical", "validated")], ["gmi"])
        self.assertEqual(list_entries("conical", "partial"), [])
        self.assertEqual([e.id for e in list_entries(tier="partial")], ["abi"])
        with self.assertRaises(ValueError):
            list_entries(tier="approved")

    def test_loaded_copies_independent(self):
        """
        Test that changing a loaded entry does not change later ones.
        """
        entry = get_entry("abi")
        entry.configurations["full_disk"].instruments[0].parameters[
            "field_of_regard"
        ] = 1
        self.assertEqual(
            get_entry("abi")
            .configurations["full_disk"]
            .instruments[0]
            .parameters["field_of_regard"],
            17.4,
        )

    def test_loads_only_requested_entry(self):
        """
        Test that getting an entry by identifier parses only its own files.
        """
        loader._load_entry.cache_clear()  # pylint: disable=protected-access
        get_entry("gmi")
        info = loader._load_entry.cache_info()  # pylint: disable=protected-access
        self.assertEqual(info.currsize, 1)


def _gmi_data() -> dict:
    with open(
        loader._ENTRIES / "gmi.yml",
        "r",
        encoding="utf-8",  # pylint: disable=protected-access
    ) as f:
        return yaml.safe_load(f)


class TestStrictFiles(unittest.TestCase):
    """
    Unit tests that errors in hand-written files are rejected.
    """

    def setUp(self):
        self.data = _gmi_data()
        self.feeds = self.data["configurations"]["feeds"]
        self.instrument = self.feeds["instruments"][0]

    def test_baseline_valid(self):
        """
        Test that the unmodified data is valid.
        """
        InstrumentEntry(**self.data)

    def _assert_invalid(self):
        with self.assertRaises(ValidationError):
            InstrumentEntry(**self.data)

    def test_misspelled_parameter(self):
        """
        Test that a misspelled instrument parameter is rejected.
        """
        self.instrument["parameters"]["cone_anglee"] = 48
        self._assert_invalid()

    def test_invalid_parameter_value(self):
        """
        Test that a parameter out of range is rejected when loaded, not when used.
        """
        self.instrument["parameters"]["cone_angle"] = 120
        self._assert_invalid()

    def test_velocity_frame_parameter(self):
        """
        Test that an instrument cannot set the velocity frame, which belongs to
        the attitude of the space system.
        """
        self.instrument["parameters"]["velocity_frame"] = "inertial"
        self.instrument["provenance"]["velocity_frame"] = {
            "kind": "nominal",
            "source": "x",
        }
        self._assert_invalid()

    def test_missing_provenance(self):
        """
        Test that each parameter needs provenance.
        """
        del self.instrument["provenance"]["cone_angle"]
        self._assert_invalid()

    def test_measured_needs_period(self):
        """
        Test that a measured value states what it was fitted to.
        """
        self.instrument["provenance"]["cone_angle"] = {
            "kind": "measured",
            "source": "x",
        }
        self._assert_invalid()

    def test_duplicate_names(self):
        """
        Test that names are unique, because instruments are found by name.
        """
        self.feeds["instruments"][1]["parameters"]["name"] = "GMI S1"
        self._assert_invalid()

    def test_unknown_field(self):
        """
        Test that an unknown field of an entry is rejected.
        """
        self.data["tier"] = "validated"
        self._assert_invalid()

    def test_unknown_platform(self):
        """
        Test that configurations refer to the entry's platforms.
        """
        self.feeds["platforms"] = ["trmm"]
        self._assert_invalid()

    def test_index_out_of_range(self):
        """
        Test that the default instrument exists.
        """
        self.feeds["instrument_index"] = 2
        self._assert_invalid()

    def test_grouping_needs_instruments(self):
        """
        Test that one instrument cannot be a group, nor several ungrouped.
        """
        self.feeds["grouping"] = "single"
        self._assert_invalid()


class TestBuild(unittest.TestCase):
    """
    Unit tests for building instruments and satellites.
    """

    def setUp(self):
        self.entry = get_entry("gmi")
        self.orbit = GeneralPerturbationsOrbit.from_tle(GPM_TLE)

    def test_payload(self):
        """
        Test that a payload has the configured instruments, attitude, and settings.
        """
        payload = self.entry.build(orbit=self.orbit, time_step=10)
        self.assertEqual([i.name for i in payload.instruments], ["GMI S1", "GMI S2"])
        self.assertIsInstance(payload.instruments[0], ConicalInstrument)
        self.assertEqual(payload.attitude.velocity_frame, VelocityFrame.INERTIAL)
        self.assertEqual(payload.instrument_index, 0)
        self.assertEqual(payload.grouping, "alternatives")
        self.assertFalse(payload.transferred)
        self.assertTrue(any("flies backward" in l for l in payload.limitations))

    def test_tailored_field_of_view(self):
        """
        Test that tailored fields of view match those of the GMI notebook, whose
        tests reported these values for this orbit.
        """
        for time_step, expected in ((1.875, 1.78), (10, 9.45), (30, 27.85)):
            payload = self.entry.build(orbit=self.orbit, time_step=time_step)
            self.assertAlmostEqual(
                payload.instruments[0].along_track_field_of_view, expected, places=2
            )

    def test_tailoring_needs_orbit_and_time_step(self):
        """
        Test that a tailored view cannot be built without what it depends on.
        """
        with self.assertRaisesRegex(ValueError, "time step"):
            self.entry.build(orbit=self.orbit)
        with self.assertRaisesRegex(ValueError, "orbit"):
            self.entry.build(time_step=10)

    def test_swath_widths(self):
        """
        Test that the swath widths match those observed at about GPM's altitude,
        within the 5 km gate and the range of altitudes of the test.
        """
        payload = self.entry.build(orbit=self.orbit, time_step=10)
        self.assertAlmostEqual(
            payload.instruments[0].get_swath_width(440e3) / 1e3, 1012.283, delta=5
        )
        self.assertAlmostEqual(
            payload.instruments[1].get_swath_width(440e3) / 1e3, 896.509, delta=5
        )

    def test_to_satellite_inertial_frame(self):
        """
        Test that the satellite has the validated attitude, and that its
        instruments pass the check of velocity frames.
        """
        satellite = self.entry.to_satellite(self.orbit, time_step=10, name="GPM")
        self.assertEqual(satellite.name, "GPM")
        self.assertEqual(satellite.attitude.velocity_frame, VelocityFrame.INERTIAL)
        for instrument in satellite.instruments:
            self.assertNotIn("velocity_frame", instrument.model_fields_set)

    def test_satellite_round_trip(self):
        """
        Test that a satellite built from an entry survives a round trip.
        """
        satellite = self.entry.to_satellite(self.orbit, time_step=10)
        copied = Satellite.model_validate_json(satellite.model_dump_json())
        self.assertEqual(copied.instruments, satellite.instruments)
        self.assertEqual(copied.attitude, satellite.attitude)

    def test_collect_observations(self):
        """
        Test each instrument of a satellite built from an entry in an analysis.
        The arc of the high-frequency feed (S2) lies inside and behind that of
        the low-frequency feed (S1), so S1 observes each point that S2 does,
        several seconds sooner.
        """
        satellite = self.entry.to_satellite(self.orbit, time_step=10)
        start = datetime(2026, 10, 3, tzinfo=timezone.utc)
        end = start + timedelta(days=1)
        delays = []
        for latitude in (-45, 0, 45):
            for longitude in (-120, 0, 120):
                point = Point(latitude=latitude, longitude=longitude)
                outer = collect_observations(point, satellite, start, end, 0)
                inner = collect_observations(point, satellite, start, end, 1)
                self.assertGreaterEqual(len(outer), len(inner))
                for epoch in inner.epoch:
                    delays.append(
                        min(
                            (epoch - other).total_seconds()
                            for other in outer.epoch
                            if other < epoch
                        )
                    )
        self.assertGreater(len(delays), 0)
        self.assertGreater(min(delays), 5)
        self.assertLess(max(delays), 60)

    def test_builds_share_no_state(self):
        """
        Test that satellites from one entry share no mutable state, and that
        changing one's instruments changes neither the other's nor the entry's.
        """
        first = self.entry.to_satellite(self.orbit, time_step=10)
        second = self.entry.to_satellite(self.orbit, time_step=10)
        self.assertIsNot(first.instruments[0], second.instruments[0])
        self.assertIsNot(first.attitude, second.attitude)
        before = copy.deepcopy(second.instruments[0].cone_angle)
        first.instruments[0].cone_angle = 1
        self.assertEqual(second.instruments[0].cone_angle, before)
        self.assertEqual(
            self.entry.to_satellite(self.orbit, time_step=10).instruments[0].cone_angle,
            before,
        )

    def test_transfer_warns(self):
        """
        Test that using a configuration on a platform on which it was not tested
        gives a warning and is reported, and that the tested platform does not.
        """
        with self.assertWarnsRegex(UserWarning, "trmm"):
            payload = self.entry.build(orbit=self.orbit, time_step=10, platform="trmm")
        self.assertTrue(payload.transferred)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            payload = self.entry.build(orbit=self.orbit, time_step=10, platform="gpm")
        self.assertFalse(payload.transferred)

    def test_replacement_attitude(self):
        """
        Test a yaw flip schedule, in the validated velocity frame, in place of
        the validated attitude.
        """
        schedule = AttitudeSchedule(
            segments=[
                AttitudeSegment(
                    start=datetime(2026, 10, 3, tzinfo=timezone.utc),
                    attitude=Attitude(
                        velocity_frame=VelocityFrame.INERTIAL, yaw_angle=180
                    ),
                )
            ]
        )
        satellite = self.entry.to_satellite(self.orbit, time_step=10, attitude=schedule)
        self.assertIs(satellite.attitude, schedule)

    def test_replacement_attitude_frame_mismatch(self):
        """
        Test that an attitude in another velocity frame is rejected, whether
        directly or in a schedule, as it would misplace the swath.
        """
        with self.assertRaisesRegex(ValueError, "inertial"):
            self.entry.build(orbit=self.orbit, time_step=10, attitude=Attitude())
        schedule = AttitudeSchedule(
            segments=[
                AttitudeSegment(
                    start=datetime(2026, 10, 3, tzinfo=timezone.utc),
                    attitude=Attitude(),
                )
            ]
        )
        with self.assertRaisesRegex(ValueError, "inertial"):
            self.entry.build(orbit=self.orbit, time_step=10, attitude=schedule)

    def test_unknown_configuration(self):
        """
        Test that an unknown configuration is an error that lists those available.
        """
        with self.assertRaisesRegex(ValueError, "feeds"):
            self.entry.build("other", orbit=self.orbit, time_step=10)


class TestEntryStatus(unittest.TestCase):
    """
    Unit tests for the status of the entries, which must not be overstated.
    """

    def test_gmi_validated(self):
        """
        Test that GMI's configurations pass with its exact values, with fitted
        values tested on held-out data, at each of the grades of its class.
        """
        validation = get_entry("gmi").validation
        self.assertEqual(validation.tier, Tier.VALIDATED)
        self.assertEqual(validation.gaps, ())
        for key in ("feeds", "sample_size"):
            configuration = validation.configurations[key]
            self.assertEqual(configuration.tier, Tier.VALIDATED, key)
            self.assertEqual(configuration.record.evidence, "rerun", key)
            self.assertEqual(configuration.record.holdout_gap_days, 75, key)
            self.assertTrue(all(r.passed for r in configuration.gates), key)
            self.assertEqual(configuration.gaps, (), key)

    def test_noaa20_entries_validated(self):
        """
        Test that the entries validated against NOAA-20 data pass with their exact
        values, with fitted values tested on held-out data, at each required grade.
        """
        for entry_id in ("atms", "viirs"):
            validation = get_entry(entry_id).validation
            self.assertEqual(validation.tier, Tier.VALIDATED, entry_id)
            self.assertEqual(validation.gaps, (), entry_id)
            for key, configuration in validation.configurations.items():
                self.assertEqual(configuration.record.evidence, "rerun", key)
                self.assertEqual(configuration.record.holdout_gap_days, 51, key)
                self.assertTrue(all(r.passed for r in configuration.gates), key)

    def test_gmi_reports_spacing_limitation(self):
        """
        Test that GMI's spacing between scans, which no gate covers, is reported
        in its record and stated as a limitation of the configuration.
        """
        entry = get_entry("gmi")
        spacing = [
            metric.value
            for metric in entry.record.configurations["sample_size"].metrics
            if metric.name == "spacing_scan_ratio_median"
        ]
        self.assertEqual(len(spacing), 2)
        for value in spacing:
            self.assertGreater(value, 1.03)
        self.assertTrue(
            any(
                "spacing between scans" in limitation
                for limitation in entry.build(
                    "sample_size",
                    orbit=GeneralPerturbationsOrbit.from_tle(GPM_TLE),
                    time_step=1.875,
                ).limitations
            )
        )

    def test_abi_partial(self):
        """
        Test that ABI's product coverage is validated, but not the entry: its
        full disk configuration has no discriminating test.
        """
        validation = get_entry("abi").validation
        self.assertEqual(validation.tier, Tier.PARTIAL)
        self.assertEqual(validation.configurations["full_disk"].tier, Tier.PARTIAL)
        for key in ("sst_good_quality", "lst_good_quality"):
            self.assertEqual(validation.configurations[key].tier, Tier.VALIDATED, key)

    def test_abi_builds_without_attitude(self):
        """
        Test that an Earth-fixed configuration without other attitude fields
        builds as it was tested, without an attitude.
        """
        payload = get_entry("abi").build()
        self.assertIsNone(payload.attitude)
        self.assertIsInstance(payload.instruments[0], Instrument)


if __name__ == "__main__":
    unittest.main()
