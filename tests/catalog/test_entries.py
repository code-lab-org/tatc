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
    COMPATIBLE_VERSIONS,
    STANDARD_VERSION,
    Grade,
    InstrumentClass,
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
            self.assertIn(entry.record.standard_version, COMPATIBLE_VERSIONS, entry.id)
            self.assertLessEqual(
                entry.record.standard_version, STANDARD_VERSION, entry.id
            )
            if entry.instrument_class == InstrumentClass.LIDAR:
                # the lidar class was added in version 2
                self.assertGreaterEqual(entry.record.standard_version, 2, entry.id)
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
            [e.id for e in list_entries(instrument_class="conical")],
            ["amsr2", "gmi", "smap_radiometer"],
        )
        self.assertEqual(
            [e.id for e in list_entries(instrument_class="geostationary")], ["abi"]
        )
        self.assertEqual([e.id for e in list_entries("conical", "validated")], ["gmi"])
        self.assertEqual(
            [e.id for e in list_entries("conical", "partial")],
            ["amsr2", "smap_radiometer"],
        )
        self.assertEqual(
            [e.id for e in list_entries(tier="partial")],
            [
                "abi",
                "airs",
                "amsr2",
                "amsua",
                "atlas",
                "avhrr",
                "cris",
                "csar",
                "dpr",
                "iasi",
                "lsar",
                "msi",
                "smap_radiometer",
                "tropomi",
            ],
        )
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

    def test_modis_entry_validated_on_terra(self):
        """
        Test that the MODIS entry passes with its exact values on held-out Terra data,
        and that its results on Aqua, a transfer, are not those of any gate.
        """
        entry = get_entry("modis")
        self.assertEqual(entry.validation.tier, Tier.VALIDATED)
        self.assertEqual(entry.validation.gaps, ())
        for key, configuration in entry.validation.configurations.items():
            record = configuration.record
            self.assertEqual(record.evidence, "rerun", key)
            self.assertEqual(record.platform, "terra", key)
            self.assertEqual(record.holdout_gap_days, 32, key)
            self.assertTrue(all(r.passed for r in configuration.gates), key)
            transfers = [m for m in record.metrics if m.scope is not None]
            self.assertTrue(transfers, key)
            gated = {r.gate.metric for r in configuration.gates}
            for metric in transfers:
                self.assertTrue(metric.name.startswith("transfer_"), metric.name)
                self.assertNotIn(metric.name, gated)

    def test_oli_entry_validated_on_landsat_9(self):
        """
        Test that the OLI entry passes on held-out Landsat 9 data, and that its results
        on Landsat 8, a transfer, are not those of any gate.
        """
        entry = get_entry("oli")
        self.assertEqual(entry.validation.tier, Tier.VALIDATED)
        configuration = entry.validation.configurations["swath"]
        record = configuration.record
        self.assertEqual(record.platform, "landsat-9")
        self.assertEqual(record.holdout_gap_days, 33)
        self.assertTrue(all(r.passed for r in configuration.gates))
        transfers = [m for m in record.metrics if m.scope is not None]
        self.assertTrue(transfers)
        for metric in transfers:
            self.assertTrue(metric.name.startswith("transfer_"), metric.name)

    def test_iasi_entry_partial_on_published_footprints(self):
        """
        Test that the IASI entry is partial because of the published footprint at its
        swath edge alone, with its swath validated on held-out data.
        """
        entry = get_entry("iasi")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        swath = entry.validation.configurations["swath"]
        self.assertEqual(swath.tier, Tier.VALIDATED)
        self.assertEqual(swath.record.holdout_gap_days, 32)
        pixels = entry.validation.configurations["pixels"]
        self.assertEqual(pixels.tier, Tier.PARTIAL)
        self.assertEqual(
            [r.gate.metric for r in pixels.gates if not r.passed],
            ["sample_size_within_3pct"],
        )

    def test_tropomi_entry_partial_on_long_time_steps(self):
        """
        Test that the TROPOMI entry, a union of views, is partial because of its swath
        area alone, which is the lowest over time steps.
        """
        entry = get_entry("tropomi")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        configuration = entry.validation.configurations["swath"]
        self.assertEqual(entry.configurations["swath"].grouping, "union")
        self.assertEqual(configuration.record.holdout_gap_days, 31)
        self.assertEqual(
            [r.gate.metric for r in configuration.gates if not r.passed],
            ["footprint_iou_median"],
        )
        self.assertTrue(
            any(
                "30 s or less" in text
                for text in entry.configurations["swath"].limitations
            )
        )

    def test_msi_entry_partial_on_long_time_steps(self):
        """
        Test that the MSI entry, tested on three satellites, is partial because of its
        swath area alone, which is the lowest over satellites and time steps.
        """
        entry = get_entry("msi")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        configuration = entry.validation.configurations["swath"]
        self.assertEqual(configuration.record.holdout_gap_days, 52)
        scopes = {m.scope for m in configuration.record.metrics}
        for platform in entry.configurations["swath"].platforms:
            self.assertIn(platform, scopes)
        self.assertEqual(
            [r.gate.metric for r in configuration.gates if not r.passed],
            ["footprint_iou_median"],
        )
        self.assertTrue(
            any(
                "30 s or less" in text
                for text in entry.configurations["swath"].limitations
            )
        )

    def test_cris_entry_partial_on_pixels(self):
        """
        Test that the CrIS entry is partial because of its pixels alone, with its swath
        validated at the timing grade on held-out data.
        """
        entry = get_entry("cris")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        swath = entry.validation.configurations["swath"]
        self.assertEqual(swath.tier, Tier.VALIDATED)
        self.assertEqual(entry.configurations["swath"].grade, Grade.TIMING)
        self.assertEqual(swath.record.holdout_gap_days, 51)
        pixels = entry.validation.configurations["pixels"]
        self.assertEqual(pixels.record.holdout_gap_days, 51)
        self.assertEqual(
            [r.gate.metric for r in pixels.gates if not r.passed],
            ["sample_size_within_3pct", "pixel_position_median"],
        )

    def test_csar_entry_partial(self):
        """
        Test that the C-SAR entry is partial: its IW mode misses points at the swath
        edges and fails the swath area, and its EW mode fails the swath area alone.
        """
        entry = get_entry("csar")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        failed = {
            key: [r.gate.metric for r in configuration.gates if not r.passed]
            for key, configuration in entry.validation.configurations.items()
        }
        self.assertEqual(
            failed,
            {
                "iw": ["pass_recall", "footprint_iou_median", "footprint_iou_p5"],
                "ew": ["footprint_iou_median", "footprint_iou_p5"],
            },
        )
        for configuration in entry.validation.configurations.values():
            self.assertEqual(configuration.record.holdout_gap_days, 51)
            self.assertEqual(configuration.record.platform, "sentinel-1c")

    def test_amsr2_entry_partial_on_published_footprints(self):
        """
        Test that the AMSR2 entry is partial because of its published footprints
        alone, with its swath validated at the timing grade on held-out data.
        """
        entry = get_entry("amsr2")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        swath = entry.validation.configurations["swath"]
        self.assertEqual(swath.tier, Tier.VALIDATED)
        self.assertEqual(entry.configurations["swath"].grade, Grade.TIMING)
        self.assertEqual(swath.record.holdout_gap_days, 36)
        sample_size = entry.validation.configurations["sample_size"]
        self.assertEqual(
            [r.gate.metric for r in sample_size.gates if not r.passed],
            ["sample_size_within_3pct"],
        )

    def test_karin_entry_validated(self):
        """
        Test that the KaRIn entry, a union of two swaths, passes with its exact values
        on held-out data.
        """
        entry = get_entry("karin")
        self.assertEqual(entry.validation.tier, Tier.VALIDATED)
        self.assertEqual(entry.configurations["swaths"].grouping, "union")
        configuration = entry.validation.configurations["swaths"]
        self.assertEqual(configuration.record.holdout_gap_days, 34)
        self.assertTrue(all(r.passed for r in configuration.gates))

    def test_dpr_entry_partial_at_long_time_steps(self):
        """
        Test that the DPR entry is partial because of its swath's area at a time step
        of 120 s alone, with its beams validated on data held out 75 days from the fit.
        """
        entry = get_entry("dpr")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        swath = entry.validation.configurations["swath"]
        self.assertEqual(
            [r.gate.metric for r in swath.gates if not r.passed],
            ["footprint_iou_median"],
        )
        self.assertEqual(swath.record.holdout_gap_days, 75)
        self.assertEqual(entry.validation.configurations["beams"].tier, Tier.VALIDATED)

    def test_mhs_entry_validated_on_published_values(self):
        """
        Test that the MHS entry, whose values are all published, is validated without
        a holdout, at the coverage grade.
        """
        entry = get_entry("mhs")
        self.assertEqual(entry.validation.tier, Tier.VALIDATED)
        self.assertEqual(entry.configurations["swath"].grade, Grade.COVERAGE)
        for configuration in entry.configurations.values():
            for instrument in configuration.instruments:
                self.assertNotIn(
                    "measured", [p.kind for p in instrument.provenance.values()]
                )

    def test_atlas_entry_partial_on_epochs(self):
        """
        Test that the ATLAS entry, of the lidar class added in version 2 of the
        standard, is partial because of its epochs alone on held-out data.
        """
        entry = get_entry("atlas")
        self.assertEqual(entry.instrument_class, InstrumentClass.LIDAR)
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        configuration = entry.validation.configurations["beam_pairs"]
        self.assertEqual(
            [r.gate.metric for r in configuration.gates if not r.passed],
            ["epoch_p95"],
        )
        self.assertEqual(configuration.record.holdout_gap_days, 34)
        self.assertEqual(entry.record.standard_version, 2)

    def test_smap_entry_partial_without_sample_size(self):
        """
        Test that the SMAP radiometer entry is partial only for want of a sample size
        configuration, with its looks validated on held-out data.
        """
        entry = get_entry("smap_radiometer")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        self.assertEqual(
            entry.validation.gaps, ("no configuration at the sample_size grade",)
        )
        looks = entry.validation.configurations["looks"]
        self.assertEqual(looks.tier, Tier.VALIDATED)
        self.assertEqual(looks.record.holdout_gap_days, 34)

    def test_lsar_entry_partial(self):
        """
        Test that the L-SAR entry is partial: its modes with a 20 or 40 MHz band miss
        the near-range edges, those with only a 5 MHz band miss a point of a small
        sample, and both fail the swath area at long time steps.
        """
        entry = get_entry("lsar")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        failed = {
            key: [r.gate.metric for r in configuration.gates if not r.passed]
            for key, configuration in entry.validation.configurations.items()
        }
        self.assertEqual(
            failed,
            {
                "swath_20_40mhz": [
                    "swath_edge_p95",
                    "footprint_iou_median",
                    "footprint_iou_p5",
                ],
                "swath_5mhz": [
                    "pass_recall",
                    "footprint_iou_median",
                    "footprint_iou_p5",
                ],
            },
        )
        for configuration in entry.validation.configurations.values():
            self.assertEqual(configuration.record.holdout_gap_days, 33)

    def test_aqua_sounder_entries_partial_on_published_footprints(self):
        """
        Test that the AIRS and AMSU-A entries, whose values are all published, are
        partial because of their published footprints alone, with their swaths
        validated, AIRS's at the timing grade and AMSU-A's at the coverage grade.
        """
        for entry_id, grade in (("airs", Grade.TIMING), ("amsua", Grade.COVERAGE)):
            entry = get_entry(entry_id)
            self.assertEqual(entry.validation.tier, Tier.PARTIAL, entry_id)
            swath = entry.validation.configurations["swath"]
            self.assertEqual(swath.tier, Tier.VALIDATED, entry_id)
            self.assertEqual(entry.configurations["swath"].grade, grade)
            footprints = entry.validation.configurations["footprints"]
            self.assertEqual(
                [r.gate.metric for r in footprints.gates if not r.passed],
                ["sample_size_within_3pct"],
                entry_id,
            )
            for configuration in entry.configurations.values():
                for instrument in configuration.instruments:
                    self.assertNotIn(
                        "measured", [p.kind for p in instrument.provenance.values()]
                    )

    def test_avhrr_entry_partial_on_published_footprint(self):
        """
        Test that the AVHRR entry is partial on data held out 33 days from the fit
        because of its published footprint alone, with its yawed swath validated at
        the timing grade.
        """
        entry = get_entry("avhrr")
        self.assertEqual(entry.validation.tier, Tier.PARTIAL)
        failed = {
            key: [r.gate.metric for r in configuration.gates if not r.passed]
            for key, configuration in entry.validation.configurations.items()
        }
        self.assertEqual(failed, {"swath": [], "pixels": ["sample_size_within_3pct"]})
        self.assertEqual(entry.configurations["swath"].grade, Grade.TIMING)
        for configuration in entry.validation.configurations.values():
            self.assertEqual(configuration.record.holdout_gap_days, 33)
        for key in ("swath", "pixels"):
            self.assertEqual(entry.configurations[key].attitude, {"yaw_angle": 0.2943})

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
