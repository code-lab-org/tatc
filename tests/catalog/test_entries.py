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
    Grade,
    InstrumentEntry,
    build_payload,
    get_entry,
    get_entry_ids,
    list_entries,
)
from tatc.catalog import loader
from tatc.catalog.standard import is_compatible
from tatc.schemas import (
    Attitude,
    AttitudeSchedule,
    AttitudeSegment,
    ConicalInstrument,
    GeneralPerturbationsOrbit,
    GlintAttitude,
    Instrument,
    Point,
    Satellite,
)
from tatc.utils import VelocityFrame

GPM_TLE = [
    "1 39574U 14009C   26277.54227627  .00004528  00000+0  10004-3 0  9997",
    "2 39574  64.9701 301.6724 0010793 268.6572  91.3332 15.45217995714704",
]
OCO2_TLE = [
    "1 40059U 14035A   26313.25381781  .00000489  00000+0  11850-3 0  9999",
    "2 40059  98.2167 252.4180 0000990  88.3344 271.7968 14.57127601657173",
]

# the validation status of each entry: its tier, and for each configuration its
# tier, the days between its fitting data and its test (`None` without fitted
# values), and its gaps (by metric, for its gates), and the gaps of the entry
EXPECTED_STATUS = {
    "abi": (
        "partial",
        {
            "full_disk": ("partial", None, ["product_limit"]),
            "sst_good_quality": ("validated", None, []),
            "lst_good_quality": ("validated", None, []),
        },
        [],
    ),
    "ahi": (
        "partial",
        {
            "full_disk": ("partial", None, ["product_limit"]),
            "cloud_top_height": ("partial", 38, ["product_limit"]),
        },
        [],
    ),
    "airs": (
        "partial",
        {
            "swath": ("validated", None, []),
            "footprints": ("partial", None, ["sample_size_within_3pct"]),
        },
        [],
    ),
    "amsr2": (
        "partial",
        {
            "swath": ("validated", 36, []),
            "sample_size": ("partial", 36, ["sample_size_within_3pct"]),
        },
        [],
    ),
    "amsua": (
        "partial",
        {
            "swath": ("validated", None, []),
            "footprints": ("partial", None, ["sample_size_within_3pct"]),
        },
        [],
    ),
    "ascat": (
        "validated",
        {
            "swaths": ("validated", 34, []),
        },
        [],
    ),
    "atlas": (
        "partial",
        {
            "beam_pairs": ("partial", 34, ["epoch_p95"]),
        },
        [],
    ),
    "atlid": (
        "partial",
        {
            "profile": ("partial", 41, ["track_offset_median"]),
        },
        [],
    ),
    "atms": (
        "validated",
        {
            "swath_noaa_20": ("validated", 51, []),
            "beams_noaa_20": ("validated", 51, []),
            "swath_noaa_21": ("validated", 34, []),
            "beams_noaa_21": ("validated", 34, []),
            "swath_suomi_npp": ("validated", 38, []),
            "beams_suomi_npp": ("validated", 38, []),
        },
        [],
    ),
    "avhrr": (
        "partial",
        {
            "swath": ("validated", 33, []),
            "pixels": ("partial", 33, ["sample_size_within_3pct"]),
        },
        [],
    ),
    "bbr": (
        "partial",
        {
            "views": ("partial", 41, ["pass_recall", "footprint_iou_median"]),
        },
        [],
    ),
    "ceres": (
        "partial",
        {
            "swath_noaa_20": ("partial", 36, ["swath_width_median"]),
            "swath_suomi_npp": ("validated", 35, []),
        },
        ["no configuration at the sample_size grade"],
    ),
    "cpr": (
        "partial",
        {
            "nadir": ("partial", None, ["track_offset_median"]),
        },
        [],
    ),
    "cris": (
        "partial",
        {
            "swath_noaa_20": ("validated", 51, []),
            "pixels_noaa_20": (
                "partial",
                51,
                ["sample_size_within_3pct", "pixel_position_median"],
            ),
            "swath_noaa_21": ("validated", 34, []),
            "pixels_noaa_21": (
                "partial",
                34,
                ["sample_size_within_3pct", "pixel_position_median"],
            ),
            "swath_suomi_npp": ("validated", 37, []),
            "pixels_suomi_npp": (
                "partial",
                37,
                ["sample_size_within_3pct", "pixel_position_median"],
            ),
        },
        [],
    ),
    "csar": (
        "partial",
        {
            "iw": ("partial", 51, ["pass_recall"]),
            "ew": ("validated", 39, []),
        },
        [],
    ),
    "dpr": (
        "validated",
        {
            "swath": ("validated", 38, []),
            "beams": ("validated", 38, []),
        },
        [],
    ),
    "earthcare_msi": (
        "validated",
        {
            "swath": ("validated", 41, []),
        },
        [],
    ),
    "ecostress": (
        "partial",
        {
            "swath": ("partial", 35, ["pass_recall"]),
        },
        ["no configuration at the sample_size grade"],
    ),
    "emit": (
        "partial",
        {
            "swath": ("partial", 35, ["pass_recall", "footprint_iou_median"]),
        },
        [],
    ),
    "gedi": (
        "partial",
        {
            "beams": ("partial", 38, ["epoch_p95"]),
        },
        [],
    ),
    "gmi": (
        "validated",
        {
            "feeds": ("validated", 75, []),
            "sample_size": ("validated", 75, []),
        },
        [],
    ),
    "harp2": (
        "validated",
        {
            "swath": ("validated", 34, []),
        },
        [],
    ),
    "iasi": (
        "partial",
        {
            "swath": ("validated", 32, []),
            "pixels": ("partial", 32, ["sample_size_within_3pct"]),
        },
        [],
    ),
    "karin": (
        "validated",
        {
            "swaths": ("validated", 34, []),
        },
        [],
    ),
    "lsar": (
        "partial",
        {
            "swath_20_40mhz": ("partial", 33, ["swath_edge_p95"]),
            "swath_5mhz": ("partial", 33, ["pass_recall"]),
            "swath_77mhz": (
                "partial",
                30,
                ["footprint_iou_median", "footprint_iou_p5"],
            ),
        },
        [],
    ),
    "mhs": (
        "validated",
        {
            "swath": ("validated", None, []),
            "beams": ("validated", None, []),
        },
        [],
    ),
    "modis": (
        "validated",
        {
            "swath": ("validated", 32, []),
            "pixels": ("validated", 32, []),
        },
        [],
    ),
    "msi": (
        "validated",
        {
            "swath": ("validated", 38, []),
        },
        [],
    ),
    "oci": (
        "partial",
        {
            "swath": ("validated", 33, []),
        },
        ["no configuration at the sample_size grade"],
    ),
    "oco2": (
        "partial",
        {
            "nadir": ("partial", 34, ["pass_recall", "footprint_iou_median"]),
            "glint": (
                "partial",
                34,
                ["pass_recall", "footprint_iou_median", "epoch_p95"],
            ),
        },
        [],
    ),
    "olci": (
        "partial",
        {
            "swath": ("partial", 31, ["swath_width_median"]),
        },
        [],
    ),
    "oli": (
        "validated",
        {
            "swath": ("validated", 33, []),
        },
        [],
    ),
    "omps_nm": (
        "partial",
        {
            "swath_noaa_21": ("partial", 35, ["swath_width_median", "epoch_p95"]),
            "swath_noaa_20": ("partial", 35, ["swath_width_median", "epoch_p95"]),
            "swath_suomi_npp": (
                "partial",
                37,
                ["pass_recall", "swath_width_median", "epoch_p95"],
            ),
        },
        [],
    ),
    "poseidon4": (
        "validated",
        {
            "nadir": ("validated", None, []),
        },
        [],
    ),
    "smap_radiometer": (
        "partial",
        {
            "looks": ("validated", 34, []),
            "footprints": ("partial", 34, ["sample_size_within_3pct"]),
        },
        [],
    ),
    "spexone": (
        "validated",
        {
            "swath": ("validated", 34, []),
        },
        [],
    ),
    "ssmis": (
        "partial",
        {
            "swath": ("validated", 30, []),
        },
        ["no configuration at the sample_size grade"],
    ),
    "tempo": (
        "validated",
        {
            "field_of_regard": ("validated", 33, []),
        },
        [],
    ),
    "tmi": (
        "validated",
        {
            "feeds": ("validated", 39, []),
            "sample_size": ("validated", 39, []),
        },
        [],
    ),
    "trmm_pr": (
        "partial",
        {
            "swath": ("validated", 39, []),
            "beams": ("partial", 39, ["sample_size_within_3pct"]),
        },
        [],
    ),
    "tropomi": (
        "validated",
        {
            "swath": ("validated", 36, []),
        },
        [],
    ),
    "viirs": (
        "validated",
        {
            "swath_noaa_20": ("validated", 51, []),
            "swath_noaa_21": ("validated", 34, []),
            "swath_suomi_npp": ("validated", 38, []),
            "zones_noaa_20": ("validated", 51, []),
            "zones_noaa_21": ("validated", 34, []),
            "zones_suomi_npp": ("validated", 38, []),
        },
        [],
    ),
    "virs": (
        "validated",
        {
            "swath": ("validated", 39, []),
            "pixels": ("validated", 39, []),
        },
        [],
    ),
}


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
            for key, record in entry.record.configurations.items():
                self.assertEqual(
                    record.digest, entry.configurations[key].digest, (entry.id, key)
                )
                self.assertTrue(
                    is_compatible(entry.instrument_class, record.standard_version),
                    (entry.id, key),
                )

    def test_platform_names(self):
        """
        Test that, in an entry with configurations fitted on more than one platform,
        each configuration's name without its platform selects it for its platform.
        """
        named = []
        for entry in list_entries():
            configurations = entry.configurations
            if len({c.platforms for c in configurations.values()}) > 1:
                for key, configuration in configurations.items():
                    (platform,) = configuration.platforms
                    name = key.removesuffix("_" + platform.replace("-", "_"))
                    self.assertEqual(
                        entry.select_configuration(name, platform=platform), key
                    )
                named.append(entry.id)
        self.assertEqual(named, ["atms", "ceres", "cris", "omps_nm", "viirs"])

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
            ["amsr2", "gmi", "smap_radiometer", "ssmis", "tmi"],
        )
        self.assertEqual(
            [e.id for e in list_entries(instrument_class="geostationary")],
            ["abi", "ahi"],
        )
        self.assertEqual(
            [e.id for e in list_entries("conical", "validated")], ["gmi", "tmi"]
        )
        self.assertEqual(
            [e.id for e in list_entries("conical", "partial")],
            ["amsr2", "smap_radiometer", "ssmis"],
        )
        self.assertEqual(
            [e.id for e in list_entries(tier="partial")],
            [k for k, (tier, _, _) in EXPECTED_STATUS.items() if tier == "partial"],
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

    def test_platform_variant_names(self):
        """
        Test that, in an entry with configurations tested on different
        platforms, each configuration's name ends with its platform.
        """
        with open(
            loader._ENTRIES / "cris.yml",  # pylint: disable=protected-access
            "r",
            encoding="utf-8",
        ) as f:
            data = yaml.safe_load(f)
        InstrumentEntry(**data)
        configurations = data["configurations"]
        configurations["swath_21"] = configurations.pop("swath_noaa_21")
        with self.assertRaisesRegex(ValidationError, "ends with"):
            InstrumentEntry(**data)

    def test_unknown_platform(self):
        """
        Test that configurations refer to the entry's platforms.
        """
        self.feeds["platforms"] = ["trmm"]
        self._assert_invalid()

    def test_tailored_pointed_instrument_is_swept(self):
        """
        Test that a pointed instrument whose view is tailored to a time step must
        be swept, so that its swaths do not depend on the time step.
        """
        with open(
            loader._ENTRIES / "msi.yml",  # pylint: disable=protected-access
            "r",
            encoding="utf-8",
        ) as f:
            data = yaml.safe_load(f)
        InstrumentEntry(**data)
        instrument = data["configurations"]["swath"]["instruments"][0]
        del instrument["parameters"]["is_swept"]
        with self.assertRaises(ValidationError):
            InstrumentEntry(**data)

    def test_default_provenance(self):
        """
        Test that the provenance of a field of view tailored to the time step and of
        `is_swept` is the same in every configuration, and that of other modeling
        choices is not.
        """
        with open(
            loader._ENTRIES / "msi.yml",  # pylint: disable=protected-access
            "r",
            encoding="utf-8",
        ) as f:
            data = yaml.safe_load(f)
        instrument = data["configurations"]["swath"]["instruments"][0]
        self.assertNotIn("along_track_field_of_view", instrument["provenance"])
        self.assertNotIn("is_swept", instrument["provenance"])
        entry = InstrumentEntry(**data)
        provenance = entry.configurations["swath"].instruments[0].provenance
        self.assertEqual(provenance["along_track_field_of_view"].kind, "derived")
        self.assertEqual(provenance["is_swept"].kind, "derived")
        del instrument["provenance"]["is_rectangular"]
        with self.assertRaisesRegex(ValidationError, "is_rectangular"):
            InstrumentEntry(**data)

    def test_unknown_instrument_type(self):
        """
        Test that an instrument's type is one of the instrument schemas.
        """
        self.instrument["type"] = "conic"
        self._assert_invalid()

    def test_index_only_of_alternatives(self):
        """
        Test that only alternatives have an instrument to analyze by default.
        """
        self.feeds["instrument_index"] = 1
        InstrumentEntry(**self.data)
        self.feeds["grouping"] = "union"
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
        del self.feeds["grouping"]
        self._assert_invalid()
        self.feeds["grouping"] = "union"
        self.feeds["instruments"] = self.feeds["instruments"][:1]
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
        self.assertEqual(payload.references, ("gmi:feeds",))
        self.assertFalse(payload.transferred)
        self.assertTrue(any("flies backward" in l for l in payload.limitations))
        satellite = payload.to_satellite(self.orbit, name="GPM")
        self.assertEqual(satellite.name, "GPM")
        self.assertEqual(satellite.instruments, payload.instruments)
        self.assertEqual(satellite.attitude, payload.attitude)

    def test_tailored_field_of_view(self):
        """
        Test that tailored fields of view match the values of 1.78, 9.45, and 27.85
        deg reported for this orbit when the configuration was validated.
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
            self.entry.build(configuration="other", orbit=self.orbit, time_step=10)


def _summarize_gaps(validation) -> list[str]:
    """
    Gets the gaps of a configuration, with those of its gates as their metrics.
    """
    metrics = [result.gate.metric for result in validation.gates]
    return [
        next((m for m in metrics if text.startswith(f"{m}:")), text)
        for text in validation.gaps
    ]


class TestEntryStatus(unittest.TestCase):
    """
    Unit tests for the status of the entries, which must not be overstated.
    """

    def test_status(self):
        """
        Test the status of each entry against the table of expected status, so
        that no change to an entry, its record, or the standard changes a status
        unnoticed.
        """
        self.assertEqual(sorted(EXPECTED_STATUS), get_entry_ids())
        for entry_id, (tier, configurations, gaps) in EXPECTED_STATUS.items():
            with self.subTest(entry_id):
                validation = get_entry(entry_id).validation
                self.assertEqual(validation.tier.value, tier)
                self.assertEqual(
                    [
                        f"no configuration at the {grade.value} grade"
                        for grade in validation.missing_grades
                    ],
                    gaps,
                )
                self.assertEqual(set(validation.configurations), set(configurations))
                for key, (c_tier, holdout, c_gaps) in configurations.items():
                    configuration = validation.configurations[key]
                    self.assertEqual(configuration.tier.value, c_tier, key)
                    self.assertEqual(configuration.holdout_gap_days, holdout, key)
                    self.assertEqual(_summarize_gaps(configuration), c_gaps, key)

    def test_atms_noaa21_models_yaw(self):
        """
        Test that NOAA-21's ATMS configurations model the measured rotation of the
        scan line with their instruments' angles, not the body's attitude, and that
        NOAA-20's do not model one.
        """
        entry = get_entry("atms")
        for key, yaw in (("swath_noaa_21", 0.9341), ("beams_noaa_21", 0.934304)):
            configuration = entry.configurations[key]
            self.assertEqual(configuration.attitude, {}, key)
            for instrument in configuration.instruments:
                self.assertEqual(instrument.parameters["yaw_angle"], yaw, key)
                self.assertEqual(instrument.parameters["pitch_angle"], -0.1362, key)
        for key in ("swath_noaa_20", "beams_noaa_20"):
            configuration = entry.configurations[key]
            self.assertEqual(configuration.attitude, {}, key)
            for instrument in configuration.instruments:
                self.assertNotIn("yaw_angle", instrument.parameters, key)

    def test_body_attitudes(self):
        """
        Test that only configurations whose body is known to turn its instruments
        (the ISS's, and OCO-2's glint mode) give the body an attitude, so that
        angles fitted to one instrument are its own.
        """
        with_attitude = {
            f"{entry_id}:{key}"
            for entry_id in get_entry_ids()
            for key, configuration in get_entry(entry_id).configurations.items()
            if configuration.attitude
        }
        self.assertEqual(
            with_attitude,
            {"ecostress:swath", "emit:swath", "gedi:beams", "oco2:glint"},
        )

    def test_entry_lists_every_gap(self):
        """
        Test that an entry's gaps include those of its configurations, after their
        identifiers, so that a partial entry states why.
        """
        validation = get_entry("lsar").validation
        self.assertEqual(validation.tier.value, "partial")
        self.assertEqual(validation.missing_grades, ())
        self.assertEqual(
            list(validation.gaps),
            [
                f"{key}: {gap}"
                for key, configuration in validation.configurations.items()
                for gap in configuration.gaps
            ],
        )
        self.assertTrue(validation.gaps)

    def test_gate_table(self):
        """
        Test that the table of gates has a row for each gate of each
        configuration, with untested gates neither passed nor failed.
        """
        for entry_id in ("viirs", "airs"):
            entry = get_entry(entry_id)
            table = entry.validation.to_frame()
            rows = [
                (key, result)
                for key, validation in entry.validation.configurations.items()
                for result in validation.gates
            ]
            self.assertEqual(len(table), len(rows))
            for (_, row), (key, result) in zip(table.iterrows(), rows):
                self.assertEqual(row["configuration"], key)
                self.assertEqual(row["grade"], entry.configurations[key].grade.value)
                self.assertEqual(row["metric"], result.gate.label)
                self.assertEqual(row["passed"], result.passed)

    def test_mhs_has_no_fitted_values(self):
        """
        Test that the MHS entry, whose values are all published, needs no holdout.
        """
        entry = get_entry("mhs")
        for configuration in entry.configurations.values():
            self.assertIsNone(configuration.get_holdout_gap("2000-01-01"))

    def test_oco2_glint_builds_glint_attitude(self):
        """
        Test that OCO-2's glint configuration builds the glint attitude it declares,
        and its nadir configuration a plain attitude.
        """
        entry = get_entry("oco2")
        orbit = GeneralPerturbationsOrbit.from_tle(OCO2_TLE)
        glint = entry.build(configuration="glint", orbit=orbit, time_step=10)
        self.assertIsInstance(glint.attitude, GlintAttitude)
        self.assertEqual(glint.attitude.max_off_nadir_angle, 51.0)
        self.assertEqual(glint.attitude.velocity_frame, VelocityFrame.EARTH_FIXED)
        nadir = entry.build(configuration="nadir", orbit=orbit, time_step=10)
        self.assertNotIsInstance(nadir.attitude, GlintAttitude)

    def test_transfers_recorded_apart(self):
        """
        Test that results on sibling platforms are recorded as transfers, which no
        gate reads, and that the configurations transferred warn on those
        platforms.
        """
        transfers = {
            "atms": "noaa-21",
            "cris": "noaa-21",
            "viirs": "noaa-21",
            "mhs": "metop-b",
            "avhrr": "metop-b",
            "csar": "sentinel-1d",
            "modis": "aqua",
            "oli": "landsat-8",
            "poseidon4": "sentinel-6b",
        }
        recorded = {
            entry.id: {
                platform
                for record in entry.record.configurations.values()
                for platform in record.transfers
            }
            for entry in list_entries()
        }
        self.assertEqual(
            {k: v for k, v in recorded.items() if v},
            {k: {v} for k, v in transfers.items()},
        )
        for entry_id, platform in transfers.items():
            entry = get_entry(entry_id)
            for key, record in entry.record.configurations.items():
                if platform not in record.transfers:
                    continue
                self.assertNotIn(platform, entry.configurations[key].platforms, key)
                with self.assertWarns(UserWarning):
                    entry.build(
                        configuration=key,
                        orbit=GeneralPerturbationsOrbit.from_tle(GPM_TLE),
                        time_step=10,
                        platform=platform,
                    )

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
                    configuration="sample_size",
                    orbit=GeneralPerturbationsOrbit.from_tle(GPM_TLE),
                    time_step=1.875,
                ).limitations
            )
        )

    def test_abi_builds_without_attitude(self):
        """
        Test that an Earth-fixed configuration without other attitude fields
        builds as it was tested, without an attitude.
        """
        payload = get_entry("abi").build()
        self.assertIsNone(payload.attitude)
        self.assertIsInstance(payload.instruments[0], Instrument)


class TestSelection(unittest.TestCase):
    """
    Unit tests for selecting configurations by platform and grade.
    """

    def setUp(self):
        self.orbit = GeneralPerturbationsOrbit.from_tle(GPM_TLE)

    def test_platform_selects_its_configuration(self):
        """
        Test that the platform modeled selects the configuration tested on it,
        without a warning, and the default configuration otherwise.
        """
        viirs = get_entry("viirs")
        self.assertEqual(viirs.select_configuration(), "swath_noaa_20")
        self.assertEqual(
            viirs.select_configuration(platform="noaa-21"), "swath_noaa_21"
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            payload = viirs.build(orbit=self.orbit, time_step=10, platform="noaa-21")
        self.assertEqual(payload.references, ("viirs:swath_noaa_21",))
        self.assertFalse(payload.transferred)

    def test_name_selects_platform_variant(self):
        """
        Test that a configuration named without its platform selects the variant
        tested on the platform, and otherwise the one of the default's platform.
        """
        cris = get_entry("cris")
        self.assertEqual(
            cris.get_variants("swath"),
            ["swath_noaa_20", "swath_noaa_21", "swath_suomi_npp"],
        )
        self.assertEqual(cris.get_variants("swath_noaa_21"), ["swath_noaa_21"])
        self.assertEqual(get_entry("modis").get_variants("swath_terra"), [])
        self.assertEqual(
            cris.select_configuration("swath", platform="noaa-21"), "swath_noaa_21"
        )
        self.assertEqual(cris.select_configuration("pixels"), "pixels_noaa_20")
        self.assertEqual(
            get_entry("omps_nm").select_configuration("swath"), "swath_noaa_21"
        )
        with self.assertWarnsRegex(UserWarning, "metop-c"):
            payload = cris.build(
                configuration="swath",
                orbit=self.orbit,
                time_step=10,
                platform="metop-c",
            )
        self.assertEqual(payload.references, ("cris:swath_noaa_20",))
        with self.assertRaisesRegex(ValueError, "unknown configuration"):
            cris.select_configuration("beams")
        with self.assertRaisesRegex(ValueError, "timing"):
            cris.select_configuration("swath", grade="sample_size")

    def test_grade_selects_a_configuration(self):
        """
        Test that a grade selects a configuration for its use, preferring the
        platform's, and that the timing grade serves coverage.
        """
        atms = get_entry("atms")
        self.assertEqual(
            atms.select_configuration(grade="sample_size", platform="noaa-21"),
            "beams_noaa_21",
        )
        self.assertEqual(
            atms.select_configuration(grade="coverage", platform="noaa-21"),
            "swath_noaa_21",
        )
        self.assertEqual(
            get_entry("gmi").select_configuration(grade="coverage"), "feeds"
        )
        self.assertEqual(
            get_entry("modis").select_configuration(grade=Grade.SAMPLE_SIZE), "pixels"
        )

    def test_transfer_when_untested(self):
        """
        Test that, without a configuration of the grade tested on the platform, one
        tested on another is selected, with a warning.
        """
        modis = get_entry("modis")
        with self.assertWarnsRegex(UserWarning, "aqua"):
            payload = modis.build(grade="sample_size", platform="aqua")
        self.assertEqual(payload.references, ("modis:pixels",))
        self.assertTrue(payload.transferred)

    def test_grade_errors(self):
        """
        Test that a grade without a configuration, or a configuration of another
        grade, is an error.
        """
        with self.assertRaisesRegex(ValueError, "sample_size"):
            get_entry("csar").select_configuration(grade="sample_size")
        with self.assertRaisesRegex(ValueError, "timing"):
            get_entry("gmi").select_configuration("feeds", grade="sample_size")

    def test_analysis_index(self):
        """
        Test the index with which a payload is analyzed: `None` for a union, the
        configured index of alternatives, and 0 for a single instrument.
        """
        orbit = self.orbit
        self.assertIsNone(
            get_entry("ascat").build(orbit=orbit, time_step=10).instrument_index
        )
        self.assertEqual(
            get_entry("viirs").build(grade="sample_size").instrument_index, 2
        )
        self.assertEqual(
            get_entry("msi").build(orbit=orbit, time_step=10).instrument_index, 0
        )

    def test_to_satellite_needs_default_index(self):
        """
        Test that a satellite is built only of a configuration analyzed with the
        default index of 0, which a satellite of the payload gives otherwise.
        """
        orbit = self.orbit
        with self.assertRaisesRegex(ValueError, "instrument_index=None"):
            get_entry("tropomi").to_satellite(orbit, time_step=10)
        with self.assertRaisesRegex(ValueError, "instrument_index=1"):
            get_entry("atms").to_satellite(orbit, grade="sample_size")
        payload = get_entry("tropomi").build(orbit=orbit, time_step=10)
        satellite = payload.to_satellite(orbit, name="Sentinel-5P")
        self.assertEqual(len(satellite.instruments), 10)
        self.assertEqual(satellite.name, "Sentinel-5P")
        self.assertEqual(
            len(get_entry("gmi").to_satellite(orbit, time_step=10).instruments), 2
        )


class TestBuildPayload(unittest.TestCase):
    """
    Unit tests for the payloads of several catalog configurations.
    """

    def setUp(self):
        self.orbit = GeneralPerturbationsOrbit.from_tle(GPM_TLE)

    def test_carried_together(self):
        """
        Test that configurations with the same attitude are carried together,
        analyzed together, with only the analyzed instrument of alternatives.
        """
        payload = build_payload(
            ["viirs", "cris", "atms:beams_noaa_20"], orbit=self.orbit, time_step=10
        )
        self.assertEqual(
            payload.references,
            ("viirs:swath_noaa_20", "cris:swath_noaa_20", "atms:beams_noaa_20"),
        )
        beams = get_entry("atms").configurations["beams_noaa_20"]
        self.assertEqual(
            [i.name for i in payload.instruments][-1],
            beams.instruments[beams.instrument_index].parameters["name"],
        )
        self.assertEqual(len(payload.instruments), 3)
        self.assertIsNone(payload.instrument_index)
        self.assertEqual(payload.attitude.velocity_frame, VelocityFrame.INERTIAL)

    def test_one_reference(self):
        """
        Test that one instrument is analyzed alone, by an alias.
        """
        payload = build_payload("GPM Microwave Imager", orbit=self.orbit, time_step=10)
        self.assertEqual(payload.references, ("gmi:feeds",))
        self.assertEqual(len(payload.instruments), 1)
        self.assertEqual(payload.instrument_index, 0)

    def test_different_attitudes(self):
        """
        Test that configurations of different attitudes cannot be carried
        together.
        """
        with self.assertRaisesRegex(ValueError, "different attitudes"):
            build_payload(["viirs", "oco2:glint"], orbit=self.orbit, time_step=10)
        with self.assertRaises(ValueError):
            build_payload([], orbit=self.orbit)

    def test_jpss_payloads(self):
        """
        Test that the configurations of the instruments of each JPSS satellite
        share the body's attitude, so that the satellite carries them together,
        each with its own angles.
        """
        for platform, references in (
            (
                "noaa-20",
                ["atms:swath", "cris:swath", "viirs:swath", "omps_nm:swath", "ceres"],
            ),
            ("noaa-21", ["atms:swath", "cris:swath", "viirs:swath", "omps_nm:swath"]),
            (
                "suomi-npp",
                ["atms:swath", "cris:swath", "viirs:swath", "omps_nm:swath", "ceres"],
            ),
        ):
            payload = build_payload(
                references, orbit=self.orbit, time_step=10, platform=platform
            )
            self.assertEqual(len(payload.instruments), len(references), platform)
            self.assertEqual(
                payload.attitude.velocity_frame, VelocityFrame.INERTIAL, platform
            )
            self.assertEqual(
                (payload.attitude.yaw_angle, payload.attitude.pitch_angle),
                (0, 0),
                platform,
            )
        omps = build_payload(
            ["atms:swath", "omps_nm:swath"],
            orbit=self.orbit,
            time_step=10,
            platform="noaa-20",
        ).instruments[1]
        self.assertEqual(
            (omps.yaw_angle, omps.roll_angle, omps.pitch_angle),
            (-0.295144, -0.2294, 0.1859),
        )

    def test_names_grade_and_attitude(self):
        """
        Test that references may name configurations without their platform,
        that a grade selects the configurations, and that an attitude replaces
        the validated one of every configuration.
        """
        payload = build_payload(
            ["atms:swath", "cris:swath"],
            orbit=self.orbit,
            time_step=10,
            platform="suomi-npp",
        )
        self.assertEqual(
            payload.references, ("atms:swath_suomi_npp", "cris:swath_suomi_npp")
        )
        payload = build_payload(
            ["atms", "cris"],
            orbit=self.orbit,
            grade="sample_size",
            platform="suomi-npp",
        )
        self.assertEqual(
            payload.references, ("atms:beams_suomi_npp", "cris:pixels_suomi_npp")
        )
        with self.assertRaisesRegex(ValueError, "timing"):
            build_payload("cris:swath", orbit=self.orbit, grade="sample_size")
        attitude = Attitude(velocity_frame=VelocityFrame.INERTIAL, roll_angle=1)
        payload = build_payload(
            ["viirs", "cris"], orbit=self.orbit, time_step=10, attitude=attitude
        )
        self.assertEqual(payload.attitude, attitude)
        with self.assertRaisesRegex(ValueError, "velocity frame"):
            build_payload("viirs", orbit=self.orbit, time_step=10, attitude=Attitude())

    def test_tilt_schedule_needs_period(self):
        """
        Test that a configuration whose tilt is scheduled by day is built for the
        period of the analysis.
        """
        orbit = GeneralPerturbationsOrbit.from_tle(OCO2_TLE)
        with self.assertRaisesRegex(ValueError, "start and end"):
            build_payload("oci", orbit=orbit, time_step=10)
        start = datetime(2026, 11, 10, tzinfo=timezone.utc)
        payload = build_payload(
            "oci", orbit=orbit, time_step=10, start=start, end=start + timedelta(days=1)
        )
        self.assertIsInstance(payload.attitude, AttitudeSchedule)


if __name__ == "__main__":
    unittest.main()
