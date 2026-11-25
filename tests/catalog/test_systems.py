"""
Unit tests for loading catalog systems, building their platforms, and the
standard that validates them.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import copy
import unittest
import warnings
from datetime import datetime, timedelta, timezone

import numpy as np
import yaml
from pydantic import ValidationError

from tatc.catalog import (
    GroundNetwork,
    RadarNetwork,
    SatelliteSystem,
    SurfaceNetwork,
    build_payload,
    clear_cache,
    get_entry,
    get_entry_ids,
    get_spacecraft,
    get_spacecraft_system,
    get_system,
    get_system_ids,
    list_entries,
    list_systems,
    loader,
)
from tatc.catalog.records import ConfigurationRecord
from tatc.catalog.system_orbits import SLOT_EPOCH, describe_orbit
from tatc.catalog.system_schemas import SystemEntry, parse_system
from tatc.catalog.system_standard import (
    COMPATIBLE_SYSTEM_VERSIONS,
    SystemClass,
    SystemGrade,
    get_orbit_type,
    get_required_system_grades,
    get_system_gates,
    serves_system,
)
from shapely.geometry import box

from tatc.schemas import (
    CircularOrbit,
    GroundStation,
    InSituInstrument,
    GeosynchronousOrbit,
    Relay,
    RepeatGroundTrackOrbit,
    SunSynchronousOrbit,
)
from tatc.analysis import collect_relay_contacts
from tatc.constants import timescale
from tatc.utils import compute_geoid_height

# the validation status of each system: its tier, and for each configuration its
# tier, the days between its fitting data and its test (`None` without fitted
# values or a test), and its gaps (by metric, for its gates), and the gaps of
# the system
EXPECTED_SYSTEM_STATUS = {
    "asos": (
        "partial",
        {
            "stations": (
                "partial",
                1,
                [
                    "membership_precision",
                    "fitted values are not tested on data at least 30 days from their fitting data",
                ],
            ),
        },
        [],
    ),
    "awos": (
        "partial",
        {
            "stations": (
                "partial",
                1,
                [
                    "horizontal_error_p95",
                    "reference_matched",
                    "membership_precision",
                    "fitted values are not tested on data at least 30 days from their fitting data",
                ],
            ),
        },
        [],
    ),
    "beidou": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "relative_phase_error_p95",
                    "relative_node_error_p95",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "cosmic_2": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "dmsp": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "local_time_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "earthcare": (
        "validated",
        {
            "plane": ("validated", 31, []),
            "ground_track": ("validated", 30, []),
        },
        [],
    ),
    "eos": (
        "partial",
        {
            "plane": (
                "partial",
                32,
                ["test_sets", "test_span_days", "not tested on aqua"],
            ),
        },
        [],
    ),
    "galileo": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "relative_phase_error_p95",
                    "relative_node_error_p95",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "gcom_w": (
        "partial",
        {
            "plane": ("partial", 35, ["test_sets"]),
            "ground_track": ("partial", 35, ["test_sets"]),
        },
        [],
    ),
    "glonass": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "relative_phase_error_p95",
                    "relative_node_error_p95",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "goes": (
        "partial",
        {
            "slot": (
                "partial",
                48,
                ["test_sets", "test_span_days", "not tested on goes-18"],
            ),
        },
        [],
    ),
    "gpm": (
        "validated",
        {
            "plane": ("validated", 81, []),
        },
        [],
    ),
    "gps": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "relative_phase_error_p95",
                    "relative_node_error_p95",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "himawari": (
        "nominal",
        {
            "slot": (
                "nominal",
                None,
                [
                    "slot_longitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "icesat_2": (
        "validated",
        {
            "plane": ("validated", 30, []),
        },
        [],
    ),
    "intelsat_40e": (
        "nominal",
        {
            "slot": (
                "nominal",
                None,
                [
                    "slot_longitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "iss": (
        "validated",
        {
            "plane": ("validated", 30, []),
        },
        [],
    ),
    "jpss": (
        "partial",
        {
            "plane": ("validated", 30, []),
            "ground_track": (
                "partial",
                33,
                ["node_longitude_error_p95", "node_time_error_p95"],
            ),
        },
        [],
    ),
    "landsat": (
        "partial",
        {
            "plane": ("partial", 34, ["test_sets"]),
        },
        [],
    ),
    "metop": (
        "validated",
        {
            "plane": ("validated", 30, []),
        },
        [],
    ),
    "metop_sg": (
        "partial",
        {
            "plane": ("partial", 36, ["test_sets", "test_span_days"]),
        },
        [],
    ),
    "mtg": (
        "nominal",
        {
            "slot": (
                "nominal",
                None,
                [
                    "slot_longitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "nasa_space_network": (
        "nominal",
        {
            "tdrs": (
                "nominal",
                None,
                [
                    "horizontal_error_p95",
                    "reference_matched",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "nexrad": (
        "validated",
        {
            "precipitation": ("validated", 30, []),
            "clear_air": ("validated", 30, []),
        },
        [],
    ),
    "nisar": (
        "validated",
        {
            "plane": ("validated", 30, []),
            "ground_track": ("validated", 30, []),
        },
        [],
    ),
    "noaa_ground_stations": (
        "nominal",
        {
            key: (
                "nominal",
                None,
                [
                    "horizontal_error_p95",
                    "reference_matched",
                    "no comparison with reference data",
                ],
            )
            for key in ("jpss", "poes", "goes")
        },
        [],
    ),
    "oco_2": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "local_time_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "pace": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "local_time_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "sentinel_1": (
        "partial",
        {
            "plane": (
                "partial",
                53,
                ["relative_phase_error_p95", "not tested on sentinel-1d"],
            ),
            "ground_track": ("validated", 52, []),
        },
        [],
    ),
    "sentinel_2": (
        "partial",
        {
            "plane": ("partial", 55, ["not tested on sentinel-2a"]),
            "ground_track": ("validated", 54, []),
        },
        [],
    ),
    "sentinel_3": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "local_time_error_p95",
                    "test_sets",
                    "test_span_days",
                    "relative_phase_error_p95",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "sentinel_5p": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "local_time_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "sentinel_6": (
        "nominal",
        {
            "plane": (
                "nominal",
                None,
                [
                    "altitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "relative_phase_error_p95",
                    "relative_node_error_p95",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "smap": (
        "partial",
        {
            "plane": ("partial", 36, ["test_sets"]),
        },
        [],
    ),
    "swot": (
        "partial",
        {
            "plane": ("partial", 35, ["test_sets"]),
        },
        [],
    ),
    "tdrs": (
        "nominal",
        {
            "slot": (
                "nominal",
                None,
                [
                    "slot_longitude_error_p95",
                    "inclination_error_p95",
                    "test_sets",
                    "test_span_days",
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "timed": (
        "validated",
        {
            "plane": ("validated", 47, []),
        },
        [],
    ),
    "trmm": (
        "validated",
        {
            "plane": ("validated", 30, []),
        },
        [],
    ),
    "tropics": (
        "partial",
        {
            "plane": ("partial", 30, ["altitude_error_p95"]),
        },
        [],
    ),
    "uscrn": (
        "partial",
        {
            "stations": (
                "partial",
                1,
                [
                    "membership_precision",
                    "fitted values are not tested on data at least 30 days from their fitting data",
                ],
            ),
        },
        [],
    ),
}

EPOCH = datetime(2026, 11, 1, tzinfo=timezone.utc)


def _system_data(system_id: str) -> dict:
    with open(
        loader._SYSTEMS / f"{system_id}.yml",  # pylint: disable=protected-access
        "r",
        encoding="utf-8",
    ) as f:
        return yaml.safe_load(f)


def _get_node(orbit, moment: datetime) -> float:
    """
    Gets the right ascension of ascending node (degrees, TEME) of an orbit at
    a time, propagated with SGP4.
    """
    satrec = orbit.to_gp_orbit().elements[0].to_satrec()
    julian = timescale.from_datetime(moment).ut1
    _, position, velocity = satrec.sgp4(int(julian), julian - int(julian))
    node = np.cross([0.0, 0.0, 1.0], np.cross(position, velocity))
    return float(np.degrees(np.arctan2(node[1], node[0])) % 360)


def _summarize_gaps(validation) -> list[str]:
    """
    Gets the gaps of a configuration, with those of its gates as their metrics.
    """
    metrics = [result.gate.metric for result in validation.gates]
    return [
        next((m for m in metrics if text.startswith(f"{m}:")), text)
        for text in validation.gaps
    ]


class TestLoading(unittest.TestCase):
    """
    Unit tests for loading systems.
    """

    def test_systems(self):
        """
        Test that the cataloged systems load, with current records.
        """
        self.assertIn("jpss", get_system_ids())
        for system in list_systems():
            self.assertEqual(get_system(system.id).id, system.id)
            if system.record is None:
                continue
            for key, record in system.record.configurations.items():
                self.assertEqual(
                    record.digest, system.get_digest(key), (system.id, key)
                )
                self.assertIn(record.standard_version, COMPATIBLE_SYSTEM_VERSIONS)

    def test_alias_lookup(self):
        """
        Test lookup by alias and name, ignoring case.
        """
        self.assertEqual(get_system("JPSS").id, "jpss")
        self.assertEqual(get_system("joint polar satellite system").id, "jpss")
        with self.assertRaisesRegex(KeyError, "jpss"):
            get_system("unknown")

    def test_list_filters(self):
        """
        Test filtering by class and tier.
        """
        self.assertEqual(
            [s.id for s in list_systems(system_class="geostationary")],
            ["goes", "himawari", "intelsat_40e", "mtg", "tdrs"],
        )
        self.assertEqual(
            [s.id for s in list_systems(tier="validated")],
            [
                k
                for k, (tier, _, _) in EXPECTED_SYSTEM_STATUS.items()
                if tier == "validated"
            ],
        )
        with self.assertRaises(ValueError):
            list_systems(system_class="balloon")

    def test_list_by_grade(self):
        """
        Test that systems are selected by their tier for a use, the lowest of
        their configurations that serve it: a system partial for the planes of its
        orbits may be validated for its ground track, a ground track configuration
        serves the plane grade, and a grade selects the systems with a
        configuration for it.
        """
        sentinel_1 = get_system("sentinel_1")
        self.assertEqual(sentinel_1.get_tier().value, "partial")
        self.assertEqual(sentinel_1.get_tier("ground_track").value, "validated")
        # the plane grade takes the plane configuration, partial, as well
        self.assertEqual(sentinel_1.get_tier("plane").value, "partial")
        self.assertIsNone(sentinel_1.get_tier("siting"))
        validated = [s.id for s in list_systems(tier="validated")]
        for_tracks = [
            s.id for s in list_systems(tier="validated", grade="ground_track")
        ]
        self.assertIn("sentinel_1", for_tracks)
        self.assertNotIn("sentinel_1", validated)
        self.assertNotIn("gpm", for_tracks)
        self.assertEqual([s.id for s in list_systems(grade="coverage")], ["nexrad"])
        self.assertEqual(
            [s.id for s in list_systems(satellite=False, grade="siting")],
            [s.id for s in list_systems(satellite=False)],
        )
        with self.assertRaises(ValueError):
            list_systems(grade="pixels")

    def test_loaded_copies_independent(self):
        """
        Test that changing a loaded system does not change later ones.
        """
        system = get_system("jpss")
        system.spacecraft[0].__dict__["name"] = "changed"
        self.assertEqual(get_system("jpss").spacecraft[0].name, "NOAA-20")

    def test_dumped_systems_parse(self):
        """
        Test that every system, dumped with its defaults (as its record is checked
        when written), parses again: values dumped as `None`, such as a
        spacecraft's unknown storage, need no provenance.
        """
        for system_id in get_system_ids():
            with self.subTest(system_id):
                parse_system(get_system(system_id).model_dump(mode="json"))

    def test_spacecraft_belong_to_one_system(self):
        """
        Test that every spacecraft belongs to one system, and that every platform
        of an instrument entry is a cataloged spacecraft.
        """
        self.assertEqual(get_spacecraft("noaa-20").norad_id, 43013)
        self.assertEqual(get_spacecraft_system("noaa-20").id, "jpss")
        with self.assertRaises(KeyError):
            get_spacecraft("unknown")
        for entry_id in get_entry_ids():
            entry = get_entry(entry_id)
            self.assertEqual(
                [s.id for s in entry.get_spacecraft()], list(entry.platforms)
            )

    def test_entries_of_spacecraft(self):
        """
        Test that the instrument entries of a spacecraft are those that list it,
        and that its payload omits the soundings of limb sounders and occultation
        receivers.
        """
        self.assertEqual(
            [e.id for e in list_entries(platform="noaa-20")],
            ["atms", "ceres", "cris", "omps_nm", "viirs"],
        )
        self.assertIn("gras", [e.id for e in list_entries(platform="metop-b")])
        self.assertNotIn("gras", get_system("metop").get_payload_references("metop-b"))

    def test_member_ids(self):
        """
        Test that members are found by identifier or name, and that stations
        are built by either.
        """
        self.assertEqual(get_system("jpss").get_member_id("NOAA-21"), "noaa-21")
        uscrn = get_system("uscrn")
        identifier, name = next(iter(uscrn.get_member_names().items()))
        self.assertEqual(uscrn.get_member_id(name.upper()), identifier)
        (station,) = uscrn.build(members=name)
        self.assertEqual(station.name, identifier)
        with self.assertRaises(KeyError):
            uscrn.get_member_id("nowhere")

    def test_clear_cache(self):
        """
        Test that clearing the cache loads the catalog again.
        """
        before = get_system("jpss")
        clear_cache()
        self.assertEqual(get_system("jpss"), before)

    def test_spacecraft_by_name(self):
        """
        Test that spacecraft are found by identifier or name, ignoring case, and
        that satellite systems or networks are listed apart.
        """
        self.assertEqual(get_spacecraft("NOAA-20").id, "noaa-20")
        self.assertEqual(get_spacecraft_system("aqua").id, "eos")
        self.assertEqual(get_system("jpss").get_spacecraft("Suomi NPP").id, "suomi-npp")
        satellites = list_systems(satellite=True)
        networks = list_systems(satellite=False)
        self.assertTrue(all(s.system_class.is_satellite() for s in satellites))
        self.assertEqual(len(satellites) + len(networks), len(get_system_ids()))

    def test_classes(self):
        """
        Test that systems are parsed as the schemas of their classes, not as
        the base class.
        """
        with self.assertRaisesRegex(TypeError, "parse_system"):
            SystemEntry(**_system_data("jpss"))
        self.assertIsInstance(get_system("jpss"), SatelliteSystem)
        self.assertIsInstance(get_system("nexrad"), RadarNetwork)
        self.assertIsInstance(get_system("uscrn"), SurfaceNetwork)
        self.assertIsInstance(get_system("noaa_ground_stations"), GroundNetwork)
        data = _system_data("nexrad")
        data["system_class"] = "surface_network"
        with self.assertRaises(ValidationError):
            parse_system(data)

    def test_summaries(self):
        """
        Test the summaries of systems, the names of stations, and the
        descriptions of orbits.
        """
        summary = get_system("jpss").summarize()
        self.assertEqual((summary["members"], summary["active"]), (3, 2))
        self.assertEqual(summary["configurations"], "plane, ground_track")
        self.assertEqual(get_system("tropics").summarize()["active"], 0)
        self.assertEqual(get_system("jpss").get_member_names()["noaa-20"], "NOAA-20")
        uscrn = get_system("uscrn")
        names = uscrn.get_member_names()
        self.assertEqual(len(names), uscrn.summarize()["members"])
        self.assertEqual(list(names), uscrn.get_members())
        orbits = get_system("jpss").configurations
        self.assertTrue(
            describe_orbit(orbits["plane"].orbits["noaa-20"]).startswith(
                "Sun-synchronous, 834.10 km"
            )
        )
        self.assertTrue(
            describe_orbit(orbits["ground_track"].orbits["noaa-20"]).startswith(
                "Repeat of 227 orbits in 16 days"
            )
        )


class TestSystemStatus(unittest.TestCase):
    """
    Unit tests for the status of the systems, which must not be overstated.
    """

    def test_status(self):
        """
        Test the status of each system against the table of expected status, so
        that no change to a system, its record, or the standard changes a status
        unnoticed.
        """
        self.assertEqual(sorted(EXPECTED_SYSTEM_STATUS), get_system_ids())
        for system_id, (tier, configurations, gaps) in EXPECTED_SYSTEM_STATUS.items():
            with self.subTest(system_id):
                validation = get_system(system_id).validation
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


class TestStrictFiles(unittest.TestCase):
    """
    Unit tests that errors in hand-written files are rejected.
    """

    def setUp(self):
        self.data = _system_data("jpss")

    def _assert_invalid(self, message: str = ""):
        with self.assertRaisesRegex(ValidationError, message):
            parse_system(self.data)

    def test_baseline_valid(self):
        """
        Test that the unmodified data is valid.
        """
        parse_system(self.data)

    def test_unknown_field(self):
        """
        Test that an unknown field is rejected.
        """
        self.data["tier"] = "validated"
        self._assert_invalid()

    def test_missing_provenance(self):
        """
        Test that each value of each orbit needs provenance.
        """
        del self.data["configurations"]["plane"]["provenance"]["noaa-20"][
            "mean_altitude"
        ]
        self._assert_invalid("provenance")

    def test_unknown_spacecraft(self):
        """
        Test that orbits are of the system's spacecraft.
        """
        orbits = self.data["configurations"]["plane"]["orbits"]
        orbits["noaa-22"] = orbits["noaa-20"]
        self.data["configurations"]["plane"]["provenance"]["noaa-22"] = self.data[
            "configurations"
        ]["plane"]["provenance"]["noaa-20"]
        self._assert_invalid("unknown spacecraft")

    def test_orbit_type_of_grade(self):
        """
        Test that an orbit at the plane grade is not a repeat ground track, and
        that a sun-synchronous system has no geostationary grade.
        """
        configurations = self.data["configurations"]
        configurations["plane"]["orbits"]["noaa-20"] = configurations["ground_track"][
            "orbits"
        ]["noaa-20"]
        configurations["plane"]["provenance"]["noaa-20"] = configurations[
            "ground_track"
        ]["provenance"]["noaa-20"]
        self._assert_invalid("of type sso")

    def test_inclined_repeat_in_sun_synchronous_system(self):
        """
        Test that a sun-synchronous system's repeat orbits have no inclination.
        """
        configuration = self.data["configurations"]["ground_track"]
        configuration["orbits"]["noaa-20"]["inclination"] = 98.7
        configuration["provenance"]["noaa-20"]["inclination"] = {
            "kind": "nominal",
            "source": "x",
        }
        self._assert_invalid("inclination")

    def test_repeat_revolutions(self):
        """
        Test that a repeat ground track's revolutions must fit its days.
        """
        self.data["configurations"]["ground_track"]["orbits"]["noaa-20"][
            "revolutions"
        ] = 1
        self._assert_invalid()

    def test_phases_of_every_orbit_or_none(self):
        """
        Test that a configuration gives the phases of every orbit or of none,
        each with its provenance.
        """
        plane = self.data["configurations"]["plane"]
        plane["orbits"]["noaa-20"]["phase"] = 90
        self._assert_invalid("provenance")
        plane["provenance"]["noaa-20"]["phase"] = {"kind": "nominal", "source": "x"}
        self._assert_invalid("every orbit or for none")

    def test_platform_periods_of_every_platform(self):
        """
        Test that a record gives the period of every platform tested, or of
        none.
        """
        record = get_system("jpss").record.configurations["plane"]
        data = record.model_dump(mode="json")
        del data["platform_periods"]["noaa-21"]
        with self.assertRaisesRegex(ValidationError, "every platform"):
            ConfigurationRecord(**data)
        data["platform_periods"] = {}
        ConfigurationRecord(**data)

    def test_snapshot_date(self):
        """
        Test that a measured table's provenance has the date of its snapshot.
        """
        self.data = _system_data("nexrad")
        self.data["stations"]["as_of"] = "2026-10-16"
        self._assert_invalid("snapshot")

    def test_patterns_in_digest(self):
        """
        Test that a radar configuration's volume coverage patterns, which select
        the volumes of its tests, are in its digest.
        """
        nexrad = get_system("nexrad")
        data = _system_data("nexrad")
        data["configurations"]["clear_air"]["volume_coverage_patterns"] = [31]
        self.assertNotEqual(
            parse_system(data).get_digest("clear_air"),
            nexrad.get_digest("clear_air"),
        )

    def test_inclined_phase_needs_node(self):
        """
        Test that an inclined plane's relative phase comes with its node.
        """
        self.data = _system_data("sentinel_6")
        configuration = self.data["configurations"]["plane"]
        for spacecraft in configuration["orbits"]:
            del configuration["orbits"][spacecraft]["right_ascension_ascending_node"]
            del configuration["provenance"][spacecraft][
                "right_ascension_ascending_node"
            ]
        self._assert_invalid("right_ascension_ascending_node")

    def test_finite_metrics(self):
        """
        Test that a metric's value is finite.
        """
        record = get_system("landsat").record.model_dump(mode="json")
        record["configurations"]["plane"]["metrics"][0]["value"] = float("nan")
        self.data = {**_system_data("landsat"), "record": record}
        self._assert_invalid("finite")

    def test_drifting_plane_needs_reference_time(self):
        """
        Test that a plane whose local time drifts has a reference time.
        """
        self.data = _system_data("eos")
        del self.data["configurations"]["plane"]["orbits"]["aqua"]["reference_time"]
        self._assert_invalid("reference_time")

    def test_networks_have_stations(self):
        """
        Test that a satellite system has no stations, and a network no
        spacecraft.
        """
        self.data["stations"] = {
            "file": "x.csv",
            "as_of": "2026-01-01",
            "provenance": {"kind": "nominal", "source": "x"},
        }
        self._assert_invalid()

    def test_record_of_other_system(self):
        """
        Test that a record names the system and its spacecraft.
        """
        record = copy.deepcopy(get_system("jpss").record.model_dump(mode="json"))
        record["entry_id"] = "metop"
        self.data["record"] = record
        self._assert_invalid("different system")
        record["entry_id"] = "jpss"
        record["configurations"]["plane"]["platforms"] = ["metop-b"]
        record["configurations"]["plane"]["platform_periods"] = {}
        self._assert_invalid("not of the system")


class TestBuild(unittest.TestCase):
    """
    Unit tests for building the platforms of systems.
    """

    def test_sun_synchronous_planes(self):
        """
        Test that planes are built at the epoch, with the cataloged altitude and
        local time, and only the operational spacecraft by default.
        """
        system = get_system("jpss")
        satellites = system.build(epoch=EPOCH, payload=False)
        self.assertEqual([s.name for s in satellites], ["NOAA-20", "NOAA-21"])
        orbit = satellites[0].orbit
        self.assertIsInstance(orbit, SunSynchronousOrbit)
        self.assertEqual(orbit.epoch, EPOCH)
        plane = system.configurations["plane"].orbits["noaa-20"]
        self.assertEqual(orbit.mean_altitude, plane.mean_altitude)
        retired = system.build(epoch=EPOCH, payload=False, members=["suomi-npp"])
        self.assertEqual([s.name for s in retired], ["Suomi NPP"])

    def test_constellation_phases(self):
        """
        Test that satellites sharing a plane keep the operator's relative phases,
        and that satellites without cataloged phases (and inclined planes
        without nodes) are spaced evenly among those built, so that none is
        built at the position of another.
        """
        phases = {
            s.name: s.orbit.true_anomaly
            for s in get_system("sentinel_2").build(epoch=EPOCH, payload=False)
        }
        self.assertEqual(
            phases, {"Sentinel-2A": 324, "Sentinel-2B": 0, "Sentinel-2C": 180}
        )
        satellites = get_system("cosmic_2").build(epoch=EPOCH, payload=False)
        spacing = 360 / len(satellites)
        for index, satellite in enumerate(satellites):
            self.assertAlmostEqual(satellite.orbit.true_anomaly, index * spacing)
            self.assertAlmostEqual(
                satellite.orbit.right_ascension_ascending_node, index * spacing
            )
        tandem = get_system("sentinel_6").build(epoch=EPOCH, payload=False)
        self.assertEqual({s.orbit.right_ascension_ascending_node for s in tandem}, {0})
        self.assertEqual([s.orbit.true_anomaly for s in tandem], [0, 358.4])
        (only,) = get_system("jpss").build(
            epoch=EPOCH, payload=False, members=["noaa-21"]
        )
        self.assertEqual(only.orbit.true_anomaly, 0)

    def test_drifting_planes(self):
        """
        Test that a plane whose local time drifts is built with its inclination,
        with the cataloged local time at its reference time, from which its node
        drifts at the J2 rate of its inclination to the epoch.
        """
        system = get_system("eos")
        plane = system.configurations["plane"].orbits["aqua"]
        sun_synchronous = plane.model_copy(
            update={"inclination": None, "reference_time": None}
        )
        nodes = {}
        for days in (0, 365):
            epoch = plane.reference_time + timedelta(days=days)
            (aqua,) = system.build(epoch=epoch, payload=False, members="aqua")
            self.assertIsInstance(aqua.orbit, CircularOrbit)
            self.assertEqual(aqua.orbit.inclination, plane.inclination)
            # the node less that of a sun-synchronous orbit of the same local time
            nodes[days] = (
                aqua.orbit.right_ascension_ascending_node
                - sun_synchronous.build(epoch).get_right_ascension_ascending_node()
            ) % 360
        self.assertAlmostEqual(nodes[0], 0, places=6)
        # Aqua's local time drifts later by more than half an hour a year
        self.assertGreater(nodes[365], 7.5)

    def test_navigation_systems(self):
        """
        Test that navigation satellites are built as transmitters, without
        cataloged instruments or warnings about them, and keep the planes of
        their constellations.
        """
        gps = get_system("gps")
        self.assertEqual(gps.role, "transmitting")
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            satellites = gps.build(epoch=EPOCH)
        self.assertEqual(len(satellites), len(gps.spacecraft))
        self.assertFalse(
            [w for w in caught if "instrument" in str(w.message)],
            [str(w.message) for w in caught],
        )
        with self.assertRaisesRegex(ValueError, "transmitters"):
            gps.build(epoch=EPOCH, payload=True)
        # GLONASS's satellites are in three planes, 120 degrees apart
        nodes = [
            s.orbit.right_ascension_ascending_node
            for s in get_system("glonass").build(epoch=EPOCH)
        ]
        offsets = [abs((n - nodes[0] + 60) % 120 - 60) for n in nodes]
        self.assertLess(max(offsets), 10)

    def test_relay_systems(self):
        """
        Test that data relay satellites are built without cataloged instruments,
        that a relay is built of those neither retired nor in storage, and that
        other systems build no relays.
        """
        tdrs = get_system("tdrs")
        self.assertEqual(tdrs.role, "relaying")
        statuses = {s.id: s.status for s in tdrs.spacecraft}
        self.assertEqual(
            tdrs.get_default_members(),
            [k for k, v in statuses.items() if v not in ("retired", "storage")],
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            relay = tdrs.build_relay(min_grazing_altitude=100e3)
        self.assertFalse([w for w in caught if "instrument" in str(w.message)])
        self.assertIsInstance(relay, Relay)
        self.assertEqual(relay.name, "tdrs")
        self.assertEqual(relay.min_grazing_altitude, 100e3)
        self.assertEqual(
            [s.name for s in relay.satellites],
            [tdrs.get_spacecraft(m).name for m in tdrs.get_default_members()],
        )
        self.assertEqual(relay.ground_stations, [])
        spares = tdrs.build_relay(members=["tdrs-6", "tdrs-7"], name="spares")
        self.assertEqual([s.name for s in spares.satellites], ["TDRS 6", "TDRS 7"])
        with self.assertRaisesRegex(ValueError, "relays"):
            tdrs.build(payload=True)
        with self.assertRaisesRegex(ValueError, "not relays"):
            get_system("goes").build_relay()

    def test_occultation_transmitters(self):
        """
        Test that the constellations of radio occultation entries are the
        catalog's navigation systems, whose satellites are their transmitters.
        """
        for entry_id in get_entry_ids():
            for configuration in get_entry(entry_id).configurations.values():
                sounding = configuration.sounding
                if sounding is None or sounding.type != "occultation":
                    continue
                for constellation in sounding.constellations:
                    self.assertEqual(get_system(constellation).role, "transmitting")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            (metop_c,) = get_system("metop").build(
                epoch=EPOCH, payload=False, members="metop-c"
            )
            payload = build_payload("gras", orbit=metop_c.orbit, platform="metop-c")
            transmitters = payload.sounding.get_transmitters(EPOCH)
        self.assertEqual(len(transmitters), len(get_system("gps").spacecraft))
        with self.assertRaisesRegex(ValueError, "alone"):
            build_payload(
                ["gras", "iasi"], orbit=metop_c.orbit, platform="metop-c", time_step=10
            )

    def test_payload_given(self):
        """
        Test that given references are carried by the satellites that list them,
        and that a reference no satellite built carries is an error.
        """
        with self.assertWarnsRegex(UserWarning, "Aura carries none"):
            satellites = get_system("eos").build(
                epoch=EPOCH, payload="modis", time_step=10
            )
        self.assertEqual(
            {s.name: [i.name for i in s.instruments] for s in satellites},
            {"Terra": ["MODIS"], "Aqua": ["MODIS"], "Aura": ["Default"]},
        )
        with self.assertRaisesRegex(ValueError, "carries"):
            get_system("eos").build(epoch=EPOCH, payload=["viirs"], time_step=10)

    def test_payload_warnings(self):
        """
        Test that a satellite without cataloged instruments other than sounders
        warns, and that instruments that need a time step say so.
        """
        with self.assertWarnsRegex(UserWarning, "default instrument"):
            get_system("cosmic_2").build(epoch=EPOCH, time_step=10)
        with self.assertRaisesRegex(ValueError, r"build\(time_step=...\)"):
            get_system("jpss").build(epoch=EPOCH, payload=True)
        # without a time step, satellites carry no cataloged instruments
        with self.assertWarnsRegex(UserWarning, "time_step"):
            satellites = get_system("jpss").build(epoch=EPOCH)
        self.assertEqual([i.name for i in satellites[0].instruments], ["Default"])

    def test_unvalidated_configuration_warns(self):
        """
        Test that building a configuration that is not validated warns.
        """
        with self.assertWarnsRegex(UserWarning, "jpss:ground_track is partial"):
            get_system("jpss").build("ground_track", payload=False)

    def test_repeat_ground_track_builds(self):
        """
        Test that a cataloged repeat ground track builds the core orbit.
        """
        orbit = get_system("jpss").configurations["ground_track"].orbits["noaa-21"]
        built = orbit.build()
        self.assertIs(type(built), RepeatGroundTrackOrbit)
        self.assertEqual(built.ascending_node, orbit.ascending_node)

    def test_drift_independent_of_epoch(self):
        """
        Test that a drifting plane built at its reference time and propagated
        has the node of the plane built at a later epoch (SGP4's rate).
        """
        plane = get_system("eos").configurations["plane"].orbits["aqua"]
        later = plane.reference_time + timedelta(days=3 * 365)
        nodes = [
            _get_node(plane.build(epoch), later)
            for epoch in (plane.reference_time, later)
        ]
        self.assertAlmostEqual(nodes[0], nodes[1], delta=0.02)

    def test_ground_tracks_at_any_epoch(self):
        """
        Test that repeat ground tracks are built without an epoch, pinned to
        their nodes, and cross them years later.
        """
        system = get_system("jpss")
        satellites = system.build("ground_track", payload=False)
        orbit = satellites[0].orbit
        self.assertIsInstance(orbit, RepeatGroundTrackOrbit)
        node = orbit.ascending_node
        later = node.time + timedelta(days=16 * 50)
        nodes = orbit.to_gp_orbit().get_ascending_nodes(
            later - timedelta(minutes=30), later + timedelta(minutes=30)
        )
        self.assertAlmostEqual(nodes.longitude.iloc[0], node.longitude, delta=1e-4)

    def test_payloads(self):
        """
        Test that satellites carry their cataloged instruments, or those given.
        """
        system = get_system("jpss")
        satellites = system.build(epoch=EPOCH, time_step=10)
        self.assertEqual(
            [i.name for i in satellites[0].instruments],
            ["ATMS", "CERES", "CrIS", "OMPS Nadir Mapper", "VIIRS"],
        )
        (satellite,) = system.build(
            epoch=EPOCH, time_step=10, members=["noaa-21"], payload=["viirs"]
        )
        self.assertEqual([i.name for i in satellite.instruments], ["VIIRS"])

    def test_geostationary_slots(self):
        """
        Test that geostationary satellites are built at their longitudes.
        """
        satellites = get_system("goes").build(epoch=EPOCH, payload=False)
        self.assertIsInstance(satellites[0].orbit, GeosynchronousOrbit)
        (default,) = get_system("goes").build(payload=False, members="goes-19")
        self.assertEqual(default.orbit.epoch, SLOT_EPOCH)
        times = [EPOCH + timedelta(days=30, hours=h) for h in range(24)]
        position = satellites[0].orbit.to_gp_orbit().get_geographic_position(times)
        self.assertAlmostEqual(
            float(np.mean(position.longitude.degrees)),
            get_system("goes").configurations["slot"].orbits["goes-19"].longitude,
            delta=0.05,
        )

    def test_radar_stations(self):
        """
        Test that radar stations are built in the order of the members, at their
        ground elevations above the ellipsoid, with a station's own lowest
        elevation angle.
        """
        stations = get_system("nexrad").build(members=["KIWA", "KFSX"])
        self.assertEqual([s.name for s in stations], ["KIWA", "KFSX"])
        kiwa, kfsx = stations
        self.assertEqual(kiwa.min_elevation_angle, 0.5)
        self.assertLess(kfsx.min_elevation_angle, 0)
        self.assertEqual(kiwa.max_elevation_angle, 19.5)
        table = get_system("nexrad").stations.load().set_index("id")
        self.assertAlmostEqual(
            kiwa.elevation,
            table.ground["KIWA"] + compute_geoid_height(kiwa.longitude, kiwa.latitude),
        )
        self.assertEqual(kiwa.height, 19)
        with self.assertRaisesRegex(ValueError, "no stations"):
            get_system("nexrad").build(members=["KXXX"])
        (ktlx,) = get_system("nexrad").build(members="KTLX")
        self.assertEqual(ktlx.name, "KTLX")
        with self.assertWarnsRegex(UserWarning, "region"):
            self.assertEqual(get_system("nexrad").build(region=box(0, 0, 1, 1)), [])

    def test_surface_stations(self):
        """
        Test that surface stations are built within a region, with their in-situ
        instrument at the cataloged or given radius.
        """
        stations = get_system("uscrn").build(region=box(-98, 35.5, -96.5, 36.5))
        self.assertGreater(len(stations), 0)
        self.assertIsInstance(stations[0].instruments[0], InSituInstrument)
        self.assertEqual(stations[0].instruments[0].radius, 5000)
        (station,) = get_system("uscrn").build(members=[stations[0].name], radius=20000)
        self.assertEqual(station.instruments[0].radius, 20000)
        with self.assertRaises(TypeError):
            get_system("nexrad").build(
                radius=1000
            )  # pylint: disable=unexpected-keyword-arg

    def test_ground_stations(self):
        """
        Test that ground stations are built for a program (its stations, with
        the configuration's minimum elevation angle, or one given), that a
        configuration's stations are of the network, and that networks of
        ground stations summarize their role as receiving.
        """
        network = get_system("noaa_ground_stations")
        stations = network.build("jpss")
        self.assertEqual(
            [s.name for s in stations], ["FCDAS", "SvalSat", "McMurdo", "TrollSat"]
        )
        self.assertTrue(all(isinstance(s, GroundStation) for s in stations))
        self.assertTrue(all(s.min_elevation_angle == 5 for s in stations))
        # (ground elevations above the ellipsoid, with the geoid height)
        fcdas = network.stations.load().set_index("id").loc["FCDAS"]
        self.assertAlmostEqual(
            stations[0].elevation,
            fcdas.ground + float(compute_geoid_height(fcdas.longitude, fcdas.latitude)),
        )
        self.assertEqual([s.name for s in network.build("goes")], ["WCDAS", "FCDAS"])
        (svalsat,) = network.build(members="SvalSat", min_elevation_angle=10)
        self.assertEqual(svalsat.min_elevation_angle, 10)
        self.assertEqual(network.summarize()["role"], "receiving")
        data = _system_data("noaa_ground_stations")
        data["configurations"]["jpss"]["members"] = ["FCDAS", "Kiruna"]
        with self.assertRaisesRegex(ValidationError, "unknown stations"):
            parse_system(data)

    def test_relay_ground_terminals(self):
        """
        Test that the Space Network's terminals build ground stations for a
        relay, through which it relays only from the relay satellites in view
        of them.
        """
        terminals = get_system("nasa_space_network").build()
        self.assertEqual([t.name for t in terminals], ["WSC"])
        relay = get_system("tdrs").build_relay(ground_stations=terminals)
        self.assertEqual(relay.ground_stations, terminals)
        satellite = get_system("jpss").build(
            members="noaa-21", epoch=EPOCH, payload=False
        )
        contacts = collect_relay_contacts(
            relay, satellite, EPOCH, EPOCH + timedelta(hours=6)
        )
        self.assertGreater(len(contacts), 0)
        self.assertEqual(set(contacts.ground_station), {"WSC"})
        # (TDRS 13, near 11.5 deg W, is not in view of White Sands)
        self.assertNotIn("TDRS 13", set(contacts.relay_satellite))

    def test_build_errors(self):
        """
        Test that orbits at the plane grade need an epoch, members must be
        cataloged, and stations' arguments do not apply to satellites.
        """
        system = get_system("jpss")
        with self.assertRaisesRegex(ValueError, "epoch"):
            system.build(payload=False)
        with self.assertRaises(KeyError):
            system.build(epoch=EPOCH, payload=False, members=["metop-b"])
        with self.assertRaises(TypeError):
            system.build(  # pylint: disable=unexpected-keyword-arg
                epoch=EPOCH, payload=False, region=object()
            )
        with self.assertRaisesRegex(ValueError, "retired"):
            get_system("tropics").build(epoch=EPOCH, payload=False)
        with self.assertRaisesRegex(ValueError, "configuration"):
            system.build("other", epoch=EPOCH)

    def test_select_configuration(self):
        """
        Test selection by name and grade.
        """
        system = get_system("jpss")
        self.assertEqual(system.select_configuration(), "plane")
        self.assertEqual(
            system.select_configuration(grade="ground_track"), "ground_track"
        )
        self.assertEqual(system.select_configuration(grade="plane"), "plane")
        with self.assertRaisesRegex(ValueError, "does not serve"):
            system.select_configuration("plane", grade="ground_track")
        with self.assertRaisesRegex(ValueError, "no configuration"):
            get_system("metop").select_configuration(grade="ground_track")


class TestSystemStandard(unittest.TestCase):
    """
    Unit tests for the system standard.
    """

    def test_grades(self):
        """
        Test the grades of each class, and those each grade serves.
        """
        self.assertEqual(
            get_required_system_grades("sun_synchronous"), (SystemGrade.PLANE,)
        )
        self.assertEqual(
            get_required_system_grades("geostationary"), (SystemGrade.GROUND_TRACK,)
        )
        self.assertTrue(serves_system("ground_track", "plane"))
        self.assertFalse(serves_system("plane", "ground_track"))
        self.assertTrue(serves_system("coverage", "siting"))
        with self.assertRaisesRegex(ValueError, "no plane grade"):
            get_system_gates(SystemClass.GEOSTATIONARY, SystemGrade.PLANE)

    def test_metrics_of_every_spacecraft(self):
        """
        Test that a gate whose metric lacks a spacecraft's value is a gap.
        """
        data = _system_data("metop")
        record = get_system("metop").record.model_dump(mode="json")
        plane = record["configurations"]["plane"]
        plane["metrics"] = [
            m
            for m in plane["metrics"]
            if not (m["name"] == "local_time_error_p95" and m["scope"] == "metop-c")
        ]
        validation = parse_system({**data, "record": record}).validation
        self.assertIn(
            "local_time_error_p95: not tested on metop-c",
            validation.configurations["plane"].gaps,
        )
        self.assertNotEqual(validation.tier.value, "validated")

    def test_orbit_types(self):
        """
        Test the orbit types of the grades of satellite classes.
        """
        self.assertEqual(get_orbit_type("sun_synchronous", "plane"), "sso")
        self.assertEqual(
            get_orbit_type("geostationary", "ground_track"), "geosynchronous"
        )
        with self.assertRaises(ValueError):
            get_orbit_type("radar_network", "coverage")

    def test_relative_gates(self):
        """
        Test that configurations with cataloged relative phases (and nodes) are
        tested on them.
        """
        plane = get_system("sentinel_2").configurations["plane"]
        self.assertEqual(plane.get_relative_values(), ("phase",))
        metrics = [
            g.metric for g in get_system_gates("sun_synchronous", "plane", ("phase",))
        ]
        self.assertIn("relative_phase_error_p95", metrics)
        self.assertNotIn(
            "relative_phase_error_p95",
            [g.metric for g in get_system_gates("sun_synchronous", "plane")],
        )
        validation = get_system("sentinel_6").validation.configurations["plane"]
        self.assertIn(
            "relative_node_error_p95", [r.gate.metric for r in validation.gates]
        )

    def test_ground_track_includes_plane(self):
        """
        Test that the gates of a ground track include those of the plane.
        """
        for system_class in ("sun_synchronous", "inclined"):
            plane = get_system_gates(system_class, "plane")
            ground_track = get_system_gates(system_class, "ground_track")
            self.assertTrue(set(plane) <= set(ground_track))

    def test_untested_spacecraft(self):
        """
        Test that a configuration tested on some of its spacecraft is not
        validated.
        """
        data = _system_data("jpss")
        record = copy.deepcopy(get_system("jpss").record.model_dump(mode="json"))
        tested = record["configurations"]["ground_track"]
        tested["platforms"] = ["noaa-20"]
        del tested["platform_periods"]["noaa-21"]
        tested["metrics"] = [m for m in tested["metrics"] if m["scope"] == "noaa-20"]
        system = parse_system({**data, "record": record})
        validation = system.validation.configurations["ground_track"]
        self.assertIn("not tested on noaa-21", validation.gaps)
        self.assertNotEqual(validation.tier.value, "validated")

    def test_holdout_by_spacecraft(self):
        """
        Test that each spacecraft's fitted values are held out from its own test
        data, or from the record's period if the record does not give each
        platform's.
        """
        system = get_system("jpss")
        configuration = system.configurations["plane"]
        fitted = {
            spacecraft: values["mean_altitude"].period
            for spacecraft, values in configuration.provenance.items()
        }
        self.assertEqual(configuration.get_holdout_gap("2027-01-01", dict(fitted)), 0)
        self.assertEqual(configuration.get_holdout_gap(fitted["noaa-20"]), 0)
        later = {spacecraft: "2030-01-01" for spacecraft in fitted}
        self.assertGreater(configuration.get_holdout_gap("2000-01-01", later), 1000)

    def test_changed_configuration(self):
        """
        Test that a configuration changed after its test is not validated.
        """
        data = _system_data("jpss")
        data["configurations"]["plane"]["orbits"]["noaa-20"]["mean_altitude"] += 1000
        data["record"] = get_system("jpss").record.model_dump(mode="json")
        validation = parse_system(data).validation.configurations["plane"]
        self.assertIn("the configuration changed after it was tested", validation.gaps)


if __name__ == "__main__":
    unittest.main()
