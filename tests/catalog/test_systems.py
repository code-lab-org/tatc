"""
Unit tests for loading catalog systems, building their platforms, and the
standard that validates them.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import copy
import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
import yaml
from pydantic import ValidationError

from tatc.catalog import (
    RadarNetwork,
    SatelliteSystem,
    SurfaceNetwork,
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
from tatc.catalog.system_schemas import describe_orbit, parse_system
from tatc.catalog.system_standard import (
    COMPATIBLE_SYSTEM_VERSIONS,
    SystemClass,
    SystemGrade,
    get_required_system_grades,
    get_system_gates,
    serves_system,
)
from shapely.geometry import box

from tatc.schemas import (
    CircularOrbit,
    InSituInstrument,
    GeosynchronousOrbit,
    RepeatGroundTrackOrbit,
    SunSynchronousOrbit,
)
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
                0,
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
                0,
                [
                    "horizontal_error_p95",
                    "membership_precision",
                    "fitted values are not tested on data at least 30 days from their fitting data",
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
            "plane": ("partial", 32, ["not tested on aqua"]),
        },
        [],
    ),
    "gcom_w": (
        "validated",
        {
            "plane": ("validated", 35, []),
            "ground_track": ("validated", 35, []),
        },
        [],
    ),
    "goes": (
        "partial",
        {
            "slot": ("partial", 48, ["not tested on goes-18"]),
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
    "himawari": (
        "nominal",
        {
            "slot": (
                "nominal",
                None,
                [
                    "slot_longitude_error_p95",
                    "inclination_error_p95",
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
        "validated",
        {
            "plane": ("validated", 34, []),
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
        "validated",
        {
            "plane": ("validated", 36, []),
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
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "nexrad": (
        "partial",
        {
            "precipitation": ("partial", 30, ["max_elevation_error_p95"]),
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
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "sentinel_1": (
        "partial",
        {
            "plane": ("partial", 53, ["not tested on sentinel-1d"]),
            "ground_track": ("validated", 52, []),
        },
        [],
    ),
    "sentinel_2": (
        "validated",
        {
            "plane": ("validated", 55, []),
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
                    "no comparison with reference data",
                ],
            ),
        },
        [],
    ),
    "smap": (
        "validated",
        {
            "plane": ("validated", 36, []),
        },
        [],
    ),
    "swot": (
        "validated",
        {
            "plane": ("validated", 35, []),
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
                0,
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
            ["goes", "himawari", "intelsat_40e", "mtg"],
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

    def test_loaded_copies_independent(self):
        """
        Test that changing a loaded system does not change later ones.
        """
        system = get_system("jpss")
        system.spacecraft[0].__dict__["name"] = "changed"
        self.assertEqual(get_system("jpss").spacecraft[0].name, "NOAA-20")

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

    def test_classes(self):
        """
        Test that systems are parsed as the schemas of their classes.
        """
        self.assertIsInstance(get_system("jpss"), SatelliteSystem)
        self.assertIsInstance(get_system("nexrad"), RadarNetwork)
        self.assertIsInstance(get_system("uscrn"), SurfaceNetwork)
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
        self.assertEqual(summary["members"], 3)
        self.assertEqual(summary["configurations"], "plane, ground_track")
        uscrn = get_system("uscrn")
        names = uscrn.get_station_names()
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
        Test that a plane whose local time drifts is built with its inclination
        and the cataloged local time at the epoch.
        """
        system = get_system("eos")
        (aqua,) = system.build(epoch=EPOCH, payload=False, members=["aqua"])
        self.assertIsInstance(aqua.orbit, CircularOrbit)
        plane = system.configurations["plane"].orbits["aqua"]
        self.assertEqual(aqua.orbit.inclination, plane.inclination)
        reference = plane.model_copy(update={"inclination": None}).build(EPOCH)
        self.assertAlmostEqual(
            aqua.orbit.right_ascension_ascending_node,
            reference.get_right_ascension_ascending_node(),
        )

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

    def test_ground_track_includes_plane(self):
        """
        Test that the gates of a ground track include those of the plane.
        """
        for system_class in ("sun_synchronous", "inclined"):
            plane = get_system_gates(system_class, "plane")
            ground_track = get_system_gates(system_class, "ground_track")
            self.assertEqual(ground_track[: len(plane)], plane)

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
