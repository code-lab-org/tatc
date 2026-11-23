"""
Unit tests for the performance and resources of catalog entries and spacecraft.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import copy
import unittest
import warnings
from datetime import datetime, timezone

import yaml
from pydantic import ValidationError

from tatc.catalog import (
    InstrumentEntry,
    get_entry,
    get_entry_ids,
    get_system,
    list_systems,
)
from tatc.catalog import loader
from tatc.catalog.schemas import PerformanceSpecification, ResourceSpecification
from tatc.schemas import SunSynchronousOrbit

EPOCH = datetime(2026, 11, 1, tzinfo=timezone.utc)
ORBIT = SunSynchronousOrbit(altitude=700e3, equator_crossing_time="13:30")
SOURCE = {"kind": "nominal", "source": "test"}


def _load(entry_id):
    with open(loader._ENTRIES / f"{entry_id}.yml", "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


class TestCatalogValues(unittest.TestCase):
    """
    Tests the cataloged values and their sources.
    """

    def test_sources(self):
        """
        Tests that every published value names where it was read, by its URL.
        """
        provenances = []
        for entry_id in get_entry_ids():
            entry = get_entry(entry_id)
            for sourced in [entry.performance, entry.resources] + [
                getattr(c, field)
                for c in entry.configurations.values()
                for field in ("performance", "resources")
            ]:
                if sourced is not None:
                    provenances += list(sourced.provenance.values())
        for system in list_systems():
            for spacecraft in getattr(system, "spacecraft", ()):
                provenances += list(spacecraft.resources.provenance.values())
        self.assertTrue(provenances)
        for provenance in provenances:
            if provenance.kind == "nominal":
                self.assertIn("(https://", provenance.source)
            else:
                self.assertEqual(provenance.kind, "derived")

    def test_values_not_in_digests(self):
        """
        Tests that performance and resources do not change the digest with which
        a validation record identifies a configuration.
        """
        entry = get_entry("viirs")
        data = _load("viirs")
        del data["performance"], data["resources"]
        stripped = InstrumentEntry(**data)
        for key, configuration in entry.configurations.items():
            self.assertEqual(configuration.digest, stripped.configurations[key].digest)


class TestBuild(unittest.TestCase):
    """
    Tests the performance and resources of built instruments.
    """

    def test_carrier(self):
        """
        Tests that every instrument has the performance, and only the one
        analyzed by default has the resources.
        """
        payload = get_entry("viirs").build(
            configuration="zones_noaa_20", orbit=ORBIT, time_step=10
        )
        self.assertEqual(payload.instrument_index, 2)
        for i, instrument in enumerate(payload.instruments):
            self.assertEqual(instrument.performance.ground_sample_distance, 750)
            self.assertEqual(instrument.resources is not None, i == 2)

    def test_union(self):
        """
        Tests that the first part of one view carries the resources.
        """
        payload = get_entry("karin").build(orbit=ORBIT, time_step=10)
        self.assertIsNone(payload.instrument_index)
        self.assertGreater(len(payload.instruments), 1)
        self.assertEqual(
            [i.resources is not None for i in payload.instruments],
            [True] + [False] * (len(payload.instruments) - 1),
        )
        satellite = payload.to_satellite(ORBIT, name="test")
        self.assertEqual(satellite.get_payload_resources().mass, 300)

    def test_override(self):
        """
        Tests that a configuration's performance replaces the entry's.
        """
        entry = get_entry("csar")
        self.assertEqual(entry.get_performance("iw").ground_sample_distance, 20)
        self.assertEqual(entry.get_performance("ew").ground_sample_distance, 40)
        self.assertEqual(entry.get_resources("ew").mass, 945)
        (instrument,) = entry.build(
            configuration="ew", orbit=ORBIT, time_step=10
        ).instruments
        self.assertEqual(instrument.performance.ground_sample_distance, 40)

    def test_built_copies(self):
        """
        Tests that built instruments do not share their values.
        """
        entry = get_entry("oli")
        first, second = (
            entry.build(orbit=ORBIT, time_step=10).instruments[0] for _ in range(2)
        )
        self.assertIsNot(first.resources, second.resources)

    def test_payload(self):
        """
        Tests the totals of a payload of several entries.
        """
        payload = loader.build_payload(
            ["viirs:zones", "atms:swath"], orbit=ORBIT, platform="noaa-20", time_step=10
        )
        total = payload.to_satellite(ORBIT, name="test").get_payload_resources()
        self.assertAlmostEqual(total.mass, 275 + 75.4)


class TestSpecifications(unittest.TestCase):
    """
    Tests the checks of sourced values.
    """

    def test_provenance(self):
        """
        Tests that each value given needs a provenance, and only those.
        """
        with self.assertRaises(ValidationError):
            ResourceSpecification(mass=10, provenance={})
        with self.assertRaises(ValidationError):
            ResourceSpecification(mass=10, provenance={"mass": SOURCE, "power": SOURCE})
        with self.assertRaises(ValidationError):
            ResourceSpecification(provenance={})
        ResourceSpecification(mass=10, provenance={"mass": SOURCE})

    def test_exponents(self):
        """
        Tests that a value that scales needs its exponent.
        """
        with self.assertRaises(ValidationError):
            PerformanceSpecification(
                reference_altitude=1,
                ground_sample_distance=1,
                provenance={
                    "reference_altitude": SOURCE,
                    "ground_sample_distance": SOURCE,
                },
            )
        with self.assertRaises(ValidationError):
            PerformanceSpecification(
                reference_altitude=1, provenance={"reference_altitude": SOURCE}
            )
        performance = PerformanceSpecification(
            reference_altitude=1,
            signal_to_noise_exponent=2,
            provenance={
                "reference_altitude": SOURCE,
                "signal_to_noise_exponent": SOURCE,
            },
        ).to_performance()
        self.assertEqual(performance.get_noise_scale(2), 4)

    def test_not_parameters(self):
        """
        Tests that an instrument's parameters cannot give its resources.
        """
        broken = copy.deepcopy(_load("oli"))
        instrument = next(iter(broken["configurations"].values()))["instruments"][0]
        instrument["parameters"]["resources"] = {"mass": 1}
        instrument["provenance"]["resources"] = SOURCE
        with self.assertRaisesRegex(ValidationError, "not to an instrument"):
            InstrumentEntry(**broken)


class TestSpacecraft(unittest.TestCase):
    """
    Tests the resources of the catalog's spacecraft.
    """

    def test_build(self):
        """
        Tests that built satellites have their spacecraft's resources.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # the configuration is nominal
            satellites = get_system("galileo").build(
                members=["galileo-gsat0101", "galileo-gsat0203"], epoch=EPOCH
            )
        self.assertEqual([s.resources.launch_mass for s in satellites], [700, 730])
        self.assertEqual(satellites[1].resources.power, 1900)
