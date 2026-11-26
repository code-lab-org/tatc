"""
Unit tests for the communications of catalog systems.
"""

import copy
import unittest
import warnings
from datetime import datetime, timezone

import yaml
from pydantic import ValidationError

from tatc.catalog import get_system, get_system_ids, loader
from tatc.catalog.schemas import (
    ReceiverSpecification,
    RelayServiceSpecification,
    SpacecraftCommunicationsSpecification,
    StationCommunicationsSpecification,
    TransmitterSpecification,
)
from tatc.catalog.system_schemas import parse_system
from tatc.schemas import Receiver, Transmitter

EPOCH = datetime(2026, 1, 1, tzinfo=timezone.utc)
SOURCE = {"kind": "nominal", "source": "a test (https://example.org)"}


def _system_data(system_id: str) -> dict:
    with open(
        loader._SYSTEMS / f"{system_id}.yml",  # pylint: disable=protected-access
        "r",
        encoding="utf-8",
    ) as f:
        return yaml.safe_load(f)


class TestSpecifications(unittest.TestCase):
    def test_transmitter(self):
        """
        Tests that a transmitter's values need provenance, but not its name
        and link.
        """
        transmitter = TransmitterSpecification(
            name="Ka-band playback",
            link="downlink",
            band="Ka",
            max_data_rate=300e6,
            provenance={"band": SOURCE, "max_data_rate": SOURCE},
        )
        core = transmitter.to_transmitter()
        self.assertIsInstance(core, Transmitter)
        self.assertNotIsInstance(core, TransmitterSpecification)
        self.assertEqual(core.max_data_rate, 300e6)
        self.assertEqual(core.name, "Ka-band playback")
        with self.assertRaises(ValidationError):
            TransmitterSpecification(
                band="Ka", max_data_rate=300e6, provenance={"band": SOURCE}
            )
        with self.assertRaises(ValidationError):
            TransmitterSpecification(
                band="Ka",
                max_data_rate=300e6,
                provenance={"band": SOURCE, "max_data_rate": SOURCE, "eirp": SOURCE},
            )

    def test_receiver(self):
        """
        Tests a receiver's provenance.
        """
        receiver = ReceiverSpecification(
            band="X", g_over_t=31.0, provenance={"band": SOURCE, "g_over_t": SOURCE}
        )
        self.assertIsInstance(receiver.to_receiver(), Receiver)
        with self.assertRaises(ValidationError):
            ReceiverSpecification(band="X", g_over_t=31.0, provenance={})

    def test_station_and_service(self):
        """
        Tests the communications of stations and relay services.
        """
        receiver = {
            "band": "X",
            "g_over_t": 31.0,
            "provenance": {"band": SOURCE, "g_over_t": SOURCE},
        }
        station = StationCommunicationsSpecification(
            antennas=2, receivers=[receiver], provenance={"antennas": SOURCE}
        )
        self.assertEqual(station.to_communications().antennas, 2)
        with self.assertRaises(ValidationError):
            StationCommunicationsSpecification(antennas=2, provenance={})
        with self.assertRaises(ValidationError):
            StationCommunicationsSpecification(
                antennas=1,
                receivers=[receiver, receiver],
                provenance={"antennas": SOURCE},
            )
        service = RelayServiceSpecification(
            antennas=2,
            receivers=[receiver],
            provenance={"antennas": SOURCE},
        )
        communications = service.to_communications()
        self.assertEqual(communications.antennas, 2)
        self.assertIsNone(communications.feeder_rate)
        with self.assertRaises(ValidationError):
            RelayServiceSpecification(
                antennas=2, feeder_rate=600e6, provenance={"antennas": SOURCE}
            )
        spacecraft = SpacecraftCommunicationsSpecification(
            transmitters=[
                {
                    "band": "X",
                    "max_data_rate": 1e6,
                    "provenance": {"band": SOURCE, "max_data_rate": SOURCE},
                }
            ]
        )
        self.assertEqual(
            spacecraft.to_communications().transmitters[0].max_data_rate, 1e6
        )


class TestSystemCommunications(unittest.TestCase):
    def test_spacecraft(self):
        """
        Tests that satellites are built with their spacecraft's communications,
        which are not part of digests.
        """
        data = _system_data("landsat")
        system = parse_system(data)
        digests = {key: system.get_digest(key) for key in system.configurations}
        changed = copy.deepcopy(data)
        for spacecraft in changed["spacecraft"]:
            spacecraft["communications"] = {
                "transmitters": [
                    {
                        "band": "X",
                        "max_data_rate": 384e6,
                        "provenance": {"band": SOURCE, "max_data_rate": SOURCE},
                    }
                ]
            }
        changed_system = parse_system(changed)
        self.assertEqual(
            {key: changed_system.get_digest(key) for key in system.configurations},
            digests,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            satellites = changed_system.build(epoch=EPOCH, payload=False)
        for satellite in satellites:
            self.assertEqual(
                satellite.communications.transmitters[0].max_data_rate, 384e6
            )

    def test_relay_service(self):
        """
        Tests that a relay is built with a service's communications, and that
        only relaying systems have services.
        """
        data = _system_data("tdrs")
        data["services"] = {
            "test": {
                "antennas": 2,
                "receivers": [
                    {
                        "band": "Ka",
                        "g_over_t": 26.5,
                        "provenance": {"band": SOURCE, "g_over_t": SOURCE},
                    }
                ],
                "provenance": {"antennas": SOURCE},
            }
        }
        system = parse_system(data)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            relay = system.build_relay(service="test")
            self.assertEqual(relay.communications.antennas, 2)
            self.assertIsNone(system.build_relay().communications)
            with self.assertRaises(ValueError):
                system.build_relay(service="unknown")
        landsat = _system_data("landsat")
        landsat["services"] = data["services"]
        with self.assertRaises(ValidationError):
            parse_system(landsat)

    def test_ground_terminals(self):
        """
        Tests that ground stations are built with their terminals, which must
        be of the network's stations.
        """
        data = _system_data("noaa_ground_stations")
        data["terminals"] = {
            "FCDAS": {"antennas": 3, "provenance": {"antennas": SOURCE}}
        }
        network = parse_system(data)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stations = {s.name: s for s in network.build("jpss")}
            self.assertEqual(stations["FCDAS"].communications.antennas, 3)
            self.assertIsNone(stations["SvalSat"].communications)
            (plain,) = network.build(members="FCDAS", communications=False)
        self.assertIsNone(plain.communications)
        data["terminals"]["Kiruna"] = data["terminals"]["FCDAS"]
        with self.assertRaises(ValidationError):
            parse_system(data)

    def test_cataloged(self):
        """
        Tests that every cataloged value of communications states its source,
        and that the cataloged relay services and terminals build.
        """
        sources = []
        for system_id in get_system_ids():
            system = get_system(system_id)
            for spacecraft in getattr(system, "spacecraft", ()):
                if spacecraft.communications is not None:
                    for transmitter in spacecraft.communications.transmitters:
                        sources += transmitter.provenance.values()
                        transmitter.to_transmitter()
            for service in getattr(system, "services", {}).values():
                sources += service.provenance.values()
                for receiver in service.receivers:
                    sources += receiver.provenance.values()
            for terminal in getattr(system, "terminals", {}).values():
                sources += terminal.provenance.values()
                for receiver in terminal.receivers:
                    sources += receiver.provenance.values()
        self.assertGreater(len(sources), 100)
        for provenance in sources:
            self.assertIn("(https://", provenance.source)
        tdrs = get_system("tdrs")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            relay = tdrs.build_relay(service="single_access")
            self.assertEqual(relay.communications.antennas, 2)
            self.assertEqual(
                {r.band for r in relay.communications.receivers}, {"S", "Ku", "Ka"}
            )
            network = get_system("near_earth_network")
            stations = network.build("ksat")
        self.assertEqual(
            [s.name for s in stations], ["SG1", "SG2", "SG3", "TR2", "TR3", "SI1"]
        )
        self.assertTrue(all(s.communications.antennas == 1 for s in stations))
