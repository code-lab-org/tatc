import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd

from tatc.analysis import collect_downlinks, collect_relay_contacts
from tatc.analysis.link_sampling import (
    compute_link_rates,
    get_generation_rate,
    sample_link_geometry,
)
from tatc.analysis.point_sampling import _get_elevation_angles
from tatc.schemas import (
    DataGeneration,
    GeosynchronousOrbit,
    GroundStation,
    Instrument,
    InstrumentResources,
    Receiver,
    Relay,
    RelayCommunications,
    Satellite,
    SpacecraftCommunications,
    SpacecraftResources,
    StationCommunications,
    SunSynchronousOrbit,
    Transmitter,
)
from tatc.utils.computation import _run

START = datetime(2025, 1, 1, tzinfo=timezone.utc)
END = START + timedelta(hours=6)


class TestLinkSampling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.satellite = Satellite(
            name="s",
            orbit=SunSynchronousOrbit(
                altitude=600e3, equator_crossing_time="10:30", epoch=START
            ),
            communications=SpacecraftCommunications(
                transmitters=[
                    Transmitter(band="X", max_data_rate=150e6),
                    Transmitter(band="Ka", link="crosslink", max_data_rate=300e6),
                ]
            ),
        )
        cls.station = GroundStation(
            name="g", latitude=78.23, longitude=15.4, min_elevation_angle=5
        )
        cls.relay = Relay(
            name="r",
            satellites=[
                Satellite(
                    name="r1", orbit=GeosynchronousOrbit(longitude=-41, epoch=START)
                )
            ],
        )
        cls.downlinks = collect_downlinks([cls.station], [cls.satellite], START, END)
        cls.contacts = pd.concat(
            [
                cls.downlinks.assign(relay=None),
                collect_relay_contacts([cls.relay], [cls.satellite], START, END),
            ],
            ignore_index=True,
        )
        cls.samples = sample_link_geometry(
            cls.contacts,
            [cls.satellite],
            [cls.station],
            [cls.relay],
            start=START,
            step=timedelta(seconds=30),
        )

    def test_cells(self):
        """
        Tests that the cells of each contact cover its duration.
        """
        seconds = self.samples.groupby("contact")["fraction"].sum() * 30
        durations = (self.contacts["end"] - self.contacts["start"]).dt.total_seconds()
        np.testing.assert_allclose(seconds.to_numpy(), durations.to_numpy(), atol=1e-6)
        self.assertTrue((self.samples["fraction"] > 0).all())
        self.assertTrue((self.samples["fraction"] <= 1).all())

    def test_ground_geometry(self):
        """
        Tests the elevations and distances of contacts with a ground station.
        """
        ground = self.samples[self.samples["relay"].isna()]
        self.assertFalse(ground.empty)
        times = [t.to_pydatetime() for t in ground["time"].dt.floor("us")]
        expected = _run(_get_elevation_angles(self.station, self.satellite, times))
        np.testing.assert_allclose(ground["elevation"], expected, atol=1e-4)
        self.assertTrue((ground["elevation"] >= 5 - 1e-6).all())
        self.assertTrue((ground["distance"] >= 600e3).all())
        self.assertTrue((ground["distance"] <= 2400e3).all())

    def test_relay_geometry(self):
        """
        Tests the distances of contacts with a geostationary relay satellite.
        """
        relayed = self.samples[self.samples["relay"].notna()]
        self.assertFalse(relayed.empty)
        self.assertTrue(relayed["elevation"].isna().all())
        self.assertTrue((relayed["distance"] > 35000e3).all())
        self.assertTrue((relayed["distance"] < 48000e3).all())

    def test_rates(self):
        """
        Tests that downlinks reach ground stations and crosslinks relays.
        """
        rates = compute_link_rates(
            self.samples, [self.satellite], [self.station], [self.relay]
        )
        transmitters = rates.groupby(rates["relay"].notna())["transmitter"].unique()
        self.assertEqual(list(transmitters[False]), ["X-band downlink"])
        self.assertEqual(list(transmitters[True]), ["Ka-band crosslink"])
        self.assertEqual(len(rates), len(self.samples))

    def test_band_compatibility(self):
        """
        Tests that a link needs a receiver in its band, if receivers are given,
        and is capped by its receiver.
        """
        station = self.station.model_copy(
            update={
                "communications": StationCommunications(
                    receivers=[Receiver(band="S", g_over_t=20)]
                )
            }
        )
        relay = self.relay.model_copy(
            update={
                "communications": RelayCommunications(
                    receivers=[Receiver(band="Ka", max_data_rate=100e6)]
                )
            }
        )
        rates = compute_link_rates(self.samples, [self.satellite], [station], [relay])
        self.assertTrue(rates["relay"].notna().all())
        self.assertTrue((rates["rate"] == 100e6).all())

    def test_downlink_rate(self):
        """
        Tests a satellite without communications that downlinks at its
        resources' rate to ground stations in any band, but not to relays.
        """
        satellite = self.satellite.model_copy(
            update={
                "communications": None,
                "resources": SpacecraftResources(downlink_rate=300e6),
            }
        )
        station = self.station.model_copy(
            update={
                "communications": StationCommunications(
                    receivers=[Receiver(band="Ka", g_over_t=30, max_data_rate=200e6)]
                )
            }
        )
        rates = compute_link_rates(self.samples, [satellite], [station], [self.relay])
        self.assertTrue(rates["relay"].isna().all())
        self.assertTrue((rates["rate"] == 200e6).all())
        self.assertEqual(set(rates["transmitter"]), {"downlink"})

    def test_generation_rate(self):
        """
        Tests the rate of generation: given, or the instruments' total.
        """
        self.assertIsNone(get_generation_rate(self.satellite))
        satellite = self.satellite.model_copy(
            update={
                "instruments": [
                    Instrument(name="a", resources=InstrumentResources(data_rate=2e6)),
                    Instrument(name="b", resources=InstrumentResources(data_rate=3e6)),
                ]
            }
        )
        self.assertEqual(get_generation_rate(satellite), 5e6)
        satellite = satellite.model_copy(
            update={
                "communications": SpacecraftCommunications(
                    generation=DataGeneration(rate=1e6)
                )
            }
        )
        self.assertEqual(get_generation_rate(satellite), 1e6)

    def test_empty(self):
        """
        Tests no contacts.
        """
        samples = sample_link_geometry(
            self.contacts.iloc[:0], [self.satellite], start=START
        )
        self.assertTrue(samples.empty)
        self.assertIn("distance", samples.columns)
        self.assertTrue(compute_link_rates(samples, [self.satellite]).empty)

    def test_unknown_member(self):
        """
        Tests that contacts must name the members given.
        """
        with self.assertRaises(ValueError):
            sample_link_geometry(self.contacts, [self.satellite], start=START)
