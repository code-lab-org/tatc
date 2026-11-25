import unittest
from datetime import timedelta

from pydantic import ValidationError

from tatc.schemas import (
    CircularOrbit,
    DataFlow,
    DataGeneration,
    GroundStation,
    LinkMode,
    MOGConstellation,
    Receiver,
    Relay,
    RelayCommunications,
    Satellite,
    ShannonModem,
    SOCConstellation,
    SpacecraftCommunications,
    StationCommunications,
    TrainConstellation,
    Transmitter,
    WalkerConstellation,
)


class TestTransmitter(unittest.TestCase):
    def test_fixed_rate(self):
        """
        Tests a transmitter given only its rate, which has no link budget.
        """
        transmitter = Transmitter(band="X", max_data_rate=150e6)
        self.assertEqual(transmitter.get_name(), "X-band downlink")
        self.assertEqual(transmitter.get_max_data_rate(), 150e6)
        self.assertFalse(transmitter.has_budget())

    def test_budget(self):
        """
        Tests transmitters with a power and antenna or EIRP and modes.
        """
        mode = LinkMode(data_rate=100e6, required_ebn0=4.5)
        self.assertTrue(Transmitter(band="X", eirp=20, modes=[mode]).has_budget())
        self.assertTrue(
            Transmitter(
                band="Ka", power=5, antenna_diameter=0.5, modes=[mode]
            ).has_budget()
        )
        self.assertFalse(Transmitter(band="X", eirp=20, max_data_rate=1e6).has_budget())

    def test_max_data_rate(self):
        """
        Tests that the highest rate is the least of the limits given.
        """
        modes = [
            LinkMode(data_rate=50e6, required_ebn0=2),
            LinkMode(data_rate=200e6, required_ebn0=6),
        ]
        self.assertEqual(
            Transmitter(band="X", eirp=20, modes=modes).get_max_data_rate(), 200e6
        )
        self.assertEqual(
            Transmitter(
                band="X", eirp=20, modes=modes, max_data_rate=100e6
            ).get_max_data_rate(),
            100e6,
        )
        modem = ShannonModem(bandwidth=100e6, max_spectral_efficiency=2)
        self.assertEqual(
            Transmitter(band="X", eirp=20, modem=modem).get_max_data_rate(), 200e6
        )
        self.assertEqual(
            Transmitter(
                band="X", eirp=20, modem=ShannonModem(bandwidth=1e6)
            ).get_max_data_rate(),
            float("inf"),
        )

    def test_invalid(self):
        """
        Tests that a transmitter rejects conflicting or missing parameters.
        """
        with self.assertRaises(ValidationError):
            Transmitter(band="X")
        with self.assertRaises(ValidationError):
            Transmitter(band="X", eirp=20, power=5, antenna_gain=10, max_data_rate=1e6)
        with self.assertRaises(ValidationError):
            Transmitter(band="X", power=5, max_data_rate=1e6)
        with self.assertRaises(ValidationError):
            Transmitter(
                band="X",
                power=5,
                antenna_gain=10,
                antenna_diameter=0.5,
                max_data_rate=1e6,
            )
        with self.assertRaises(ValidationError):
            Transmitter(band="L", max_data_rate=1e6)


class TestReceiver(unittest.TestCase):
    def test_merit(self):
        """
        Tests receivers given a G/T or an antenna and noise temperature.
        """
        self.assertTrue(Receiver(band="X", g_over_t=31).has_budget())
        self.assertTrue(
            Receiver(
                band="X", antenna_diameter=11, system_noise_temperature=150
            ).has_budget()
        )
        self.assertFalse(Receiver(band="X", max_data_rate=1e6).has_budget())

    def test_invalid(self):
        """
        Tests that a receiver rejects conflicting or partial parameters.
        """
        with self.assertRaises(ValidationError):
            Receiver(band="X", g_over_t=31, antenna_diameter=11)
        with self.assertRaises(ValidationError):
            Receiver(band="X", antenna_diameter=11)

    def test_one_per_band(self):
        """
        Tests that stations and relays take at most one receiver per band.
        """
        receivers = [Receiver(band="X", g_over_t=31), Receiver(band="X", g_over_t=30)]
        with self.assertRaises(ValidationError):
            StationCommunications(receivers=receivers)
        with self.assertRaises(ValidationError):
            RelayCommunications(receivers=receivers)


class TestCommunications(unittest.TestCase):
    def setUp(self):
        self.communications = SpacecraftCommunications(
            transmitters=[Transmitter(band="X", max_data_rate=150e6)],
            generation=DataGeneration(rate=10e6, activity="sunlit"),
        )

    def test_defaults(self):
        """
        Tests the defaults: no communications, and generation always at the
        instruments' rate.
        """
        orbit = CircularOrbit(altitude=500e3, inclination=50)
        self.assertIsNone(Satellite(name="s", orbit=orbit).communications)
        self.assertIsNone(
            GroundStation(name="g", latitude=0, longitude=0).communications
        )
        generation = SpacecraftCommunications().generation
        self.assertIsNone(generation.rate)
        self.assertEqual(generation.activity, "always")
        self.assertEqual(DataFlow().step, timedelta(seconds=30))
        with self.assertRaises(ValidationError):
            DataFlow(step=timedelta(0))

    def test_constellation_members(self):
        """
        Tests that the members of every kind of constellation have its
        communications.
        """
        orbit = CircularOrbit(altitude=500e3, inclination=50)
        constellations = [
            WalkerConstellation(
                name="w",
                orbit=orbit,
                number_satellites=4,
                number_planes=2,
                communications=self.communications,
            ),
            TrainConstellation(
                name="t",
                orbit=orbit,
                number_satellites=2,
                interval=timedelta(minutes=5),
                communications=self.communications,
            ),
            MOGConstellation(
                name="m",
                orbit=orbit,
                parallel_axis=1000,
                transverse_axis=500,
                number_satellites=3,
                communications=self.communications,
            ),
            SOCConstellation(
                name="s",
                orbit=orbit,
                polar=False,
                swath_width=3000e3,
                packing_distance=1,
                communications=self.communications,
            ),
            SOCConstellation(
                name="sp",
                orbit=CircularOrbit(altitude=780e3, inclination=86.4),
                swath_width=3000e3,
                packing_distance=1,
                communications=self.communications,
            ),
        ]
        for constellation in constellations:
            members = constellation.generate_members()
            self.assertTrue(members)
            for member in members:
                self.assertEqual(member.communications, self.communications)
            self.assertIsNot(members[0].communications, constellation.communications)

    def test_relay(self):
        """
        Tests a relay's communications.
        """
        relay = Relay(
            name="r",
            satellites=[Satellite(name="r1", orbit=CircularOrbit(altitude=35786e3))],
            communications=RelayCommunications(
                antennas=2,
                receivers=[Receiver(band="Ka", g_over_t=26.5)],
                feeder_rate=600e6,
            ),
        )
        self.assertEqual(relay.communications.antennas, 2)

    def test_serialization(self):
        """
        Tests that communications survive serialization.
        """
        satellite = Satellite(
            name="s",
            orbit=CircularOrbit(altitude=500e3, inclination=50),
            communications=self.communications,
        )
        loaded = Satellite.model_validate_json(satellite.model_dump_json())
        self.assertEqual(loaded, satellite)
