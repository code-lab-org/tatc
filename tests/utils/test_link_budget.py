import unittest

import numpy as np

from tatc.schemas import LinkMode, Receiver, ShannonModem, Transmitter
from tatc.utils import (
    compute_antenna_gain,
    compute_atmospheric_loss,
    compute_carrier_to_noise_density,
    compute_free_space_path_loss,
    compute_link_rate,
    compute_mode_rate,
    compute_shannon_rate,
    get_band,
    get_eirp,
    get_frequency,
    get_g_over_t,
    get_link_modes,
)
from tatc.utils.link_budget import _load_table


class TestLinkBudget(unittest.TestCase):
    def test_path_loss(self):
        """
        Tests free-space path loss against 92.45 + 20 log10(km) + 20 log10(GHz).
        """
        self.assertAlmostEqual(
            compute_free_space_path_loss(2574e3, 8.16e9),
            92.45 + 20 * np.log10(2574) + 20 * np.log10(8.16),
            delta=0.01,
        )
        losses = compute_free_space_path_loss([1e6, 2e6], 8e9)
        self.assertAlmostEqual(losses[1] - losses[0], 20 * np.log10(2), places=9)

    def test_antenna_gain(self):
        """
        Tests the gain of a parabolic antenna.
        """
        # efficiency 1: (pi D / lambda) ** 2
        wavelength = 299792458.0 / 8e9
        self.assertAlmostEqual(
            compute_antenna_gain(1, 8e9, 1),
            20 * np.log10(np.pi / wavelength),
            places=9,
        )
        self.assertAlmostEqual(
            compute_antenna_gain(1, 8e9, 0.5) - compute_antenna_gain(1, 8e9, 1),
            10 * np.log10(0.5),
            places=9,
        )

    def test_atmospheric_loss(self):
        """
        Tests the scaling of atmospheric loss by elevation.
        """
        self.assertAlmostEqual(compute_atmospheric_loss(90, 0.1), 0.1)
        self.assertAlmostEqual(compute_atmospheric_loss(30, 0.1), 0.2)
        # held at the minimum elevation
        np.testing.assert_allclose(
            compute_atmospheric_loss([0, 5], 0.1),
            0.1 / np.sin(np.radians(5)),
        )

    def test_aqua_direct_broadcast(self):
        """
        Tests a published link budget: the 15 Mbps X-band direct broadcast of
        EOS Aqua (GSFC 422-11-19-11, Appendix A, Table A-2), whose received
        Eb/N0 is 15.8 dB at 5 degrees elevation (2574 km) and 18.9 dB at 90
        degrees (705 km), with a nominal margin of 5.7 and 8.8 dB.
        """
        cases = [
            # EIRP, range, atmospheric and polarization losses, G/T, Eb/N0, margin
            (17.8, 2574e3, 0.5, 0.4, 21.0, 15.8, 5.7),
            (6.9, 705e3, 0.05, 0.6, 23.4, 18.9, 8.8),
        ]
        for eirp, distance, atmosphere, polarization, g_over_t, ebn0, margin in cases:
            cn0 = compute_carrier_to_noise_density(
                eirp,
                compute_free_space_path_loss(distance, 8.16e9),
                g_over_t,
                atmosphere + polarization,
            )
            self.assertAlmostEqual(cn0 - 10 * np.log10(15e6), ebn0, delta=0.1)
            # required 6.8 dB (BER 1e-3, uncoded) with 3.3 dB implementation loss
            mode = LinkMode(data_rate=15e6, required_ebn0=6.8 + 3.3)
            self.assertEqual(compute_mode_rate(cn0, [mode], margin - 0.15), 15e6)
            self.assertEqual(compute_mode_rate(cn0, [mode], margin + 0.15), 0)

    def test_aqua_link_rate(self):
        """
        Tests the rate of the Aqua link at 5 degrees elevation with schemas,
        with the atmospheric loss scaled from zenith, within 0.5 dB of the
        published margin.
        """
        receiver = Receiver(band="X", g_over_t=21.0, zenith_attenuation=0.05)
        for margin, rate in [(5.2, 15e6), (6.2, 0)]:
            transmitter = Transmitter(
                band="X",
                frequency=8.16e9,
                power=10**1.4,
                antenna_gain=6.8 - 3.0,
                losses=0.4,
                margin=margin,
                modes=[LinkMode(data_rate=15e6, required_ebn0=6.8 + 3.3)],
            )
            self.assertAlmostEqual(get_eirp(transmitter), 17.8, places=6)
            self.assertEqual(compute_link_rate(transmitter, receiver, 2574e3, 5), rate)

    def test_mode_rate(self):
        """
        Tests that the fastest mode that closes is used.
        """
        modes = [
            LinkMode(data_rate=10e6, required_ebn0=2),
            LinkMode(data_rate=100e6, required_ebn0=6),
        ]
        # 100 Mbps needs 80 + 6 = 86 dB-Hz; 10 Mbps needs 72
        np.testing.assert_array_equal(
            compute_mode_rate([70, 75, 86, 90], modes), [0, 10e6, 100e6, 100e6]
        )
        self.assertEqual(compute_mode_rate(86, modes, margin=1), 10e6)

    def test_shannon_rate(self):
        """
        Tests a modem following the capacity of its channel less a gap.
        """
        modem = ShannonModem(bandwidth=1e6, gap=0)
        # S/N = 1 (0 dB): 1 bit/s/Hz
        self.assertAlmostEqual(compute_shannon_rate(60, modem), 1e6)
        capped = ShannonModem(bandwidth=1e6, gap=0, max_spectral_efficiency=0.5)
        self.assertAlmostEqual(compute_shannon_rate(60, capped), 0.5e6)
        self.assertLess(compute_shannon_rate(60, modem, margin=3), 1e6)

    def test_link_rate(self):
        """
        Tests link rates with and without what a budget needs.
        """
        transmitter = Transmitter(
            band="X",
            eirp=10,
            modes=get_link_modes("dvb-s2", 50e6),
            max_data_rate=150e6,
        )
        receiver = Receiver(band="X", g_over_t=30)
        rates = compute_link_rate(transmitter, receiver, [800e3, 2500e3], [80, 5])
        # decreasing with range, capped by the transmitter
        self.assertGreater(rates[0], rates[1])
        self.assertLessEqual(rates[0], 150e6)
        self.assertGreater(rates[1], 0)
        # without a receiver: the transmitter's highest rate
        self.assertEqual(compute_link_rate(transmitter, None, 1e9), 150e6)
        # capped by the receiver
        capped = Receiver(band="X", g_over_t=30, max_data_rate=20e6)
        self.assertEqual(compute_link_rate(transmitter, capped, 800e3, 80), 20e6)
        # between satellites, without an atmosphere
        self.assertGreaterEqual(
            compute_link_rate(transmitter, receiver, 2500e3),
            compute_link_rate(transmitter, receiver, 2500e3, 5),
        )

    def test_power_and_antenna(self):
        """
        Tests the EIRP of a power into a parabolic antenna and the G/T of an
        antenna and noise temperature.
        """
        transmitter = Transmitter(
            band="Ka", power=10, antenna_diameter=0.5, max_data_rate=1e6
        )
        frequency = get_frequency(transmitter)
        self.assertEqual(frequency, get_band("Ka").frequency)
        self.assertAlmostEqual(
            get_eirp(transmitter), 10 + compute_antenna_gain(0.5, frequency)
        )
        receiver = Receiver(
            band="X", antenna_diameter=11.3, system_noise_temperature=200
        )
        self.assertAlmostEqual(
            get_g_over_t(receiver, 8.2e9),
            compute_antenna_gain(11.3, 8.2e9) - 10 * np.log10(200),
        )


class TestTables(unittest.TestCase):
    def test_bands(self):
        """
        Tests the nominal values of bands and their sources.
        """
        for band in ["S", "X", "Ku", "Ka"]:
            defaults = get_band(band)
            self.assertGreater(defaults.frequency, 1e9)
            self.assertGreater(defaults.zenith_attenuation, 0)
            for source in defaults.sources.values():
                self.assertIn("(https://", source)
        self.assertLess(
            get_band("X").zenith_attenuation, get_band("Ka").zenith_attenuation
        )
        with self.assertRaises(ValueError):
            get_band("L")

    def test_modes(self):
        """
        Tests the modes of standards.
        """
        modes = get_link_modes("dvb-s2", 100e6)
        self.assertEqual(len(modes), 28)
        rates = [mode.data_rate for mode in modes]
        self.assertEqual(rates, sorted(rates))
        # QPSK 1/2: Es/N0 1.00 dB at 0.988858 bits per symbol
        (half,) = get_link_modes("dvb-s2", 100e6, ["QPSK 1/2"])
        self.assertAlmostEqual(half.data_rate, 98.8858e6)
        self.assertAlmostEqual(half.required_ebn0, 1.0 - 10 * np.log10(0.988858))
        (loss,) = get_link_modes("ccsds", 1e6, ["LDPC 7/8 QPSK"], 1.5)
        self.assertAlmostEqual(loss.required_ebn0, 5.5)
        for standard in _load_table("modcods.yml")["standards"].values():
            self.assertIn("(https://", standard["source"])
        with self.assertRaises(ValueError):
            get_link_modes("dvb-s3", 1e6)
        with self.assertRaises(ValueError):
            get_link_modes("dvb-s2", 1e6, ["QPSK 7/8"])
