"""
Unit tests for the performance and resources of instruments and spacecraft.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import timedelta

import numpy as np
from pydantic import ValidationError

from tatc.schemas import (
    CircularOrbit,
    Instrument,
    InstrumentPerformance,
    InstrumentResources,
    PointedInstrument,
    Satellite,
    SpacecraftResources,
    WalkerConstellation,
)


class TestInstrumentPerformance(unittest.TestCase):
    """
    Tests the scaling of performance with altitude.
    """

    def test_ground_sample_distance(self):
        """
        Tests that the ground sample distance scales by its exponent.
        """
        performance = InstrumentPerformance(
            reference_altitude=800e3, ground_sample_distance=30
        )
        self.assertAlmostEqual(performance.get_ground_sample_distance(400e3), 15)
        sar = InstrumentPerformance(
            reference_altitude=800e3,
            ground_sample_distance=20,
            ground_sample_distance_exponent=0,
        )
        self.assertAlmostEqual(sar.get_ground_sample_distance(400e3), 20)
        altimeter = InstrumentPerformance(
            reference_altitude=1000e3,
            ground_sample_distance=20e3,
            ground_sample_distance_exponent=0.5,
        )
        self.assertAlmostEqual(altimeter.get_ground_sample_distance(4000e3), 40e3)

    def test_signal_to_noise_ratio(self):
        """
        Tests that the signal-to-noise ratio scales by its exponent, inversely to
        the noise, and that arrays are scaled elementwise.
        """
        lidar = InstrumentPerformance(
            reference_altitude=500e3,
            signal_to_noise_ratio=100,
            signal_to_noise_exponent=2,
        )
        np.testing.assert_allclose(
            lidar.get_signal_to_noise_ratio([250e3, 1000e3]), [400, 25]
        )
        self.assertAlmostEqual(lidar.get_noise_scale(1000e3), 4)
        passive = InstrumentPerformance(
            reference_altitude=500e3, signal_to_noise_ratio=100
        )
        self.assertAlmostEqual(passive.get_signal_to_noise_ratio(1000e3), 100)
        self.assertIsInstance(passive.get_noise_scale(1000e3), float)

    def test_missing_values(self):
        """
        Tests that values not given, or an altitude that is not positive, raise.
        """
        performance = InstrumentPerformance(reference_altitude=500e3)
        with self.assertRaises(ValueError):
            performance.get_ground_sample_distance(500e3)
        with self.assertRaises(ValueError):
            performance.get_signal_to_noise_ratio(500e3)
        with self.assertRaises(ValueError):
            performance.get_noise_scale(0)
        with self.assertRaises(ValidationError):
            InstrumentPerformance(reference_altitude=0)


class TestResources(unittest.TestCase):
    """
    Tests the resources of instruments and spacecraft.
    """

    def test_total(self):
        """
        Tests that totals sum the values given, and are `None` if none is given.
        """
        total = InstrumentResources.total(
            [
                InstrumentResources(mass=10, data_rate=1e6),
                None,
                InstrumentResources(mass=5),
            ]
        )
        self.assertEqual(total.mass, 15)
        self.assertEqual(total.data_rate, 1e6)
        self.assertIsNone(total.power)

    def test_data_volume(self):
        """
        Tests the volume of data produced while observing.
        """
        resources = InstrumentResources(data_rate=1e6)
        self.assertEqual(resources.get_data_volume(timedelta(minutes=1)), 60e6)
        with self.assertRaises(ValueError):
            InstrumentResources(mass=1).get_data_volume(timedelta(minutes=1))

    def test_spacecraft_masses(self):
        """
        Tests that the dry mass cannot exceed the launch mass.
        """
        with self.assertRaises(ValidationError):
            SpacecraftResources(launch_mass=100, dry_mass=200)
        SpacecraftResources(launch_mass=200, dry_mass=100)

    def test_ballistic_coefficient(self):
        """
        Tests that the ballistic coefficient is the mass (by default, dry,
        or else at launch) over the drag coefficient (by default, 2.2) times
        the drag area.
        """
        resources = SpacecraftResources(launch_mass=300, dry_mass=220, drag_area=2)
        self.assertAlmostEqual(resources.get_ballistic_coefficient(), 50)
        self.assertAlmostEqual(resources.get_ballistic_coefficient(440), 100)
        launch = SpacecraftResources(launch_mass=330, drag_area=1.5, drag_coefficient=2)
        self.assertAlmostEqual(launch.get_ballistic_coefficient(), 110)
        with self.assertRaises(ValueError):
            SpacecraftResources(dry_mass=100).get_ballistic_coefficient()
        with self.assertRaises(ValueError):
            SpacecraftResources(drag_area=1).get_ballistic_coefficient()
        with self.assertRaises(ValidationError):
            SpacecraftResources(specific_impulse=0)

    def test_unknown_fields(self):
        """
        Tests that unknown fields are rejected.
        """
        with self.assertRaises(ValidationError):
            InstrumentResources(weight=10)
        with self.assertRaises(ValidationError):
            SpacecraftResources(mass=10)


class TestPlatformResources(unittest.TestCase):
    """
    Tests the resources of instruments carried by space systems.
    """

    def setUp(self):
        self.instrument = PointedInstrument(
            cross_track_field_of_view=10,
            along_track_field_of_view=1,
            performance={"reference_altitude": 500e3, "ground_sample_distance": 30},
            resources={"mass": 10, "power": 20},
        )

    def test_payload_resources(self):
        """
        Tests the totals of a satellite's instruments, counting those without
        resources as nothing.
        """
        satellite = Satellite(
            name="test",
            orbit=CircularOrbit(altitude=500e3, inclination=50),
            instruments=[self.instrument, self.instrument, Instrument()],
            resources={"launch_mass": 150},
        )
        total = satellite.get_payload_resources()
        self.assertEqual((total.mass, total.power, total.data_rate), (20, 40, None))
        self.assertIsNone(
            Satellite(
                name="default", orbit=CircularOrbit(altitude=500e3, inclination=50)
            )
            .get_payload_resources()
            .mass
        )

    def test_constellation_members(self):
        """
        Tests that the members of a constellation have its resources.
        """
        constellation = WalkerConstellation(
            name="test",
            orbit=CircularOrbit(altitude=500e3, inclination=50),
            number_satellites=4,
            number_planes=2,
            instruments=[self.instrument],
            resources={"launch_mass": 150, "design_life": 5},
        )
        for member in constellation.generate_members():
            self.assertEqual(member.resources.launch_mass, 150)
            self.assertEqual(member.instruments[0].resources.mass, 10)

    def test_serialization(self):
        """
        Tests that performance and resources survive serialization.
        """
        satellite = Satellite(
            name="test",
            orbit=CircularOrbit(altitude=500e3, inclination=50),
            instruments=[self.instrument],
            resources={"storage": 1e12},
        )
        loaded = Satellite.model_validate_json(satellite.model_dump_json())
        self.assertEqual(loaded, satellite)
        self.assertEqual(loaded.instruments[0].performance.ground_sample_distance, 30)
