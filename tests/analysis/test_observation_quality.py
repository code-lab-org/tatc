"""
Unit tests for the quality of observations in tatc.analysis.
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
from shapely.geometry import Point

from tatc.analysis import (
    collect_footprint_observations,
    collect_footprints,
    collect_limb_observations,
    collect_observations,
    compute_observation_quality,
)
from tatc.schemas import (
    Aircraft,
    Instrument,
    InstrumentPerformance,
    PointedInstrument,
    RadarStation,
    Satellite,
    SunSynchronousOrbit,
    Trajectory,
)
from tatc.utils import VelocityFrame

START = datetime(2026, 3, 1, tzinfo=timezone.utc)
END = START + timedelta(days=1)
POINTS = [Point(-100, 35), Point(-60, -10)]
PERFORMANCE = InstrumentPerformance(
    reference_altitude=830e3,
    ground_sample_distance=16e3,
    signal_to_noise_ratio=100,
    signal_to_noise_exponent=2,
)


class TestObservationQuality(unittest.TestCase):
    """
    Unit tests for `compute_observation_quality`.
    """

    def setUp(self):
        self.scanner = PointedInstrument(
            name="scanner",
            field_of_regard=105.6,
            cross_track_field_of_view=105.6,
            along_track_field_of_view=1.1,
            cross_track_pixels=96,
            is_rectangular=True,
            velocity_frame=VelocityFrame.INERTIAL,
            performance=PERFORMANCE,
        )
        self.simple = Instrument(
            name="simple",
            field_of_regard=100,
            performance=InstrumentPerformance(
                reference_altitude=830e3, ground_sample_distance=1e3
            ),
        )
        self.satellite = Satellite(
            name="sat",
            orbit=SunSynchronousOrbit(
                altitude=830e3, equator_crossing_time="13:30", epoch=START
            ),
            instruments=[self.scanner, self.simple],
        )

    def test_scaled_performance(self):
        """Performance scales by its power laws to the slant range."""
        observations = collect_observations(
            POINTS, self.satellite, START, END, instrument_index=0
        )
        quality = compute_observation_quality(observations, self.satellite)
        self.assertGreater(len(quality), 0)
        ratio = quality.platform_range / 830e3
        np.testing.assert_allclose(quality.noise_scale, ratio**2)
        np.testing.assert_allclose(quality.signal_to_noise_ratio, 100 / ratio**2)
        np.testing.assert_allclose(quality.ground_sample_distance, 16e3 * ratio)
        np.testing.assert_allclose(
            quality.resolution_along,
            np.maximum(quality.sample_along, quality.spacing_along),
        )
        np.testing.assert_allclose(
            quality.resolution_cross,
            np.maximum(quality.sample_cross, quality.spacing_cross),
        )
        self.assertTrue(quality.sample_azimuth.notna().all())

    def test_without_samples(self):
        """An instrument without pixels resolves its nominal ground sample distance."""
        observations = collect_observations(
            POINTS, self.satellite, START, END, instrument_index=1
        )
        quality = compute_observation_quality(observations, self.satellite)
        self.assertGreater(len(quality), 0)
        self.assertTrue(quality.sample_along.isna().all())
        np.testing.assert_allclose(quality.noise_scale, 1)
        self.assertTrue(quality.signal_to_noise_ratio.isna().all())
        expected = 1e3 * quality.platform_range / 830e3
        np.testing.assert_allclose(quality.resolution_along, expected)
        np.testing.assert_allclose(quality.resolution_cross, expected)

    def test_without_performance(self):
        """An instrument without performance has no scaled performance."""
        satellite = self.satellite.model_copy(
            update={
                "instruments": [self.scanner.model_copy(update={"performance": None})]
            }
        )
        observations = collect_observations(POINTS, satellite, START, END)
        quality = compute_observation_quality(observations, satellite)
        self.assertTrue(
            quality[["ground_sample_distance", "noise_scale", "signal_to_noise_ratio"]]
            .isna()
            .all(axis=None)
        )
        self.assertTrue(quality.resolution_along.notna().all())

    def test_instrument_override(self):
        """A given instrument's samples and performance replace the observing one's."""
        observations = collect_observations(
            POINTS, self.satellite, START, END, instrument_index=1
        )
        quality = compute_observation_quality(
            observations, self.satellite, instrument={"simple": self.scanner}
        )
        ratio = quality.platform_range / 830e3
        np.testing.assert_allclose(quality.noise_scale, ratio**2)
        self.assertTrue(quality.sample_along.notna().all())

    def test_other_platforms(self):
        """A radar station's observations have a slant range only."""
        aircraft = Aircraft(
            name="P-3",
            trajectory=Trajectory.from_waypoints(
                [(-100.2, 35.0), (-99.8, 35.0)], START, 150.0, 8000.0
            ),
            instruments=[
                self.scanner.model_copy(
                    update={
                        "velocity_frame": VelocityFrame.EARTH_FIXED,
                        "performance": InstrumentPerformance(
                            reference_altitude=8000, signal_to_noise_ratio=50
                        ),
                    }
                )
            ],
        )
        radar = RadarStation(name="KXXX", latitude=35.0, longitude=-100.3)
        times = list(aircraft.trajectory.time)
        footprints = collect_footprints([aircraft, radar], times, elevation=3000)
        observations = collect_footprint_observations(Point(-100, 35, 3000), footprints)
        quality = compute_observation_quality(observations, [aircraft, radar])
        flown = quality[quality.aircraft.notna()]
        self.assertEqual(len(flown), 1)
        self.assertAlmostEqual(flown.signal_to_noise_ratio.iloc[0], 50)
        station = quality[quality.radar_station.notna()]
        self.assertGreater(station.platform_range.iloc[0], 20e3)
        self.assertTrue(
            station[["noise_scale", "resolution_along"]].isna().all(axis=None)
        )

    def test_limb_observations(self):
        """Observations without instruments (limb soundings) have a slant range only."""
        observations = collect_limb_observations(
            self.satellite,
            [START + timedelta(minutes=m) for m in range(3)],
            scan_azimuth=0,
            scan_elevations=[10e3, 30e3, 50e3],
            scan_duration=timedelta(seconds=20),
        )
        quality = compute_observation_quality(observations, self.satellite)
        self.assertEqual(len(quality), len(observations))
        self.assertTrue((quality.platform_range > 1000e3).all())
        self.assertTrue(quality.noise_scale.isna().all())

    def test_empty(self):
        """Empty observations get empty quality columns."""
        observations = collect_observations(
            [Point(0, 89.9)], self.satellite, START, START + timedelta(minutes=1)
        )
        quality = compute_observation_quality(observations, self.satellite)
        self.assertEqual(len(quality), 0)
        self.assertIn("noise_scale", quality.columns)
