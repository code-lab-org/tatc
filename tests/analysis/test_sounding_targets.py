"""
Unit tests for the soundings of limb sounders and occultation receivers as
observations of targets.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone
from functools import partial

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point as ShapelyPoint
from shapely.geometry import box

from tatc.analysis import (
    collect_limb_soundings,
    collect_occultations,
    collect_sounding_observations,
)
from tatc.schemas import (
    CircularOrbit,
    Instrument,
    LimbSounder,
    OccultationReceiver,
    Satellite,
    ScanStart,
    SunSynchronousOrbit,
    SurfacePlatform,
)
from tatc.tradespace import Observations
from tatc.utils.geometry import hash_geometry

START = datetime(2026, 1, 1, tzinfo=timezone.utc)


class TestLimbSoundingObservations(unittest.TestCase):
    """
    Unit tests for the soundings of limb sounders as observations.
    """

    def setUp(self):
        self.end = START + timedelta(minutes=10)
        self.satellite = Satellite(
            name="S",
            orbit=SunSynchronousOrbit(
                altitude=700e3, equator_crossing_time="13:30", epoch=START
            ),
            instruments=[
                Instrument(),
                LimbSounder(
                    name="Limb",
                    scan_azimuth=90,
                    scan_elevations=[10e3, 50e3, 90e3],
                    scan_duration=timedelta(seconds=15),
                    scans=[
                        ScanStart(direction="upward"),
                        ScanStart(offset=timedelta(seconds=40), direction="downward"),
                    ],
                    cycle_period=timedelta(seconds=80),
                ),
            ],
        )
        self.soundings = collect_limb_soundings(
            self.satellite, START, self.end, instrument_index=1
        )

    def test_points_near_soundings(self):
        """
        Test that a point is observed by the soundings within the radius, with
        the point's hash, as `collect_observations` gives it.
        """
        sounding = self.soundings.iloc[5]
        near = ShapelyPoint(sounding.position.x, sounding.position.y)
        far = ShapelyPoint(sounding.position.x + 90, 0)
        observations = collect_sounding_observations(
            [near, far], self.satellite, START, self.end, radius=50e3
        )
        self.assertEqual(len(observations.index), 1)
        observation = observations.iloc[0]
        self.assertEqual(
            observation.target_hash, hash_geometry(ShapelyPoint(*near.coords[0], 0))
        )
        self.assertTrue(observation.geometry.equals(ShapelyPoint(*near.coords[0], 0)))
        self.assertEqual(observation.satellite, "S")
        self.assertEqual(observation.instrument, "Limb")
        self.assertEqual(observation.epoch, sounding.time)
        self.assertEqual(observation.start, sounding.start)
        self.assertTrue(np.isnan(observation.platform_alt))

    def test_regions_contain_soundings(self):
        """
        Test that a region is observed by the soundings within it, at their
        representative points.
        """
        positions = gpd.GeoSeries(self.soundings.position)
        longitude, latitude = positions.x, positions.y
        region = box(
            longitude[:4].min() - 0.1,
            latitude[:4].min() - 0.1,
            longitude[:4].max() + 0.1,
            latitude[:4].max() + 0.1,
        )
        inside = positions.within(region)
        observations = collect_sounding_observations(
            gpd.GeoDataFrame(geometry=[region]),
            self.satellite,
            START,
            self.end,
            radius=1,
        )
        self.assertEqual(len(observations.index), int(inside.sum()))
        self.assertTrue((observations.target_hash == hash_geometry(region)).all())
        self.assertTrue(
            all(
                observed.equals(position)
                for observed, position in zip(
                    observations.geometry, self.soundings.position[inside]
                )
            )
        )

    def test_without_sounders(self):
        """
        Test that platforms without sounders have no observations, that an
        index of another instrument raises, and that the radius is positive.
        """
        point = ShapelyPoint(0, 0)
        for platform in (
            Satellite(name="T", orbit=self.satellite.orbit),
            SurfacePlatform(name="P", latitude=0, longitude=0, elevation=0, height=10),
        ):
            observations = collect_sounding_observations(
                point, platform, START, self.end, radius=1e6
            )
            self.assertTrue(observations.empty)
            self.assertIn("target_hash", observations.columns)
        with self.assertRaisesRegex(ValueError, "not a sounder"):
            collect_sounding_observations(
                point, self.satellite, START, self.end, radius=1e6, instrument_index=0
            )
        with self.assertRaisesRegex(ValueError, "positive"):
            collect_sounding_observations(
                point, self.satellite, START, self.end, radius=0
            )

    def test_analysis_of_observations(self):
        """
        Test that the observations of a trade study take the soundings as their
        analysis.
        """
        targets = gpd.GeoDataFrame(
            geometry=[
                ShapelyPoint(p.x, p.y) for p in self.soundings.position.iloc[::3]
            ],
            crs="EPSG:4326",
        )
        analysis = partial(collect_sounding_observations, radius=100e3)
        observations = Observations(targets, START, self.end - START, analysis=analysis)
        pd.testing.assert_frame_equal(
            observations([self.satellite]).reset_index(drop=True),
            analysis(targets, self.satellite, START, self.end).reset_index(drop=True),
            check_like=True,
        )


class TestOccultationObservations(unittest.TestCase):
    """
    Unit tests for the soundings of occultation receivers as observations.
    """

    def test_points_near_occultations(self):
        """
        Test that a point is observed by the occultations within the radius,
        which need transmitters.
        """
        end = START + timedelta(hours=2)
        receiver = Satellite(
            name="RX",
            orbit=CircularOrbit(altitude=520e3, inclination=24, epoch=START),
            instruments=[OccultationReceiver(name="RO")],
        )
        transmitters = [
            Satellite(
                name=f"TX{i}",
                orbit=CircularOrbit(
                    altitude=20200e3,
                    inclination=55,
                    right_ascension_ascending_node=60 * i,
                    epoch=START,
                ),
            )
            for i in range(3)
        ]
        occultations = collect_occultations(receiver, transmitters, START, end)
        self.assertGreater(len(occultations.index), 0)
        position = occultations.position.iloc[0]
        observations = collect_sounding_observations(
            ShapelyPoint(position.x, position.y),
            receiver,
            START,
            end,
            radius=10e3,
            transmitters=transmitters,
        )
        self.assertEqual(observations.satellite.tolist()[0], "RX")
        self.assertIn(occultations.time.iloc[0], observations.epoch.tolist())
        with self.assertRaisesRegex(ValueError, "transmitters"):
            collect_sounding_observations(
                ShapelyPoint(0, 0), receiver, START, end, radius=10e3
            )


if __name__ == "__main__":
    unittest.main()
