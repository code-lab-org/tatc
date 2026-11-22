"""
Unit tests for the InSituInstrument schema and its footprints on surface
platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
from pydantic import TypeAdapter, ValidationError
from pyproj import Geod

from tatc.analysis import (
    collect_ground_track,
    collect_surface_track,
    find_coincident_ground_track,
)
from tatc.schemas import (
    Aircraft,
    InSituInstrument,
    Instrument,
    Satellite,
    SunSynchronousOrbit,
    SurfacePlatform,
)
from tatc.schemas.air import Trajectory
from tatc.schemas.instrument import AllInstruments

GEOD = Geod(ellps="WGS84")


class TestInSituInstrument(unittest.TestCase):
    """
    Unit tests for the InSituInstrument schema.
    """

    def setUp(self):
        self.instrument = InSituInstrument(name="Thermometer", radius=10000)
        self.station = SurfacePlatform(
            name="Station",
            latitude=35.0,
            longitude=-100.0,
            elevation=800,
            height=1.5,
            instruments=[self.instrument],
        )

    def test_good_data(self):
        """
        Test that the instrument initializes and parses from data by its type or
        its radius.
        """
        self.assertEqual(self.instrument.type, "in_situ")
        adapter = TypeAdapter(AllInstruments)
        self.assertIsInstance(
            adapter.validate_python({"name": "T", "radius": 100}), InSituInstrument
        )
        self.assertIsInstance(
            adapter.validate_python(self.instrument.model_dump()), InSituInstrument
        )
        with self.assertRaises(ValidationError):
            InSituInstrument(name="T", radius=0)
        with self.assertRaises(ValidationError):
            InSituInstrument(name="T")

    def test_only_on_surface_platforms(self):
        """
        Test that satellites and aircraft reject in-situ instruments.
        """
        with self.assertRaisesRegex(ValidationError, "surface platform"):
            Satellite(
                name="S",
                orbit=SunSynchronousOrbit(
                    altitude=700e3, equator_crossing_time="13:30"
                ),
                instruments=[self.instrument],
            )
        start = datetime(2026, 1, 1, tzinfo=timezone.utc)
        trajectory = Trajectory.from_waypoints(
            [(-100, 35), (-99, 35)], start=start, speed=200, altitude=8000
        )
        with self.assertRaisesRegex(ValidationError, "surface platform"):
            Aircraft(name="A", trajectory=trajectory, instruments=[self.instrument])

    def test_footprint_is_disk(self):
        """
        Test that the footprint is a geodesic disk of the radius about the
        station, at the target elevation.
        """
        track = collect_surface_track(self.station)
        footprint = track.geometry.iloc[0]
        self.assertTrue(track.valid_obs.iloc[0])
        area, _ = GEOD.geometry_area_perimeter(footprint)
        self.assertAlmostEqual(abs(area), np.pi * 10000**2, delta=0.002 * np.pi * 1e8)
        longitudes, latitudes = footprint.exterior.coords.xy
        _, _, distances = GEOD.inv(
            np.full(len(longitudes), -100.0),
            np.full(len(latitudes), 35.0),
            longitudes,
            latitudes,
        )
        np.testing.assert_allclose(distances, 10000, rtol=1e-6)
        self.assertTrue(
            all(z == self.station.elevation for _, _, z in footprint.exterior.coords)
        )
        self.assertEqual(self.instrument.get_swath_width(1.5), 20000)

    def test_footprint_across_antimeridian(self):
        """
        Test that a disk across the anti-meridian is split rather than spanning
        the globe.
        """
        station = self.station.model_copy(update={"longitude": 179.95})
        footprint = collect_surface_track(station).geometry.iloc[0]
        self.assertEqual(footprint.geom_type, "MultiPolygon")
        area, _ = GEOD.geometry_area_perimeter(footprint)
        self.assertAlmostEqual(abs(area), np.pi * 10000**2, delta=0.002 * np.pi * 1e8)

    def test_coincident_with_satellite(self):
        """
        Test that a satellite's footprints pair with a station's disk when they
        overlap in time.
        """
        start = datetime(2026, 1, 1, tzinfo=timezone.utc)
        times = [start + timedelta(minutes=m) for m in range(0, 24 * 60)]
        satellite = Satellite(
            name="S",
            orbit=SunSynchronousOrbit(
                altitude=700e3, equator_crossing_time="13:30", epoch=start
            ),
            instruments=[Instrument(name="Imager", field_of_regard=100)],
        )
        stations = collect_surface_track(self.station, times)
        footprints = collect_ground_track(satellite, times)
        pairs = find_coincident_ground_track(
            pd.concat([footprints, stations], ignore_index=True), timedelta(minutes=1)
        )
        self.assertGreater(len(pairs), 0)


if __name__ == "__main__":
    unittest.main()
