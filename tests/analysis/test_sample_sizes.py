"""
Unit tests for the sample size functions in tatc.analysis.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point as ShapelyPoint

from tatc.analysis import (
    collect_feature_observations,
    collect_observations,
    compute_sample_sizes,
)
from tatc.schemas import (
    ConicalInstrument,
    Instrument,
    PointedInstrument,
    Satellite,
    SunSynchronousOrbit,
)
from tatc.utils import VelocityFrame, compute_along_track_field_of_view


class TestSampleSizes(unittest.TestCase):
    """
    Unit tests for `compute_sample_sizes`.
    """

    def setUp(self):
        self.start = datetime(2026, 3, 1, tzinfo=timezone.utc)
        self.end = self.start + timedelta(days=1)
        # a cross-track scanner like ATMS: 2.2 deg beams at 1.1 deg spacing
        self.scanner = PointedInstrument(
            name="scanner",
            field_of_regard=105.6,
            cross_track_field_of_view=105.6,
            along_track_field_of_view=2.2,
            cross_track_pixels=96,
            cross_track_oversampling=0.5,
            is_rectangular=True,
            velocity_frame=VelocityFrame.INERTIAL,
        )
        self.simple = Instrument(name="simple", field_of_regard=100)
        self.satellite = Satellite(
            name="sat",
            orbit=SunSynchronousOrbit(
                altitude=830e3, equator_crossing_time="13:30", epoch=self.start
            ),
            instruments=[self.scanner, self.simple],
        )
        self.points = [ShapelyPoint(-100, 35), ShapelyPoint(-60, -10)]

    def test_point_observations(self):
        """
        Test that the pixels of a pointed instrument are larger across track
        than along track off nadir, that the spacing of overlapping pixels is
        half their size across track, and that the scan angle lies within
        the field of regard.
        """
        observations = collect_observations(
            self.points, self.satellite, self.start, self.end, instrument_index=0
        )
        sizes = compute_sample_sizes(observations, self.satellite)
        self.assertGreater(len(sizes), 0)
        self.assertEqual(list(sizes.index), list(observations.index))
        self.assertTrue((sizes.sample_along > 0).all())
        self.assertTrue((sizes.sample_cross >= sizes.sample_along).all())
        np.testing.assert_allclose(
            sizes.spacing_cross / sizes.sample_cross, 0.5, rtol=0.02
        )
        np.testing.assert_allclose(sizes.spacing_along, sizes.sample_along)
        self.assertTrue((sizes.scan_angle.abs() <= 105.6 / 2 + 0.1).all())
        # the slant range is at least the satellite's altitude
        self.assertTrue((sizes.platform_range > 800e3).all())

    def test_instrument_without_pixels(self):
        """
        Test that an instrument without pixels has undefined sample sizes, but
        a slant range and scan angle.
        """
        observations = collect_observations(
            self.points, self.satellite, self.start, self.end, instrument_index=1
        )
        sizes = compute_sample_sizes(observations, self.satellite)
        self.assertGreater(len(sizes), 0)
        self.assertTrue(sizes.sample_along.isna().all())
        self.assertTrue(sizes.spacing_cross.isna().all())
        self.assertTrue(sizes.platform_range.notna().all())
        self.assertTrue(sizes.scan_angle.notna().all())

    def test_feature_observations(self):
        """
        Test that the samples of feature observations are computed at the
        centroids of the viewed parts, at their times.
        """
        satellite = self.satellite.model_copy(
            update={
                "instruments": [
                    self.scanner.model_copy(
                        update={
                            "along_track_field_of_view": compute_along_track_field_of_view(
                                830e3, 30, 98.7
                            )
                        }
                    )
                ]
            }
        )
        times = [self.start + timedelta(seconds=30 * (i + 0.5)) for i in range(2880)]
        features = gpd.GeoDataFrame(
            {
                "start": pd.DatetimeIndex([self.start], dtype="datetime64[ns, UTC]"),
                "end": pd.DatetimeIndex([self.end], dtype="datetime64[ns, UTC]"),
            },
            geometry=[ShapelyPoint(-100, 35).buffer(0.5)],
            crs="EPSG:4326",
        )
        observations = collect_feature_observations(features, satellite, times)
        sizes = compute_sample_sizes(observations, satellite)
        self.assertGreater(len(sizes), 0)
        np.testing.assert_allclose(
            sizes.platform_range, observations.platform_range, rtol=1e-9
        )
        self.assertTrue((sizes.sample_cross > 0).all())

    def test_unknown_platform(self):
        """
        Test that observations by a platform not given raise an error.
        """
        observations = collect_observations(
            self.points, self.satellite, self.start, self.end
        )
        other = self.satellite.model_copy(update={"name": "other"})
        with self.assertRaises(ValueError):
            compute_sample_sizes(observations, other)

    def test_empty_observations(self):
        """
        Test that no observations give no sample sizes, with their columns.
        """
        observations = collect_observations(
            self.points[0],
            self.satellite,
            self.start,
            self.start + timedelta(minutes=1),
        )
        sizes = compute_sample_sizes(observations, self.satellite)
        self.assertEqual(len(sizes), 0)
        self.assertIn("sample_cross", sizes.columns)

    def test_instrument_override(self):
        """
        Test that the samples of another instrument are sized at the same
        observations: one instrument for every observation, or one for the
        observations of a named instrument (the others sized as their own),
        equal to those of the instrument carried by the satellite.
        """
        observations = collect_observations(
            self.points, self.satellite, self.start, self.end, instrument_index=None
        )
        finer = self.scanner.model_copy(
            update={"name": "finer", "along_track_field_of_view": 1.1}
        )
        carried = compute_sample_sizes(
            observations[observations.instrument == "scanner"],
            self.satellite.model_copy(
                update={"instruments": [finer.model_copy(update={"name": "scanner"})]}
            ),
        )
        sized = compute_sample_sizes(observations, self.satellite, {"scanner": finer})
        scanner = (sized.instrument == "scanner").to_numpy()
        np.testing.assert_allclose(
            sized.sample_along[scanner], carried.sample_along, rtol=1e-12
        )
        self.assertTrue(sized.sample_along[~scanner].isna().all())
        everything = compute_sample_sizes(observations, self.satellite, finer)
        self.assertTrue((everything.sample_along > 0).all())
        np.testing.assert_allclose(
            everything.sample_along[scanner], carried.sample_along, rtol=1e-12
        )

    def test_other_platforms(self):
        """
        Test that the observations of an aircraft at a layer are sized from
        its height above the layer, and that those of a radar station have a
        slant range only.
        """
        from tatc.analysis import collect_footprint_observations, collect_footprints
        from tatc.schemas import Aircraft, RadarStation, Trajectory

        aircraft = Aircraft(
            name="P-3",
            trajectory=Trajectory.from_waypoints(
                [(-100.2, 35.0), (-99.8, 35.0)], self.start, 150.0, 8000.0
            ),
            instruments=[
                self.scanner.model_copy(
                    update={"velocity_frame": VelocityFrame.EARTH_FIXED}
                )
            ],
        )
        radar = RadarStation(name="KXXX", latitude=35.0, longitude=-100.3)
        times = list(aircraft.trajectory.time)
        footprints = collect_footprints([aircraft, radar], times, elevation=3000)
        observations = collect_footprint_observations(
            ShapelyPoint(-100, 35, 3000), footprints
        )
        sizes = compute_sample_sizes(observations, [aircraft, radar])
        flown = sizes[sizes.aircraft.notna()]
        self.assertEqual(len(flown), 1)
        self.assertAlmostEqual(flown.platform_range.iloc[0], 5000, delta=200)
        # a 2.2 degree beam from 5 km above the layer is about 190 m wide
        self.assertAlmostEqual(flown.sample_along.iloc[0], 192, delta=10)
        station = sizes[sizes.radar_station.notna()]
        self.assertEqual(len(station), 1)
        self.assertGreater(station.platform_range.iloc[0], 20e3)
        self.assertTrue(station[["scan_angle", "sample_along"]].isna().all(axis=None))
        with self.assertRaisesRegex(ValueError, "aircraft named 'P-3'"):
            compute_sample_sizes(observations, radar)

    def test_feature_elevation(self):
        """
        Test that the observations of features are sized at the elevation of
        the footprints that viewed them (a surface platform's ground).
        """
        from shapely.geometry import box

        from tatc.schemas import SurfacePlatform

        camera = PointedInstrument(
            name="camera",
            field_of_regard=60,
            cross_track_field_of_view=60,
            along_track_field_of_view=60,
            cross_track_pixels=600,
            along_track_pixels=600,
            is_rectangular=True,
        )
        tower = SurfacePlatform(
            name="Tower",
            latitude=35.0,
            longitude=-100.0,
            elevation=1500,
            height=100,
            instruments=[camera],
        )
        features = gpd.GeoDataFrame(
            {
                "start": pd.DatetimeIndex([self.start], dtype="datetime64[ns, UTC]"),
                "end": pd.DatetimeIndex([self.end], dtype="datetime64[ns, UTC]"),
            },
            geometry=[box(-100.0002, 34.9998, -99.9998, 35.0002)],
            crs="EPSG:4326",
        )
        observations = collect_feature_observations(features, tower, [self.start])
        sizes = compute_sample_sizes(observations, tower)
        self.assertAlmostEqual(sizes.platform_range.iloc[0], 100, delta=1)
        # a 0.1 degree pixel from 100 m is about 17 cm wide
        self.assertAlmostEqual(sizes.sample_cross.iloc[0], 0.17, delta=0.02)

    def test_conical_observations(self):
        """
        Test that the samples of a conical scanner with a beam width are
        computed, and longer along track (the oblique direction of a forward
        view) than across track.
        """
        conical = ConicalInstrument(
            name="conical",
            cone_angle=47.5,
            along_track_field_of_view=0.82,
            scan_half_width=75,
            beam_width=1.8,
            scan_sample_spacing=0.8,
        )
        satellite = self.satellite.model_copy(update={"instruments": [conical]})
        observations = collect_observations(
            self.points, satellite, self.start, self.end
        )
        sizes = compute_sample_sizes(observations, satellite)
        self.assertGreater(len(sizes), 0)
        self.assertTrue(sizes[["sample_along", "sample_cross"]].notna().all().all())
        self.assertTrue((sizes.spacing_along > 0).all())
