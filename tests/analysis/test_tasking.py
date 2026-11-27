"""
Unit tests for the tatc.analysis.tasking module.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
from shapely.geometry import Point

from tatc.analysis import collect_imaging_opportunities
from tatc.analysis.tasking import _get_direction, _Pointing
from tatc.schemas import (
    CircularOrbit,
    PointedInstrument,
    Satellite,
    WalkerConstellation,
)

START = datetime(2026, 6, 1, tzinfo=timezone.utc)
END = START + timedelta(days=1)
IMAGER = PointedInstrument(
    name="imager",
    field_of_regard=60,
    cross_track_field_of_view=1,
    along_track_field_of_view=1,
)
SATELLITE = Satellite(
    name="agile",
    orbit=CircularOrbit(altitude=500e3, inclination=97.4, epoch=START),
    instruments=[IMAGER],
)


def _requests(count=30, seed=1, **columns):
    """Requests at random points over the western United States."""
    rng = np.random.default_rng(seed)
    points = [
        Point(x, y)
        for x, y in zip(rng.uniform(-120, -80, count), rng.uniform(30, 48, count))
    ]
    return gpd.GeoDataFrame(columns, geometry=points, crs="EPSG:4326")


def _off_nadir(roll, pitch):
    return np.degrees(np.arccos(np.cos(np.radians(roll)) * np.cos(np.radians(pitch))))


class TestOpportunities(unittest.TestCase):
    """
    Tests the opportunities to image requests.
    """

    def test_within_field_of_regard(self):
        """
        Tests that opportunities are the periods in which the target is within
        the field of regard: at its edge at their ends, and within it between.
        """
        opportunities = collect_imaging_opportunities(
            SATELLITE, _requests(), START, END
        )
        self.assertGreater(len(opportunities), 0)
        self.assertTrue(opportunities["start"].is_monotonic_increasing)
        for row in opportunities.itertuples():
            pointing = _Pointing(
                SATELLITE,
                IMAGER,
                row.geometry,
                row.start,
                row.end,
                timedelta(seconds=1),
            )
            off_nadir = _off_nadir(pointing.roll, pointing.pitch)
            self.assertLessEqual(off_nadir.max(), 30 + 0.01)
            self.assertAlmostEqual(off_nadir[0], 30, delta=0.01)
            self.assertAlmostEqual(off_nadir[-1], 30, delta=0.01)

    def test_clipped_to_requests(self):
        """
        Tests that opportunities lie between the onsets and expiry of their
        requests.
        """
        rng = np.random.default_rng(2)
        onset = [START + timedelta(hours=h) for h in rng.uniform(0, 12, 30)]
        expires = [time + timedelta(hours=6) for time in onset]
        requests = _requests(onset=onset, expires=expires)
        everything = collect_imaging_opportunities(SATELLITE, _requests(), START, END)
        clipped = collect_imaging_opportunities(SATELLITE, requests, START, END)
        self.assertLess(len(clipped), len(everything))
        for row in clipped.itertuples():
            self.assertGreaterEqual(row.start, requests.onset[row.request])
            self.assertLessEqual(row.end, requests.expires[row.request])

    def test_constellations(self):
        """
        Tests that constellations are their members, and that requests
        without opportunities give an empty frame.
        """
        walker = WalkerConstellation(
            name="walker",
            orbit=SATELLITE.orbit,
            instruments=[IMAGER],
            number_satellites=2,
            number_planes=2,
        )
        opportunities = collect_imaging_opportunities(walker, _requests(), START, END)
        self.assertEqual(
            set(opportunities["satellite"]),
            {member.name for member in walker.generate_members()},
        )
        empty = collect_imaging_opportunities(
            SATELLITE, _requests(), START, START + timedelta(minutes=1)
        )
        self.assertEqual(len(empty), 0)
        self.assertEqual(list(empty.columns), list(opportunities.columns))


class TestPointing(unittest.TestCase):
    """
    Tests the pointing of views.
    """

    def test_directions(self):
        """
        Tests that the directions of views are unit vectors whose angle from
        nadir is the off-nadir angle of their roll and pitch.
        """
        roll, pitch = np.array([0, 20, -15, 10]), np.array([0, 0, 25, -30])
        directions = _get_direction(roll, pitch)
        np.testing.assert_allclose(np.linalg.norm(directions, axis=1), 1)
        np.testing.assert_allclose(
            np.degrees(np.arccos(directions[:, 2])), _off_nadir(roll, pitch)
        )
