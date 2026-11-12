"""
Regression tests that recompute the result of the TEMPO entry from a reduced fixture
written by `docs/validation/validate_catalog_tempo.py`: grid points near and away
from the edges of the area a held-out scan imaged, with whether each was imaged.

These detect changes in TAT-C's results without network access. They do not themselves
establish agreement with the reference data, which the script does.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import unittest
from datetime import datetime, timedelta

import numpy as np
from skyfield.api import wgs84
from skyfield.framelib import itrs

from tatc.analysis import collect_observations
from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit, Point

from .common import recorded_earth_orientation

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "tempo.npz"
MAX_FIXTURE_BYTES = 200e3
CHECKED = 40


_EARTH_ORIENTATION = recorded_earth_orientation()


def setUpModule():  # pylint: disable=invalid-name
    """Computes with the UT1 table the catalog records were made with."""
    _EARTH_ORIENTATION.start()


def tearDownModule():  # pylint: disable=invalid-name
    """Restores the installed skyfield's UT1 table."""
    _EARTH_ORIENTATION.stop()


class TestTempoFixture(unittest.TestCase):
    """
    Regression tests of the TEMPO entry's field of regard.
    """

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("tempo")
        cls.time = datetime.fromisoformat(str(cls.data["time"]))
        cls.satellite = cls.entry.to_satellite(
            GeneralPerturbationsOrbit.from_tle(cls.data["tle"].tolist()),
            platform="intelsat-40e",
        )
        cls.track = cls.satellite.orbit.to_gp_orbit().get_orbit_track([cls.time])
        cls.recorded = cls.entry.record.configurations[
            "field_of_regard"
        ].fixture_metrics

    def predict(self, latitude, longitude) -> np.ndarray:
        """Gets whether points are in the view and above the horizon."""
        target = wgs84.latlon(latitude, longitude)
        in_view = self.satellite.get_instrument(0).is_in_field_of_view(
            self.track, target
        )
        position = np.asarray(self.track.frame_xyz(itrs).m).reshape(3, 1)
        points = np.asarray(target.itrs_xyz.m)
        lat, lon = np.radians(latitude), np.radians(longitude)
        up = np.stack(
            [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)]
        )
        return in_view & (np.sum((position - points) * up, axis=0) > 0)

    def test_fixture_size(self):
        """
        Test that the fixture stays small, as it is part of the repository.
        """
        self.assertLess(FIXTURE.stat().st_size, MAX_FIXTURE_BYTES)

    def test_points_agree(self):
        """
        Test that the share of points whose prediction agrees with the imaged area
        does not move.
        """
        latitude, longitude, imaged = self.data["points"].T
        predicted = self.predict(latitude, longitude)
        self.assertAlmostEqual(
            100 * float(np.mean(predicted == imaged.astype(bool))),
            self.recorded["points_agreeing"],
            places=6,
        )

    def test_collect_observations_agrees(self):
        """
        Test that `collect_observations` observes the points that the view does,
        for points near the edges.
        """
        latitude, longitude, _ = self.data["points"][:CHECKED].T
        predicted = self.predict(latitude, longitude)
        for lat, lon, expected in zip(latitude, longitude, predicted):
            found = collect_observations(
                Point(latitude=float(lat), longitude=float(lon)),
                self.satellite,
                self.time - timedelta(minutes=10),
                self.time + timedelta(minutes=10),
            )
            self.assertEqual(not found.empty, bool(expected), (lat, lon))


if __name__ == "__main__":
    unittest.main()
