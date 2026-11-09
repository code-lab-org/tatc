"""
Regression tests that recompute the results of the ICESat-2 ATLAS entry from a reduced
fixture, written by `docs/validation/validate_catalog_atlas.py`.

These detect changes in TAT-C's results without network access. They do not themselves
establish agreement with the reference data, which the script does.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import unittest
from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from tatc.analysis import collect_observations
from tatc.catalog import get_entry
from tatc.catalog.metrics import summarize
from tatc.schemas import GeneralPerturbationsOrbit, Point

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "atlas.npz"
MAX_FIXTURE_BYTES = 200e3


class TestAtlasFixture(unittest.TestCase):
    """
    Regression tests of the ATLAS entry: beam segments sampled over Antarctica, each
    tested against the view of its pair of beams.
    """

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("atlas")
        cls.orbit = GeneralPerturbationsOrbit.from_tle(cls.data["tle"].tolist())
        cls.day = datetime.fromisoformat(str(cls.data["day"]))

    def test_fixture_size(self):
        """
        Test that the fixture stays small, as it is part of the repository.
        """
        self.assertLess(FIXTURE.stat().st_size, MAX_FIXTURE_BYTES)

    def test_segments(self):
        """
        Test that the recorded segments are observed by their pairs' views, with the
        recorded epochs.
        """
        satellite = self.entry.to_satellite(
            self.orbit, "beam_pairs", time_step=1, platform="icesat-2"
        )
        differences = []
        for latitude, longitude, seconds, pair in self.data["passes"]:
            time = self.day + timedelta(seconds=float(seconds))
            found = collect_observations(
                Point(latitude=float(latitude), longitude=float(longitude)),
                satellite,
                time - timedelta(minutes=2),
                time + timedelta(minutes=2),
                int(pair),
            )
            differences.append(
                min(
                    ((e - pd.Timestamp(time)).total_seconds() for e in found.epoch),
                    key=abs,
                )
                if not found.empty
                else np.nan
            )
        differences = np.array(differences, dtype=float)
        found = np.isfinite(differences)
        recorded = self.entry.record.configurations["beam_pairs"].fixture_metrics
        self.assertEqual(100 * float(found.mean()), recorded["passes_found"])
        self.assertAlmostEqual(
            float(summarize(differences[found])["p95"]),
            recorded["epoch_p95"],
            places=4,
        )

    def test_pairs_ordered_left_to_right(self):
        """
        Test that the views of the pairs are ordered from left to right, without gaps
        between them.
        """
        satellite = self.entry.to_satellite(
            self.orbit, "beam_pairs", time_step=1, platform="icesat-2"
        )
        rolls = [i.roll_angle for i in satellite.instruments]
        self.assertEqual(rolls, sorted(rolls, reverse=True))
        for left, right in zip(satellite.instruments, satellite.instruments[1:]):
            self.assertLess(
                left.roll_angle - left.cross_track_field_of_view / 2,
                right.roll_angle + right.cross_track_field_of_view / 2,
            )


if __name__ == "__main__":
    unittest.main()
