"""
Regression tests that recompute the GMI entry's results from a reduced fixture of
observed passes, written by `docs/validation/validate_catalog_gmi.py`.

These detect changes in TAT-C's results without network access. They do not
themselves establish agreement with the reference data, which the script does.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import unittest
from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from tatc.analysis import collect_observations
from tatc.catalog import get_entry, get_gates
from tatc.catalog.metrics import summarize
from tatc.schemas import GeneralPerturbationsOrbit, Point

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "gmi.npz"
MAX_FIXTURE_BYTES = 200e3
SCAN_PERIOD = 1.875


class TestGmiFixture(unittest.TestCase):
    """
    Regression tests of the GMI entry.
    """

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("gmi")
        cls.orbit = GeneralPerturbationsOrbit.from_tle(cls.data["tle"].tolist())
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        satellite = cls.entry.to_satellite(
            cls.orbit, time_step=SCAN_PERIOD, platform="gpm"
        )
        differences = []
        for (latitude, longitude, seconds), index in zip(
            cls.data["passes"], cls.data["feeds"]
        ):
            time = cls.day + timedelta(seconds=float(seconds))
            found = collect_observations(
                Point(latitude=float(latitude), longitude=float(longitude)),
                satellite,
                time - timedelta(minutes=10),
                time + timedelta(minutes=10),
                int(index),
            )
            differences.append(
                min(
                    ((e - pd.Timestamp(time)).total_seconds() for e in found.epoch),
                    key=abs,
                )
                if not found.empty
                else np.nan
            )
        cls.differences = np.array(differences, dtype=float)
        cls.recorded = cls.entry.record.configurations["feeds"].fixture_metrics

    def test_fixture_size(self):
        """
        Test that the fixture stays small, as it is part of the repository.
        """
        self.assertLess(FIXTURE.stat().st_size, MAX_FIXTURE_BYTES)

    def test_fixture_covers_both_feeds(self):
        """
        Test that the fixture samples passes of each feed.
        """
        self.assertEqual(set(self.data["feeds"].tolist()), {0, 1})

    def test_passes_found(self):
        """
        Test that the recorded passes are predicted, and that this does not move.
        """
        found = 100 * float(np.isfinite(self.differences).mean())
        self.assertEqual(found, self.recorded["passes_found"])
        self.assertGreaterEqual(found, 99)

    def test_epochs(self):
        """
        Test that epochs meet the gate of the configuration and do not move.
        """
        entry = self.entry
        gates = get_gates(entry.instrument_class, entry.configurations["feeds"].grade)
        limit = next(g.limit for g in gates if g.metric == "epoch_p95")
        p95 = summarize(self.differences)["p95"]
        self.assertLessEqual(p95, limit)
        self.assertAlmostEqual(p95, self.recorded["epoch_p95"], places=4)


if __name__ == "__main__":
    unittest.main()
