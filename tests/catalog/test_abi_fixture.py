"""
Regression tests that recompute the ABI entry's results from a reduced fixture of
GOES-19 pixels, written by `docs/validation/validate_catalog_abi.py`.

These detect changes in TAT-C's results without network access. They do not
themselves establish agreement with the reference data, which the script does.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import unittest
from datetime import datetime, timedelta

import numpy as np

from tatc.analysis import collect_observations
from tatc.catalog import get_entry, get_gates
from tatc.catalog.metrics import compute_product_limit_error
from tatc.schemas import GeneralPerturbationsOrbit, Point

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "abi.npz"
MAX_FIXTURE_BYTES = 200e3


class TestAbiFixture(unittest.TestCase):
    """
    Regression tests of the ABI entry.
    """

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("abi")
        cls.orbit = GeneralPerturbationsOrbit.from_tle(cls.data["tle"].tolist())
        cls.time = datetime.fromisoformat(str(cls.data["time"]))

    def _predict(self, key):
        rows = self.data[key]
        satellite = self.entry.to_satellite(self.orbit, key, platform="goes-19")
        results = [
            collect_observations(
                Point(latitude=float(latitude), longitude=float(longitude)),
                satellite,
                self.time,
                self.time + timedelta(minutes=10),
            )
            for latitude, longitude in rows[:, :2]
        ]
        observed = np.array([not r.empty for r in results])
        elevation = np.array(
            [float(r.sat_alt.iloc[0]) if not r.empty else np.nan for r in results]
        )
        return rows[:, 2], observed, elevation

    def _recorded(self, key):
        return self.entry.record.configurations[key].fixture_metrics

    def test_fixture_size(self):
        """
        Test that the fixture stays small, as it is part of the repository.
        """
        self.assertLess(FIXTURE.stat().st_size, MAX_FIXTURE_BYTES)

    def test_product_limits(self):
        """
        Test that the product limits meet the gate and do not move.
        """
        limit = get_gates(
            self.entry.instrument_class, self.entry.configurations["full_disk"].grade
        )[0].limit
        for key in ("sst_good_quality", "lst_good_quality"):
            zenith, observed, _ = self._predict(key)
            error = compute_product_limit_error(
                zenith, observed, float(self.data[f"{key}_bound"])
            )
            self.assertLessEqual(error, limit, key)
            self.assertAlmostEqual(
                error, self._recorded(key)["product_limit"], places=6, msg=key
            )

    def test_full_disk(self):
        """
        Test that all sampled pixels of the Earth's disk are observed, with the
        same elevation angles.
        """
        zenith, observed, elevation = self._predict("full_disk")
        recorded = self._recorded("full_disk")
        self.assertEqual(100 * observed.mean(), recorded["imaged_pixels_observed"])
        self.assertAlmostEqual(
            float(np.nanmax(np.abs(elevation - (90 - zenith)))),
            recorded["elevation_error_max"],
            places=6,
        )


if __name__ == "__main__":
    unittest.main()
