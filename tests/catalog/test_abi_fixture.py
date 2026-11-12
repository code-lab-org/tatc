"""
Regression tests that recompute the results of the geostationary entries from
reduced fixtures of pixels: GOES-19 ABI pixels, written by
`docs/validation/validate_catalog_abi.py`, and Himawari-9 AHI pixels, written by
`docs/validation/validate_catalog_ahi.py`.

These detect changes in TAT-C's results without network access. They do not
themselves establish agreement with the reference data, which the scripts do.

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

from .common import recorded_earth_orientation

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
MAX_FIXTURE_BYTES = 200e3


_EARTH_ORIENTATION = recorded_earth_orientation()


def setUpModule():  # pylint: disable=invalid-name
    """Computes with the UT1 table the catalog records were made with."""
    _EARTH_ORIENTATION.start()


def tearDownModule():  # pylint: disable=invalid-name
    """Restores the installed skyfield's UT1 table."""
    _EARTH_ORIENTATION.stop()


class _GeostationaryFixture:
    """
    Regression tests shared by geostationary entries.
    """

    entry_id: str
    platform: str
    limit_keys: tuple[str, ...]
    # the keys whose product limit meets the gate
    passing_keys: tuple[str, ...]

    @classmethod
    def setUpClass(cls):
        cls.fixture = FIXTURES / f"{cls.entry_id}.npz"
        with np.load(cls.fixture, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry(cls.entry_id)
        cls.orbit = GeneralPerturbationsOrbit.from_tle(cls.data["tle"].tolist())
        cls.time = datetime.fromisoformat(str(cls.data["time"]))

    def _predict(self, key):
        rows = self.data[key]
        satellite = self.entry.to_satellite(self.orbit, key, platform=self.platform)
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
        self.assertLess(self.fixture.stat().st_size, MAX_FIXTURE_BYTES)

    def test_product_limits(self):
        """
        Test that the product limits do not move, and that those that met the gate
        still do.
        """
        limit = get_gates(
            self.entry.instrument_class, self.entry.configurations["full_disk"].grade
        )[0].limit
        for key in self.limit_keys:
            zenith, observed, _ = self._predict(key)
            error = compute_product_limit_error(
                zenith, observed, float(self.data[f"{key}_bound"])
            )
            if key in self.passing_keys:
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


class TestAbiFixture(_GeostationaryFixture, unittest.TestCase):
    """
    Regression tests of the ABI entry.
    """

    entry_id = "abi"
    platform = "goes-19"
    limit_keys = ("sst_good_quality", "lst_good_quality")
    passing_keys = limit_keys


class TestAhiFixture(_GeostationaryFixture, unittest.TestCase):
    """
    Regression tests of the AHI entry, whose cloud top height's limit, near the
    limb, does not meet the gate.
    """

    entry_id = "ahi"
    platform = "himawari-9"
    limit_keys = ("cloud_top_height",)
    passing_keys = ()


if __name__ == "__main__":
    unittest.main()
