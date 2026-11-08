"""
Regression tests that recompute the results of the entries of pushbroom imagers and
synthetic aperture radars from reduced fixtures, written by the scripts
`docs/validation/validate_catalog_<entry>.py`.

These detect changes in TAT-C's results without network access. They do not themselves
establish agreement with the reference data, which the scripts do.

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

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
MAX_FIXTURE_BYTES = 200e3
TIME_STEP = 10


class _PushbroomFixture:
    """
    Tests shared by the entries of swath instruments, for one configuration.
    """

    entry_id: str
    platform: str
    configuration: str = "swath"
    passes_key: str = "passes"
    instrument_index: int | None = 0

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURES / f"{cls.entry_id}.npz", allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry(cls.entry_id)
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.satellite = cls.entry.to_satellite(
            GeneralPerturbationsOrbit.from_tle(cls.data["tle"].tolist()),
            cls.configuration,
            time_step=TIME_STEP,
            platform=cls.platform,
        )
        cls.recorded = cls.entry.record.configurations[
            cls.configuration
        ].fixture_metrics

    def test_fixture_size(self):
        """
        Test that the fixture stays small, as it is part of the repository.
        """
        size = (FIXTURES / f"{self.entry_id}.npz").stat().st_size
        self.assertLess(size, MAX_FIXTURE_BYTES)

    def test_passes(self):
        """
        Test that the recorded passes are predicted, and that this does not move.
        """
        differences = []
        for latitude, longitude, seconds in self.data[self.passes_key]:
            time = self.day + timedelta(seconds=float(seconds))
            found = collect_observations(
                Point(latitude=float(latitude), longitude=float(longitude)),
                self.satellite,
                time - timedelta(minutes=10),
                time + timedelta(minutes=10),
                self.instrument_index,
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
        self.assertEqual(
            100 * float(np.isfinite(differences).mean()), self.recorded["passes_found"]
        )
        self.assertAlmostEqual(
            float(summarize(differences)["p95"]), self.recorded["epoch_p95"], places=4
        )


class TestOliFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the OLI entry.
    """

    entry_id = "oli"
    platform = "landsat-9"


class TestTropomiFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the TROPOMI entry, whose views are analyzed as a union.
    """

    entry_id = "tropomi"
    platform = "sentinel-5p"
    instrument_index = None


class TestMsiFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the MSI entry, on Sentinel-2B.
    """

    entry_id = "msi"
    platform = "sentinel-2b"


class TestCsarIwFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the IW mode of the C-SAR entry, with its roll profile.
    """

    entry_id = "csar"
    platform = "sentinel-1c"
    configuration = "iw"
    passes_key = "passes_iw"


class TestCsarEwFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the EW mode of the C-SAR entry, with its roll profile.
    """

    entry_id = "csar"
    platform = "sentinel-1c"
    configuration = "ew"
    passes_key = "passes_ew"


if __name__ == "__main__":
    unittest.main()
