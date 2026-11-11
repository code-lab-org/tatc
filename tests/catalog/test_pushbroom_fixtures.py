"""
Regression tests that recompute the results of the entries of pushbroom imagers,
synthetic aperture radars, ASCAT, the SMAP radiometer, SSMIS, and OCI's swath from reduced
fixtures, written by
the scripts
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
    tle_key: str = "tle"
    day_key: str = "day"
    instrument_index: int | None = 0
    time_step: float = TIME_STEP

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURES / f"{cls.entry_id}.npz", allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry(cls.entry_id)
        cls.day = datetime.fromisoformat(str(cls.data[cls.day_key]))
        # the day is the period of the analysis, for a schedule of the attitude
        cls.satellite = cls.entry.to_satellite(
            GeneralPerturbationsOrbit.from_tle(cls.data[cls.tle_key].tolist()),
            cls.configuration,
            time_step=cls.time_step,
            platform=cls.platform,
            start=cls.day,
            end=cls.day + timedelta(days=1),
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
    tle_key = "tle_iw"
    day_key = "day_iw"


class TestCsarEwFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the EW mode of the C-SAR entry, with its roll profile.
    """

    entry_id = "csar"
    platform = "sentinel-1c"
    configuration = "ew"
    passes_key = "passes_ew"
    tle_key = "tle_ew"
    day_key = "day_ew"


class TestKarinFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the KaRIn entry, whose two swaths are analyzed as a union.
    """

    entry_id = "karin"
    platform = "swot"
    configuration = "swaths"
    instrument_index = None


class TestAscatFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the ASCAT entry, a union of two swaths beside a nadir gap.
    """

    entry_id = "ascat"
    platform = "metop-c"
    configuration = "swaths"
    instrument_index = None


class TestSmapRadiometerFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the SMAP radiometer entry, a full-circle conical scanner,
    tailored to its scan period.
    """

    entry_id = "smap_radiometer"
    platform = "smap"
    configuration = "looks"
    time_step = 4.11


class TestSsmisFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the SSMIS entry, a conical scanner, tailored to its scan
    period.
    """

    entry_id = "ssmis"
    platform = "dmsp-f18"
    time_step = 1.899


class TestOlciFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the OLCI entry, a pushbroom rolled away from the sun.
    """

    entry_id = "olci"
    platform = "sentinel-3a"


class TestOciFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the OCI entry's swath, a cross-track scanner whose tilt
    follows a pitch profile scheduled by day.
    """

    entry_id = "oci"
    platform = "pace"


class TestLsarWideBandFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the L-SAR entry's modes with a 20 or 40 MHz band.
    """

    entry_id = "lsar"
    platform = "nisar"
    configuration = "swath_20_40mhz"
    passes_key = "passes_swath_20_40mhz"


class TestLsarNarrowBandFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the L-SAR entry's modes with only a 5 MHz band.
    """

    entry_id = "lsar"
    platform = "nisar"
    configuration = "swath_5mhz"
    passes_key = "passes_swath_5mhz"


class TestLsar77MhzFixture(_PushbroomFixture, unittest.TestCase):
    """
    Regression tests of the L-SAR entry's 77 MHz mode, tested on several days, whose
    fixture keeps the passes of the first with its element set.
    """

    entry_id = "lsar"
    platform = "nisar"
    configuration = "swath_77mhz"
    passes_key = "passes_swath_77mhz"
    tle_key = "tle_swath_77mhz"
    day_key = "day_swath_77mhz"


if __name__ == "__main__":
    unittest.main()
