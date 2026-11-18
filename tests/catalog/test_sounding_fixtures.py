"""
Regression tests that recompute the results of the sounder entries from reduced
fixtures: the limb sounders MLS and SABER (scans and reference tangent points, written
by `validation/catalog/validate_catalog_mls.py` and
`validation/catalog/validate_catalog_saber.py`, see `validation/catalog/sounding.py`),
and the COSMIC-2 radio occultation receiver TGRS (occultations and the tangent points
of the straight lines between recorded positions, written by
`validation/catalog/validate_catalog_tgrs.py`).

These detect changes in TAT-C's results without network access. They do not themselves
establish agreement with the reference data, which the scripts do.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import unittest
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
from pyproj import Geod

from tatc.analysis import collect_limb_observations
from tatc.analysis.limb_sampling import ScanDirection
from tatc.catalog import get_entry
from tatc.schemas import GeneralPerturbationsOrbit, Satellite

from .common import use_fixture_earth_orientation

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
MAX_FIXTURE_BYTES = 200e3
GEOD = Geod(ellps="WGS84")
# the window either side of an occultation in which it is predicted, as the script
# predicts the fixture's occultations
OCCULTATION_WINDOW = timedelta(minutes=10)


class _LimbFixture:
    """
    Tests shared by the limb sounders: the representative points of scans started
    at the reference's times, the period of the cycle, and the extent of latitudes
    of a day of scans.
    """

    entry_id: str
    configuration: str

    @classmethod
    def setUpClass(cls):
        cls.fixture = FIXTURES / f"{cls.entry_id}.npz"
        with np.load(cls.fixture, allow_pickle=False) as data:
            use_fixture_earth_orientation(cls, data)
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry(cls.entry_id)
        key = cls.configuration
        cls.sounding = cls.entry.configurations[key].sounding
        cls.orbit = GeneralPerturbationsOrbit.from_tle(cls.data[f"tle_{key}"].tolist())
        cls.satellite = cls.entry.to_satellite(cls.orbit, configuration=key)
        cls.day = datetime.fromisoformat(str(cls.data[f"day_{key}"]))
        cls.recorded = cls.entry.record.configurations[key].fixture_metrics

    def _get(self, name):
        return self.data[f"{self.configuration}_{name}"]

    def test_fixture_size(self):
        """
        Test that the fixture stays small, as it is part of the repository.
        """
        self.assertLess(self.fixture.stat().st_size, MAX_FIXTURE_BYTES)

    def test_matched_scans(self):
        """
        Test that scans started at the reference's times reach the sample elevation
        at the recorded distances and times from the reference's tangent points.
        """
        starts, upward = self._get("start"), self._get("upward")
        elapsed = np.full(len(starts), np.nan)
        latitude = np.full(len(starts), np.nan)
        longitude = np.full(len(starts), np.nan)
        for direction in ScanDirection:
            index = np.flatnonzero(upward == (direction == ScanDirection.UPWARD))
            if len(index) == 0:
                continue
            times = [self.day + timedelta(seconds=float(starts[i])) for i in index]
            found = collect_limb_observations(
                self.satellite,
                times,
                self.sounding.scan_azimuth,
                list(self.sounding.scan_elevations),
                timedelta(seconds=self.sounding.scan_duration),
                scan_direction=direction,
                sample_elevation=self.sounding.sample_elevation,
            )
            self.assertEqual(len(found), len(index))
            seconds = np.array([(t - self.day).total_seconds() for t in times])
            for row in found.itertuples():
                start = (row.start.to_pydatetime() - self.day).total_seconds()
                i = index[int(np.argmin(np.abs(seconds - start)))]
                elapsed[i] = (row.time.to_pydatetime() - self.day).total_seconds() - (
                    starts[i]
                )
                latitude[i] = row.position.y
                longitude[i] = row.position.x
        distance = (
            GEOD.inv(
                longitude, latitude, self._get("longitude"), self._get("latitude")
            )[2]
            / 1e3
        )
        error = elapsed - self._get("elapsed")
        self.assertAlmostEqual(
            float(np.median(distance)),
            self.recorded["tangent_position_median"],
            places=4,
        )
        self.assertAlmostEqual(
            float(np.percentile(np.abs(error), 95)),
            self.recorded["epoch_p95"],
            places=4,
        )

    def test_scan_period(self):
        """
        Test that the period of the cycle of scans has the recorded error.
        """
        measured = float(self._get("measured_period"))
        self.assertAlmostEqual(
            100 * abs(self.sounding.cycle_period - measured) / measured,
            self.recorded["scan_period_error"],
            places=9,
        )

    def test_latitude_extent(self):
        """
        Test that a day of scans repeating from midnight spans the reference's
        latitudes with the recorded error.
        """
        found = self.sounding.collect_observations(
            self.satellite, self.day, self.day + timedelta(days=1)
        )
        latitude = np.array([p.y for p in found.position])
        least, greatest = self._get("latitude_extent")
        error = max(abs(latitude.min() - least), abs(latitude.max() - greatest))
        self.assertAlmostEqual(error, self.recorded["latitude_extent_error"], places=6)


class TestMLSFixture(_LimbFixture, unittest.TestCase):
    """
    Regression tests of MLS's scans on Aura.
    """

    entry_id = "mls"
    configuration = "scans"


class TestSABERFixture(_LimbFixture, unittest.TestCase):
    """
    Regression tests of SABER's south-viewing scans on TIMED.
    """

    entry_id = "saber"
    configuration = "south_viewing"


class TestSABERNorthFixture(_LimbFixture, unittest.TestCase):
    """
    Regression tests of SABER's north-viewing scans on TIMED.
    """

    entry_id = "saber"
    configuration = "north_viewing"


class TestTGRSFixture(unittest.TestCase):
    """
    Regression tests of COSMIC-2's occultations: each sampled occultation is
    predicted near its time, at the recorded distances and times from the tangent
    points of the straight lines between the recorded positions.
    """

    @classmethod
    def setUpClass(cls):
        cls.fixture = FIXTURES / "tgrs.npz"
        with np.load(cls.fixture, allow_pickle=False) as data:
            use_fixture_earth_orientation(cls, data)
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("tgrs")
        cls.sounding = cls.entry.configurations["occultations"].sounding
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.satellites = {
            name: Satellite(
                name=name,
                orbit=GeneralPerturbationsOrbit.from_tle(lines.split("\n")),
            )
            for name, lines in zip(cls.data["satellites"], cls.data["element_sets"])
        }
        cls.recorded = cls.entry.record.configurations["occultations"].fixture_metrics

    def test_fixture_size(self):
        """
        Test that the fixture stays small, as it is part of the repository.
        """
        self.assertLess(self.fixture.stat().st_size, MAX_FIXTURE_BYTES)

    def test_occultations(self):
        """
        Test that every sampled occultation is predicted, at the recorded distances
        and times from the true tangent points.
        """
        distances, errors = [], []
        for receiver, transmitter, is_rising, seconds, longitude, latitude in zip(
            self.data["receiver"],
            self.data["transmitter"],
            self.data["is_rising"],
            self.data["true_seconds"],
            self.data["true_longitude"],
            self.data["true_latitude"],
        ):
            time = self.day + timedelta(seconds=float(seconds))
            found = self.sounding.collect_observations(
                self.satellites[str(receiver)],
                self.satellites[str(transmitter)],
                time - OCCULTATION_WINDOW,
                time + OCCULTATION_WINDOW,
                time_step=timedelta(seconds=1),
            )
            found = found[found.is_rising == bool(is_rising)]
            self.assertFalse(found.empty)
            k = int(
                np.argmin(np.abs((found.time - pd.Timestamp(time)).dt.total_seconds()))
            )
            errors.append((found.time.iloc[k].to_pydatetime() - time).total_seconds())
            position = found.position.iloc[k]
            distances.append(
                GEOD.inv(position.x, position.y, longitude, latitude)[2] / 1e3
            )
        self.assertAlmostEqual(
            float(np.median(distances)),
            self.recorded["tangent_position_median"],
            places=4,
        )
        self.assertAlmostEqual(
            float(np.percentile(np.abs(errors), 95)),
            self.recorded["epoch_p95"],
            places=4,
        )


if __name__ == "__main__":
    unittest.main()
