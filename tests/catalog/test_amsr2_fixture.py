"""
Regression tests that recompute the AMSR2 entry's results from a reduced fixture of
observed passes, written by `docs/validation/validate_catalog_amsr2.py`.

These detect changes in TAT-C's results without network access. They do not
themselves establish agreement with the reference data, which the script does.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd

from tatc.analysis import collect_observations, compute_sample_sizes
from tatc.catalog import get_entry, get_gates
from tatc.catalog.metrics import summarize
from tatc.schemas import CircularOrbit, GeneralPerturbationsOrbit, Point

from .common import recorded_earth_orientation

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "amsr2.npz"
MAX_FIXTURE_BYTES = 200e3
SCAN_PERIOD = 1.5
PLATFORM = "gcom-w"


_EARTH_ORIENTATION = recorded_earth_orientation()


def setUpModule():  # pylint: disable=invalid-name
    """Computes with the UT1 table the catalog records were made with."""
    _EARTH_ORIENTATION.start()


def tearDownModule():  # pylint: disable=invalid-name
    """Restores the installed skyfield's UT1 table."""
    _EARTH_ORIENTATION.stop()


def _observe(satellite, index, time, longitude, latitude):
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": satellite.instruments[index].name,
            "epoch": pd.DatetimeIndex(time),
        },
        geometry=gpd.points_from_xy(longitude, latitude),
        crs="EPSG:4326",
    )


class TestAmsr2Fixture(unittest.TestCase):
    """
    Regression tests of the AMSR2 entry's swath configuration.
    """

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("amsr2")
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.recorded = cls.entry.record.configurations["swath"].fixture_metrics

    def test_fixture_size(self):
        """
        Test that the fixture stays small, as it is part of the repository.
        """
        self.assertLess(FIXTURE.stat().st_size, MAX_FIXTURE_BYTES)

    def test_passes(self):
        """
        Test that the recorded passes are predicted with epochs that meet the gate of
        the configuration, and that this does not move.
        """
        satellite = self.entry.to_satellite(
            GeneralPerturbationsOrbit.from_tle(self.data["tle"].tolist()),
            time_step=SCAN_PERIOD,
            platform=PLATFORM,
        )
        differences = []
        for latitude, longitude, seconds in self.data["passes"]:
            time = self.day + timedelta(seconds=float(seconds))
            found = collect_observations(
                Point(latitude=float(latitude), longitude=float(longitude)),
                satellite,
                time - timedelta(minutes=10),
                time + timedelta(minutes=10),
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
        found = 100 * float(np.isfinite(differences).mean())
        p95 = float(summarize(differences)["p95"])
        self.assertEqual(found, self.recorded["passes_found"])
        gates = get_gates(
            self.entry.instrument_class, self.entry.configurations["swath"].grade
        )
        self.assertLessEqual(
            p95, next(g.limit for g in gates if g.metric == "epoch_p95")
        )
        self.assertAlmostEqual(p95, self.recorded["epoch_p95"], places=4)


class TestAmsr2SampleSizeFixture(unittest.TestCase):
    """
    Regression tests of the AMSR2 entry's sample size configuration.
    """

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("amsr2")
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.recorded = cls.entry.record.configurations["sample_size"].fixture_metrics

    def _build(self, orbit):
        return self.entry.to_satellite(
            orbit, "sample_size", time_step=SCAN_PERIOD, platform=PLATFORM
        )

    def test_published_dimensions(self):
        """
        Test that the share of published footprint dimensions within 3% does not move,
        and that it does not meet the gate, as the record states.
        """
        epoch = datetime(2026, 10, 1, tzinfo=timezone.utc)
        time = epoch + timedelta(minutes=10)
        errors = []
        for index, altitude, along, cross in self.data["footprints"]:
            orbit = CircularOrbit(
                altitude=altitude * 1e3, inclination=98.2, epoch=epoch
            )
            satellite = self._build(orbit)
            target = satellite.get_instrument(int(index)).compute_footprint_center(
                orbit.to_gp_orbit().get_orbit_track(time)
            )
            sizes = compute_sample_sizes(
                _observe(
                    satellite,
                    int(index),
                    [time],
                    [target.longitude.degrees],
                    [target.latitude.degrees],
                ),
                satellite,
            )
            errors += [
                abs(sizes.sample_along.iloc[0] / 1e3 / along - 1),
                abs(sizes.sample_cross.iloc[0] / 1e3 / cross - 1),
            ]
        within = 100 * float(np.mean(np.array(errors) <= 0.03))
        self.assertEqual(len(errors), 12)
        self.assertEqual(within, self.recorded["sample_size_within_3pct"])
        gates = get_gates(
            self.entry.instrument_class,
            self.entry.configurations["sample_size"].grade,
        )
        self.assertLess(within, gates[0].limit)

    def test_spacing(self):
        """
        Test that the reported spacing of the samples of each group does not move.
        """
        satellite = self._build(
            GeneralPerturbationsOrbit.from_tle(self.data["tle"].tolist())
        )
        names = [instrument.name for instrument in satellite.instruments]
        for group, name in (("S1", "AMSR2 10.65 GHz"), ("S5", "AMSR2 89.0 GHz")):
            index = names.index(name)
            rows = self.data["spacing"][self.data["spacing_instruments"] == index]
            sizes = compute_sample_sizes(
                _observe(
                    satellite,
                    index,
                    [self.day + timedelta(seconds=float(s)) for s in rows[:, 2]],
                    rows[:, 1],
                    rows[:, 0],
                ),
                satellite,
            )
            self.assertAlmostEqual(
                float(np.median(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
                self.recorded[f"spacing_arc_ratio_median_{group}"],
                places=4,
            )
            self.assertAlmostEqual(
                float(np.median(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 4])),
                self.recorded[f"spacing_scan_ratio_median_{group}"],
                places=4,
            )


if __name__ == "__main__":
    unittest.main()
