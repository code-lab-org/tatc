"""
Regression tests that recompute the GMI entry's results from a reduced fixture of
observed passes, written by `docs/validation/validate_catalog_gmi.py`.

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

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "gmi.npz"
MAX_FIXTURE_BYTES = 200e3
SCAN_PERIOD = 1.875


_EARTH_ORIENTATION = recorded_earth_orientation()


def setUpModule():  # pylint: disable=invalid-name
    """Computes with the UT1 table the catalog records were made with."""
    _EARTH_ORIENTATION.start()


def tearDownModule():  # pylint: disable=invalid-name
    """Restores the installed skyfield's UT1 table."""
    _EARTH_ORIENTATION.stop()


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


class TestGmiSampleSizeFixture(unittest.TestCase):
    """
    Regression tests of the GMI entry's sample size configuration.
    """

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("gmi")
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.recorded = cls.entry.record.configurations["sample_size"].fixture_metrics

    def _build(self, orbit):
        return self.entry.to_satellite(
            orbit, "sample_size", time_step=SCAN_PERIOD, platform="gpm"
        )

    def test_published_dimensions(self):
        """
        Test that published footprint dimensions meet the gate and do not move.
        """
        epoch = datetime(2026, 10, 1, tzinfo=timezone.utc)
        time = epoch + timedelta(minutes=10)
        errors = []
        for index, altitude, along, cross in self.data["footprints"]:
            orbit = CircularOrbit(altitude=altitude * 1e3, inclination=65, epoch=epoch)
            satellite = self._build(orbit)
            target = satellite.get_instrument(int(index)).compute_footprint_center(
                orbit.to_gp_orbit().get_orbit_track(time)
            )
            sizes = compute_sample_sizes(
                gpd.GeoDataFrame(
                    {
                        "satellite": [satellite.name],
                        "instrument": [satellite.instruments[int(index)].name],
                        "epoch": pd.DatetimeIndex([time]),
                    },
                    geometry=gpd.points_from_xy(
                        [target.longitude.degrees], [target.latitude.degrees]
                    ),
                    crs="EPSG:4326",
                ),
                satellite,
            )
            errors += [
                abs(sizes.sample_along.iloc[0] / 1e3 / along - 1),
                abs(sizes.sample_cross.iloc[0] / 1e3 / cross - 1),
            ]
        within = 100 * float(np.mean(np.array(errors) <= 0.03))
        self.assertEqual(len(errors), 22)
        gates = get_gates(
            self.entry.instrument_class,
            self.entry.configurations["sample_size"].grade,
        )
        self.assertGreaterEqual(within, gates[0].limit)
        self.assertEqual(within, self.recorded["sample_size_within_3pct"])

    def test_spacing(self):
        """
        Test that the reported spacing of the samples of each feed does not move.
        """
        satellite = self._build(
            GeneralPerturbationsOrbit.from_tle(self.data["tle"].tolist())
        )
        for feed, index in (("S1", 4), ("S2", 5)):
            rows = self.data["spacing"][self.data["spacing_instruments"] == index]
            sizes = compute_sample_sizes(
                gpd.GeoDataFrame(
                    {
                        "satellite": satellite.name,
                        "instrument": satellite.instruments[index].name,
                        "epoch": pd.DatetimeIndex(
                            [self.day + timedelta(seconds=float(s)) for s in rows[:, 2]]
                        ),
                    },
                    geometry=gpd.points_from_xy(rows[:, 1], rows[:, 0]),
                    crs="EPSG:4326",
                ),
                satellite,
            )
            self.assertAlmostEqual(
                float(np.median(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
                self.recorded[f"spacing_arc_ratio_median_{feed}"],
                places=4,
            )
            self.assertAlmostEqual(
                float(np.median(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 4])),
                self.recorded[f"spacing_scan_ratio_median_{feed}"],
                places=4,
            )


if __name__ == "__main__":
    unittest.main()
