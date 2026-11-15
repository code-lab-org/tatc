"""
Regression tests that recompute the TRMM TMI entry's results from a reduced fixture,
written by `validation/catalog/validate_catalog_tmi.py`.

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

from .common import use_fixture_earth_orientation

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "tmi.npz"
MAX_FIXTURE_BYTES = 200e3
SCAN_PERIOD = 1.899
PLATFORM = "trmm"
# published footprints at the center of the scan: instrument, altitude (km), and
# sizes along and across track (km)
FOOTPRINTS = (
    ("TMI 10.65 GHz", 403, 68.1, 41.3),
    ("TMI 19.35 GHz", 403, 34.6, 20.9),
    ("TMI 21.3 GHz", 403, 31.1, 18.9),
    ("TMI 37.0 GHz", 403, 18.3, 11.1),
    ("TMI 85.5 GHz", 403, 7.8, 4.8),
)
SPACING_INSTRUMENTS = {
    "S1": "TMI 10.65 GHz",
    "S2": "TMI 19.35 GHz",
    "S3": "TMI 85.5 GHz",
}


def _observe(satellite, names, times, longitude, latitude):
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": names,
            "epoch": pd.DatetimeIndex(times),
        },
        geometry=gpd.points_from_xy(longitude, latitude),
        crs="EPSG:4326",
    )


class _TmiFixture:
    """
    Loads the fixture and the entry.
    """

    key: str

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            use_fixture_earth_orientation(cls, data)
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("tmi")
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.recorded = cls.entry.record.configurations[cls.key].fixture_metrics

    def build(self, orbit):
        """
        Gets the configuration's satellite on an orbit.
        """
        return self.entry.to_satellite(
            orbit, self.key, time_step=SCAN_PERIOD, platform=PLATFORM
        )


class TestTmiFixture(_TmiFixture, unittest.TestCase):
    """
    Regression tests of the TMI entry's feeds, from passes of the first feed (S1).
    """

    key = "feeds"

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
        satellite = self.build(
            GeneralPerturbationsOrbit.from_tle(self.data["tle"].tolist())
        )
        differences = []
        for latitude, longitude, seconds in self.data["passes"]:
            time = self.day + timedelta(seconds=float(seconds))
            found = collect_observations(
                Point(latitude=float(latitude), longitude=float(longitude)),
                satellite,
                time - timedelta(minutes=10),
                time + timedelta(minutes=10),
                0,
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
        p95 = float(summarize(differences)["p95"])
        self.assertEqual(
            100 * float(np.isfinite(differences).mean()), self.recorded["passes_found"]
        )
        gates = get_gates(
            self.entry.instrument_class, self.entry.configurations[self.key].grade
        )
        self.assertLessEqual(
            p95, next(g.limit for g in gates if g.metric == "epoch_p95")
        )
        self.assertAlmostEqual(p95, self.recorded["epoch_p95"], places=4)


class TestTmiSampleSizeFixture(_TmiFixture, unittest.TestCase):
    """
    Regression tests of the TMI entry's sample size configuration.
    """

    key = "sample_size"

    def test_published_dimensions(self):
        """
        Test that the share of published footprint dimensions within 3% does not move,
        and that it meets the gate, as the record states.
        """
        epoch = datetime(2026, 10, 1, tzinfo=timezone.utc)
        time = epoch + timedelta(minutes=10)
        errors = []
        for name, altitude, along, cross in FOOTPRINTS:
            orbit = CircularOrbit(altitude=altitude * 1e3, inclination=35, epoch=epoch)
            satellite = self.build(orbit)
            names = [instrument.name for instrument in satellite.instruments]
            target = satellite.get_instrument(
                names.index(name)
            ).compute_footprint_center(orbit.to_gp_orbit().get_orbit_track(time))
            sizes = compute_sample_sizes(
                _observe(
                    satellite,
                    [name],
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
        self.assertEqual(within, self.recorded["sample_size_within_3pct"])
        gates = get_gates(
            self.entry.instrument_class, self.entry.configurations[self.key].grade
        )
        self.assertGreaterEqual(within, gates[0].limit)

    def test_spacing(self):
        """
        Test that the reported spacing of the samples of each group, between scans and
        along the scan's arc, does not move.
        """
        satellite = self.build(
            GeneralPerturbationsOrbit.from_tle(self.data["tle"].tolist())
        )
        for group, name in SPACING_INSTRUMENTS.items():
            rows = self.data[f"spacing_{group}"]
            seconds = self.data[f"spacing_seconds_{group}"]
            sizes = compute_sample_sizes(
                _observe(
                    satellite,
                    [name] * len(rows),
                    [self.day + timedelta(seconds=float(s)) for s in seconds],
                    rows[:, 1],
                    rows[:, 0],
                ),
                satellite,
            )
            self.assertAlmostEqual(
                float(np.nanmedian(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 2])),
                self.recorded[f"spacing_scan_ratio_median_{group}"],
                places=4,
            )
            self.assertAlmostEqual(
                float(np.nanmedian(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
                self.recorded[f"spacing_arc_ratio_median_{group}"],
                places=4,
            )


if __name__ == "__main__":
    unittest.main()
