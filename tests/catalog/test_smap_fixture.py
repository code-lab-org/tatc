"""
Regression tests that recompute the results of the SMAP radiometer entry's sample size
configuration from a reduced fixture, written by
`docs/validation/validate_catalog_smap.py`.

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

from tatc.analysis import compute_sample_sizes
from tatc.catalog import get_entry, get_gates
from tatc.schemas import CircularOrbit, GeneralPerturbationsOrbit

FIXTURE = pathlib.Path(__file__).parent / "fixtures" / "smap_radiometer.npz"
SCAN_PERIOD = 4.11
PLATFORM = "smap"
CONFIGURATION = "footprints"


def _observe(satellite, time, longitude, latitude):
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": satellite.instruments[0].name,
            "epoch": pd.DatetimeIndex(time),
        },
        geometry=gpd.points_from_xy(longitude, latitude),
        crs="EPSG:4326",
    )


class TestSmapFootprintsFixture(unittest.TestCase):
    """
    Regression tests of the SMAP radiometer entry's sample size configuration.
    """

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURE, allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry("smap_radiometer")
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.recorded = cls.entry.record.configurations[CONFIGURATION].fixture_metrics

    def _build(self, orbit):
        return self.entry.to_satellite(
            orbit, CONFIGURATION, time_step=SCAN_PERIOD, platform=PLATFORM
        )

    def test_published_dimensions(self):
        """
        Test that the share of published footprint dimensions within 3% does not move,
        and that it does not meet the gate, as the record states.
        """
        epoch = datetime(2026, 10, 1, tzinfo=timezone.utc)
        time = epoch + timedelta(minutes=10)
        errors = []
        for altitude, along, cross in self.data["footprints"]:
            orbit = CircularOrbit(
                altitude=altitude * 1e3, inclination=98.1, epoch=epoch
            )
            satellite = self._build(orbit)
            target = satellite.get_instrument(0).compute_footprint_center(
                orbit.to_gp_orbit().get_orbit_track(time)
            )
            sizes = compute_sample_sizes(
                _observe(
                    satellite,
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
        self.assertEqual(len(errors), 4)
        self.assertEqual(within, self.recorded["sample_size_within_3pct"])
        gates = get_gates(
            self.entry.instrument_class,
            self.entry.configurations[CONFIGURATION].grade,
        )
        self.assertLess(within, gates[0].limit)

    def test_spacing(self):
        """
        Test that the reported spacing of the samples does not move.
        """
        satellite = self._build(
            GeneralPerturbationsOrbit.from_tle(self.data["tle"].tolist())
        )
        rows = self.data["spacing"]
        sizes = compute_sample_sizes(
            _observe(
                satellite,
                [self.day + timedelta(seconds=float(s)) for s in rows[:, 2]],
                rows[:, 1],
                rows[:, 0],
            ),
            satellite,
        )
        self.assertAlmostEqual(
            float(np.median(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
            self.recorded["spacing_arc_ratio_median"],
            places=4,
        )
        self.assertAlmostEqual(
            float(np.median(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 4])),
            self.recorded["spacing_scan_ratio_median"],
            places=4,
        )


if __name__ == "__main__":
    unittest.main()
