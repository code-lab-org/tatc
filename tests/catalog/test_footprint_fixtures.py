"""
Regression tests that recompute the results of the footprints configurations of beam
instruments (GEDI, ATLID, and CPR) from reduced fixtures, written by their validation
scripts in `validation/catalog/` (see `validation_io.validate_footprints`): the sizes
of their footprints against published ones, and the spacing of their samples along
track.

These detect changes in TAT-C's results without network access. They do not themselves
establish agreement with the reference data, which the scripts do.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import pathlib
import unittest
from datetime import datetime, timedelta, timezone

import geopandas as gpd
import numpy as np
import pandas as pd

from tatc.analysis import compute_sample_sizes
from tatc.catalog import get_entry
from tatc.catalog.standard import get_gates
from tatc.schemas import CircularOrbit, GeneralPerturbationsOrbit
from tatc.utils import VelocityFrame, compute_projected_ray_position

from .common import use_fixture_earth_orientation

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
PUBLISHED_EPOCH = datetime(2026, 10, 1, tzinfo=timezone.utc)
KEY = "footprints"


def _observations(satellite, names, epochs, lon, lat) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": list(names),
            "epoch": pd.DatetimeIndex(epochs),
        },
        geometry=gpd.points_from_xy(lon, lat),
        crs="EPSG:4326",
    )


class _FootprintsFixture:
    """
    Tests shared by the footprints configurations of beam instruments.
    """

    entry_id: str
    platform: str
    # published footprints: instrument, altitude (km), inclination (deg), scan angle
    # (deg), and sizes along and across track (km)
    footprints: tuple

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURES / f"{cls.entry_id}.npz", allow_pickle=False) as data:
            use_fixture_earth_orientation(cls, data)
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry(cls.entry_id)
        cls.orbit = GeneralPerturbationsOrbit.from_tle(cls.data["tle"].tolist())
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.satellite = cls.entry.build(
            configuration=KEY, orbit=cls.orbit, platform=cls.platform
        ).to_satellite(cls.orbit, name=cls.entry.name)
        cls.recorded = cls.entry.record.configurations[KEY].fixture_metrics

    def test_published_dimensions(self):
        """
        Test that the published footprints, compared without the gate, do not move.
        """
        errors = []
        for name, altitude, inclination, angle, along, cross in self.footprints:
            orbit = CircularOrbit(
                altitude=altitude * 1e3, inclination=inclination, epoch=PUBLISHED_EPOCH
            )
            satellite = self.entry.build(configuration=KEY, orbit=orbit).to_satellite(
                orbit, name=self.entry.name
            )
            time = PUBLISHED_EPOCH + timedelta(minutes=10)
            target = compute_projected_ray_position(
                orbit.to_gp_orbit().get_orbit_track(time),
                0,
                0,
                angle,
                0,
                velocity_frame=VelocityFrame.INERTIAL,
            )
            sizes = compute_sample_sizes(
                _observations(
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
        self.assertNotIn("sample_size_within_3pct", self.recorded)
        self.assertAlmostEqual(
            100 * max(errors), self.recorded["sample_size_largest_error"], places=6
        )

    def test_spacing(self):
        """
        Test that the spacing of samples along track does not move, and meets the
        gate as recorded.
        """
        rows = self.data["spacing"]
        sizes = compute_sample_sizes(
            _observations(
                self.satellite,
                self.data["spacing_instruments"],
                [
                    self.day + timedelta(seconds=float(s))
                    for s in self.data["spacing_seconds"]
                ],
                rows[:, 1],
                rows[:, 0],
            ),
            self.satellite,
        )
        ratio = float(np.median(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 2]))
        self.assertAlmostEqual(
            ratio, self.recorded["spacing_along_ratio_median"], places=4
        )
        gates = get_gates(self.entry.instrument_class, "sample_size")
        limit = next(g.limit for g in gates if g.metric == "spacing_along_error")
        self.assertLessEqual(100 * abs(ratio - 1), limit)


class TestGediFootprints(_FootprintsFixture, unittest.TestCase):
    """
    Regression tests of the GEDI entry's footprints, on the International Space
    Station.
    """

    entry_id = "gedi"
    platform = "iss"
    footprints = (("GEDI BEAM0101", 410, 51.6, 0, 0.023, 0.023),)


class TestAtlidFootprints(_FootprintsFixture, unittest.TestCase):
    """
    Regression tests of the ATLID entry's footprint, on EarthCARE.
    """

    entry_id = "atlid"
    platform = "earthcare"
    footprints = (("ATLID", 393, 97.05, 0, 0.015, 0.015),)


class TestCprFootprints(_FootprintsFixture, unittest.TestCase):
    """
    Regression tests of the CPR entry's footprint, on EarthCARE.
    """

    entry_id = "cpr"
    platform = "earthcare"
    footprints = (("CPR", 393.14, 97.05, 0, 1.0, 0.65),)


if __name__ == "__main__":
    unittest.main()
