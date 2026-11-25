"""
Regression tests that recompute the results of the pixels configurations of
geostationary imagers (ABI, AHI, and FCI) from reduced fixtures, written by their
validation scripts in `validation/catalog/` (see `validation_io.validate_fixed_grid`):
their pixels projected from the nominal position of the fixed grid, the spacing of
their pixels, and the published sizes at the sub-satellite point.

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
from pyproj import Transformer

from tatc.analysis import compute_sample_sizes
from tatc.catalog import get_entry
from tatc.catalog.standard import get_gates
from tatc.schemas import CircularOrbit, GeosynchronousOrbit
from tatc.utils import VelocityFrame, compute_projected_ray_position

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
PUBLISHED_EPOCH = datetime(2026, 10, 1, tzinfo=timezone.utc)
KEY = "pixels"
TO_ECEF = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)


def _ecef(lon, lat) -> np.ndarray:
    return np.stack(TO_ECEF.transform(lon, lat, 0.0), axis=-1)


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


class _FixedGridFixture:
    """
    Tests shared by the pixels configurations of geostationary imagers.
    """

    entry_id: str
    instrument: str
    # the published size at the sub-satellite point: instrument, altitude (km),
    # inclination (deg), scan angle (deg), and sizes east-west and north-south (km)
    footprints: tuple

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURES / f"{cls.entry_id}_pixels.npz") as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry(cls.entry_id)
        cls.time = datetime.fromisoformat(str(cls.data["time"]))
        orbit = GeosynchronousOrbit(
            longitude=float(cls.data["longitude"]), epoch=cls.time, inclination=0
        )
        cls.satellite = cls.entry.build(configuration=KEY, orbit=orbit).to_satellite(
            orbit, name=cls.entry.name
        )
        cls.recorded = cls.entry.record.configurations[KEY].fixture_metrics
        cls.gates = {
            g.metric: g.limit
            for g in get_gates(cls.entry.instrument_class, "sample_size")
        }

    def test_published_dimensions(self):
        """
        Test that the published size at the sub-satellite point, compared without
        the gate, does not move.
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

    def test_pixels(self):
        """
        Test that pixels projected from the nominal position lie as recorded from the
        grid's, within the gate.
        """
        size = int(self.data["size"])
        view = self.satellite.get_instrument(0)
        track = self.satellite.orbit.to_gp_orbit().get_orbit_track(self.time)
        differences = []
        for row, column, lat, lon in self.data["pixels"]:
            # TAT-C numbers pixels from the south and from the east
            position = view.compute_projected_pixel_position(
                track, size - 1 - int(row), size - 1 - int(column)
            )
            predicted = _ecef(position.longitude.degrees, position.latitude.degrees)
            differences.append(np.linalg.norm(predicted - _ecef(lon, lat)) / 1e3)
        median = float(np.median(differences))
        self.assertAlmostEqual(
            median, self.recorded["pixel_position_median"], delta=0.01
        )
        self.assertLessEqual(median, self.gates["pixel_position_median"])

    def test_spacing(self):
        """
        Test that the spacing of pixels, east-west and north-south, does not move.
        """
        rows = self.data["spacing"]
        sizes = compute_sample_sizes(
            _observations(
                self.satellite,
                [self.instrument] * len(rows),
                [self.time] * len(rows),
                rows[:, 1],
                rows[:, 0],
            ),
            self.satellite,
        )
        self.assertAlmostEqual(
            float(np.median(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 2])),
            self.recorded["spacing_along_ratio_median"],
            places=4,
        )
        self.assertAlmostEqual(
            float(np.median(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
            self.recorded["spacing_cross_ratio_median"],
            places=4,
        )


class TestAbiPixels(_FixedGridFixture, unittest.TestCase):
    """
    Regression tests of the ABI entry's fixed grid (sweep about the x axis).
    """

    entry_id = "abi"
    instrument = "ABI"
    footprints = (("ABI", 35786.023, 0, 0, 2.0, 2.0),)


class TestAhiPixels(_FixedGridFixture, unittest.TestCase):
    """
    Regression tests of the AHI entry's fixed grid (sweep about the y axis).
    """

    entry_id = "ahi"
    instrument = "AHI"
    footprints = (("AHI", 35785.863, 0, 0, 2.0, 2.0),)


class TestFciPixels(_FixedGridFixture, unittest.TestCase):
    """
    Regression tests of the FCI entry's fixed grid (sweep about the y axis).
    """

    entry_id = "fci"
    instrument = "FCI"
    footprints = (("FCI", 35786.4, 0, 0, 2.0, 2.0),)


if __name__ == "__main__":
    unittest.main()
