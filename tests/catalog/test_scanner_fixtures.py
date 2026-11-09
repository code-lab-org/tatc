"""
Regression tests that recompute the results of the entries of cross-track scanners
from reduced fixtures, written by `docs/validation/validate_catalog_atms.py`,
`docs/validation/validate_catalog_viirs.py`,
`docs/validation/validate_catalog_modis.py`,
`docs/validation/validate_catalog_iasi.py`,
`docs/validation/validate_catalog_cris.py`,
`docs/validation/validate_catalog_dpr.py`, and
`docs/validation/validate_catalog_mhs.py`.

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

from tatc.analysis import collect_observations, compute_sample_sizes
from tatc.catalog import get_entry, get_gates
from tatc.catalog.metrics import summarize
from tatc.schemas import CircularOrbit, GeneralPerturbationsOrbit, Point
from tatc.utils import VelocityFrame, compute_projected_ray_position

FIXTURES = pathlib.Path(__file__).parent / "fixtures"
MAX_FIXTURE_BYTES = 200e3
TO_ECEF = Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
PUBLISHED_EPOCH = datetime(2026, 10, 1, tzinfo=timezone.utc)
SIZE_TOLERANCE = 0.03


def _ecef(lon, lat) -> np.ndarray:
    return np.stack(TO_ECEF.transform(lon, lat, 0.0), axis=-1)


def _observations(satellite, names, epochs, lon, lat) -> gpd.GeoDataFrame:
    return gpd.GeoDataFrame(
        {
            "satellite": satellite.name,
            "instrument": names,
            "epoch": pd.DatetimeIndex(epochs),
        },
        geometry=gpd.points_from_xy(lon, lat),
        crs="EPSG:4326",
    )


class _ScannerFixture:
    """
    Tests shared by the entries of cross-track scanners.
    """

    entry_id: str
    sizes_key: str
    footprints: tuple
    platform = "noaa-20"
    meets_published_gate = True
    meets_pixel_gate = True

    @classmethod
    def setUpClass(cls):
        with np.load(FIXTURES / f"{cls.entry_id}.npz", allow_pickle=False) as data:
            cls.data = {key: data[key] for key in data.files}
        cls.entry = get_entry(cls.entry_id)
        cls.orbit = GeneralPerturbationsOrbit.from_tle(cls.data["tle"].tolist())
        cls.day = datetime.fromisoformat(str(cls.data["day"]))
        cls.satellite = cls.entry.to_satellite(
            cls.orbit, cls.sizes_key, platform=cls.platform
        )
        cls.names = np.array([i.name for i in cls.satellite.instruments])
        cls.recorded = cls.entry.record.configurations[cls.sizes_key].fixture_metrics

    def time_step(self) -> float:
        """
        Gets the analysis time step of the fixture's passes.
        """
        raise NotImplementedError

    def spacing_names(self, rows: np.ndarray) -> np.ndarray:
        """
        Gets the names of the instruments of rows of measured spacing.
        """
        raise NotImplementedError

    def recompute_pixel_differences(self) -> np.ndarray:
        """
        Gets the differences (km) between recomputed pixels and the fixture's samples.
        """
        raise NotImplementedError

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
        satellite = self.entry.to_satellite(
            self.orbit, "swath", time_step=self.time_step(), platform=self.platform
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
        recorded = self.entry.record.configurations["swath"].fixture_metrics
        self.assertEqual(
            100 * float(np.isfinite(differences).mean()), recorded["passes_found"]
        )
        self.assertAlmostEqual(
            float(summarize(differences)["p95"]), recorded["epoch_p95"], places=4
        )

    def test_published_dimensions(self):
        """
        Test that published footprint dimensions meet the gate and do not move.
        """
        errors = []
        for name, altitude, inclination, angle, along, cross in self.footprints:
            orbit = CircularOrbit(
                altitude=altitude * 1e3, inclination=inclination, epoch=PUBLISHED_EPOCH
            )
            satellite = self.entry.to_satellite(orbit, self.sizes_key)
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
        within = 100 * float(np.mean(np.array(errors) <= SIZE_TOLERANCE))
        configuration = self.entry.configurations[self.sizes_key]
        gates = get_gates(self.entry.instrument_class, configuration.grade)
        limit = next(g.limit for g in gates if g.metric == "sample_size_within_3pct")
        self.assertEqual(within >= limit, self.meets_published_gate)
        self.assertEqual(within, self.recorded["sample_size_within_3pct"])

    def test_pixels(self):
        """
        Test that pixels are positioned as recorded, and meet the gate if the entry
        does.
        """
        differences = self.recompute_pixel_differences()
        median = float(np.median(differences))
        configuration = self.entry.configurations[self.sizes_key]
        gates = get_gates(self.entry.instrument_class, configuration.grade)
        limit = next(g.limit for g in gates if g.metric == "pixel_position_median")
        self.assertEqual(median <= limit, self.meets_pixel_gate)
        self.assertAlmostEqual(
            median, self.recorded["pixel_position_median"], delta=0.01
        )

    def test_spacing(self):
        """
        Test that the reported spacing of samples does not move.
        """
        rows = self.data["spacing"]
        sizes = compute_sample_sizes(
            _observations(
                self.satellite,
                self.spacing_names(rows),
                [
                    self.day + timedelta(seconds=float(s))
                    for s in self.data["spacing_seconds"]
                ],
                rows[:, 1],
                rows[:, 0],
            ),
            self.satellite,
        )
        self.assertAlmostEqual(
            float(np.nanmedian(sizes.spacing_along.to_numpy() / 1e3 / rows[:, 2])),
            self.recorded["spacing_along_ratio_median"],
            places=4,
        )
        self.assertAlmostEqual(
            float(np.nanmedian(sizes.spacing_cross.to_numpy() / 1e3 / rows[:, 3])),
            self.recorded["spacing_cross_ratio_median"],
            places=4,
        )


class TestAtmsFixture(_ScannerFixture, unittest.TestCase):
    """
    Regression tests of the ATMS entry.
    """

    entry_id = "atms"
    sizes_key = "beams"
    footprints = (
        ("ATMS channels 1-2", 824, 98.7, 0, 74.8, 74.8),
        ("ATMS channels 1-2", 824, 98.7, 52.725, 141.8, 323.1),
        ("ATMS channels 3-16", 824, 98.7, 0, 31.6, 31.6),
        ("ATMS channels 3-16", 824, 98.7, 52.725, 60.0, 136.7),
        ("ATMS channels 17-22", 824, 98.7, 0, 15.8, 15.8),
        ("ATMS channels 17-22", 824, 98.7, 52.725, 30.0, 68.4),
    )

    def time_step(self) -> float:
        return 8 / 3

    def spacing_names(self, rows: np.ndarray) -> np.ndarray:
        return np.full(len(rows), "ATMS channels 3-16")

    def recompute_pixel_differences(self) -> np.ndarray:
        instrument = self.satellite.get_instrument(1)
        orbit = self.satellite.orbit.to_gp_orbit()
        differences = []
        for seconds, pixel, _, lat, lon in self.data["pixels"]:
            position = instrument.compute_projected_pixel_position(
                orbit.get_orbit_track(self.day + timedelta(seconds=float(seconds))),
                int(pixel),
                0,
            )
            predicted = _ecef(position.longitude.degrees, position.latitude.degrees)
            differences.append(np.linalg.norm(predicted - _ecef(lon, lat)) / 1e3)
        return np.array(differences)


class TestViirsFixture(_ScannerFixture, unittest.TestCase):
    """
    Regression tests of the VIIRS entry.
    """

    entry_id = "viirs"
    sizes_key = "zones"
    footprints = (
        ("VIIRS zone 2", 824, 98.7, 0, 0.742, 0.776),
        ("VIIRS zone 4", 824, 98.7, 56.0, 1.60, 1.58),
    )

    def time_step(self) -> float:
        return float(self.data["scan_period"])

    def spacing_names(self, rows: np.ndarray) -> np.ndarray:
        return self.names[rows[:, 4].astype(int)]

    def recompute_pixel_differences(self) -> np.ndarray:
        orbit = self.satellite.orbit.to_gp_orbit()
        differences = []
        for seconds, zone, pixel, _, lat, lon in self.data["pixels"]:
            track = orbit.get_orbit_track(self.day + timedelta(seconds=float(seconds)))
            instrument = self.satellite.get_instrument(int(zone))
            fore, aft = (
                instrument.compute_projected_pixel_position(track, int(pixel), detector)
                for detector in (7, 8)
            )
            predicted = (
                _ecef(fore.longitude.degrees, fore.latitude.degrees)
                + _ecef(aft.longitude.degrees, aft.latitude.degrees)
            ) / 2
            differences.append(np.linalg.norm(predicted - _ecef(lon, lat)) / 1e3)
        return np.array(differences)


class TestModisFixture(_ScannerFixture, unittest.TestCase):
    """
    Regression tests of the MODIS entry.
    """

    entry_id = "modis"
    sizes_key = "pixels"
    platform = "terra"
    footprints = (
        ("MODIS", 705, 98.2, 0, 1.0, 1.0),
        ("MODIS", 705, 98.2, 55.0, 2.0, 4.8),
    )

    def time_step(self) -> float:
        return 1.477

    def spacing_names(self, rows: np.ndarray) -> np.ndarray:
        return np.full(len(rows), "MODIS")

    def recompute_pixel_differences(self) -> np.ndarray:
        instrument = self.satellite.get_instrument(0)
        orbit = self.satellite.orbit.to_gp_orbit()
        differences = []
        for seconds, frame, pixel, height, _, lat, lon in self.data["pixels"]:
            position = instrument.compute_projected_pixel_position(
                orbit.get_orbit_track(self.day + timedelta(seconds=float(seconds))),
                int(frame),
                int(pixel),
                elevation=float(height),
            )
            predicted = _ecef(position.longitude.degrees, position.latitude.degrees)
            differences.append(np.linalg.norm(predicted - _ecef(lon, lat)) / 1e3)
        return np.array(differences)


class TestIasiFixture(_ScannerFixture, unittest.TestCase):
    """
    Regression tests of the IASI entry, whose published footprint at the swath edge
    does not meet the gate, so that the entry is partial.
    """

    entry_id = "iasi"
    sizes_key = "pixels"
    platform = "metop-b"
    meets_published_gate = False
    footprints = (
        ("IASI", 817, 98.7, 0, 12.0, 12.0),
        ("IASI", 817, 98.7, 48.3, 20.0, 39.0),
    )

    def time_step(self) -> float:
        return 8

    def spacing_names(self, rows: np.ndarray) -> np.ndarray:
        return np.full(len(rows), "IASI")

    def recompute_pixel_differences(self) -> np.ndarray:
        instrument = self.satellite.get_instrument(0)
        orbit = self.satellite.orbit.to_gp_orbit()
        differences = []
        for seconds, pixel, _, lat, lon in self.data["pixels"]:
            position = instrument.compute_projected_pixel_position(
                orbit.get_orbit_track(self.day + timedelta(seconds=float(seconds))),
                int(pixel),
                0,
            )
            predicted = _ecef(position.longitude.degrees, position.latitude.degrees)
            differences.append(np.linalg.norm(predicted - _ecef(lon, lat)) / 1e3)
        return np.array(differences)


class _BeamScannerFixture(_ScannerFixture):
    """
    Tests shared by entries whose fixtures give a pixel of one instrument per row.
    """

    size_instrument: str

    def spacing_names(self, rows: np.ndarray) -> np.ndarray:
        return np.full(len(rows), self.size_instrument)

    def recompute_pixel_differences(self) -> np.ndarray:
        instrument = self.satellite.get_instrument(0)
        orbit = self.satellite.orbit.to_gp_orbit()
        differences = []
        for seconds, pixel, _, lat, lon in self.data["pixels"]:
            position = instrument.compute_projected_pixel_position(
                orbit.get_orbit_track(self.day + timedelta(seconds=float(seconds))),
                int(pixel),
                0,
            )
            predicted = _ecef(position.longitude.degrees, position.latitude.degrees)
            differences.append(np.linalg.norm(predicted - _ecef(lon, lat)) / 1e3)
        return np.array(differences)


class TestDprFixture(_BeamScannerFixture, unittest.TestCase):
    """
    Regression tests of the DPR entry.
    """

    entry_id = "dpr"
    sizes_key = "beams"
    size_instrument = "DPR"
    platform = "gpm"
    footprints = (
        ("DPR", 400, 65, 0, 5.04, 5.04),
        ("DPR", 400, 65, 17.0, 5.04, 5.57),
        ("DPR", 435, 65, 0, 5.48, 5.48),
        ("DPR", 435, 65, 17.0, 5.48, 6.05),
    )

    def time_step(self) -> float:
        return 10


class TestMhsFixture(_BeamScannerFixture, unittest.TestCase):
    """
    Regression tests of the MHS entry.
    """

    entry_id = "mhs"
    sizes_key = "beams"
    size_instrument = "MHS"
    platform = "metop-c"
    footprints = (
        ("MHS", 817, 98.7, 0, 15.9, 15.9),
        ("MHS", 817, 98.7, 49.444, 27.1, 52.8),
    )

    def time_step(self) -> float:
        return 8 / 3


class TestCrisFixture(_ScannerFixture, unittest.TestCase):
    """
    Regression tests of the CrIS entry, whose pixels do not meet the gates, so that
    the entry is partial.
    """

    entry_id = "cris"
    sizes_key = "pixels"
    meets_published_gate = False
    meets_pixel_gate = False
    footprints = (
        ("CrIS", 824, 98.7, 0, 14.0, 14.0),
        ("CrIS", 824, 98.7, 49.43, 24.0, 48.0),
    )

    def time_step(self) -> float:
        return 8

    def recompute_pixel_differences(self) -> np.ndarray:
        instrument = self.satellite.get_instrument(0)
        orbit = self.satellite.orbit.to_gp_orbit()
        differences = []
        for seconds, cross, along, _, lat, lon in self.data["pixels"]:
            position = instrument.compute_projected_pixel_position(
                orbit.get_orbit_track(self.day + timedelta(seconds=float(seconds))),
                int(cross),
                int(along),
            )
            predicted = _ecef(position.longitude.degrees, position.latitude.degrees)
            differences.append(np.linalg.norm(predicted - _ecef(lon, lat)) / 1e3)
        return np.array(differences)

    def test_spacing(self):
        """
        The spacing of CrIS's samples is not measured.
        """
        self.skipTest("the spacing of CrIS's samples is not measured")


if __name__ == "__main__":
    unittest.main()
