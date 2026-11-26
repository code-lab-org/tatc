"""
Unit tests for exporting observations to the files of other tools.
"""

import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pytest
from shapely.geometry import Point, box

from tatc.export import to_cf_netcdf, to_ioda

xr = pytest.importorskip("xarray")
netCDF4 = pytest.importorskip("netCDF4")

START = datetime(2026, 1, 1, 6, tzinfo=timezone.utc)


def _observations():
    """Observations of points (one at a height) and a region, with values."""
    return gpd.GeoDataFrame(
        {
            "epoch": [START + timedelta(seconds=s) for s in (0, 90.4, 3600)],
            "satellite": ["a", "a", "b"],
            "platform_alt": [60.0, 80.0, np.nan],
            "platform_az": [10.0, 200.0, 30.0],
            "solar_alt": [30.0, 40.0, 50.0],
            "solar_az": [100.0, 110.0, 120.0],
            "temperature": [280.0, 281.5, np.nan],
            "temperature_observed": [280.4, 281.1, np.nan],
            "temperature_error_std": [0.5, 0.6, 0.7],
            "flag": [0, 1, 0],
            "valid_obs": [True, True, False],
        },
        geometry=[Point(10, 20, 1500), Point(-30, -5), box(0, 0, 2, 2)],
        crs="EPSG:4326",
    )


class TestCFNetCDF(unittest.TestCase):
    """
    Unit tests for `to_cf_netcdf`.
    """

    def test_round_trip(self):
        """A CF point file holds the observations' places, times, and values."""
        observations = _observations()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "observations.nc"
            dataset = to_cf_netcdf(
                observations,
                path,
                attributes={"temperature": {"units": "K"}},
                global_attributes={"title": "test"},
            )
            with xr.open_dataset(path) as read:
                self.assertEqual(read.attrs["featureType"], "point")
                self.assertEqual(read.attrs["title"], "test")
                np.testing.assert_allclose(read.lat[:2], [20, -5])
                self.assertTrue(0 <= float(read.lat[2]) <= 2)
                np.testing.assert_allclose(read.height[0], 1500)
                self.assertTrue(np.isnan(read.height[1]))
                self.assertEqual(read.temperature.attrs["units"], "K")
                self.assertTrue(np.isnan(read.temperature[2]))
                np.testing.assert_array_equal(read.time.values, dataset.time.values)
                self.assertEqual(
                    read.height.attrs["standard_name"],
                    "height_above_reference_ellipsoid",
                )
        # numeric columns by default, not booleans or names
        self.assertIn("flag", dataset)
        self.assertNotIn("valid_obs", dataset)
        self.assertNotIn("satellite", dataset)

    def test_selected_variables(self):
        """Only selected columns are exported, which must exist."""
        dataset = to_cf_netcdf(_observations(), variables=["temperature"])
        self.assertEqual(list(dataset.data_vars), ["temperature"])
        self.assertNotIn("height", to_cf_netcdf(_observations().iloc[1:]).coords)
        with self.assertRaisesRegex(ValueError, "no columns"):
            to_cf_netcdf(_observations(), variables=["missing"])


class TestIODA(unittest.TestCase):
    """
    Unit tests for `to_ioda`.
    """

    def test_layout(self):
        """An IODA file holds the groups, names, types, and fill values of JEDI."""
        observations = _observations()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "observations.nc4"
            to_ioda(
                observations,
                path,
                variables={"temperature_observed": "airTemperature"},
                qc={"temperature_observed": "flag"},
                units={"temperature_observed": "K"},
            )
            with netCDF4.Dataset(path) as dataset:
                self.assertEqual(len(dataset.dimensions["Location"]), 3)
                self.assertEqual(
                    set(dataset.groups), {"MetaData", "ObsValue", "ObsError", "PreQC"}
                )
                meta = dataset["MetaData"]
                self.assertEqual(meta["dateTime"].dtype, np.int64)
                self.assertEqual(
                    meta["dateTime"].units, "seconds since 1970-01-01T00:00:00Z"
                )
                epoch = int(START.timestamp())
                np.testing.assert_array_equal(
                    meta["dateTime"][:], [epoch, epoch + 90, epoch + 3600]
                )
                np.testing.assert_allclose(meta["sensorZenithAngle"][:2], [30, 10])
                self.assertTrue(np.ma.is_masked(meta["sensorZenithAngle"][2]))
                np.testing.assert_allclose(meta["solarZenithAngle"][:], [60, 50, 40])
                np.testing.assert_allclose(meta["height"][0], 1500)
                value = dataset["ObsValue"]["airTemperature"]
                self.assertEqual(value.dtype, np.float32)
                self.assertEqual(value.units, "K")
                np.testing.assert_allclose(value[:2], [280.4, 281.1], rtol=1e-6)
                self.assertAlmostEqual(
                    float(value._FillValue), 9.969209968386869e36, delta=1e30
                )
                self.assertTrue(np.ma.is_masked(value[2]))
                np.testing.assert_allclose(
                    dataset["ObsError"]["airTemperature"][:], [0.5, 0.6, 0.7], rtol=1e-6
                )
                qc = dataset["PreQC"]["airTemperature"]
                self.assertEqual(qc.dtype, np.int32)
                np.testing.assert_array_equal(qc[:], [0, 1, 0])
                np.testing.assert_array_equal(dataset["Location"][:], [0, 1, 2])

    def test_profiles(self):
        """The observations of profiles are numbered by their observation."""
        observations = _observations().assign(observation=[4, 4, 7])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "profiles.nc4"
            to_ioda(
                observations,
                path,
                variables={"temperature": "airTemperature"},
                metadata={"platform_az": "platformAzimuth"},
            )
            with netCDF4.Dataset(path) as dataset:
                np.testing.assert_array_equal(
                    dataset["MetaData"]["sequenceNumber"][:], [0, 0, 1]
                )
                self.assertIn("platformAzimuth", dataset["MetaData"].variables)
                # the error of "temperature" is found by its name
                self.assertIn("airTemperature", dataset["ObsError"].variables)

    def test_missing_columns(self):
        """Exported columns must exist."""
        with self.assertRaisesRegex(ValueError, "no columns"):
            to_ioda(_observations(), "unused.nc4", variables={"missing": "x"})
