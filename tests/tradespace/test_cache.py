"""
Unit tests for caches of analyses by satellite.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import importlib.util
import pickle
import tempfile
import unittest
import warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

from tatc.analysis import collect_observations
from tatc.schemas import Instrument, Point, Satellite, SunSynchronousOrbit
from tatc.tradespace import (
    DesignSpace,
    Integer,
    SatelliteCache,
    Study,
)
from tatc.tradespace.cache import get_satellite_key

HAS_JOBLIB = importlib.util.find_spec("joblib") is not None

EPOCH = datetime(2026, 1, 1, tzinfo=timezone.utc)


def _satellite(name: str, anomaly: float = 0, altitude: float = 700e3) -> Satellite:
    return Satellite(
        name=name,
        orbit=SunSynchronousOrbit(
            altitude=altitude,
            equator_crossing_time="13:30",
            epoch=EPOCH,
            true_anomaly=anomaly,
        ),
        instruments=[Instrument(name="imager", field_of_regard=100)],
    )


def _name(satellite: Satellite) -> str:
    return satellite.name


def _frame(satellite: Satellite) -> pd.DataFrame:
    return pd.DataFrame({"satellite": [satellite.name] * 2})


class TestSatelliteKey(unittest.TestCase):
    """
    Unit tests for satellite keys.
    """

    def test_name_is_ignored(self):
        """Satellites that differ only in name share a key."""
        self.assertEqual(
            get_satellite_key(_satellite("a")), get_satellite_key(_satellite("b"))
        )

    def test_content_is_keyed(self):
        """Satellites with other orbits or instruments have other keys."""
        base = _satellite("a")
        self.assertNotEqual(
            get_satellite_key(base), get_satellite_key(_satellite("a", anomaly=90))
        )
        other = base.model_copy(
            update={"instruments": [Instrument(name="imager", field_of_regard=90)]}
        )
        self.assertNotEqual(get_satellite_key(base), get_satellite_key(other))


class TestSatelliteCache(unittest.TestCase):
    """
    Unit tests for caches of analyses by satellite.
    """

    def setUp(self):
        self.calls = []

        def analysis(satellite):
            self.calls.append(satellite.name)
            return pd.DataFrame(
                {"satellite": [satellite.name] * 2, "value": [1.0, 2.0]}
            )

        self.analysis = analysis

    def test_hits_and_relabeling(self):
        """Results are computed once per content and relabeled with the name."""
        cache = SatelliteCache(self.analysis)
        first = cache([_satellite("a"), _satellite("b", anomaly=90)])
        second = cache([_satellite("c"), _satellite("d", anomaly=90)])
        self.assertEqual(self.calls, ["a", "b"])
        self.assertEqual((cache.hits, cache.misses), (2, 2))
        self.assertEqual(list(first.satellite), ["a", "a", "b", "b"])
        self.assertEqual(list(second.satellite), ["c", "c", "d", "d"])
        # the cached result is not changed by relabeling
        self.assertEqual(list(cache([_satellite("a")]).satellite), ["a", "a"])

    def test_disk_cache(self):
        """Results on disk are shared by caches of the same namespace."""
        with tempfile.TemporaryDirectory() as directory:
            SatelliteCache(self.analysis, directory, namespace="x")(_satellite("a"))
            again = SatelliteCache(self.analysis, directory, namespace="x")
            self.assertEqual(list(again(_satellite("b")).satellite), ["b", "b"])
            self.assertEqual((again.hits, again.misses), (1, 0))
            other = SatelliteCache(self.analysis, directory, namespace="y")
            other(_satellite("c"))
            self.assertEqual(other.misses, 1)
        self.assertEqual(self.calls, ["a", "c"])

    def test_disk_cache_of_other_functions(self):
        """Results on disk are not shared by functions of other definitions."""

        def make(rows):
            return SatelliteCache(
                lambda satellite, rows=rows: pd.DataFrame(
                    {"satellite": [satellite.name] * rows}
                ),
                directory,
                namespace="x",
            )

        with tempfile.TemporaryDirectory() as directory:
            for rows in (1, 1):
                self.assertEqual(len(make(rows)(_satellite("a"))), rows)
            # a namespace with results only of other definitions warns
            with self.assertWarnsRegex(
                UserWarning, r"1 results.*defaults\[0\]: 1 -> 2"
            ):
                cache = make(2)
            # once, of the caches of the definition (such as of parallel workers)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                make(2)
            self.assertEqual(caught, [])
            self.assertEqual(len(cache(_satellite("a"))), 2)
            namespaces = list(Path(directory).iterdir())
            self.assertEqual(len(namespaces), 1)
            self.assertEqual((namespaces[0] / "namespace.txt").read_text(), "x")
            folders = [p for p in namespaces[0].iterdir() if p.is_dir()]
            self.assertEqual(len(folders), 2)
            self.assertTrue(all((p / "definition.json").exists() for p in folders))
            # results of other definitions are removed
            self.assertEqual(cache.prune(), 1)
            self.assertEqual(len([p for p in namespaces[0].iterdir() if p.is_dir()]), 1)
            self.assertEqual(SatelliteCache(_name).prune(), 0)  # in memory

    def test_disk_cache_of_other_data(self):
        """Results on disk are keyed by the data the function closes over."""

        def make(frame):
            return SatelliteCache(
                lambda satellite: frame.assign(satellite=satellite.name), directory
            )

        with tempfile.TemporaryDirectory() as directory:
            make(pd.DataFrame({"value": [1.0]}))(_satellite("a"))
            same = make(pd.DataFrame({"value": [1.0]}))
            same(_satellite("b"))
            self.assertEqual((same.hits, same.misses), (1, 0))
            with self.assertWarns(UserWarning):
                other = make(pd.DataFrame({"value": [2.0]}))
            self.assertEqual(list(other(_satellite("c")).value), [2.0])
            self.assertEqual(other.misses, 1)

    def test_unreadable_file(self):
        """A result that cannot be read from disk is computed again and replaced."""
        with tempfile.TemporaryDirectory() as directory:
            cache = SatelliteCache(self.analysis, directory, namespace="x")
            cache(_satellite("a"))
            for path in Path(directory).rglob("*.pkl"):
                path.write_bytes(b"not a pickle")
            again = SatelliteCache(self.analysis, directory, namespace="x")
            self.assertEqual(list(again(_satellite("b")).satellite), ["b", "b"])
            self.assertEqual(again.misses, 1)
            # the file is replaced
            third = SatelliteCache(self.analysis, directory, namespace="x")
            third(_satellite("c"))
            self.assertEqual((third.hits, third.misses), (1, 0))

    @unittest.skipUnless(HAS_JOBLIB, "requires joblib")
    def test_disk_cache_in_parallel(self):
        """Workers of a study share a cache on disk, storing each result once."""
        satellites = [_satellite(f"s{i}", anomaly=30 * (i % 4)) for i in range(40)]
        with tempfile.TemporaryDirectory() as directory:
            cache = SatelliteCache(_frame, directory, namespace="parallel")
            study = Study(
                DesignSpace(i=Integer(0, 39)),
                lambda d: satellites[d.i],
                lambda satellite: {"rows": len(cache(satellite))},
                n_jobs=4,
                chunks_per_job=5,
            )
            results = study.run()
            self.assertTrue((results.rows == 2).all())
            files = sorted(p.suffix for p in Path(directory).rglob("*.*"))
            self.assertEqual(files, [".json"] + [".pkl"] * 4 + [".txt"])

    def test_pickled_cache_is_empty(self):
        """Copies of a cache (as sent to workers) start without results."""
        cache = SatelliteCache(_name)
        cache.get(_satellite("a"))
        copy = pickle.loads(pickle.dumps(cache))
        self.assertEqual(copy._memory, {})  # pylint: disable=protected-access
        cache.clear()
        self.assertEqual((cache.hits, cache.misses), (0, 0))

    def test_no_satellites(self):
        """No satellites give an empty frame."""
        result = SatelliteCache(self.analysis)([])
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(len(result), 0)

    def test_non_frame_results(self):
        """Results that are not data frames are listed."""
        cache = SatelliteCache(lambda satellite: satellite.name)
        self.assertEqual(cache([_satellite("a"), _satellite("b")]), ["a", "a"])

    def test_matches_analysis(self):
        """Cached observations of satellites match those of the constellation."""
        points = [Point(id=i, latitude=10 * i, longitude=0) for i in range(5)]
        end = EPOCH + timedelta(days=1)
        satellites = [_satellite(f"s{i}", anomaly=120 * i) for i in range(3)]
        cache = SatelliteCache(
            lambda satellite: collect_observations(points, satellite, EPOCH, end)
        )
        cached = cache(satellites).sort_values(["target_hash", "epoch"])
        direct = collect_observations(points, satellites, EPOCH, end).sort_values(
            ["target_hash", "epoch"]
        )
        pd.testing.assert_frame_equal(
            cached.reset_index(drop=True)[direct.columns],
            direct.reset_index(drop=True),
        )
