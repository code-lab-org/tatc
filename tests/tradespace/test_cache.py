"""
Unit tests for caches of analyses by satellite.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import pickle
import tempfile
import unittest
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
    get_satellite_key,
)

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
            with self.assertRaises(ValueError):
                SatelliteCache(self.analysis, directory=directory)
            SatelliteCache(self.analysis, directory, namespace="x")(_satellite("a"))
            again = SatelliteCache(self.analysis, directory, namespace="x")
            self.assertEqual(list(again(_satellite("b")).satellite), ["b", "b"])
            self.assertEqual((again.hits, again.misses), (1, 0))
            other = SatelliteCache(self.analysis, directory, namespace="y")
            other(_satellite("c"))
            self.assertEqual(other.misses, 1)
        self.assertEqual(self.calls, ["a", "c"])

    def test_unreadable_file(self):
        """A result that cannot be read from disk is computed again."""
        with tempfile.TemporaryDirectory() as directory:
            cache = SatelliteCache(self.analysis, directory, namespace="x")
            cache(_satellite("a"))
            for path in Path(directory).iterdir():
                path.write_bytes(b"not a pickle")
            again = SatelliteCache(self.analysis, directory, namespace="x")
            self.assertEqual(list(again(_satellite("b")).satellite), ["b", "b"])
            self.assertEqual(again.misses, 1)

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
            files = sorted(p.suffix for p in Path(directory).iterdir())
            self.assertEqual(files, [".pkl"] * 4)

    def test_pickled_cache_is_empty(self):
        """Copies of a cache (as sent to workers) start without results."""
        cache = SatelliteCache(_name)
        cache.get(_satellite("a"))
        copy = pickle.loads(pickle.dumps(cache))
        self.assertEqual(copy._memory, {})  # pylint: disable=protected-access
        cache.clear()
        self.assertEqual((cache.hits, cache.misses), (0, 0))

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
