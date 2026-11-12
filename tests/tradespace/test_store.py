"""
Unit tests for persistent stores of results.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import math
import sqlite3
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from datetime import date, datetime, time, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from tatc.tradespace import DesignSpace, Integer, ResultStore, Study, Uncertain


class _Interrupt(Exception):
    """An interruption of a study."""


def _stop_at_four(d):
    if d.a == 4:
        raise ValueError("stop")
    return {"value": d.a}


class TestResultStore(unittest.TestCase):
    """
    Unit tests for stores of results.
    """

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / "results.sqlite"
        self.space = DesignSpace(a=Integer(1, 4), noise=Uncertain.uniform(0, 1))
        self.calls = []
        self.fail_at = None

    def tearDown(self):
        self.directory.cleanup()

    def _study(self, store, **kwargs):
        def evaluate(d):
            if self.fail_at is not None and len(self.calls) == self.fail_at:
                raise _Interrupt()
            self.calls.append((d.a, d.noise))
            return {"value": d.a + d.noise}

        return Study(
            self.space, lambda d: d, evaluate, replicates=3, store=store, **kwargs
        )

    def test_store_and_get(self):
        """Records are stored and retrieved by design and replicate."""
        store = ResultStore(self.path, "x")
        store.put([("d1", 0, [{"value": 1}]), ("d1", 1, [{"value": 2}])])
        self.assertEqual(len(store), 2)
        self.assertEqual(store.get([("d1", 1), ("d2", 0)]), {("d1", 1): [{"value": 2}]})
        # records of another context are separate
        store.put([("d1", 1, [{"value": 3}])], context="c")
        self.assertEqual(store.get([("d1", 1)], "c"), {("d1", 1): [{"value": 3}]})
        self.assertEqual(store.get([("d1", 1)])[("d1", 1)], [{"value": 2}])
        self.assertEqual(store.contexts(), {"": 2, "c": 1})
        other = ResultStore(self.path, "y")
        self.assertEqual(len(other), 0)
        store.clear()
        self.assertEqual(len(store), 0)
        # a store without a namespace has its own records
        self.assertEqual(len(ResultStore(self.path)), 0)

    def test_resume(self):
        """An interrupted study resumes, evaluating only what is not stored."""
        store = ResultStore(self.path, "study")
        self.fail_at = 5
        with self.assertRaises(_Interrupt):
            self._study(store).run()
        self.assertEqual(len(store), 5)
        self.fail_at = None
        resumed = self._study(store).run()
        self.assertEqual(len(self.calls), 12)  # 5 before, 7 after
        complete = self._study(None).run()
        columns = [c for c in complete.columns if c != "seconds"]
        pd.testing.assert_frame_equal(resumed[columns], complete[columns])
        # a later run takes every evaluation from the store
        self.calls.clear()
        self._study(store).run()
        self.assertEqual(self.calls, [])

    def test_more_replicates(self):
        """Extending a study evaluates only the new replicates."""
        store = ResultStore(self.path, "study")
        self._study(store).run()
        self.calls.clear()
        results = self._study(store).run(replicates=range(5))
        self.assertEqual(len(self.calls), 4 * 2)
        self.assertEqual(len(results), 4 * 5)

    def test_parallel(self):
        """Parallel studies store their evaluations as chunks complete."""
        store = ResultStore(self.path, "study")
        study = Study(
            self.space,
            lambda d: d,
            lambda d: {"value": d.a + d.noise},
            replicates=3,
            store=store,
            n_jobs=2,
        )
        results = study.run()
        self.assertEqual(len(store), 12)
        again = study.run()
        pd.testing.assert_frame_equal(results, again)

    def test_other_definitions(self):
        """Records of studies of other definitions are not reused."""
        store = ResultStore(self.path, "study")
        self._study(store).run()
        self.calls.clear()
        self._study(store).run()
        self.assertEqual(self.calls, [])
        # another seed draws other uncertain values
        with self.assertWarns(UserWarning):
            reseeded = self._study(store, seed=1).run()
        self.assertEqual(len(self.calls), 12)
        self.assertNotEqual(
            reseeded.noise.tolist(), self._study(None).run().noise.tolist()
        )
        # another fixed value of the design space
        self.calls.clear()
        self.space = DesignSpace(a=Integer(1, 4), b=2.0, noise=Uncertain.uniform(0, 1))
        with self.assertWarns(UserWarning):
            self._study(store).run()
        self.assertEqual(len(self.calls), 12)
        # another value that the evaluation closes over
        self.calls.clear()

        def study(offset):
            return Study(
                self.space,
                lambda d: d,
                lambda d, offset=offset: self.calls.append(1) or {"value": offset},
                store=ResultStore(self.path, "closure"),
            )

        study(0.0).run()
        with self.assertWarnsRegex(UserWarning, r"defaults\[0\]: 0.0 -> 1.0"):
            study(1.0).run()
        self.assertEqual(len(self.calls), 8)

    def test_records_as_json(self):
        """Records are kept as JSON, and read as they were written."""
        record = {
            "float": 1.5,
            "nan": float("nan"),
            "infinite": -math.inf,
            "numpy": np.float64(2.5),
            "integer": np.int64(3),
            "flag": np.bool_(True),
            "none": None,
            "when": datetime(2026, 1, 1, 12, tzinfo=timezone.utc),
            "naive": datetime(2026, 1, 1, 12, 30),
            "day": date(2026, 1, 2),
            "time": time(13, 30),
            "duration": timedelta(days=1, microseconds=5),
            "pair": (1, "a"),
            "list": [1, (2, 3)],
            "nested": {"a": {"b": 1}},
            "keys": {1: "a", ("x", 2): "b"},
            "tagged": {"$float": "x"},
            "array": np.arange(3),
        }
        store = ResultStore(self.path, "x")
        store.put([("d1", 0, [record])])
        read = store.get([("d1", 0)])[("d1", 0)][0]
        self.assertEqual(set(read), set(record))
        self.assertTrue(math.isnan(read.pop("nan")))
        np.testing.assert_array_equal(read.pop("array"), record["array"])
        for key, value in read.items():
            # numpy scalars are read as Python's
            written = record[key]
            written = written.item() if isinstance(written, np.generic) else written
            self.assertEqual(value, written, key)
            self.assertIs(type(value), type(written), key)
        # the records are JSON, which SQLite reads
        with closing(sqlite3.connect(self.path)) as connection:
            text, value = connection.execute(
                "SELECT records, json_extract(records, '$[0].float') FROM results"
            ).fetchone()
        self.assertEqual(json.loads(text)[0]["integer"], 3)
        self.assertEqual(value, 1.5)

    def test_definitions(self):
        """A store keeps the definitions of studies, which explain changes."""
        store = ResultStore(self.path, "study")
        study = self._study(store)
        study.run()
        self.assertEqual(store.definition(study.fingerprint), study.definition)
        self.assertIsNone(store.definition("other"))
        reseeded = self._study(store, seed=1)
        with self.assertWarnsRegex(UserWarning, "latest differs in: seed: 0 -> 1"):
            reseeded.run()
        # contexts in the order they were defined, the latest last
        self.assertEqual(
            list(store.contexts()), [study.fingerprint, reseeded.fingerprint]
        )
        # the records of other definitions are removed
        self.assertEqual(store.prune(reseeded.fingerprint), 12)
        self.assertEqual(store.contexts(), {reseeded.fingerprint: 12})
        self.assertIsNone(store.definition(study.fingerprint))
        store.clear()
        self.assertIsNone(store.definition(reseeded.fingerprint))

    def test_parallel_interruption(self):
        """Workers store each evaluation as it completes, before an error."""
        store = ResultStore(self.path, "interrupted")
        # an executor that waits for the chunks in progress as it shuts down
        with ThreadPoolExecutor(2) as executor:
            study = Study(
                DesignSpace(a=Integer(1, 8)),
                lambda d: d,
                _stop_at_four,
                store=store,
                executor=executor,
                chunks_per_job=1,
            )
            with self.assertRaises(ValueError):
                study.run()
        # the designs evaluated before the error in its chunk (1 to 4) are kept
        designs = study.space.enumerate()
        stored = store.get(((d.key, 0) for d in designs), study.fingerprint)
        values = sorted(records[0]["value"] for records in stored.values())
        self.assertEqual(values[:3], [1, 2, 3])
        self.assertNotIn(4, values)

    def test_variants(self):
        """Variants of replicates are kept apart from plain replicates."""
        store = ResultStore(self.path, "x")
        store.put([("d1", 0, [{"value": 1}])])
        store.put([("d1", 0, [{"value": 2}])], variant='{"u": 5}')
        self.assertEqual(store.get([("d1", 0)]), {("d1", 0): [{"value": 1}]})
        self.assertEqual(
            store.get([("d1", 0)], variant='{"u": 5}'), {("d1", 0): [{"value": 2}]}
        )
        self.assertEqual(store.get([("d1", 0)], variant='{"u": 6}'), {})
        self.assertEqual(len(store), 2)

    def test_shifted_runs(self):
        """Shifted replicates are stored and reused apart from plain ones."""
        store = ResultStore(self.path, "study")
        plain = self._study(store).run()
        shifted = self._study(store).run(shift={"noise": 3})
        self.assertEqual(len(self.calls), 2 * 4 * 3)
        self.assertEqual(list(shifted.replicate), list(plain.replicate))
        later = self._study(store).run(replicates=range(3, 6))
        np.testing.assert_array_equal(shifted.noise, later.noise)
        self.calls.clear()
        again = self._study(store).run(shift={"noise": 3})
        self.assertEqual(self.calls, [])
        pd.testing.assert_frame_equal(
            again.drop(columns="seconds"), shifted.drop(columns="seconds")
        )
        unshifted = self._study(store).run(shift={"noise": 0})
        self.assertEqual(self.calls, [])
        pd.testing.assert_frame_equal(
            unshifted.drop(columns="seconds"), plain.drop(columns="seconds")
        )
        with self.assertRaises(ValueError):
            self._study(store).run(shift={"other": 1})

    def test_table_without_variants(self):
        """Evaluations stored before replicates had variants are plain."""
        with closing(sqlite3.connect(self.path)) as connection:
            with connection:
                connection.execute(
                    "CREATE TABLE results (namespace TEXT, context TEXT, design TEXT, "
                    + "replicate INTEGER, records TEXT, "
                    + "PRIMARY KEY (namespace, context, design, replicate))"
                )
                connection.execute(
                    "INSERT INTO results VALUES ('x', '', 'd1', 0, '[{\"value\": 1}]')"
                )
        store = ResultStore(self.path, "x")
        self.assertEqual(store.get([("d1", 0)]), {("d1", 0): [{"value": 1}]})
        store.put([("d1", 0, [{"value": 2}])], variant="v")
        self.assertEqual(len(store), 2)
        self.assertEqual(len(ResultStore(self.path, "x")), 2)

    def test_earlier_format(self):
        """Files with evaluations of an earlier format are read as empty."""
        with closing(sqlite3.connect(self.path)) as connection:
            with connection:
                connection.execute(
                    "CREATE TABLE evaluations (namespace TEXT, context TEXT, "
                    + "design TEXT, replicate INTEGER, records BLOB)"
                )
                connection.execute(
                    "INSERT INTO evaluations VALUES ('x', '', 'd1', 0, x'00')"
                )
        store = ResultStore(self.path, "x")
        self.assertEqual(len(store), 0)
        store.put([("d1", 0, [{"value": 1}])])
        self.assertEqual(store.get([("d1", 0)]), {("d1", 0): [{"value": 1}]})
