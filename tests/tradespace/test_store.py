"""
Unit tests for persistent stores of results.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from tatc.tradespace import DesignSpace, Integer, ResultStore, Study, Uncertain


class _Interrupt(Exception):
    """An interruption of a study."""


class TestResultStore(unittest.TestCase):
    """
    Unit tests for stores of results.
    """

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name) / "results.sqlite"
        self.space = DesignSpace(a=Integer(1, 4), noise=Uncertain.uniform(0, 1))
        self.calls = []

    def tearDown(self):
        self.directory.cleanup()

    def _study(self, store, fail_at=None, **kwargs):
        def evaluate(d):
            if fail_at is not None and len(self.calls) == fail_at:
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
        other = ResultStore(self.path, "y")
        self.assertEqual(len(other), 0)
        store.clear()
        self.assertEqual(len(store), 0)
        with self.assertRaises(ValueError):
            ResultStore(self.path, "")

    def test_resume(self):
        """An interrupted study resumes, evaluating only what is not stored."""
        store = ResultStore(self.path, "study")
        with self.assertRaises(_Interrupt):
            self._study(store, fail_at=5).run()
        self.assertEqual(len(store), 5)
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
