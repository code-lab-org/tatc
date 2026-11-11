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
        with self.assertRaises(ValueError):
            ResultStore(self.path, "")

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
        for offset in (0.0, 1.0):
            Study(
                self.space,
                lambda d: d,
                lambda d, offset=offset: self.calls.append(1) or {"value": offset},
                store=ResultStore(self.path, "closure"),
            ).run()
        self.assertEqual(len(self.calls), 8)
