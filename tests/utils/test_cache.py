"""
Unit tests for the cache utility functions.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import copy
import pickle
import unittest

from pydantic import BaseModel

from tatc.utils.cache import get_cached


class Unpicklable:
    """A value that cannot be pickled (as a Skyfield satellite)."""

    def __reduce__(self):
        raise TypeError("cannot pickle")


class Model(BaseModel):
    """A model on which values are cached."""

    value: float


class TestGetCached(unittest.TestCase):
    """
    Unit tests for get_cached.
    """

    def test_cached_until_key_changes(self):
        """
        Test that a value is computed once per key, and recomputed if the
        key changes.
        """
        model = Model(value=1)
        calls = []

        def compute():
            calls.append(model.value)
            return model.value * 2

        self.assertEqual(get_cached(model, "double", model.value, compute), 2)
        self.assertEqual(get_cached(model, "double", model.value, compute), 2)
        self.assertEqual(len(calls), 1)
        model.value = 3
        self.assertEqual(get_cached(model, "double", model.value, compute), 6)
        self.assertEqual(len(calls), 2)

    def test_cache_not_copied_or_pickled(self):
        """
        Test that cached values (here, unpicklable) are neither pickled nor
        copied with their object, nor shared with a shallow copy.
        """
        model = Model(value=1)
        cached = get_cached(model, "value", None, Unpicklable)
        for other in (
            pickle.loads(pickle.dumps(model)),
            copy.deepcopy(model),
            model.model_copy(deep=True),
            model.model_copy(),
        ):
            self.assertEqual(other, model)
            self.assertIsNot(get_cached(other, "value", None, Unpicklable), cached)
        self.assertIs(get_cached(model, "value", None, Unpicklable), cached)
