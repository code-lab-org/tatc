"""
Unit tests for fingerprints of the definitions of trade studies.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import os
import subprocess
import sys
import unittest
from unittest import mock
from datetime import datetime, timedelta, timezone
from functools import partial

import numpy as np
import pandas as pd

from tatc.schemas import Instrument
from tatc.tradespace import (
    CostModel,
    DesignSpace,
    Integer,
    Real,
    SatelliteCost,
    Study,
    Uncertain,
    get_fingerprint,
    get_source_hash,
)

THRESHOLD = 4


def _helper(x):
    return x * THRESHOLD


def _uses_helper(d):
    return {"value": _helper(d.a)}


MODEL = None
CACHES = {}


def _uses_cache(d):
    if d not in CACHES:
        CACHES[d] = object()
    return CACHES[d]


def _uses_model(d):
    return {"cost": MODEL.estimate(d.layout).total}


def _space(high=4, fixed=1.0, spread=1.0):
    return DesignSpace(
        a=Integer(1, high),
        b=Real(0, lambda d: d.a * 2, step=0.5),
        fixed=fixed,
        noise=Uncertain.normal(0, spread),
    )


_SCRIPT = """
from tatc.tradespace import DesignSpace, Integer, Study, Uncertain
space = DesignSpace(a=Integer(1, 4), names={"x", "y", "z"}, noise=Uncertain.uniform(0, 1))
study = Study(space, lambda d: d, lambda d: {"v": d.a in {"p", "q"}})
print(study.fingerprint)
"""


class TestFingerprint(unittest.TestCase):
    """
    Unit tests for fingerprints.
    """

    def test_same_definitions(self):
        """Separately made definitions that are the same have the same fingerprint."""
        self.assertEqual(get_fingerprint(_space()), get_fingerprint(_space()))
        self.assertEqual(
            get_fingerprint(lambda d: d.a + 1), get_fingerprint(lambda d: d.a + 1)
        )

    def test_space_definitions(self):
        """Bounds, fixed values, and distributions change the fingerprint."""
        base = get_fingerprint(_space())
        self.assertNotEqual(base, get_fingerprint(_space(high=5)))
        self.assertNotEqual(base, get_fingerprint(_space(fixed=2.0)))
        self.assertNotEqual(base, get_fingerprint(_space(spread=2.0)))

    def test_functions(self):
        """Code, closures, defaults, and plain globals change the fingerprint."""
        self.assertNotEqual(
            get_fingerprint(lambda d: d.a + 1), get_fingerprint(lambda d: d.a + 2)
        )
        fingerprints = set()
        for scale in (1, 2):
            fingerprints.add(get_fingerprint(lambda d, s=scale: d.a * s))
            fingerprints.add(get_fingerprint(partial(max, scale)))
        self.assertEqual(len(fingerprints), 4)
        global THRESHOLD  # pylint: disable=global-statement
        before = get_fingerprint(_uses_helper)
        THRESHOLD = 5
        try:
            # a plain global of a function of the same module that it calls
            self.assertNotEqual(before, get_fingerprint(_uses_helper))
        finally:
            THRESHOLD = 4

    def test_values(self):
        """Plain values, models, and arrays are described by content."""
        start = datetime(2026, 1, 1, tzinfo=timezone.utc)
        self.assertNotEqual(
            get_fingerprint(start), get_fingerprint(start + timedelta(seconds=1))
        )
        self.assertNotEqual(
            get_fingerprint(Instrument(name="a", field_of_regard=30)),
            get_fingerprint(Instrument(name="a", field_of_regard=40)),
        )
        self.assertNotEqual(
            get_fingerprint(np.arange(3)), get_fingerprint(np.arange(3) + 1)
        )
        self.assertEqual(get_fingerprint({1, 2, 3}), get_fingerprint({3, 2, 1}))
        # other objects are described by their type only
        self.assertEqual(
            get_fingerprint(pd.DataFrame({"a": [1]})),
            get_fingerprint(pd.DataFrame({"a": [2]})),
        )

    def test_global_state(self):
        """Mutable globals (such as caches) do not change the fingerprint."""
        before = get_fingerprint(_uses_cache)
        _uses_cache(1)
        self.assertEqual(before, get_fingerprint(_uses_cache))

    def test_dataclasses(self):
        """Dataclasses, such as cost models, are described by their fields."""
        cheap = CostModel(SatelliteCost(first_unit=10))
        dear = CostModel(SatelliteCost(first_unit=12))
        self.assertNotEqual(get_fingerprint(cheap), get_fingerprint(dear))
        global MODEL  # pylint: disable=global-statement
        MODEL = cheap
        before = get_fingerprint(_uses_model)
        MODEL = dear
        # a dataclass that a function refers to as a global
        self.assertNotEqual(before, get_fingerprint(_uses_model))

    def test_study(self):
        """A study's fingerprint changes with its seed and random numbers."""

        def study(**kwargs):
            return Study(_space(), lambda d: d, lambda d: {"value": d.a}, **kwargs)

        fingerprints = {
            study().fingerprint,
            study(seed=1).fingerprint,
            study(common_random_numbers=False).fingerprint,
        }
        self.assertEqual(len(fingerprints), 3)
        # but not with the settings of its execution
        self.assertEqual(study().fingerprint, study(n_jobs=2).fingerprint)

    def test_source(self):
        """The source of TAT-C (its code and resources) changes the fingerprint."""
        self.assertEqual(len(get_source_hash()), 16)
        self.assertIs(get_source_hash(), get_source_hash())  # computed once
        before = get_fingerprint(_space())
        with mock.patch(
            "tatc.tradespace.fingerprint.get_source_hash", return_value="changed"
        ):
            self.assertNotEqual(before, get_fingerprint(_space()))
        self.assertEqual(before, get_fingerprint(_space()))

    def test_stable_between_sessions(self):
        """Fingerprints do not change between sessions (with other hash seeds)."""
        prints = set()
        for seed in ("1", "2"):
            output = subprocess.run(
                [sys.executable, "-c", _SCRIPT],
                capture_output=True,
                text=True,
                check=True,
                env={**os.environ, "PYTHONHASHSEED": seed},
            )
            prints.add(output.stdout.strip())
        self.assertEqual(len(prints), 1)
