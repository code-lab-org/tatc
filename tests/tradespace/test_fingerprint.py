"""
Unit tests for fingerprints of the definitions of trade studies.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import ast
import os
import subprocess
import sys
import unittest
from unittest import mock
from datetime import datetime, timedelta, timezone
from functools import partial

import numpy as np
import pandas as pd

from shapely.geometry import Point

from tatc.analysis import collect_observations
from tatc.generation import generate_points_rotated_lattice
from tatc.schemas import Instrument
from tatc.tradespace import fingerprint
from tatc.tradespace import (
    CostModel,
    DesignSpace,
    Integer,
    Real,
    SatelliteCost,
    Study,
    Uncertain,
)
from tatc.tradespace.fingerprint import (
    compare_definitions,
    get_definition,
    get_fingerprint,
)

THRESHOLD = 4


def _helper(x):
    return x * THRESHOLD


def _uses_helper(d):
    return {"value": _helper(d.a)}


MODEL = None
TABLE = pd.DataFrame({"a": [1.0]})


class _Opaque:
    """An object described by its type only."""

    def __init__(self, value):
        self.value = value


def _uses_table(d):
    return {"value": TABLE.a.sum() * d.a}


def _uses_attribute(d):
    # an attribute named as the global `values`, which the function does not use
    return {"value": d.values}


values = pd.DataFrame({"a": [1.0]})  # pylint: disable=invalid-name


CACHES = {}


def _uses_cache(d):
    if d not in CACHES:
        CACHES[d] = object()
    return CACHES[d]


def _uses_model(d):
    return {"cost": MODEL.estimate(d.layout).total}


START = datetime(2026, 1, 1, tzinfo=timezone.utc)


def _observe(satellite):
    return collect_observations([], satellite, START, START + timedelta(days=1))


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
        # data frames and series by their content (and geometries by theirs)
        self.assertEqual(
            get_fingerprint(pd.DataFrame({"a": [1]})),
            get_fingerprint(pd.DataFrame({"a": [1]})),
        )
        self.assertNotEqual(
            get_fingerprint(pd.DataFrame({"a": [1]})),
            get_fingerprint(pd.DataFrame({"a": [2]})),
        )
        self.assertNotEqual(
            get_fingerprint(pd.Series([1.0], name="a")),
            get_fingerprint(pd.Series([1.0], name="b")),
        )
        points = generate_points_rotated_lattice(5, seed=1)
        moved = points.copy()
        moved.loc[0, "geometry"] = Point(0, 1e-9)
        self.assertEqual(
            get_fingerprint(points),
            get_fingerprint(generate_points_rotated_lattice(5, seed=1)),
        )
        self.assertNotEqual(get_fingerprint(points), get_fingerprint(moved))
        # other objects are described by their type only
        self.assertEqual(get_fingerprint(_Opaque(1)), get_fingerprint(_Opaque(2)))

    def test_attributes_are_not_globals(self):
        """Attributes named as globals do not refer to the globals."""
        global values  # pylint: disable=global-statement,invalid-name
        before = get_fingerprint(_uses_attribute)
        values = pd.DataFrame({"a": [2.0]})
        try:
            self.assertEqual(before, get_fingerprint(_uses_attribute))
        finally:
            values = pd.DataFrame({"a": [1.0]})

    def test_data_globals(self):
        """Data that functions refer to as globals is described by content."""
        global TABLE  # pylint: disable=global-statement
        before = get_fingerprint(_uses_table)
        TABLE = pd.DataFrame({"a": [2.0]})
        try:
            self.assertNotEqual(before, get_fingerprint(_uses_table))
        finally:
            TABLE = pd.DataFrame({"a": [1.0]})

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

    def test_modules(self):
        """A definition has the modules of TAT-C (and resources) values depend on."""
        modules = get_definition(_observe)["modules"]
        # the module that defines a function imported from a package, not
        # the package, and the modules that it imports
        self.assertIn("tatc.analysis.point_sampling", modules)
        self.assertIn("tatc.schemas.space.satellite", modules)
        self.assertNotIn("tatc.analysis", modules)
        self.assertNotIn("tatc.analysis.ro_sampling", modules)
        self.assertFalse(any(m.startswith("tatc.tradespace") for m in modules))
        self.assertIn("resources/defaults.yml", get_definition(_observe)["resources"])
        # studies depend on the modules of studies and their design spaces
        study = Study(_space(), lambda d: d, lambda d: {"value": d.a}).definition
        self.assertIn("tatc.tradespace.study", study["modules"])
        self.assertIn("tatc.tradespace.nodes", study["modules"])
        self.assertNotIn("tatc.tradespace.pareto", study["modules"])
        # plain values depend on no modules
        self.assertEqual(get_definition(1, a="x")["modules"], {})

    def test_source(self):
        """The source of the modules values depend on changes the fingerprint."""
        original = fingerprint._hash  # pylint: disable=protected-access

        def edited(name):
            return lambda path: "edited" if path.name == name else original(path)

        before = get_fingerprint(_observe)
        with mock.patch.object(fingerprint, "_hash", edited("point_sampling.py")):
            self.assertNotEqual(before, get_fingerprint(_observe))
        with mock.patch.object(fingerprint, "_hash", edited("pareto.py")):
            self.assertEqual(before, get_fingerprint(_observe))
        self.assertEqual(before, get_fingerprint(_observe))

    def test_imports_resolve(self):
        """Every import of a module of TAT-C resolves to modules of TAT-C."""
        # pylint: disable=protected-access
        files = fingerprint._module_files()
        self.assertIn("tatc.tradespace.fingerprint", files)
        for module in files:
            for _, other in fingerprint._imports(module):
                self.assertIn(other, files, f"{other} (imported by {module})")
        modules, resources = fingerprint._get_dependencies(files)
        self.assertEqual(set(modules), set(files))
        self.assertTrue(resources)
        # a constant of a package depends on its source only
        self.assertEqual(
            fingerprint._resolve("tatc", "__version__"), {("file", "tatc")}
        )

    def test_static_imports(self):
        """Modules of TAT-C import modules by statements only (not by name)."""
        # pylint: disable=protected-access
        for module, path in fingerprint._module_files().items():
            tree = ast.parse(path.read_bytes())
            for node in ast.walk(tree):
                if isinstance(node, ast.Call):
                    name = getattr(node.func, "id", getattr(node.func, "attr", None))
                    self.assertNotIn(name, ("import_module", "__import__"), module)

    def test_compare(self):
        """Comparisons of definitions name what differs."""

        def study(high=4, offset=0.0, seed=0):
            return Study(
                _space(high=high),
                lambda d: d,
                lambda d, offset=offset: {"value": d.a + offset},
                seed=seed,
            ).definition

        base = study()
        self.assertEqual(compare_definitions(base, study()), [])
        self.assertEqual(
            compare_definitions(base, study(high=5)),
            ["space.children.a.high: 4 -> 5"],
        )
        self.assertEqual(
            compare_definitions(base, study(offset=1.0)),
            ["evaluate.defaults[0]: 0.0 -> 1.0"],
        )
        self.assertEqual(
            compare_definitions(base, study(high=5, seed=1), limit=1),
            ["seed: 0 -> 1", "1 more"],
        )
        other = {
            **base,
            "python": "cpython 0.0",
            "modules": {**base["modules"], "tatc.tradespace.nodes": "edited"},
        }
        self.assertEqual(
            compare_definitions(base, other),
            [
                f"python: {base['python']} -> cpython 0.0",
                "module tatc.tradespace.nodes",
            ],
        )
        self.assertTrue(base["python"].endswith(f"{sys.version_info[1]}"))

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
