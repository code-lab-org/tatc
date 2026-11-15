"""
Unit tests for fingerprints of the definitions of trade studies.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

import ast
import functools
import importlib
import os
import tempfile
import subprocess
import sys
import unittest
import warnings
from types import SimpleNamespace
from unittest import mock
from datetime import datetime, timedelta, timezone
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

from shapely.geometry import Point, box

from tatc.analysis import collect_observations
from tatc.generation import generate_points_fibonacci_lattice
from tatc.schemas import Instrument
from tatc.tradespace import fingerprint
from tatc.tradespace import (
    CostModel,
    DesignSpace,
    Integer,
    Real,
    SatelliteCache,
    SatelliteCost,
    Study,
    Uncertain,
)
from tatc.tradespace.fingerprint import (
    compare_definitions,
    describe,
    get_definition,
    get_fingerprint,
    warn_stale,
)

THRESHOLD = 4


def _helper(x):
    return x * THRESHOLD


def _uses_helper(d):
    return {"value": _helper(d.a)}


MODEL = None
TABLE = pd.DataFrame({"a": [1.0]})


class _Own:
    """An object of a class of the study's own, described by its attributes."""

    def __init__(self, value):
        self.value = value


def _uses_table(d):
    return {"value": TABLE.a.sum() * d.a}


def _uses_attribute(d):
    # an attribute named as the global `values`, which the function does not use
    return {"value": d.values}


values = pd.DataFrame({"a": [1.0]})  # pylint: disable=invalid-name


class _Cache:
    """State of the study's own, kept in an attribute (not described)."""

    def __init__(self):
        self.results = {}


CACHE = _Cache()


def _uses_cache(d):
    if d not in CACHE.results:
        CACHE.results[d] = object()
    return CACHE.results[d]


PAYLOAD = Instrument(name="imager", field_of_regard=30)
INSTRUMENTS = {30: Instrument(name="a", field_of_regard=30)}
SCORER = _Own(3)
METRIC = partial(max, 3)


def _uses_payload(d):
    return PAYLOAD.field_of_regard * d


def _uses_instruments(d):
    return INSTRUMENTS[d]


def _uses_scorer(d):
    return SCORER.value * d


def _uses_metric(d):
    return METRIC(d)


@functools.lru_cache(maxsize=None)
def _cached(x):
    return x * THRESHOLD


def _uses_cached(d):
    return _cached(d)


def _uses_model(d):
    return {"cost": MODEL.estimate(d.layout).total}


def _offset_cache(offset):
    """A cache of an analysis of satellites that gives an offset."""
    return SatelliteCache(lambda satellite: offset)


SATELLITES = _offset_cache(0.0)
BASE = DesignSpace(n=Integer(1, 2)).design({"n": 1})


def _uses_satellites(satellites):
    return {"m": SATELLITES(satellites)}


def _uses_base(d):
    return {"m": BASE.n * d}


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
        points = generate_points_fibonacci_lattice(count=5, rotate=True, seed=1)
        moved = points.copy()
        moved.loc[0, "geometry"] = Point(0, 1e-9)
        self.assertEqual(
            get_fingerprint(points),
            get_fingerprint(
                generate_points_fibonacci_lattice(count=5, rotate=True, seed=1)
            ),
        )
        self.assertNotEqual(get_fingerprint(points), get_fingerprint(moved))
        # other objects are described by their type only
        # objects of installed classes by their type, and of the study's own
        # classes by their attributes (other than mutable containers)
        self.assertEqual(
            get_fingerprint(SimpleNamespace(value=1)),
            get_fingerprint(SimpleNamespace(value=2)),
        )
        self.assertNotEqual(get_fingerprint(_Own(1)), get_fingerprint(_Own(2)))

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
        """State in the attributes of objects (such as caches) does not change
        the fingerprint, but global containers (such as of instruments) do."""
        before = get_fingerprint(_uses_cache)
        _uses_cache(1)
        self.assertEqual(before, get_fingerprint(_uses_cache))
        global INSTRUMENTS  # pylint: disable=global-statement
        before = get_fingerprint(_uses_instruments)
        INSTRUMENTS = {30: Instrument(name="a", field_of_regard=45)}
        try:
            self.assertNotEqual(before, get_fingerprint(_uses_instruments))
        finally:
            INSTRUMENTS = {30: Instrument(name="a", field_of_regard=30)}

    def test_globals_by_content(self):
        """Globals are described as the values they are: models, objects of the
        study's own classes, partial functions, and decorated functions."""
        global THRESHOLD  # pylint: disable=global-statement
        cases = [
            ("PAYLOAD", _uses_payload, Instrument(name="imager", field_of_regard=60)),
            ("SCORER", _uses_scorer, _Own(5)),
            ("METRIC", _uses_metric, partial(max, 5)),
        ]
        for name, function, changed in cases:
            with self.subTest(name=name):
                original = globals()[name]
                before = get_fingerprint(function)
                globals()[name] = changed
                try:
                    self.assertNotEqual(before, get_fingerprint(function))
                finally:
                    globals()[name] = original
        # a function that a decorator wraps, by the wrapped function's code
        before = get_fingerprint(_uses_cached)
        THRESHOLD = 5
        try:
            self.assertNotEqual(before, get_fingerprint(_uses_cached))
        finally:
            THRESHOLD = 4

    def test_one_rule(self):
        """A global is described as the same value given directly is."""
        for name, function in (
            ("PAYLOAD", _uses_payload),
            ("INSTRUMENTS", _uses_instruments),
            ("SCORER", _uses_scorer),
            ("METRIC", _uses_metric),
            ("TABLE", _uses_table),
        ):
            with self.subTest(name=name):
                described = get_definition(f=function)["values"]["f"]["globals"]
                given = get_definition(globals()[name])["values"]["0"]
                self.assertEqual(described[name], given)

    def test_versions(self):
        """The same study with other code is found across versions of Python,
        whose code orders, groups, and repeats constants differently, and a
        study of other literals is not."""
        # pylint: disable=protected-access
        # the constants of the same function in Python 3.10 and 3.14
        older = {
            "code": "aa",
            "constants": [
                "target_hash",
                0,
                "sufficient",
                "0.0",
                100,
                {"code": "bb", "constants": [1], "names": []},
                "evaluate.<locals>.<listcomp>",
                ["sufficient", "n"],
                None,
            ],
            "names": ["groupby"],
        }
        newer = {
            "code": "cc",
            "constants": ["target_hash", "sufficient", "0.0", "n", 0, 100, 100, 1],
            "names": ["groupby", "size"],
        }
        self.assertEqual(
            fingerprint._without_code({"f": older}),
            fingerprint._without_code({"f": newer}),
        )
        other = {**newer, "constants": [*newer["constants"], 4]}
        self.assertNotEqual(
            fingerprint._without_code({"f": older}),
            fingerprint._without_code({"f": other}),
        )
        # of definitions
        definition = {"values": {"f": newer}}
        self.assertEqual(
            fingerprint._find_stale(
                definition,
                {"a": {"values": {"f": older}}, "b": None, "c": {"values": {}}},
                lambda key: 1,
            ),
            {"a": 1},
        )

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
            fingerprint._resolve_import("tatc", "__version__"), {("file", "tatc")}
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


_BOX = box(0, 0, 1, 1)


def _uses_box(d):
    return {"area": _BOX.area * d.a}


def _function(source, **names):
    """A function `f` of source code, with globals."""
    namespace = dict(names)
    exec(source, namespace)  # pylint: disable=exec-used
    return namespace["f"]


class TestDescriptions(unittest.TestCase):
    """
    Unit tests for the descriptions of values that fingerprints hash.
    """

    def test_geometries(self):
        """Geometries are described by their content, also as globals."""

        def closing(geometry):
            return lambda d: geometry.area * d.a

        self.assertNotEqual(
            get_fingerprint(f=closing(box(0, 0, 1, 1))),
            get_fingerprint(f=closing(box(0, 0, 2, 2))),
        )
        self.assertEqual(
            get_fingerprint(f=closing(box(0, 0, 1, 1))),
            get_fingerprint(f=closing(box(0, 0, 1, 1))),
        )
        described = get_definition(f=_uses_box)["values"]["f"]["globals"]["_BOX"]
        self.assertIn("data", described)

    def test_docstrings(self):
        """The docstrings of functions do not change fingerprints."""
        one = _function('def f(x):\n    """One."""\n    return x\n')
        other = _function('def f(x):\n    """Another."""\n    return x\n')
        changed = _function('def f(x):\n    """One."""\n    return x + 1\n')
        self.assertEqual(get_fingerprint(f=one), get_fingerprint(f=other))
        self.assertNotEqual(get_fingerprint(f=one), get_fingerprint(f=changed))

    def test_own_modules(self):
        """Functions of modules that are not installed are described by code."""
        with tempfile.TemporaryDirectory() as directory:
            Path(directory, "tatc_own_helpers.py").write_text(
                "def helper(x):\n    return x + 1\n", encoding="utf-8"
            )
            sys.path.insert(0, directory)
            try:
                module = importlib.import_module("tatc_own_helpers")
                calling = _function(
                    "def f(x):\n    return helper(x)\n", helper=module.helper
                )
                described = get_definition(f=calling)["values"]["f"]["globals"]
                self.assertIn("code", described["helper"])
                # functions of installed modules by their names
                installed = _function("def f(x):\n    return isna(x)\n", isna=pd.isna)
                described = get_definition(f=installed)["values"]["f"]["globals"]
                self.assertEqual(set(described["isna"]), {"name"})
            finally:
                sys.path.remove(directory)
                sys.modules.pop("tatc_own_helpers", None)

    def test_own_module_attributes(self):
        """Functions and values of a module of the study's own that a function
        reaches through the module are described as globals."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory, "tatc_own_attributes.py")
            path.write_text(
                "LIMIT = 2\n\ndef score(x):\n    return x * 1\n", encoding="utf-8"
            )
            sys.path.insert(0, directory)
            try:
                module = importlib.import_module("tatc_own_attributes")
                calling = _function(
                    "def f(x):\n    return helpers.score(x) + helpers.LIMIT\n",
                    helpers=module,
                )
                before = get_fingerprint(f=calling)
                described = get_definition(f=calling)["values"]["f"]["globals"]
                self.assertEqual(
                    set(described["helpers"]["attributes"]), {"score", "LIMIT"}
                )
                # (of another size, so that its compiled code is not reused)
                path.write_text(
                    "LIMIT = 2\n\ndef score(x):\n    return x * 20\n", encoding="utf-8"
                )
                importlib.reload(module)
                self.assertNotEqual(before, get_fingerprint(f=calling))
            finally:
                sys.path.remove(directory)
                sys.modules.pop("tatc_own_attributes", None)


class TestOwnClasses(unittest.TestCase):
    """
    Unit tests for the descriptions of the classes and objects of a study's
    own code.
    """

    def test_objects(self):
        """Callable objects of the study's own classes are described by their
        attributes and the code of their classes' methods."""

        def classes(body):
            namespace = {"__name__": "__main__"}
            exec(  # pylint: disable=exec-used
                "class Score:\n"
                "    def __init__(self, k):\n"
                "        self.k = k\n"
                "        self.calls = []\n"
                "    def __call__(self, d):\n"
                f"        return {body}\n",
                namespace,
            )
            return namespace["Score"]

        score = classes("{'v': d.a * self.k}")
        self.assertNotEqual(get_fingerprint(f=score(1)), get_fingerprint(f=score(3)))
        self.assertEqual(get_fingerprint(f=score(1)), get_fingerprint(f=score(1)))
        # mutable attributes (such as records of calls) are state, not described
        counted = score(1)
        counted.calls.append(1)
        self.assertEqual(get_fingerprint(f=counted), get_fingerprint(f=score(1)))
        # the code of the methods of the class
        edited = classes("{'v': d.a + self.k}")
        self.assertNotEqual(get_fingerprint(f=score(1)), get_fingerprint(f=edited(1)))


class _Described:
    """An object that describes itself by what changes results."""

    def __init__(self, value, calls=0):
        self.value, self.calls = value, calls

    def __fingerprint__(self):
        return {"value": self.value}


class TestSelfDescribed(unittest.TestCase):
    """
    Unit tests for objects that describe themselves (by `__fingerprint__`),
    such as designs and caches of satellites.
    """

    def test_hook(self):
        """An object is described by what its `__fingerprint__` gives, not by
        its other attributes."""
        self.assertEqual(
            get_fingerprint(_Described(1)), get_fingerprint(_Described(1, calls=5))
        )
        self.assertNotEqual(
            get_fingerprint(_Described(1)), get_fingerprint(_Described(2))
        )
        # the class is not described by the hook (it is not an object of it)
        self.assertNotIn("fingerprint", get_definition(_Described)["values"]["0"])

    def test_satellite_cache(self):
        """A study whose function refers to a cache of satellites changes its
        fingerprint with the cache's function, not with its use."""
        global SATELLITES  # pylint: disable=global-statement
        space = DesignSpace(n=Integer(1, 2))

        def fingerprint_of():
            return Study(space, lambda d: [], _uses_satellites).fingerprint

        before = fingerprint_of()
        SATELLITES([])
        SATELLITES.hits += 3
        self.assertEqual(fingerprint_of(), before)
        try:
            SATELLITES = _offset_cache(1.0)
            self.assertNotEqual(fingerprint_of(), before)
        finally:
            SATELLITES = _offset_cache(0.0)
        self.assertEqual(fingerprint_of(), before)

    def test_design(self):
        """A design is described by its design space, decisions, replicate,
        and scenario."""
        global BASE  # pylint: disable=global-statement
        space = DesignSpace(n=Integer(1, 2), u=Uncertain.uniform(0, 1))
        before = get_fingerprint(_uses_base)
        try:
            BASE = BASE.space.design({"n": 2})
            self.assertNotEqual(get_fingerprint(_uses_base), before)
        finally:
            BASE = BASE.space.design({"n": 1})
        self.assertEqual(get_fingerprint(_uses_base), before)
        design = space.design({"n": 1})
        prints = {
            get_fingerprint(design),
            get_fingerprint(space.design({"n": 2})),
            get_fingerprint(design.realize(0)),
            get_fingerprint(design.realize(1)),
            get_fingerprint(design.realize(0, overrides={"u": 0.5}, scenario="half")),
            get_fingerprint(DesignSpace(n=Integer(1, 3)).design({"n": 1})),
        }
        self.assertEqual(len(prints), 6)
        self.assertEqual(
            get_fingerprint(design), get_fingerprint(space.design({"n": 1}))
        )

    def test_describe(self):
        """A value is described as a definition describes it."""
        values = {"limit": 2.0, "names": ("a", "b")}
        self.assertEqual(describe(values), get_definition(v=values)["values"]["v"])
        self.assertEqual(describe(values), {"limit": "2.0", "names": ["a", "b"]})


class TestStale(unittest.TestCase):
    """
    Unit tests for the results of the same study with other code.
    """

    def setUp(self):
        self.older = {"values": {"f": {"code": "aa", "constants": [1], "names": []}}}
        self.newer = {"values": {"f": {"code": "bb", "constants": [1], "names": []}}}
        self.other = {"values": {"f": {"code": "bb", "constants": [2], "names": []}}}

    def test_stale(self):
        """Definitions of the same values but other code are stale, in order,
        with their numbers of results."""
        others = {"b": self.older, "c": self.other, "d": None, "a": self.older}
        counts = {"a": 0, "b": 3, "c": 5, "d": 1}
        # pylint: disable-next=protected-access
        found = fingerprint._find_stale(self.newer, others, counts.__getitem__)
        self.assertEqual(list(found.items()), [("b", 3), ("a", 0)])

    def test_warn_stale(self):
        """A warning names the stale results and the differences of the latest
        of them, once if asked."""
        definitions = {"b": self.older}
        with self.assertWarnsRegex(
            UserWarning, r"^7 things with other code.*\(remove them\).*f\.code"
        ):
            warned = warn_stale(
                {"a": 0, "b": 7},
                self.newer,
                definitions.get,
                lambda count: f"{count} things",
                "remove them",
            )
        self.assertTrue(warned)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            self.assertFalse(
                warn_stale({"a": 0}, self.newer, definitions.get, str, "remove them")
            )
            self.assertFalse(
                warn_stale(
                    {"b": 7},
                    self.newer,
                    definitions.get,
                    str,
                    "remove them",
                    once=lambda: False,
                )
            )
