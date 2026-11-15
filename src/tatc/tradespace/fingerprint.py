"""
Fingerprints of the definitions of trade studies (design spaces, functions,
and settings), which key persistent results so that results computed with
other definitions are not reused.

A definition (see `get_definition`) describes values as follows (the rules
of `_Describer.describe`), whether they are given (such as the design space
and functions of a study), or are referred to by functions (as globals, the
values they close over, and their defaults), or are within other values:

- By content:

  - plain values (`None`, booleans, numbers, strings, bytes, dates, times,
    and durations) and their containers (lists, tuples, sets, and
    dictionaries, also as globals, such as a dictionary of instruments);
  - pydantic models (such as satellites and instruments);
  - arrays and data frames (by a hash of their data, with geometries by
    their well-known binary), and geometries;
  - dataclasses (such as cost models), by their fields, but those whose
    metadata has `"fingerprint": False` (such as directories of caches,
    which do not change results);
  - the nodes of design spaces, by their attributes;
  - objects that describe themselves by a method `__fingerprint__()`, by
    the value that it returns, which is described in their place: designs
    (by their design space, decisions, replicate, and scenario), caches
    of satellites (by their function and namespace, not their counts of
    hits and misses, nor what they hold), and statistics (by their
    attributes, such as their column and quantile);
  - partial functions (by their function and arguments), bound methods (by
    their function and object), and functions that decorators wrap (such
    as by `functools.lru_cache`, by the wrapped function).

- By code:

  - functions, by their compiled code (which ties a definition to the
    version of Python, which it records, but not their docstrings or line
    numbers, which do not change results), the values they close over,
    their defaults, and the globals they refer to: those of functions of
    their own module or of a module of the study's own (one that is not
    installed, such as `__main__` or a module next to a notebook) are
    described by code, also as attributes of a module of the study's own
    (`helpers.score`), so that a notebook's function that calls a function
    of such a module is described as the module's copy of it is;
  - the classes of the study's own modules, by the code of their methods
    (and their constants), and their objects, by their classes and their
    attributes, but those that are mutable containers (lists, dictionaries,
    and sets), which are often state (such as counters and records of
    calls) rather than definitions.

- By name: the functions of installed modules (as globals of functions),
  classes and built-in functions (but the study's own classes), and modules
  (but the attributes of the study's own modules that functions refer to).

- By type: other objects, such as those of installed classes (a connection,
  or a lock) and those of TAT-C's classes that none of the rules above
  describe, whose state a definition takes not to change results.

- Not at all: data that functions read (such as from files) and the code of
  installed modules other than TAT-C's. Identify them, and the state of
  objects described by type that changes results, by a namespace (see
  `ResultStore` and `SatelliteCache`).

TAT-C is described by the source of the modules that the values depend on:
the modules that define the functions, classes, and objects of TAT-C they
refer to, the modules of TAT-C that those import (anywhere in their code,
including imports within functions), and so on, with the resource files
these modules name (such as the instrument catalog). A name imported from a
package depends on the module that defines it, not on every module of the
package, so that `from tatc.analysis import collect_observations` depends
on the module of point sampling (and those it imports), not on the other
analyses. Results computed by a study are then reused after changes to TAT-C
that cannot change them (such as to modules the study does not use), and
not after changes that can.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import ast
import dataclasses
import dis
import functools
import hashlib
import json
import sys
import sysconfig
import types
from collections.abc import Callable, Iterable, Mapping
from datetime import date, datetime, time, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import shapely
from pydantic import BaseModel

from .alerts import _warn
from .hashing import canonical, digest
from .nodes import Node

_MAX_DEPTH = 24
"""The deepest nesting described, beyond which values are described by type."""

_MAX_ITEMS = 1000
"""The most items of a container (or array) described one by one."""

_TATC = "tatc"
"""The name of the package (TAT-C) whose modules a definition describes by
source."""

_ROOT = Path(__file__).resolve().parents[1]
"""The directory of the package."""

_RESOURCES = "resources"
"""The directory of the package's resource files (data, not modules)."""


def _type(value: Any) -> dict:
    kind = type(value)
    return {"type": f"{kind.__module__}.{kind.__qualname__}"}


def _name(value: Any) -> str:
    return f"{getattr(value, '__module__', None)}.{value.__qualname__}"


def _is_tatc_module(module: Any) -> bool:
    """Whether a module name is of TAT-C (or is TAT-C's package)."""
    return isinstance(module, str) and (
        module == _TATC or module.startswith(_TATC + ".")
    )


def _is_dataclass(value: Any) -> bool:
    """Whether a value is an instance of a dataclass (such as a cost model)."""
    return dataclasses.is_dataclass(value) and not isinstance(value, type)


def _has_fingerprint(value: Any) -> bool:
    """
    Whether a value is an object that describes itself (by a method
    `__fingerprint__`, see `_Describer.fingerprinted`), not a class.
    """
    return callable(getattr(type(value), "__fingerprint__", None))


@functools.lru_cache(maxsize=1)
def _installed_directories() -> tuple[Path, ...]:
    """Gets the directories of installed modules (and the standard library)."""
    paths = sysconfig.get_paths()
    names = ("stdlib", "platstdlib", "purelib", "platlib")
    found = {Path(paths[name]).resolve() for name in names if name in paths}
    found.add(Path(sys.base_prefix).resolve())
    return tuple(sorted(found))


def _is_own_module(module: Any) -> bool:
    """
    Whether a module is of the study's own: neither of TAT-C nor installed
    (such as `__main__`, or a module next to a notebook), so that its
    functions are described by their code (see the module).
    """
    if not isinstance(module, str) or _is_tatc_module(module):
        return False
    if module == "__main__":
        return True
    file = getattr(sys.modules.get(module), "__file__", None)
    if file is None:
        return False  # built in, or not loaded
    path = Path(file).resolve()
    return not any(path.is_relative_to(d) for d in _installed_directories())


def _is_own_class(value: Any) -> bool:
    """Whether a value is a class of the study's own (see `_is_own_module`)."""
    return isinstance(value, type) and _is_own_module(value.__module__)


def _is_own_object(value: Any) -> bool:
    """Whether a value is an object of a class of the study's own."""
    return _is_own_module(type(value).__module__)


def _wrapped(value: Any) -> types.FunctionType | None:
    """
    Gets the function that a decorator wraps (such as `functools.lru_cache`,
    or another of `functools.wraps`), if a value is such a wrapper.
    """
    if isinstance(value, (types.FunctionType, type)) or not callable(value):
        return None
    function = getattr(value, "__wrapped__", None)
    return function if isinstance(function, types.FunctionType) else None


def _hash_column(column: pd.Series) -> bytes:
    """
    Gets the hashes of the values of a column (geometries by their
    well-known binary, which is exact).
    """
    if str(column.dtype) == "geometry":
        return b"".join(
            hashlib.sha1(b"" if wkb is None else wkb).digest()
            for wkb in column.to_wkb()
        )
    try:
        return pd.util.hash_pandas_object(column, index=False).to_numpy().tobytes()
    except TypeError:  # such as of unhashable objects
        return "\x1f".join(repr(item) for item in column.tolist()).encode()


def _frame(value: pd.DataFrame | pd.Series) -> dict:
    """
    Describes a data frame (or series) by its columns, types, and shape, and
    a hash of its index and values.
    """
    frame = value.to_frame() if isinstance(value, pd.Series) else value
    hasher = hashlib.sha1()
    hasher.update(pd.util.hash_pandas_object(frame.index).to_numpy().tobytes())
    for i in range(frame.shape[1]):
        hasher.update(_hash_column(frame.iloc[:, i]))
    return {
        **_type(value),
        "columns": [str(column) for column in frame.columns],
        "types": [str(dtype) for dtype in frame.dtypes],
        "index": str(frame.index.dtype),
        "shape": list(value.shape),
        "data": hasher.hexdigest(),
    }


def _has_contents(cell: Any) -> bool:
    try:
        cell.cell_contents  # pylint: disable=pointless-statement
    except ValueError:  # an empty cell
        return False
    return True


_GLOBAL_LOADS = ("LOAD_GLOBAL", "LOAD_NAME")
"""The instructions that load a global (rather than an attribute, whose names
code also lists)."""

_ATTRIBUTE_LOADS = ("LOAD_ATTR", "LOAD_METHOD")
"""The instructions that load an attribute of the value loaded before."""


def _global_names(code: types.CodeType) -> dict[str, set[tuple[str, ...]]]:
    """
    Gets the global names a code (and its nested code) loads (not the names of
    attributes, such as `values` of `design.values`, which may be the names
    of unrelated globals), each with the chains of attributes that it loads
    of the global (such as `("score",)` of `helpers.score`).
    """
    found: dict[str, set[tuple[str, ...]]] = {}
    instructions = list(dis.get_instructions(code))
    for i, instruction in enumerate(instructions):
        if instruction.opname not in _GLOBAL_LOADS:
            continue
        chain = []
        for following in instructions[i + 1 :]:
            if following.opname not in _ATTRIBUTE_LOADS:
                break
            chain.append(following.argval)
        chains = found.setdefault(instruction.argval, set())
        if chain:
            chains.add(tuple(chain))
    for constant in code.co_consts:
        if isinstance(constant, types.CodeType):
            for name, chains in _global_names(constant).items():
                found.setdefault(name, set()).update(chains)
    return found


def _attribute(module: types.ModuleType, chain: tuple[str, ...]) -> tuple[str, Any]:
    """
    Gets the value that a chain of attributes of a module reaches, through
    its submodules (stopping at the first value that is not a module), and
    the path to it, or `("", None)` if there is none.
    """
    value: Any = module
    path = []
    for name in chain:
        if not isinstance(value, types.ModuleType):
            break
        if not hasattr(value, name):
            return "", None
        value = getattr(value, name)
        path.append(name)
    return ".".join(path), value


_PLAIN = (type(None), bool, int, float, str, bytes, date, time, timedelta, np.generic)


def _is_constant(value: Any, depth: int = 0) -> bool:
    """
    Whether a value is a constant: a number, string, date or time, or a
    (small) tuple or frozen set of them. Mutable containers (lists,
    dictionaries, and sets) are not constants: as attributes of classes,
    they are often state (such as registries) rather than definitions.
    """
    if isinstance(value, _PLAIN):
        return True
    if depth > 3:
        return False
    if isinstance(value, (tuple, frozenset)):
        return len(value) <= _MAX_ITEMS and all(
            _is_constant(v, depth + 1) for v in value
        )
    return False


class _Describer:
    """
    Describes values as plain data (see the module, and `_RULES`), noting
    the modules of the package that they depend on.
    """

    # pylint: disable=too-many-public-methods

    def __init__(self):
        self.modules: set[str] = set()
        """The modules of the package that the described values refer to."""

    def note(self, module: Any) -> None:
        """Notes a module that a value refers to (if it is of the package)."""
        if _is_tatc_module(module):
            self.modules.add(module)

    def describe(self, value: Any, depth: int = 0, seen: set | None = None) -> Any:
        """
        Describes a value as plain data, by the first rule that it matches
        (see `_RULES`), or by its type if none does.
        """
        seen = set() if seen is None else seen
        if depth <= _MAX_DEPTH:
            for test, rule in _RULES:
                if isinstance(value, test) if isinstance(test, tuple) else test(value):
                    return rule(self, value, depth, seen)
        self.note(type(value).__module__)
        return _type(value)

    def noted(self, module: Any, description: Any) -> Any:
        """Notes a module (see `note`), and gives a description."""
        self.note(module)
        return description

    # the rules of descriptions that are more than one expression (see `_RULES`)

    def dataclass(self, value: Any, depth: int, seen: set) -> dict:
        """
        Describes a dataclass by its fields (but those whose metadata has
        `"fingerprint": False`), and its class, if of the study's own.
        """
        self.note(type(value).__module__)
        described = {
            **_type(value),
            **{
                item.name: self.describe(getattr(value, item.name), depth + 1, seen)
                for item in dataclasses.fields(value)
                if item.metadata.get("fingerprint", True)
            },
        }
        if _is_own_module(type(value).__module__):
            described["class"] = self.klass(type(value), depth + 1, seen)
        return described

    def fingerprinted(self, value: Any, depth: int, seen: set) -> dict:
        """
        Describes an object by the value that its method `__fingerprint__()`
        returns (such as a design by its design space and decisions), which
        identifies what of it changes results.
        """
        self.note(type(value).__module__)
        if id(value) in seen:
            return {**_type(value), "recursive": True}
        seen = seen | {id(value)}
        return {
            **_type(value),
            "fingerprint": self.describe(value.__fingerprint__(), depth + 1, seen),
        }

    def sequence(self, value: list | tuple, depth: int, seen: set) -> Any:
        """Describes a list or tuple by its items (or its length, if long)."""
        if len(value) > _MAX_ITEMS:
            return {**_type(value), "length": len(value)}
        return [self.describe(v, depth + 1, seen) for v in value]

    def collection(self, value: set | frozenset, depth: int, seen: set) -> Any:
        """Describes a set by its sorted items (or its length, if long)."""
        if len(value) > _MAX_ITEMS:
            return {**_type(value), "length": len(value)}
        return sorted((self.describe(v, depth + 1, seen) for v in value), key=canonical)

    def mapping(self, value: dict, depth: int, seen: set) -> Any:
        """
        Describes a dictionary by its items (sorted pairs, if its keys are
        not strings), or its length, if long.
        """
        if len(value) > _MAX_ITEMS:
            return {**_type(value), "length": len(value)}
        if all(isinstance(key, str) for key in value):
            return {key: self.describe(v, depth + 1, seen) for key, v in value.items()}
        pairs = (
            [self.describe(k, depth + 1, seen), self.describe(v, depth + 1, seen)]
            for k, v in value.items()
        )
        return sorted(pairs, key=canonical)

    def node(self, value: Node, depth: int, seen: set) -> dict:
        """Describes a node of a design space by its attributes."""
        self.note(type(value).__module__)
        if id(value) in seen:
            return {**_type(value), "recursive": True}
        seen = seen | {id(value)}
        return {
            **_type(value),
            **{
                name: self.describe(attribute, depth + 1, seen)
                for name, attribute in sorted(vars(value).items())
            },
        }

    def code(
        self, code: types.CodeType, depth: int, seen: set, doc: str | None = None
    ) -> dict:
        """
        Describes compiled code by its instructions, constants, and names (not
        its file or line numbers, which change as code moves, nor the
        docstring `doc` of its function, which does not change results). The
        constants include the small integers that some versions of Python
        load by an instruction of their own (`LOAD_SMALL_INT`, from 3.14)
        rather than from the constants of the code, so that the constants of
        code are its literals in every version.
        """
        constants = list(code.co_consts)
        if doc is not None and constants and constants[0] == doc:
            constants = constants[1:]
        constants.extend(
            instruction.argval
            for instruction in dis.get_instructions(code)
            if instruction.opname == "LOAD_SMALL_INT"
        )
        return {
            "code": code.co_code.hex(),
            "constants": [self.describe(c, depth + 1, seen) for c in constants],
            "names": list(code.co_names),
        }

    def function(self, function: types.FunctionType, depth: int, seen: set) -> dict:
        """
        Describes a function: its code, the values it closes over, its
        defaults, and the globals it refers to (see `reference`).
        """
        self.note(function.__module__)
        if id(function) in seen:
            return {"function": function.__qualname__, "recursive": True}
        seen = seen | {id(function)}
        closure = {
            name: self.describe(cell.cell_contents, depth + 1, seen)
            for name, cell in zip(
                function.__code__.co_freevars, function.__closure__ or ()
            )
            if _has_contents(cell)
        }
        referenced = {
            name: self.reference(
                function.__globals__[name], function.__module__, chains, depth, seen
            )
            for name, chains in _global_names(function.__code__).items()
            if name in function.__globals__
        }
        return {
            "function": function.__qualname__,
            **self.code(function.__code__, depth, seen, function.__doc__),
            "closure": closure,
            "defaults": self.describe(function.__defaults__, depth + 1, seen),
            "keyword defaults": self.describe(function.__kwdefaults__, depth + 1, seen),
            "globals": referenced,
        }

    def reference(
        self,
        value: Any,
        module: str,
        chains: Iterable[tuple[str, ...]],
        depth: int,
        seen: set,
    ) -> Any:
        """
        Describes a global that a function of a module refers to: a function
        of another module (or one that a decorator wraps) by its name, unless
        it is of a module of the study's own, a module of the study's own by
        the attributes of it that the function refers to (by their `chains`,
        see `_global_names`), and other values as any other (see `describe`).
        """
        # pylint: disable=too-many-arguments,too-many-positional-arguments
        function = _wrapped(value) or value
        if (
            isinstance(function, types.FunctionType)
            and function.__module__ != module
            and not _is_own_module(function.__module__)
        ):
            self.note(function.__module__)
            return {"name": _name(function)}
        if isinstance(value, types.ModuleType) and _is_own_module(value.__name__):
            attributes = {}
            for chain in sorted(chains):
                path, reached = _attribute(value, chain)
                if path:
                    attributes[path] = self.reference(
                        reached, value.__name__, (), depth + 1, seen
                    )
            return {"module": value.__name__, "attributes": attributes}
        return self.describe(value, depth + 1, seen)

    def method(self, value: types.MethodType, depth: int, seen: set) -> dict:
        """Describes a bound method by its function and its object."""
        return {
            "method": self.describe(value.__func__, depth + 1, seen),
            "self": self.describe(value.__self__, depth + 1, seen),
        }

    def partial(self, value: functools.partial, depth: int, seen: set) -> dict:
        """Describes a partial function by its function and arguments."""
        return {
            "partial": self.describe(value.func, depth + 1, seen),
            "arguments": self.describe(value.args, depth + 1, seen),
            "keywords": self.describe(value.keywords, depth + 1, seen),
        }

    def wrapper(self, value: Any, depth: int, seen: set) -> dict:
        """Describes a decorator's wrapper by its type and wrapped function."""
        self.note(type(value).__module__)
        return {
            **_type(value),
            "wrapped": self.describe(_wrapped(value), depth + 1, seen),
        }

    def klass(self, cls: type, depth: int, seen: set) -> dict:
        """
        Describes a class of the study's own (see `_is_own_module`): its name,
        bases (other than `object`), and members: its functions (methods,
        static and class methods, and the functions of properties) by their
        code, and its constants.
        """
        if id(cls) in seen:
            return {"class": cls.__qualname__, "recursive": True}
        seen = seen | {id(cls)}
        members: dict[str, Any] = {}
        for name, member in sorted(vars(cls).items()):
            if isinstance(member, (staticmethod, classmethod)):
                member = member.__func__
            if isinstance(member, property):
                members[name] = {
                    part: self.function(accessor, depth + 1, seen)
                    for part, accessor in (("get", member.fget), ("set", member.fset))
                    if isinstance(accessor, types.FunctionType)
                }
            elif isinstance(member, types.FunctionType):
                members[name] = self.function(member, depth + 1, seen)
            elif not name.startswith("__") and _is_constant(member):
                members[name] = self.describe(member, depth + 1, seen)
        return {
            "class": cls.__qualname__,
            "bases": [
                self.describe(base, depth + 1, seen)
                for base in cls.__bases__
                if base is not object
            ],
            "members": members,
        }

    def instance(self, value: Any, depth: int, seen: set) -> dict:
        """
        Describes an object of a class of the study's own: its class (see
        `klass`) and its attributes, but those that are mutable containers
        (lists, dictionaries, and sets), which are often state (such as
        counters and records of calls).
        """
        if id(value) in seen:
            return {**_type(value), "recursive": True}
        seen = seen | {id(value)}
        if hasattr(value, "__dict__"):
            attributes = dict(vars(value))
        else:
            slots = getattr(type(value), "__slots__", ())
            slots = [slots] if isinstance(slots, str) else list(slots)
            attributes = {n: getattr(value, n) for n in slots if hasattr(value, n)}
        attributes = {
            name: attribute
            for name, attribute in attributes.items()
            if not isinstance(attribute, (list, dict, set, bytearray))
        }
        return {
            **_type(value),
            "class": self.klass(type(value), depth + 1, seen),
            "attributes": self.describe(attributes, depth + 1, seen),
        }


def _array(value: np.ndarray) -> dict:
    """Describes an array by its type, shape, and a hash of its data."""
    if value.dtype == object:
        return {**_type(value), "shape": list(value.shape)}
    return {
        "array": str(value.dtype),
        "shape": list(value.shape),
        "data": hashlib.sha1(np.ascontiguousarray(value).tobytes()).hexdigest(),
    }


_RULES: tuple[tuple[tuple[type, ...] | Callable[[Any], bool], Callable], ...] = (
    # by content: plain values as themselves (floats by their representation,
    # exactly, and numpy scalars as Python's), bytes by their hash, dates and
    # times by their ISO format, and durations by their seconds
    ((type(None), bool, int, str), lambda _, value, *__: value),
    ((float,), lambda _, value, *__: repr(value)),
    (
        (np.generic,),
        lambda describer, value, *rest: describer.describe(value.item(), *rest),
    ),
    ((bytes,), lambda _, value, *__: hashlib.sha1(value).hexdigest()),
    ((datetime, date, time), lambda _, value, *__: value.isoformat()),
    ((timedelta,), lambda _, value, *__: {"timedelta": value.total_seconds()}),
    (_has_fingerprint, _Describer.fingerprinted),
    (_is_dataclass, _Describer.dataclass),
    # pydantic models by their data, arrays and data frames by hashes of their
    # data (see `_array` and `_frame`), and geometries by a hash of their
    # well-known binary
    (
        (BaseModel,),
        lambda describer, value, *_: describer.noted(
            type(value).__module__, {**_type(value), **value.model_dump(mode="json")}
        ),
    ),
    ((np.ndarray,), lambda _, value, *__: _array(value)),
    ((pd.DataFrame, pd.Series), lambda _, value, *__: _frame(value)),
    (
        (shapely.Geometry,),
        lambda _, value, *__: {
            **_type(value),
            "data": hashlib.sha1(value.wkb).hexdigest(),
        },
    ),
    ((list, tuple), _Describer.sequence),
    ((set, frozenset), _Describer.collection),
    ((dict,), _Describer.mapping),
    ((Node,), _Describer.node),
    ((functools.partial,), _Describer.partial),
    ((types.MethodType,), _Describer.method),
    (lambda value: _wrapped(value) is not None, _Describer.wrapper),
    # by code
    ((types.CodeType,), _Describer.code),
    ((types.FunctionType,), _Describer.function),
    (_is_own_class, _Describer.klass),
    (_is_own_object, _Describer.instance),
    # by name: classes and built-in functions, and modules
    (
        (types.BuiltinFunctionType, type),
        lambda describer, value, *_: describer.noted(
            getattr(value, "__module__", None), {"name": _name(value)}
        ),
    ),
    (
        (types.ModuleType,),
        lambda describer, value, *_: describer.noted(
            value.__name__, {"module": value.__name__}
        ),
    ),
)
"""What values are described by (see the module): the first rule whose test
(types, or a predicate) a value matches describes it, and a value that
matches none is described by its type (as is one nested too deep)."""


@functools.lru_cache(maxsize=1)
def _module_files() -> dict[str, Path]:
    """
    Gets the source files of the modules of the package, by module name
    (packages by their `__init__.py`).
    """
    found = {}
    for path in sorted(_ROOT.rglob("*.py")):
        if "__pycache__" in path.parts:
            continue
        parts = path.relative_to(_ROOT).with_suffix("").parts
        if parts[-1] == "__init__":
            parts = parts[:-1]
        found[".".join((_TATC, *parts))] = path
    return found


def _is_init_module(module: str) -> bool:
    """Whether a module of TAT-C is a package (an `__init__.py`)."""
    path = _module_files().get(module)
    return path is not None and path.name == "__init__.py"


@functools.lru_cache(maxsize=None)
def _tree(module: str) -> ast.Module:
    """Gets the syntax tree of a module of the package."""
    return ast.parse(_module_files()[module].read_bytes())


def _absolute(module: str, node: ast.ImportFrom) -> str:
    """Gets the absolute name of the module of an import from a module."""
    if node.level == 0:
        return node.module or ""
    package = module if _is_init_module(module) else module.rpartition(".")[0]
    parts = package.split(".")
    parts = parts[: len(parts) - (node.level - 1)]
    return ".".join(parts + ([node.module] if node.module else []))


def _is_constant_node(node: ast.AST) -> bool:
    """Whether an expression is a constant (or a tuple or list of them)."""
    if isinstance(node, ast.Constant):
        return True
    return isinstance(node, (ast.Tuple, ast.List)) and all(
        _is_constant_node(item) for item in node.elts
    )


def _assigned(node: ast.stmt) -> tuple[list[str], ast.AST | None]:
    """Gets the names a statement assigns, and the value it assigns them."""
    if isinstance(node, ast.Assign):
        return [t.id for t in node.targets if isinstance(t, ast.Name)], node.value
    if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
        return [node.target.id], node.value
    return [], None


@functools.lru_cache(maxsize=None)
def _resolve_import(module: str, name: str) -> frozenset[tuple[str, str]]:
    """
    Gets what a name imported from a module of the package depends on: a
    submodule of that name, or (from a package) the module that defines it,
    following the package's imports. Items are `("module", m)`, for a module
    and those it imports, or `("file", m)`, for the source of a module only
    (a package through which a name is imported, or that defines it as a
    constant).
    """
    files = _module_files()
    if f"{module}.{name}" in files:
        return frozenset({("module", f"{module}.{name}")})
    if module not in files:
        return frozenset()
    if not _is_init_module(module):
        return frozenset({("module", module)})
    for node in _tree(module).body:
        if isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if (alias.asname or alias.name) == name:
                    source = _absolute(module, node)
                    if not _is_tatc_module(source) or source == module:
                        return frozenset({("file", module)})
                    found = _resolve_import(source, alias.name)
                    return frozenset({("file", module)}) | found
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.asname == name and _is_tatc_module(alias.name):
                    return frozenset({("file", module), ("module", alias.name)})
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if node.name == name:
                return frozenset({("module", module)})
        else:
            names, value = _assigned(node)
            if name in names:
                if value is not None and _is_constant_node(value):
                    return frozenset({("file", module)})
                return frozenset({("module", module)})
    return frozenset({("module", module)})


@functools.lru_cache(maxsize=None)
def _imports(module: str) -> frozenset[tuple[str, str]]:
    """
    Gets what the imports of a module of the package (anywhere in its code)
    depend on (see `_resolve_import`).
    """
    found: set[tuple[str, str]] = set()
    for node in ast.walk(_tree(module)):
        if isinstance(node, ast.Import):
            found.update(
                ("module", alias.name)
                for alias in node.names
                if _is_tatc_module(alias.name)
            )
        elif isinstance(node, ast.ImportFrom):
            source = _absolute(module, node)
            if not _is_tatc_module(source):
                continue
            for alias in node.names:
                if alias.name == "*":
                    found.add(("module", source))
                else:
                    found.update(_resolve_import(source, alias.name))
    return frozenset(found)


@functools.lru_cache(maxsize=None)
def _resources_of(module: str) -> frozenset[str]:
    """
    Gets the resource files a module of the package names: if it refers to
    the directory of resources, those whose names, or the names of whose
    directories (within it), are strings of its code.
    """
    strings = {
        node.value
        for node in ast.walk(_tree(module))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }
    if _RESOURCES not in strings:
        return frozenset()
    directory = _ROOT / _RESOURCES
    found = set()
    for path in directory.rglob("*"):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        folders = path.relative_to(directory).parts[:-1]
        if path.name in strings or any(folder in strings for folder in folders):
            found.add(path.relative_to(_ROOT).as_posix())
    return frozenset(found)


def _get_dependencies(modules: Iterable[str]) -> tuple[list[str], list[str]]:
    """
    Gets the modules of TAT-C (and resource files) that modules depend on:
    the modules and the modules of TAT-C they import (anywhere in their
    code), and so on, where a name imported from a package depends on the
    module that defines it (and the source of the package's `__init__.py`),
    not on every module of the package, and the resource files that these
    modules name.

    Args:
        modules (Iterable[str]): The names of the modules (others than of
            TAT-C are ignored).

    Returns:
        tuple[list[str], list[str]]: the names of the modules, and the paths
            of the resource files (within the package)
    """
    files = _module_files()
    whole: set[str] = set()
    sources: set[str] = set()
    pending = [module for module in modules if _is_tatc_module(module)]
    while pending:
        module = pending.pop()
        if module in whole or module not in files:
            continue
        whole.add(module)
        for kind, other in _imports(module):
            if kind == "module":
                pending.append(other)
            else:
                sources.add(other)
    resources = set().union(*(_resources_of(module) for module in whole))
    # the trees are no longer needed: what they give is kept
    _tree.cache_clear()
    return sorted(whole | sources), sorted(resources)


@functools.lru_cache(maxsize=None)
def _hash_file(path: Path, modified: int, size: int) -> str:
    del modified, size  # keys of the cache, so that a changed file is hashed again
    return digest(path.read_bytes(), 16)


def _hash(path: Path) -> str:
    """Gets a hash of the content of a file (computed again if it changes)."""
    status = path.stat()
    return _hash_file(path, status.st_mtime_ns, status.st_size)


def describe(value: Any) -> Any:
    """
    Describes a value as plain data, as a definition describes the values it
    is given (see the module and `get_definition`), such as to identify the
    overrides of a scenario by their content.

    Args:
        value (Any): The value.

    Returns:
        Any: the description, as plain data
    """
    return _Describer().describe(value)


def get_definition(*values: Any, **named: Any) -> dict:
    """
    Gets the definition of values (such as a design space and the functions
    of a study) that their fingerprint hashes (see `get_fingerprint`): a
    description of each value (see the module), the version of Python (as
    functions are described by their compiled code), and a hash of the
    source of each module of TAT-C (and of each resource file) that the
    values depend on (see the module).

    Args:
        *values (Any): The values, described by their positions.
        **named (Any): The values, described by their names.

    Returns:
        dict: the definition, as plain data (`python`, `values`, `modules`,
            and `resources`), which `compare_definitions` compares
    """
    describer = _Describer()
    described = {str(i): describer.describe(value) for i, value in enumerate(values)}
    described.update({name: describer.describe(value) for name, value in named.items()})
    modules, resources = _get_dependencies(describer.modules)
    files = _module_files()
    return {
        "python": f"{sys.implementation.name} {sys.version_info[0]}.{sys.version_info[1]}",
        "values": described,
        "modules": {module: _hash(files[module]) for module in modules},
        "resources": {path: _hash(_ROOT / path) for path in resources},
    }


def _fingerprint_of(definition: dict) -> str:
    """Gets the fingerprint of a definition (see `get_definition`)."""
    return digest(canonical(definition), 16)


def get_fingerprint(*values: Any, **named: Any) -> str:
    """
    Gets a fingerprint of values (such as a design space and the functions
    of a study): a hash of their definition (see `get_definition`), which
    changes if they do, or if the source of TAT-C that they depend on does.

    Args:
        *values (Any): The values, described by their positions.
        **named (Any): The values, described by their names.

    Returns:
        str: the fingerprint
    """
    return _fingerprint_of(get_definition(*values, **named))


_SHORT = 40
"""The longest value shown in a difference between definitions."""


def _show(value: Any) -> str:
    return value if isinstance(value, str) else json.dumps(value, default=str)


def _differences(before: Any, after: Any, path: str) -> list[str]:
    """
    Gets the differences between descriptions, by the path to each value
    that differs (and the values, if they are short).
    """
    if before == after:
        return []
    if isinstance(before, dict) and isinstance(after, dict):
        found = []
        for key in sorted(set(before) | set(after), key=str):
            inner = f"{path}.{key}" if path else str(key)
            if key not in before:
                found.append(f"{inner} (added)")
            elif key not in after:
                found.append(f"{inner} (removed)")
            else:
                found.extend(_differences(before[key], after[key], inner))
        return found
    if (
        isinstance(before, list)
        and isinstance(after, list)
        and len(before) == len(after)
    ):
        found = []
        for i, (old, new) in enumerate(zip(before, after)):
            found.extend(_differences(old, new, f"{path}[{i}]"))
        return found
    path = path or "values"
    if all(
        not isinstance(v, (dict, list)) and len(_show(v)) <= _SHORT
        for v in (before, after)
    ):
        return [f"{path}: {_show(before)} -> {_show(after)}"]
    return [path]


def _is_code(description: Any) -> bool:
    """Whether a description is of compiled code (see `_Describer.code`)."""
    return isinstance(description, dict) and {"code", "constants", "names"} <= set(
        description
    )


_LOCALS = ".<locals>."
"""The part of the qualified names of nested code, which some versions of
Python (3.10) keep as constants of the code that makes their functions."""


def _literals(constants: list) -> list:
    """
    Gets the literals of the described constants of code, as the same code
    gives them in every version of Python: the items of tuples and frozen
    sets, and the constants of nested code (whose code is inlined by some
    versions), each once, sorted, without `None` (which some versions keep
    as a placeholder of a docstring) and the qualified names of nested code.
    """
    found: dict[str, Any] = {}

    def visit(constant: Any) -> None:
        if constant is None or (isinstance(constant, str) and _LOCALS in constant):
            return
        if isinstance(constant, list):
            for item in constant:
                visit(item)
        elif _is_code(constant):
            for item in constant["constants"]:
                visit(item)
        else:
            found.setdefault(canonical(constant), constant)

    for constant in constants:
        visit(constant)
    return [found[key] for key in sorted(found)]


def _without_code(description: Any) -> Any:
    """
    Gets a description of values without the code of their functions (the
    instructions and names of compiled code, see `_Describer.code`), which is
    the same for the same values with other code (such as a function after
    its expressions are edited, or compiled by another version of Python).
    The constants of the code (such as the numbers of a constraint) are
    values, as its literals (see `_literals`, which every version of Python
    gives alike): functions of other literals are of other studies, whose
    results are not those of the same study with other code.
    """
    if _is_code(description):
        return {
            key: (_literals(value) if key == "constants" else _without_code(value))
            for key, value in description.items()
            if key not in ("code", "names")
        }
    if isinstance(description, dict):
        return {key: _without_code(value) for key, value in description.items()}
    if isinstance(description, list):
        return [_without_code(value) for value in description]
    return description


def _find_stale(
    definition: dict,
    others: Mapping[Any, dict | None],
    count: Callable[[Any], int],
) -> dict[Any, int]:
    """
    Gets the results of the same study (or function) as a definition with
    other code (of its functions, of TAT-C, or of Python), which are not
    reused: the keys of the other definitions (such as the contexts of a
    store, or the folders of a cache) that have the same values but other
    code (see `_without_code`), in their order, with their numbers of
    results (by `count`, a function of a key, 0 for a definition kept
    without results). The results of other studies (such as variants of
    other settings) are not stale, as studies share stores and caches
    freely, and definitions that are `None` (unknown) are of no study.
    """
    values = _without_code(definition.get("values", {}))
    return {
        key: count(key)
        for key, other in others.items()
        if other is not None and _without_code(other.get("values", {})) == values
    }


def warn_stale(
    found: Mapping[Any, int],
    definition: dict,
    definition_of: Callable[[Any], dict | None],
    subject: Callable[[int], str],
    remedy: str,
    once: Callable[[], bool] | None = None,
) -> bool:
    """
    Warns (at the user's call) that the results of the same study (or
    function) as a definition with other code (see `_find_stale`) are not
    reused, if any are found, with the differences of the latest of their
    definitions from the definition (see `compare_definitions`).

    Args:
        found (Mapping[Any, int]): The numbers of stale results, by the key
            of their definition, the latest first (see `_find_stale`).
        definition (dict): The definition (see `get_definition`).
        definition_of (Callable[[Any], dict | None]): A function of a key
            that gives its definition.
        subject (Callable[[int], str]): A function of the number of stale
            results that says which they are (such as `"the store has 12
            evaluations of namespace 'a' by this study"`).
        remedy (str): How to remove them.
        once (Callable[[], bool] | None): A function called before warning
            that gives whether to warn (such as only for the first of the
            parallel workers of a study), or `None` to warn.

    Returns:
        bool: whether it warned
    """
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    found = {key: number for key, number in found.items() if number}
    if not found or (once is not None and not once()):
        return False
    latest = definition_of(next(iter(found))) or {}
    _warn(
        f"{subject(sum(found.values()))} with other code, which are not reused "
        + f"({remedy}); the latest differs in: "
        + "; ".join(compare_definitions(latest, definition))
    )
    return True


def compare_definitions(before: dict, after: dict, limit: int = 5) -> list[str]:
    """
    Compares definitions (see `get_definition`), as to explain why results
    of one are not reused by the other: the version of Python, the values
    that differ (by their path in the description, such as
    `evaluate.closure.offset: 0.0 -> 1.0`, or `evaluate.code` for a change
    to the code of a function), and the modules of TAT-C (and resource
    files) whose source differs, or that only one depends on.

    Args:
        before (dict): A definition.
        after (dict): Another definition.
        limit (int): The most differences listed (and a count of others).

    Returns:
        list[str]: the differences (none if the definitions are the same)
    """
    found = []
    if before.get("python") != after.get("python"):
        found.append(f"python: {before.get('python')} -> {after.get('python')}")
    found.extend(_differences(before.get("values", {}), after.get("values", {}), ""))
    for kind, label in (("modules", "module"), ("resources", "resource")):
        old, new = before.get(kind, {}), after.get(kind, {})
        for name in sorted(set(old) | set(new)):
            if name not in old:
                found.append(f"{label} {name} (added)")
            elif name not in new:
                found.append(f"{label} {name} (removed)")
            elif old[name] != new[name]:
                found.append(f"{label} {name}")
    if len(found) > limit:
        return found[:limit] + [f"{len(found) - limit} more"]
    return found
