"""
Fingerprints of the definitions of trade studies (design spaces, functions,
and settings), which key persistent results so that results computed with
other definitions are not reused.

A definition (see `get_definition`) describes values by their content where
it can: plain values (numbers, strings, dates and times, and their
containers), pydantic models (such as satellites and instruments), arrays
and data frames (by a hash of their data, with geometries by their
well-known binary), the nodes of design spaces, dataclasses (such as cost
models) by their fields (but those whose metadata has `"fingerprint":
False`, such as directories of caches, which do not change results), and
functions, by their compiled code (which ties a definition to the version of
Python, which it records), the values they close over and their defaults,
and the constants (numbers, strings, dates and times, and tuples of them),
data (arrays and data frames), dataclasses, and functions of the same module
that they refer to as globals (not mutable containers, which are often
state, such as caches). Functions and classes of other modules are described
by their names, and other objects by their type only. Data that functions
read (such as from files) and the code of other modules they call (other
than TAT-C's) are not described: identify them by a namespace (see
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
import types
from collections.abc import Iterable
from datetime import date, datetime, time, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from pydantic import BaseModel

_MAX_DEPTH = 24
"""The deepest nesting described, beyond which values are described by type."""

_MAX_ITEMS = 1000
"""The most items of a container (or array) described one by one."""

_PACKAGE = "tatc"
"""The name of the package whose modules a definition describes by source."""

_ROOT = Path(__file__).resolve().parents[1]
"""The directory of the package."""

_RESOURCES = "resources"
"""The directory of the package's resource files (data, not modules)."""


def _type(value: Any) -> dict:
    kind = type(value)
    return {"type": f"{kind.__module__}.{kind.__qualname__}"}


def _name(value: Any) -> str:
    return f"{getattr(value, '__module__', None)}.{value.__qualname__}"


def _sort_key(value: Any) -> str:
    return json.dumps(value, sort_keys=True, default=str)


def _is_package_module(module: Any) -> bool:
    """Whether a module name is of the package (TAT-C)."""
    return isinstance(module, str) and (
        module == _PACKAGE or module.startswith(_PACKAGE + ".")
    )


def _is_dataclass(value: Any) -> bool:
    """Whether a value is an instance of a dataclass (such as a cost model)."""
    return dataclasses.is_dataclass(value) and not isinstance(value, type)


def _is_data(value: Any) -> bool:
    """Whether a value is data (an array or a data frame or series)."""
    return isinstance(value, (np.ndarray, pd.DataFrame, pd.Series))


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
    digest = hashlib.sha1()
    digest.update(pd.util.hash_pandas_object(frame.index).to_numpy().tobytes())
    for i in range(frame.shape[1]):
        digest.update(_hash_column(frame.iloc[:, i]))
    return {
        **_type(value),
        "columns": [str(column) for column in frame.columns],
        "types": [str(dtype) for dtype in frame.dtypes],
        "index": str(frame.index.dtype),
        "shape": list(value.shape),
        "data": digest.hexdigest(),
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


def _global_names(code: types.CodeType) -> list[str]:
    """
    Gets the global names a code (and its nested code) loads (not the names of
    attributes, such as `values` of `design.values`, which may be the names
    of unrelated globals).
    """
    names = [
        instruction.argval
        for instruction in dis.get_instructions(code)
        if instruction.opname in _GLOBAL_LOADS
    ]
    for constant in code.co_consts:
        if isinstance(constant, types.CodeType):
            names.extend(_global_names(constant))
    return list(dict.fromkeys(names))


_PLAIN = (type(None), bool, int, float, str, bytes, date, time, timedelta, np.generic)


def _is_constant(value: Any, depth: int = 0) -> bool:
    """
    Whether a value is a constant: a number, string, date or time, or a
    (small) tuple or frozen set of them. Mutable containers (lists,
    dictionaries, and sets) are not constants: as globals, they are often
    state (such as caches) rather than definitions, whose contents change
    as a study runs.
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
    Describes values as plain data (see the module), noting the modules of
    the package that they depend on.
    """

    def __init__(self):
        self.modules: set[str] = set()
        """The modules of the package that the described values refer to."""

    def note(self, module: Any) -> None:
        """Notes a module that a value refers to (if it is of the package)."""
        if _is_package_module(module):
            self.modules.add(module)

    def code(self, code: types.CodeType, depth: int, seen: set) -> dict:
        """
        Describes compiled code by its instructions, constants, and names (not
        its file or line numbers, which change as code moves).
        """
        return {
            "code": code.co_code.hex(),
            "constants": [self.describe(c, depth + 1, seen) for c in code.co_consts],
            "names": list(code.co_names),
        }

    def function(self, function: types.FunctionType, depth: int, seen: set) -> dict:
        """
        Describes a function: its code, the values it closes over, its
        defaults, the constants, dataclasses, and functions of its module it
        refers to as globals, and the names of the functions, classes, and
        modules of other modules it refers to as globals.
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
        referenced: dict[str, Any] = {}
        module = function.__module__
        for name in _global_names(function.__code__):
            if name not in function.__globals__:
                continue
            value = function.__globals__[name]
            if isinstance(value, types.FunctionType):
                if value.__module__ == module:
                    referenced[name] = self.function(value, depth + 1, seen)
                else:
                    self.note(value.__module__)
                    referenced[name] = {"name": _name(value)}
            elif (
                isinstance(value, (type, types.BuiltinFunctionType, types.ModuleType))
                or _is_constant(value)
                or _is_data(value)
                or _is_dataclass(value)
            ):
                referenced[name] = self.describe(value, depth + 1, seen)
            else:
                # state (such as a cache), whose code may be of the package
                self.note(type(value).__module__)
        return {
            "function": function.__qualname__,
            **self.code(function.__code__, depth, seen),
            "closure": closure,
            "defaults": self.describe(function.__defaults__, depth + 1, seen),
            "keyword defaults": self.describe(function.__kwdefaults__, depth + 1, seen),
            "globals": referenced,
        }

    def describe(self, value: Any, depth: int = 0, seen: set | None = None) -> Any:
        """
        Describes a value as plain data (see the module).
        """
        # pylint: disable=too-many-return-statements,too-many-branches,import-outside-toplevel
        seen = set() if seen is None else seen
        if depth > _MAX_DEPTH:
            self.note(type(value).__module__)
            return _type(value)
        if value is None or isinstance(value, (bool, int, str)):
            return value
        if isinstance(value, float):
            return repr(value)
        if isinstance(value, np.generic):
            return self.describe(value.item(), depth, seen)
        if isinstance(value, bytes):
            return hashlib.sha1(value).hexdigest()
        if isinstance(value, (datetime, date, time)):
            return value.isoformat()
        if isinstance(value, timedelta):
            return {"timedelta": value.total_seconds()}
        if _is_dataclass(value):
            self.note(type(value).__module__)
            return {
                **_type(value),
                **{
                    item.name: self.describe(getattr(value, item.name), depth + 1, seen)
                    for item in dataclasses.fields(value)
                    if item.metadata.get("fingerprint", True)
                },
            }
        if isinstance(value, BaseModel):
            self.note(type(value).__module__)
            return {**_type(value), **value.model_dump(mode="json")}
        if isinstance(value, np.ndarray):
            if value.dtype == object:
                return {**_type(value), "shape": list(value.shape)}
            return {
                "array": str(value.dtype),
                "shape": list(value.shape),
                "data": hashlib.sha1(np.ascontiguousarray(value).tobytes()).hexdigest(),
            }
        if isinstance(value, (pd.DataFrame, pd.Series)):
            return _frame(value)
        if isinstance(value, (list, tuple)):
            if len(value) > _MAX_ITEMS:
                return {**_type(value), "length": len(value)}
            return [self.describe(v, depth + 1, seen) for v in value]
        if isinstance(value, (set, frozenset)):
            if len(value) > _MAX_ITEMS:
                return {**_type(value), "length": len(value)}
            return sorted(
                (self.describe(v, depth + 1, seen) for v in value), key=_sort_key
            )
        if isinstance(value, dict):
            if len(value) > _MAX_ITEMS:
                return {**_type(value), "length": len(value)}
            if all(isinstance(key, str) for key in value):
                return {
                    key: self.describe(item, depth + 1, seen)
                    for key, item in value.items()
                }
            return sorted(
                (
                    [
                        self.describe(k, depth + 1, seen),
                        self.describe(v, depth + 1, seen),
                    ]
                    for k, v in value.items()
                ),
                key=_sort_key,
            )
        if isinstance(value, types.CodeType):
            return self.code(value, depth, seen)
        if isinstance(value, types.FunctionType):
            return self.function(value, depth, seen)
        if isinstance(value, types.MethodType):
            return {
                "method": self.describe(value.__func__, depth + 1, seen),
                "self": self.describe(value.__self__, depth + 1, seen),
            }
        if isinstance(value, functools.partial):
            return {
                "partial": self.describe(value.func, depth + 1, seen),
                "arguments": self.describe(value.args, depth + 1, seen),
                "keywords": self.describe(value.keywords, depth + 1, seen),
            }
        if isinstance(value, (types.BuiltinFunctionType, type)):
            self.note(getattr(value, "__module__", None))
            return {"name": _name(value)}
        if isinstance(value, types.ModuleType):
            self.note(value.__name__)
            return {"module": value.__name__}
        self.note(type(value).__module__)
        from .nodes import Node

        if isinstance(value, Node):
            if id(value) in seen:
                return {**_type(value), "recursive": True}
            return {
                **_type(value),
                **{
                    name: self.describe(attribute, depth + 1, seen | {id(value)})
                    for name, attribute in sorted(vars(value).items())
                },
            }
        return _type(value)


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
        found[".".join((_PACKAGE, *parts))] = path
    return found


def _is_package(module: str) -> bool:
    """Whether a module of the package is a package (an `__init__.py`)."""
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
    package = module if _is_package(module) else module.rpartition(".")[0]
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
def _resolve(module: str, name: str) -> frozenset[tuple[str, str]]:
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
    if not _is_package(module):
        return frozenset({("module", module)})
    for node in _tree(module).body:
        if isinstance(node, ast.ImportFrom):
            for alias in node.names:
                if (alias.asname or alias.name) == name:
                    source = _absolute(module, node)
                    if not _is_package_module(source) or source == module:
                        return frozenset({("file", module)})
                    return frozenset({("file", module)}) | _resolve(source, alias.name)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.asname == name and _is_package_module(alias.name):
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
    depend on (see `_resolve`).
    """
    found: set[tuple[str, str]] = set()
    for node in ast.walk(_tree(module)):
        if isinstance(node, ast.Import):
            found.update(
                ("module", alias.name)
                for alias in node.names
                if _is_package_module(alias.name)
            )
        elif isinstance(node, ast.ImportFrom):
            source = _absolute(module, node)
            if not _is_package_module(source):
                continue
            for alias in node.names:
                if alias.name == "*":
                    found.add(("module", source))
                else:
                    found.update(_resolve(source, alias.name))
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
    pending = [module for module in modules if _is_package_module(module)]
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
    return hashlib.sha1(path.read_bytes()).hexdigest()[:16]


def _hash(path: Path) -> str:
    """Gets a hash of the content of a file (computed again if it changes)."""
    status = path.stat()
    return _hash_file(path, status.st_mtime_ns, status.st_size)


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
    text = json.dumps(definition, sort_keys=True, default=str)
    return hashlib.sha1(text.encode()).hexdigest()[:16]


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


def _without_code(description: Any) -> Any:
    """
    Gets a description of values without the code of their functions (the
    instructions, constants, and names of compiled code, see
    `_Describer.code`), which is the same for the same values with other
    code (such as a function after it is edited, or compiled by another
    version of Python).
    """
    if isinstance(description, dict):
        code = {"code", "constants", "names"} <= set(description)
        return {
            key: _without_code(value)
            for key, value in description.items()
            if not (code and key in ("code", "constants", "names"))
        }
    if isinstance(description, list):
        return [_without_code(value) for value in description]
    return description


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
