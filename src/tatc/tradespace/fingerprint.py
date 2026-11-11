"""
Fingerprints of the definitions of trade studies (design spaces, functions,
and settings), which key persistent results so that results computed with
other definitions are not reused.

A fingerprint describes values by their content where it can: plain values
(numbers, strings, dates and times, and their containers), pydantic models
(such as satellites and instruments), small arrays, the nodes of design
spaces, dataclasses (such as cost models) by their fields, and functions, by
their code, the values they close over and their defaults, and the
constants (numbers, strings, dates and times, and tuples of them),
dataclasses, and functions of the same module that they refer to as globals
(not mutable containers, which are often state, such as caches). Other
objects (such as data frames) are described by their type only: identify
them by a namespace (see `ResultStore` and `SatelliteCache`). Functions of
other modules (such as TAT-C's) are described by their names, as the version
of TAT-C is part of a fingerprint.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
import functools
import hashlib
import json
import types
from datetime import date, datetime, time, timedelta
from typing import Any

import numpy as np
from pydantic import BaseModel

_MAX_DEPTH = 24
"""The deepest nesting described, beyond which values are described by type."""

_MAX_ITEMS = 1000
"""The most items of a container (or array) described one by one."""


def _type(value: Any) -> dict:
    kind = type(value)
    return {"type": f"{kind.__module__}.{kind.__qualname__}"}


def _sort_key(value: Any) -> str:
    return json.dumps(value, sort_keys=True, default=str)


def _code(code: types.CodeType, depth: int, seen: set) -> dict:
    """
    Describes compiled code by its instructions, constants, and names (not
    its file or line numbers, which change as code moves).
    """
    return {
        "code": code.co_code.hex(),
        "constants": [_describe(c, depth + 1, seen) for c in code.co_consts],
        "names": list(code.co_names),
    }


def _function(function: types.FunctionType, depth: int, seen: set) -> dict:
    """
    Describes a function: its code, the values it closes over, its defaults,
    and the constants, dataclasses, and functions of its module it refers to
    as globals.
    """
    if id(function) in seen:
        return {"function": function.__qualname__, "recursive": True}
    seen = seen | {id(function)}
    closure = [
        _describe(cell.cell_contents, depth + 1, seen)
        for cell in (function.__closure__ or ())
        if _has_contents(cell)
    ]
    referenced = {}
    module = function.__module__
    for name in _global_names(function.__code__):
        if name not in function.__globals__:
            continue
        value = function.__globals__[name]
        if isinstance(value, types.FunctionType):
            if value.__module__ == module:
                referenced[name] = _function(value, depth + 1, seen)
        elif _is_constant(value) or _is_dataclass(value):
            referenced[name] = _describe(value, depth + 1, seen)
    return {
        "function": function.__qualname__,
        **_code(function.__code__, depth, seen),
        "closure": closure,
        "defaults": _describe(function.__defaults__, depth + 1, seen),
        "keyword defaults": _describe(function.__kwdefaults__, depth + 1, seen),
        "globals": referenced,
    }


def _is_dataclass(value: Any) -> bool:
    """Whether a value is an instance of a dataclass (such as a cost model)."""
    return dataclasses.is_dataclass(value) and not isinstance(value, type)


def _has_contents(cell: Any) -> bool:
    try:
        cell.cell_contents  # pylint: disable=pointless-statement
    except ValueError:  # an empty cell
        return False
    return True


def _global_names(code: types.CodeType) -> list[str]:
    """
    Gets the names a code (and its nested code) refers to.
    """
    names = list(code.co_names)
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


def _describe(value: Any, depth: int = 0, seen: set | None = None) -> Any:
    """
    Describes a value as plain data (see the module).
    """
    # pylint: disable=too-many-return-statements,too-many-branches,import-outside-toplevel
    seen = set() if seen is None else seen
    if depth > _MAX_DEPTH:
        return _type(value)
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, np.generic):
        return _describe(value.item(), depth, seen)
    if isinstance(value, bytes):
        return hashlib.sha1(value).hexdigest()
    if isinstance(value, (datetime, date, time)):
        return value.isoformat()
    if isinstance(value, timedelta):
        return {"timedelta": value.total_seconds()}
    if _is_dataclass(value):
        return {
            **_type(value),
            **{
                item.name: _describe(getattr(value, item.name), depth + 1, seen)
                for item in dataclasses.fields(value)
            },
        }
    if isinstance(value, BaseModel):
        return {**_type(value), **value.model_dump(mode="json")}
    if isinstance(value, np.ndarray):
        if value.dtype == object:
            return {**_type(value), "shape": list(value.shape)}
        return {
            "array": str(value.dtype),
            "shape": list(value.shape),
            "data": hashlib.sha1(np.ascontiguousarray(value).tobytes()).hexdigest(),
        }
    if isinstance(value, (list, tuple)):
        if len(value) > _MAX_ITEMS:
            return {**_type(value), "length": len(value)}
        return [_describe(v, depth + 1, seen) for v in value]
    if isinstance(value, (set, frozenset)):
        if len(value) > _MAX_ITEMS:
            return {**_type(value), "length": len(value)}
        return sorted((_describe(v, depth + 1, seen) for v in value), key=_sort_key)
    if isinstance(value, dict):
        if len(value) > _MAX_ITEMS:
            return {**_type(value), "length": len(value)}
        return sorted(
            (
                [_describe(k, depth + 1, seen), _describe(v, depth + 1, seen)]
                for k, v in value.items()
            ),
            key=_sort_key,
        )
    if isinstance(value, types.CodeType):
        return _code(value, depth, seen)
    if isinstance(value, types.FunctionType):
        return _function(value, depth, seen)
    if isinstance(value, types.MethodType):
        return {
            "method": _describe(value.__func__, depth + 1, seen),
            "self": _describe(value.__self__, depth + 1, seen),
        }
    if isinstance(value, functools.partial):
        return {
            "partial": _describe(value.func, depth + 1, seen),
            "arguments": _describe(value.args, depth + 1, seen),
            "keywords": _describe(value.keywords, depth + 1, seen),
        }
    if isinstance(value, (types.BuiltinFunctionType, type)):
        return {"name": f"{value.__module__}.{value.__qualname__}"}
    from .nodes import Node

    if isinstance(value, Node):
        if id(value) in seen:
            return {**_type(value), "recursive": True}
        return {
            **_type(value),
            **{
                name: _describe(attribute, depth + 1, seen | {id(value)})
                for name, attribute in sorted(vars(value).items())
            },
        }
    return _type(value)


def get_fingerprint(*values: Any) -> str:
    """
    Gets a fingerprint of values (such as a design space and the functions
    of a study) and of the version of TAT-C: a key that changes if their
    definitions do (see the module for what it describes).

    Args:
        *values (Any): The values.

    Returns:
        str: the fingerprint
    """
    # pylint: disable=import-outside-toplevel
    from .. import __version__

    text = json.dumps(
        {"tatc": __version__, "values": [_describe(value) for value in values]},
        sort_keys=True,
        default=str,
    )
    return hashlib.sha1(text.encode()).hexdigest()[:16]
