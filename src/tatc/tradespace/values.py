"""
Views of the values of designs, the scopes given to the callables of design
spaces, and the paths of values (see `Values`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import difflib
import re
from collections.abc import Iterable, Iterator, Mapping
from typing import Any

import numpy as np


class _ChoiceValue:
    """The value of a choice: the chosen branch and its value."""

    __slots__ = ("choice", "value")

    def __init__(self, choice: str, value: Any = None):
        self.choice, self.value = choice, value


class _RepeatValue:
    """The value of a repetition: its count and its items."""

    __slots__ = ("count", "items")

    def __init__(self, count: int):
        self.count, self.items = count, []


def _view(value: Any) -> Any:
    """
    Wraps a decided value for attribute access.
    """
    if isinstance(value, dict):
        return Values(value)
    if isinstance(value, _ChoiceValue):
        return ChoiceValues(value)
    if isinstance(value, _RepeatValue):
        return RepeatValues(value)
    return value


_SEGMENT = re.compile(r"\[(\d+)\]|([^.\[\]]+)")
"""A segment of a path: the index of an item (`[k]`) or a name."""


def _segments(path: str) -> list[str | int]:
    """
    Gets the segments of a path of values (see `Values`), such as
    `["layout", "sso", "planes", 0, "crossing"]` of
    `"layout.sso.planes[0].crossing"`.
    """
    segments: list[str | int] = []
    for index, name in _SEGMENT.findall(path):
        segments.append(int(index) if index else name)
    rebuilt = "".join(f"[{s}]" if isinstance(s, int) else f".{s}" for s in segments)
    if not segments or rebuilt.lstrip(".") != path.lstrip("."):
        raise KeyError(f"{path!r} is not a path of values")
    return segments


def _join(path: str, segment: str | int) -> str:
    """Gets the path of a segment within a path (see `_segments`)."""
    if isinstance(segment, int):
        return f"{path}[{segment}]"
    return f"{path}.{segment}" if path else segment


_ITEM_INDEX = re.compile(r"\[\d+\]")
"""The index of an item of a repetition in a path (`[k]`)."""


def generic_path(path: str) -> str:
    """
    Gets the path of the node of a value: its path (see `Design.flat`) with
    `[*]` for the index of any item of a repetition, as the paths of
    uncertain parameters (see `DesignSpace.uncertain_paths`), shifts, and
    the overrides of scenarios name them (such as `planes[*].phase` of
    `planes[2].phase`).

    Args:
        path (str): The path of a value.

    Returns:
        str: the path, with `[*]` for the indices of items
    """
    return _ITEM_INDEX.sub("[*]", path)


def suggestion(name: Any, names: Iterable[Any]) -> str:
    """
    Gets a suggestion of the name that a misspelled name may be: the closest
    of the names, as `": did you mean 'fov'?"`, or `""` if none is close.

    Args:
        name (Any): The misspelled name.
        names (Iterable[Any]): The names.

    Returns:
        str: the suggestion, or `""`
    """
    close = difflib.get_close_matches(str(name), [str(n) for n in names], n=1)
    return f": did you mean {close[0]!r}?" if close else ""


def named(name: Any, names: Iterable[Any]) -> str:
    """
    Gets a name with a suggestion of the name it may be a misspelling of (see
    `suggestion`), such as `"'fvo' (did you mean 'fov'?)"`.

    Args:
        name (Any): The name.
        names (Iterable[Any]): The names it may be a misspelling of.

    Returns:
        str: the name (quoted), and the suggestion, if any
    """
    hint = suggestion(name, names)
    return f"{name!r} ({hint[2:]})" if hint else repr(name)


def missing_value(name: Any, names: Iterable[Any]) -> str:
    """
    Gets the message of a value of a name that is not among the names of
    values, such as `"no value 'fvo' (it has ['fov', 'n']: did you mean
    'fov'?)"`.

    Args:
        name (Any): The name.
        names (Iterable[Any]): The names of the values.

    Returns:
        str: the message
    """
    names = list(names)
    return f"no value {name!r} (it has {names}{suggestion(name, names)})"


def is_number(value: Any) -> bool:
    """
    Whether a value is a number: an integer or a float (also of numpy), but
    not a boolean.

    Args:
        value (Any): The value.

    Returns:
        bool: `True` if the value is a number
    """
    return isinstance(value, (int, float, np.number)) and not isinstance(
        value, (bool, np.bool_)
    )


def _step(value: Any, segment: str | int, path: str) -> Any:
    """
    Gets the value of a segment of a path within a value (see `Values`): a
    child of a group, an item of a repetition (`[k]`) or its `count`, or, of
    a choice, its chosen branch (by its name) or a value of that branch.
    """
    where = path or "the design"
    if isinstance(value, RepeatValues):
        if isinstance(segment, int):
            if not 0 <= segment < len(value):
                raise KeyError(f"{where} has {len(value)} items, not an item {segment}")
            return value[segment]
        if segment == "count":
            return value.count
        raise KeyError(f"{where} is a repetition: give an item ([k]) or its count")
    if isinstance(value, ChoiceValues):
        if segment == value.choice:
            return value.value
        try:
            return _step(value.value, segment, path)
        except KeyError:
            raise KeyError(
                f"{where} chose {value.choice!r}, which has no value {segment!r}"
                + suggestion(segment, _names_of(value.value))
            ) from None
    if isinstance(value, Values) and isinstance(segment, str):
        names = list(value)
        if segment in names:
            return value[segment]
        raise KeyError(f"{where} has {missing_value(segment, names)}")
    raise KeyError(f"{where} is a value ({value!r}), which has no {segment!r}")


def _names_of(value: Any) -> list[str]:
    """Gets the names of the values of a value (of a group), if any."""
    return list(value) if isinstance(value, Values) else []


def _resolve_path(value: Any, path: str) -> Any:
    """Gets the value at a path within a value (see `Values`)."""
    walked = ""
    for segment in _segments(path):
        value = _step(value, segment, walked)
        walked = _join(walked, segment)
    return value


class Values(Mapping):
    """
    A read-only view of the values of a design space (a group of nodes), by
    attribute or by key. A key may also be a path of values, as the columns
    of the results of a study name them (see `Design.flat`): `"a.b"` for the
    value `b` of group `a`, `"c.<branch>.d"` for the value `d` of the chosen
    branch of choice `c` (an error if it chose another), `"e.count"` for the
    count of repetition `e`, and `"e[k].f"` for the value `f` of its item
    `k`. Fixed values, which the columns of results omit, are also given::

        d["layout.walker.planes"]  # d.layout.walker.planes, if it chose walker
        d["layout.sso.planes[0].crossing"]
    """

    def __init__(self, values: dict, index: int | None = None):
        self._values, self._index = values, index

    def __getitem__(self, name: str) -> Any:
        if name in self._values:
            return _view(self._values[name])
        if isinstance(name, str) and ("." in name or "[" in name):
            return _resolve_path(self, name)
        if not isinstance(name, str):
            raise KeyError(name)
        raise KeyError(missing_value(name, self._values))

    def __contains__(self, name: object) -> bool:
        if name in self._values:
            return True
        if not isinstance(name, str) or not ("." in name or "[" in name):
            return False
        try:
            _resolve_path(self, name)
        except KeyError:
            return False
        return True

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        if name == "index" and self._index is not None:
            return self._index
        try:
            return self[name]
        except KeyError:
            raise AttributeError(
                missing_value(name, self._values), name=name, obj=self
            ) from None

    def __dir__(self) -> list[str]:
        # the values, for completion (and suggestions of misspelled names)
        index = ["index"] if self._index is not None else []
        return sorted({*object.__dir__(self), *self._values, *index})

    def __iter__(self) -> Iterator[str]:
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)

    def __repr__(self) -> str:
        return f"Values({ {k: _view(v) for k, v in self._values.items()} })"


class ChoiceValues:
    """
    A read-only view of the value of a choice: the name of the chosen branch
    (`choice`), which the view equals (as the columns of the results of a
    study give it, see `Design.flat`), and, as attributes (or keys), the
    values of the branch, which are also those of the branch's name::

        d.layout == "walker"  # d.layout.choice == "walker"
        d.layout.planes  # d.layout.walker.planes, if it chose walker

    An attribute is a value of the branch (of a branch that is a group, a
    choice, or a repetition) if the branch has one of that name, and
    otherwise the branch, if it is the branch's name; a path (see `Values`)
    names the branch first, as the columns of results do.
    """

    def __init__(self, value: _ChoiceValue):
        self._value = value

    @property
    def choice(self) -> str:
        """The name of the chosen branch."""
        return self._value.choice

    @property
    def value(self) -> Any:
        """The value of the chosen branch."""
        return _view(self._value.value)

    def __getattr__(self, name: str) -> Any:
        # the values of the branch first (as before branches were attributes),
        # then the branch by its name (which paths give first, see `Values`)
        if name.startswith("_"):
            raise AttributeError(name)
        value = self.value
        if isinstance(value, (Values, ChoiceValues, RepeatValues)):
            # not the attributes of a branch that is a plain value (a variable)
            try:
                return getattr(value, name)
            except AttributeError:
                pass
        if name == self._value.choice:
            return value
        raise AttributeError(
            f"the choice chose {self.choice!r}, which has no value {name!r}"
            + suggestion(name, _names_of(value)),
            name=name,
            obj=self,
        )

    def __dir__(self) -> list[str]:
        # the values of the branch and its name, for completion (and
        # suggestions of misspelled names)
        value = self.value
        branch = isinstance(value, (Values, ChoiceValues, RepeatValues))
        names = dir(value) if branch else []
        return sorted({*object.__dir__(self), *names, self._value.choice})

    def __getitem__(self, name: str) -> Any:
        return _resolve_path(self, name)

    def __contains__(self, name: object) -> bool:
        try:
            self[name]  # type: ignore[index]
        except (KeyError, TypeError):
            return False
        return True

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.choice == other
        if isinstance(other, ChoiceValues):
            return self.choice == other.choice and self.value == other.value
        return NotImplemented

    def __hash__(self) -> int:
        return hash(self.choice)

    __iter__ = None  # type: ignore[assignment]

    def __str__(self) -> str:
        return self.choice

    def __repr__(self) -> str:
        return f"ChoiceValues({self.choice!r}, {self.value!r})"


class RepeatValues:
    """
    A read-only view of the value of a repetition: a sequence of the values
    of its items, whose length is also given as `count`.
    """

    def __init__(self, value: _RepeatValue):
        self._value = value

    @property
    def count(self) -> int:
        """The number of items."""
        return self._value.count

    def __getitem__(self, index: int | slice) -> Any:
        if isinstance(index, slice):
            return [self[k] for k in range(len(self))[index]]
        position = index + len(self) if index < 0 else index
        if not 0 <= position < len(self):
            raise IndexError(f"the repetition has {len(self)} items, not item {index}")
        index = position
        item = self._value.items[index]
        return Values(item, index) if isinstance(item, dict) else _view(item)

    def __len__(self) -> int:
        return len(self._value.items)

    def __iter__(self) -> Iterator[Any]:
        return (self[k] for k in range(len(self)))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, (RepeatValues, list, tuple)):
            return list(self) == list(other)
        return NotImplemented

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f"RepeatValues({list(self)!r})"


class Scope:
    """
    The values decided so far, as given to the callables of a design space.
    A name resolves to the innermost enclosing group that declares it (and
    has decided it), and `index` to the position of the innermost enclosing
    repeated item.
    """

    def __init__(self, frames: tuple = ()):
        self._frames = frames

    def push(self, values: dict, index: int | None = None) -> Scope:
        """
        Gets the scope of a group within this scope (or of a repeated item,
        at `index`).
        """
        return Scope(self._frames + ((values, index),))

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        for values, index in reversed(self._frames):
            if name in values:
                return _view(values[name])
            if name == "index" and index is not None:
                return index
        raise AttributeError(f"{name!r} is not decided before it is used")

    def get(self, name: str, default: Any = None) -> Any:
        """
        Gets a value by its name, or by a path within it (see `__getitem__`),
        if it is decided.

        Args:
            name (str): The name or path.
            default (Any): The value if it is not decided.

        Returns:
            Any: the value
        """
        try:
            return self[name]
        except KeyError:
            return default

    def __contains__(self, name: object) -> bool:
        if not isinstance(name, str):
            return False
        try:
            self[name]
        except KeyError:
            return False
        return True

    def __getitem__(self, name: str) -> Any:
        """
        Gets a value by its name, or by a path within it (see `Values`), such
        as `"layout.walker.planes"`.
        """
        segments = _segments(name)
        if not isinstance(segments[0], str):
            raise KeyError(f"{name!r} does not start with a name")
        try:
            value = getattr(self, segments[0])
        except AttributeError as error:
            raise KeyError(str(error)) from None
        walked = segments[0]
        for segment in segments[1:]:
            value = _step(value, segment, walked)
            walked = _join(walked, segment)
        return value
