"""
Views of the values of designs, and the scopes given to the callables of
design spaces.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterator, Mapping
from typing import Any


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


class Values(Mapping):
    """
    A read-only view of the values of a design space (a group of nodes), by
    attribute or by key.
    """

    def __init__(self, values: dict, index: int | None = None):
        self._values, self._index = values, index

    def __getitem__(self, name: str) -> Any:
        return _view(self._values[name])

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_"):
            raise AttributeError(name)
        if name == "index" and self._index is not None:
            return self._index
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name) from None

    def __iter__(self) -> Iterator[str]:
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)

    def __repr__(self) -> str:
        return f"Values({ {k: _view(v) for k, v in self._values.items()} })"


class ChoiceValues:
    """
    A read-only view of the value of a choice: the name of the chosen branch
    (`choice`) and, as attributes, the values of the branch.
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
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self.value, name)

    def __getitem__(self, name: str) -> Any:
        return self.value[name]

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

    def __getitem__(self, index: int) -> Any:
        item = self._value.items[index]
        return Values(item, index) if isinstance(item, dict) else _view(item)

    def __len__(self) -> int:
        return len(self._value.items)

    def __iter__(self) -> Iterator[Any]:
        return (self[k] for k in range(len(self)))

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

    def __getitem__(self, name: str) -> Any:
        return getattr(self, name)


def configure(base: Any, values: Any = None, **overrides: Any) -> Any:
    """
    Configures a dataclass (such as a `RestorationPolicy` or a
    `DeploymentSchedule`) from the values of a design: the fields of `base`
    that the values name (such as those of `policy_space`, or design
    variables named after fields) are replaced by them, as are those given
    as keyword arguments. Values that name no field are ignored.

    Args:
        base (Any): The dataclass instance, whose other fields are kept.
        values (Any): The values of a design (a `Values` or `ChoiceValues`
            view, or a mapping), or `None`.
        **overrides (Any): Other values of fields.

    Returns:
        Any: a copy of `base` with the values of its fields
    """
    if not dataclasses.is_dataclass(base) or isinstance(base, type):
        raise TypeError("configure requires a dataclass instance")
    if isinstance(values, ChoiceValues):
        values = values.value
    names = {field.name for field in dataclasses.fields(base)}
    changes = {name: values[name] for name in (values or {}) if name in names}
    unknown = set(overrides) - names
    if unknown:
        raise ValueError(f"no fields named {sorted(unknown)}")
    return dataclasses.replace(base, **changes, **overrides)
