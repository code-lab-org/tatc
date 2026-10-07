"""
Methods to cache computed values on objects.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from collections.abc import Callable, Hashable
from typing import Any, TypeVar

T = TypeVar("T")


def get_cached(
    obj: Any,
    name: str,
    key: Hashable,
    compute: Callable[[], T],
    lazy_load: bool = True,
) -> T:
    """
    Gets a value cached on an object, computing (and caching) it if it is not
    cached with the same key or if `lazy_load` is False.

    The value is stored in the object's `__dict__` alongside the key, which
    identifies the inputs from which it is computed (such as the object's
    field values). A copy of a pydantic model (for example, from
    `model_copy(update=...)`) copies its `__dict__`, and so its cache, which
    the key then invalidates if the copy's inputs differ.

    Args:
        obj (Any): The object on which the value is cached.
        name (str): The name of the cached value.
        key (Hashable): The key of the inputs from which the value is computed.
        compute (Callable[[], T]): Computes the value.
        lazy_load (bool): True, if a value cached with the same key should be used.

    Returns:
        T: the value
    """
    cached = obj.__dict__.get(name) if lazy_load else None
    if cached is None or cached[0] != key:
        cached = (key, compute())
        obj.__dict__[name] = cached
    return cached[1]
