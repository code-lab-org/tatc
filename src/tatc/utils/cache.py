"""
Methods to cache computed values on objects.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from collections.abc import Callable, Hashable
from typing import Any, TypeVar

T = TypeVar("T")

_CACHE = "_cache"
"""Name of the attribute in an object's `__dict__` that holds its cache."""


class _Cache(dict):
    """
    Values cached on an object, which are neither copied nor pickled with
    it: some (such as Skyfield satellites) cannot be pickled, as needed to
    send objects to other processes. A copy of an object (for example, from
    pydantic's `model_copy`) or an unpickled object starts with an empty
    cache.
    """

    def __init__(self, owner: int | None = None):
        """
        Initializes a cache.

        Args:
            owner (int | None): The identity (`id`) of the object that owns the cache.
        """
        super().__init__()
        self.owner = owner

    def __copy__(self) -> "_Cache":
        return _Cache()

    def __deepcopy__(self, memo: dict) -> "_Cache":
        return _Cache()

    def __reduce__(self) -> tuple:
        return (_Cache, ())


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

    The value is stored alongside the key, which identifies the inputs from
    which it is computed (such as the object's field values), in a cache in
    the object's `__dict__` that is neither copied nor pickled with the
    object. A shallow copy of a pydantic model (from `model_copy`) shares its
    `__dict__` values, so a cache owned by another object is replaced.

    Args:
        obj (Any): The object on which the value is cached.
        name (str): The name of the cached value.
        key (Hashable): The key of the inputs from which the value is computed.
        compute (Callable[[], T]): Computes the value.
        lazy_load (bool): True, if a value cached with the same key should be used.

    Returns:
        T: the value
    """
    cache = obj.__dict__.get(_CACHE)
    if not isinstance(cache, _Cache) or cache.owner != id(obj):
        cache = _Cache(id(obj))
        obj.__dict__[_CACHE] = cache
    cached = cache.get(name) if lazy_load else None
    if cached is None or cached[0] != key:
        cached = (key, compute())
        cache[name] = cached
    return cached[1]
