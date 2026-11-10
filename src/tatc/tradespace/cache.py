"""
Caches of analyses by satellite, which designs of a trade study share.

Most analyses of a constellation (such as `tatc.analysis.collect_observations`
or `tatc.analysis.collect_feature_observations`) compute the result of each
satellite independently of the others, then combine them. Designs in a
trade study often share satellites (the same orbit and instruments, if not
the same name), so caching the result of each satellite avoids computing it
again for every design that flies it.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import tempfile
import time
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from ..schemas import Satellite


def get_satellite_key(satellite: Satellite) -> str:
    """
    Gets a key of a satellite's content: its orbit, instruments, and
    attitude, but not its name, so that satellites that differ only in their
    names share a key.

    Args:
        satellite (Satellite): The satellite.

    Returns:
        str: the key
    """
    data = satellite.model_dump(mode="json", exclude={"name"})
    text = json.dumps(
        {"type": type(satellite).__name__, **data}, sort_keys=True, default=str
    )
    return hashlib.sha1(text.encode()).hexdigest()


class SatelliteCache:
    """
    A cache of an analysis of single satellites, keyed by their content (see
    `get_satellite_key`), in memory and optionally on disk.

    The function computes the result of one satellite, with every other
    argument fixed (for example, by a closure over the targets and the
    analysis period): the cache assumes that its result depends on nothing
    else. A result that is a data frame with a `satellite` column is
    relabeled with the name of the satellite it is retrieved for.

    A cache on disk (`directory`) persists between sessions and is shared by
    the parallel workers of a trade study (see `Study`), each of which
    otherwise has its own cache in memory. Its `namespace` must identify the
    analysis (the function and its fixed arguments), as results in the
    directory are only keyed by the namespace and the satellite: change it
    whenever the analysis changes.

    Args:
        function (Callable): The analysis of one satellite.
        directory (str | os.PathLike | None): A directory in which to store
            results, or `None` to keep them in memory only.
        namespace (str | None): An identifier of the analysis, required with
            a `directory`.
    """

    def __init__(
        self,
        function: Callable[[Satellite], Any],
        directory: str | os.PathLike | None = None,
        namespace: str | None = None,
    ):
        if directory is not None and not namespace:
            raise ValueError("a cache on disk requires a namespace")
        self.function = function
        self.directory = None if directory is None else Path(directory)
        self.namespace = namespace
        self.hits = 0
        """int: The number of results retrieved from the cache."""
        self.misses = 0
        """int: The number of results computed."""
        self._memory: dict[str, Any] = {}
        if self.directory is not None:
            self.directory.mkdir(parents=True, exist_ok=True)

    def __getstate__(self) -> dict:
        # copies sent to parallel workers start with empty memories
        state = self.__dict__.copy()
        state["_memory"] = {}
        return state

    def _path(self, key: str) -> Path:
        name = hashlib.sha1(f"{self.namespace}/{key}".encode()).hexdigest()
        return self.directory / f"{name}.pkl"

    def _load(self, key: str) -> Any:
        if key in self._memory:
            return self._memory[key]
        if self.directory is not None:
            path = self._path(key)
            for attempt in range(5):
                if not path.exists():
                    break
                try:
                    with open(path, "rb") as file:
                        result = pickle.load(file)
                except (OSError, EOFError, pickle.UnpicklingError):
                    # on Windows, another worker may be replacing the file:
                    # try again, then compute the result instead
                    time.sleep(0.02 * (attempt + 1))
                    continue
                self._memory[key] = result
                return result
        raise KeyError(key)

    def _store(self, key: str, result: Any) -> None:
        self._memory[key] = result
        if self.directory is None:
            return
        path = self._path(key)
        if path.exists():
            # another worker stored the (same) result
            return
        # write to a temporary file, then rename it, so that a parallel
        # worker never reads a partial file
        descriptor, temporary = tempfile.mkstemp(dir=self.directory, suffix=".tmp")
        with os.fdopen(descriptor, "wb") as file:
            pickle.dump(result, file)
        try:
            os.replace(temporary, path)
        except PermissionError:
            # on Windows, another worker is replacing or reading the file,
            # which then holds the same result
            os.remove(temporary)
            if not path.exists():
                raise

    def get(self, satellite: Satellite) -> Any:
        """
        Gets the result of a satellite, computing it if it is not cached.

        Args:
            satellite (Satellite): The satellite.

        Returns:
            Any: the result
        """
        key = get_satellite_key(satellite)
        try:
            result = self._load(key)
            self.hits += 1
        except KeyError:
            result = self.function(satellite)
            self.misses += 1
            self._store(key, result)
        if isinstance(result, pd.DataFrame) and "satellite" in result.columns:
            if len(result) > 0 and (result["satellite"] != satellite.name).any():
                result = result.copy()
                result["satellite"] = satellite.name
        return result

    def __call__(self, satellites: Satellite | Iterable[Satellite]) -> Any:
        """
        Gets the results of satellites, concatenated if they are data frames
        (or the list of results otherwise).

        Args:
            satellites (Satellite | Iterable[Satellite]): The satellites.

        Returns:
            Any: the concatenated data frame, or the list of results
        """
        if isinstance(satellites, Satellite):
            satellites = [satellites]
        results = [self.get(satellite) for satellite in satellites]
        if results and all(isinstance(result, pd.DataFrame) for result in results):
            non_empty = [result for result in results if len(result) > 0]
            if not non_empty:
                return results[0]
            return pd.concat(non_empty, ignore_index=True)
        return results

    def clear(self) -> None:
        """
        Clears the cache in memory (not on disk) and its statistics.
        """
        self._memory.clear()
        self.hits = self.misses = 0

    def __repr__(self) -> str:
        return (
            f"SatelliteCache({len(self._memory)} results, "
            + f"{self.hits} hits, {self.misses} misses)"
        )
