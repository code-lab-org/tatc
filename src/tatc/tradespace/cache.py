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
import shutil
import tempfile
import time
import warnings
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from ..schemas import Satellite
from .fingerprint import _fingerprint_of, compare_definitions, get_definition


def _write(path: Path, data: bytes, replace: bool = False) -> None:
    """
    Writes a file (if it does not exist, or to `replace` it) through a
    temporary file that is renamed, so that a parallel worker never reads a
    partial file.
    """
    if path.exists() and not replace:
        return
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    with os.fdopen(descriptor, "wb") as file:
        file.write(data)
    try:
        os.replace(temporary, path)
    except PermissionError:
        # on Windows, another worker is replacing or reading the file, which
        # then holds the same content
        os.remove(temporary)
        if not path.exists():
            raise


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
    otherwise has its own cache in memory. Results in the directory are keyed
    by the fingerprint of the function (its code, the values it closes over
    or refers to, such as data frames of targets, by their content, and the
    source of the modules of TAT-C it depends on, see `get_definition`), by
    a `namespace`, which must identify what the fingerprint cannot see (data
    the function reads, such as from files, and the code of other modules it
    calls, other than TAT-C's: change it whenever these change), and by the
    satellite, in a folder of
    each namespace (by a hash of it, with the namespace in `namespace.txt`),
    and in it a folder of each fingerprint (with the definition in
    `definition.json`). A cache warns if its namespace has results only of
    other definitions of its function, which it does not reuse, with the
    differences from the latest (see `compare_definitions`), and removes
    them with `prune`. A result that cannot be read (such as one written
    with other versions of the libraries of its data) is computed again.

    Args:
        function (Callable): The analysis of one satellite.
        directory (str | os.PathLike | None): A directory in which to store
            results, or `None` to keep them in memory only.
        namespace (str): An identifier of the analysis (of what the
            fingerprint of the function cannot see).
    """

    def __init__(
        self,
        function: Callable[[Satellite], Any],
        directory: str | os.PathLike | None = None,
        namespace: str = "",
    ):
        self.function = function
        self.directory = None if directory is None else Path(directory)
        self.namespace = namespace
        self.fingerprint = ""
        """str: The fingerprint of the function, which keys results on disk."""
        self.hits = 0
        """int: The number of results retrieved from the cache."""
        self.misses = 0
        """int: The number of results computed."""
        self._memory: dict[str, Any] = {}
        self._unreadable: set[str] = set()
        self._folder: Path | None = None
        if self.directory is not None:
            definition = get_definition(function=function)
            self.fingerprint = _fingerprint_of(definition)
            space = hashlib.sha1(namespace.encode()).hexdigest()[:16]
            self._folder = self.directory / space / self.fingerprint
            self._folder.mkdir(parents=True, exist_ok=True)
            _write(self._folder.parent / "namespace.txt", namespace.encode())
            text = json.dumps(definition, sort_keys=True, default=str)
            _write(self._folder / "definition.json", text.encode())
            self._warn_other_definitions(definition)

    def __getstate__(self) -> dict:
        # copies sent to parallel workers start with empty memories
        state = self.__dict__.copy()
        state["_memory"], state["_unreadable"] = {}, set()
        return state

    def _others(self) -> list[Path]:
        """
        Gets the folders of the results of the namespace by other definitions
        of the function (that hold results), the latest last.
        """
        folders = [
            folder
            for folder in self._folder.parent.iterdir()
            if folder.is_dir()
            and folder != self._folder
            and next(folder.glob("*.pkl"), None) is not None
        ]
        return sorted(folders, key=lambda folder: folder.stat().st_mtime)

    def _warn_other_definitions(self, definition: dict) -> None:
        """
        Warns if the namespace has results only of other definitions of the
        function, with the differences from the latest, once (of the caches
        of the definition, such as those of parallel workers).
        """
        if next(self._folder.glob("*.pkl"), None) is not None:
            return
        others = self._others()
        if not others:
            return
        try:
            # the first cache of the definition to warn marks its folder
            with open(self._folder / "warned", "x", encoding="utf-8"):
                pass
        except FileExistsError:
            return
        count = sum(1 for folder in others for _ in folder.glob("*.pkl"))
        message = (
            f"the cache has {count} results of namespace {self.namespace!r} by "
            + "other definitions of its function, which are not reused (remove "
            + "them with SatelliteCache.prune)"
        )
        try:
            before = json.loads((others[-1] / "definition.json").read_text())
        except (OSError, ValueError):
            before = None
        if before is not None:
            differences = compare_definitions(before, definition)
            message += "; the latest differs in: " + "; ".join(differences)
        warnings.warn(message, stacklevel=3)

    def prune(self) -> int:
        """
        Removes the results on disk of the namespace by other definitions of
        the function (see `SatelliteCache`).

        Returns:
            int: the number of results removed
        """
        if self._folder is None:
            return 0
        removed = 0
        for folder in self._folder.parent.iterdir():
            if folder.is_dir() and folder != self._folder:
                removed += sum(1 for _ in folder.glob("*.pkl"))
                shutil.rmtree(folder, ignore_errors=True)
        return removed

    def _path(self, key: str) -> Path:
        return self._folder / f"{key}.pkl"

    def _load(self, key: str) -> Any:
        if key in self._memory:
            return self._memory[key]
        if self._folder is not None:
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
                except Exception:  # pylint: disable=broad-exception-caught
                    # written with other versions of the libraries of its data
                    break
                self._memory[key] = result
                return result
            if path.exists():
                self._unreadable.add(key)
        raise KeyError(key)

    def _store(self, key: str, result: Any) -> None:
        self._memory[key] = result
        if self._folder is None:
            return
        path = self._path(key)
        if path.exists() and key not in self._unreadable:
            # another worker stored the (same) result
            return
        self._unreadable.discard(key)
        _write(path, pickle.dumps(result), replace=True)

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
        (or the list of results otherwise), or an empty data frame if there are
        no satellites.

        Args:
            satellites (Satellite | Iterable[Satellite]): The satellites.

        Returns:
            Any: the concatenated data frame, or the list of results
        """
        if isinstance(satellites, Satellite):
            satellites = [satellites]
        results = [self.get(satellite) for satellite in satellites]
        if not results:
            # no satellites (such as if all fail): an empty frame
            return pd.DataFrame()
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
