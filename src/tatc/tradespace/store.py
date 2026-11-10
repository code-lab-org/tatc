"""
Persistent stores of the results of trade studies, so that a study that is
interrupted resumes where it stopped, and a study that is extended (with
more designs or replicates) evaluates only what is new.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import os
import pickle
import sqlite3
from collections.abc import Iterable
from contextlib import closing
from pathlib import Path


class ResultStore:
    """
    A store of the records of the evaluations of a study (see `Study.run`),
    by the key of each design and its replicate, in a SQLite database file.

    The `namespace` identifies the study: its design space's build and
    evaluation functions and their fixed inputs (such as the targets and the
    period of the analysis), as the records are only keyed by the namespace,
    the design, and the replicate. Change it whenever the study changes.
    Several studies may share a file, in their own namespaces.

    Args:
        path (str | os.PathLike): The database file (created if needed).
        namespace (str): An identifier of the study.
    """

    def __init__(self, path: str | os.PathLike, namespace: str):
        if not namespace:
            raise ValueError("a store requires a namespace")
        self.path, self.namespace = Path(path), namespace
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with closing(self._connect()) as connection:
            with connection:
                connection.execute(
                    "CREATE TABLE IF NOT EXISTS results (namespace TEXT, design TEXT, "
                    + "replicate INTEGER, records BLOB, "
                    + "PRIMARY KEY (namespace, design, replicate))"
                )

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path, timeout=60)

    def get(self, keys: Iterable[tuple[str, int]]) -> dict[tuple[str, int], list[dict]]:
        """
        Gets the stored records of evaluations.

        Args:
            keys (Iterable[tuple[str, int]]): The evaluations: the key of a
                design and its replicate.

        Returns:
            dict: the records of each stored evaluation (others are omitted)
        """
        keys = list(keys)
        found = {}
        with closing(self._connect()) as connection:
            for start in range(0, len(keys), 500):
                part = keys[start : start + 500]
                designs = sorted({design for design, _ in part})
                marks = ",".join("?" * len(designs))
                rows = connection.execute(
                    "SELECT design, replicate, records FROM results "
                    + f"WHERE namespace = ? AND design IN ({marks})",
                    [self.namespace, *designs],
                ).fetchall()
                wanted = set(part)
                for design, replicate, blob in rows:
                    if (design, replicate) in wanted:
                        found[design, replicate] = pickle.loads(blob)
        return found

    def put(self, items: Iterable[tuple[str, int, list[dict]]]) -> None:
        """
        Stores the records of evaluations, replacing any stored before.

        Args:
            items (Iterable[tuple[str, int, list[dict]]]): The evaluations: the
                key of a design, its replicate, and its records.
        """
        rows = [
            (self.namespace, design, int(replicate), pickle.dumps(records))
            for design, replicate, records in items
        ]
        if rows:
            with closing(self._connect()) as connection:
                with connection:
                    connection.executemany(
                        "INSERT OR REPLACE INTO results VALUES (?, ?, ?, ?)", rows
                    )

    def __len__(self) -> int:
        with closing(self._connect()) as connection:
            (count,) = connection.execute(
                "SELECT COUNT(*) FROM results WHERE namespace = ?", [self.namespace]
            ).fetchone()
        return count

    def clear(self) -> None:
        """
        Removes the stored records of the namespace.
        """
        with closing(self._connect()) as connection:
            with connection:
                connection.execute(
                    "DELETE FROM results WHERE namespace = ?", [self.namespace]
                )

    def __repr__(self) -> str:
        return f"ResultStore({str(self.path)!r}, {self.namespace!r}: {len(self)} evaluations)"
