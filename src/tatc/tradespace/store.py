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

    Records are also keyed by a `context`: a study gives its fingerprint
    (see `Study.fingerprint`), which changes with its seed, its use of
    common random numbers, the definitions of its design space and its build
    and evaluation functions (their code, and the plain values they refer
    to), and the version of TAT-C, so that records computed with other
    definitions are not reused. The `namespace` identifies what a
    fingerprint cannot see: the data the functions use (such as the targets
    of the analysis) and the code of other modules they call. Change it
    whenever these change. Several studies may share a file, in their own
    namespaces.

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
                    "CREATE TABLE IF NOT EXISTS evaluations (namespace TEXT, "
                    + "context TEXT, design TEXT, replicate INTEGER, records BLOB, "
                    + "PRIMARY KEY (namespace, context, design, replicate))"
                )

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path, timeout=60)

    def get(
        self, keys: Iterable[tuple[str, int]], context: str = ""
    ) -> dict[tuple[str, int], list[dict]]:
        """
        Gets the stored records of evaluations.

        Args:
            keys (Iterable[tuple[str, int]]): The evaluations: the key of a
                design and its replicate.
            context (str): The context of the evaluations (such as the
                fingerprint of a study).

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
                    "SELECT design, replicate, records FROM evaluations "
                    + f"WHERE namespace = ? AND context = ? AND design IN ({marks})",
                    [self.namespace, context, *designs],
                ).fetchall()
                wanted = set(part)
                for design, replicate, blob in rows:
                    if (design, replicate) in wanted:
                        found[design, replicate] = pickle.loads(blob)
        return found

    def put(
        self, items: Iterable[tuple[str, int, list[dict]]], context: str = ""
    ) -> None:
        """
        Stores the records of evaluations, replacing any stored before.

        Args:
            items (Iterable[tuple[str, int, list[dict]]]): The evaluations: the
                key of a design, its replicate, and its records.
            context (str): The context of the evaluations (such as the
                fingerprint of a study).
        """
        rows = [
            (self.namespace, context, design, int(replicate), pickle.dumps(records))
            for design, replicate, records in items
        ]
        if rows:
            with closing(self._connect()) as connection:
                with connection:
                    connection.executemany(
                        "INSERT OR REPLACE INTO evaluations VALUES (?, ?, ?, ?, ?)",
                        rows,
                    )

    def contexts(self) -> dict[str, int]:
        """
        Gets the contexts of the namespace's stored evaluations.

        Returns:
            dict[str, int]: the number of evaluations of each context
        """
        with closing(self._connect()) as connection:
            rows = connection.execute(
                "SELECT context, COUNT(*) FROM evaluations WHERE namespace = ? "
                + "GROUP BY context",
                [self.namespace],
            ).fetchall()
        return dict(rows)

    def __len__(self) -> int:
        return sum(self.contexts().values())

    def clear(self) -> None:
        """
        Removes the stored records of the namespace (of every context).
        """
        with closing(self._connect()) as connection:
            with connection:
                connection.execute(
                    "DELETE FROM evaluations WHERE namespace = ?", [self.namespace]
                )

    def __repr__(self) -> str:
        return f"ResultStore({str(self.path)!r}, {self.namespace!r}: {len(self)} evaluations)"
