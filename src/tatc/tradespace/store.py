"""
Persistent stores of the results of trade studies, so that a study that is
interrupted resumes where it stopped, and a study that is extended (with
more designs or replicates) evaluates only what is new.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import base64
import json
import math
import os
import pickle
import sqlite3
from collections.abc import Callable, Iterable
from contextlib import closing
from datetime import date, datetime, time, timedelta, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .hashing import _canonical

_RESULTS_TABLE = (
    "CREATE TABLE IF NOT EXISTS results (namespace TEXT, context TEXT, "
    + "design TEXT, replicate INTEGER, variant TEXT, records TEXT, "
    + "PRIMARY KEY (namespace, context, design, replicate, variant))"
)
"""The table of the records of evaluations."""

_JSON = (bool, int, str)
"""Types whose values are kept as they are in JSON (exactly these types)."""


def _encode(value: Any) -> Any:
    """
    Encodes a value as JSON data: plain values as they are, and others
    tagged (as `{"$tag": data}`), so that they are decoded exactly (numpy
    scalars as Python's), or pickled if JSON cannot hold them.
    """
    # pylint: disable=too-many-return-statements
    kind = type(value)
    if value is None or kind in _JSON:
        return value
    if kind is float:
        return value if math.isfinite(value) else {"$float": repr(value)}
    if isinstance(value, np.generic):
        item = value.item()
        if type(item) in (*_JSON, float):
            return _encode(item)
    elif kind is list:
        return [_encode(item) for item in value]
    elif kind is tuple:
        return {"$tuple": [_encode(item) for item in value]}
    elif kind is dict:
        if all(type(key) is str and not key.startswith("$") for key in value):
            return {key: _encode(item) for key, item in value.items()}
        return {"$dict": [[_encode(k), _encode(v)] for k, v in value.items()]}
    elif kind is datetime and value.tzinfo in (None, timezone.utc):
        return {"$datetime": value.isoformat()}
    elif kind is date:
        return {"$date": value.isoformat()}
    elif kind is time and value.tzinfo is None:
        return {"$time": value.isoformat()}
    elif kind is timedelta:
        return {"$timedelta": [value.days, value.seconds, value.microseconds]}
    return {"$pickle": base64.b64encode(pickle.dumps(value)).decode("ascii")}


_DECODERS: dict[str, Callable[[Any], Any]] = {
    "$float": float,
    "$tuple": lambda items: tuple(_decode(item) for item in items),
    "$dict": lambda pairs: {_decode(k): _decode(v) for k, v in pairs},
    "$datetime": datetime.fromisoformat,
    "$date": date.fromisoformat,
    "$time": time.fromisoformat,
    "$timedelta": lambda parts: timedelta(*parts),
    "$pickle": lambda text: pickle.loads(base64.b64decode(text)),
}


def _decode(value: Any) -> Any:
    """Decodes a value from JSON data (see `_encode`)."""
    if isinstance(value, list):
        return [_decode(item) for item in value]
    if isinstance(value, dict):
        if len(value) == 1:
            tag, data = next(iter(value.items()))
            if tag in _DECODERS:
                return _DECODERS[tag](data)
        return {key: _decode(item) for key, item in value.items()}
    return value


class ResultStore:
    """
    A store of the records of the evaluations of a study (see `Study.run`),
    by the key of each design, its replicate, and the variant of the
    replicate, in a SQLite database file. The variant of a plain replicate
    is `""`; that of a replicate whose parameters are drawn as in other
    replicates (see `Study.run`), or that is evaluated in a scenario (see
    `Study`), names the shift and a hash of the overrides of the scenario
    (by their content, not the scenario's name), so that the evaluations of
    a scenario are those of any scenario of the same overrides, and a
    scenario without overrides shares the evaluations of plain replicates.

    Records are also keyed by a `context`: a study gives its fingerprint
    (see `Study.fingerprint`), which changes with its seed, its use of
    common random numbers, the definitions of its design space and its build
    and evaluation functions (their code, and the plain values they refer
    to), and the source of the modules of TAT-C they depend on, so that
    records computed with other definitions are not reused. The store keeps
    the definition of each context (see `define`), from which a study
    explains why it does not reuse records (see `compare_definitions`). A
    fingerprint describes the data that the functions are given (such as
    data frames of targets they close over or refer to, by their content),
    but not data that they read (such as from files) nor the code of other
    modules they call (other than TAT-C's): a `namespace` identifies these,
    and must change whenever they do. Studies of other definitions (such as
    variants of a study) may share a namespace, as their fingerprints keep
    their records apart, and several namespaces may share a file.

    Records are kept as JSON (in the column `records` of the table `results`,
    a list of records per evaluation), which other tools can read: numbers
    (numpy's as Python's), strings, and booleans as they are, and other
    values (such as infinite numbers, dates and times, and tuples) tagged,
    as `{"$tag": data}`, so that they are read exactly as they were written
    (or pickled, as `{"$pickle": data}`, if JSON cannot hold them).

    Args:
        path (str | os.PathLike): The database file (created if needed).
        namespace (str): An identifier of what the fingerprints of the
            studies cannot see (such as data they read from files), if any.
    """

    def __init__(self, path: str | os.PathLike, namespace: str = ""):
        self.path, self.namespace = Path(path), namespace
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with closing(self._connect()) as connection:
            connection.isolation_level = None
            connection.execute("BEGIN IMMEDIATE")
            try:
                columns = {
                    row[1] for row in connection.execute("PRAGMA table_info(results)")
                }
                if columns and "variant" not in columns:
                    # evaluations stored before replicates had variants are of
                    # plain replicates
                    connection.execute("ALTER TABLE results RENAME TO _invariant")
                    connection.execute(_RESULTS_TABLE)
                    connection.execute(
                        "INSERT INTO results SELECT namespace, context, design, "
                        + "replicate, '', records FROM _invariant"
                    )
                    connection.execute("DROP TABLE _invariant")
                connection.execute(_RESULTS_TABLE)
                connection.execute(
                    "CREATE TABLE IF NOT EXISTS definitions (namespace TEXT, "
                    + "context TEXT, definition TEXT, defined TEXT, "
                    + "PRIMARY KEY (namespace, context))"
                )
                connection.execute("COMMIT")
            except BaseException:
                connection.execute("ROLLBACK")
                raise

    def _connect(self) -> sqlite3.Connection:
        return sqlite3.connect(self.path, timeout=60)

    def get(
        self, keys: Iterable[tuple[str, int]], context: str = "", variant: str = ""
    ) -> dict[tuple[str, int], list[dict]]:
        """
        Gets the stored records of evaluations.

        Args:
            keys (Iterable[tuple[str, int]]): The evaluations: the key of a
                design and its replicate.
            context (str): The context of the evaluations (such as the
                fingerprint of a study).
            variant (str): The variant of the replicates (`""` for plain
                replicates).

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
                    "SELECT design, replicate, records FROM results WHERE "
                    + "namespace = ? AND context = ? AND variant = ? "
                    + f"AND design IN ({marks})",
                    [self.namespace, context, variant, *designs],
                ).fetchall()
                wanted = set(part)
                for design, replicate, text in rows:
                    if (design, replicate) in wanted:
                        found[design, replicate] = _decode(json.loads(text))
        return found

    def put(
        self,
        items: Iterable[tuple[str, int, list[dict]]],
        context: str = "",
        variant: str = "",
    ) -> None:
        """
        Stores the records of evaluations, replacing any stored before.

        Args:
            items (Iterable[tuple[str, int, list[dict]]]): The evaluations: the
                key of a design, its replicate, and its records.
            context (str): The context of the evaluations (such as the
                fingerprint of a study).
            variant (str): The variant of the replicates (`""` for plain
                replicates).
        """
        rows = [
            (
                self.namespace,
                context,
                design,
                int(replicate),
                variant,
                json.dumps(_encode(records), allow_nan=False),
            )
            for design, replicate, records in items
        ]
        if rows:
            with closing(self._connect()) as connection:
                with connection:
                    connection.executemany(
                        "INSERT OR REPLACE INTO results VALUES (?, ?, ?, ?, ?, ?)",
                        rows,
                    )

    def define(self, context: str, definition: dict) -> None:
        """
        Keeps the definition of a context (such as that of a study, see
        `Study.definition`), if it is not kept already.

        Args:
            context (str): The context (such as the fingerprint of a study).
            definition (dict): Its definition (see `get_definition`).
        """
        defined = datetime.now(timezone.utc).isoformat(timespec="seconds")
        with closing(self._connect()) as connection:
            with connection:
                connection.execute(
                    "INSERT OR IGNORE INTO definitions VALUES (?, ?, ?, ?)",
                    [
                        self.namespace,
                        context,
                        _canonical(definition),
                        defined,
                    ],
                )

    def definition(self, context: str) -> dict | None:
        """
        Gets the definition of a context (see `define`).

        Args:
            context (str): The context.

        Returns:
            dict | None: the definition, or `None` if it is not kept
        """
        with closing(self._connect()) as connection:
            row = connection.execute(
                "SELECT definition FROM definitions WHERE namespace = ? AND context = ?",
                [self.namespace, context],
            ).fetchone()
        return None if row is None else json.loads(row[0])

    def contexts(self) -> dict[str, int]:
        """
        Gets the contexts of the namespace: those of its stored evaluations
        and of its kept definitions (see `define`, such as of a study whose
        evaluations all failed), in the order they were defined (those not
        defined first).

        Returns:
            dict[str, int]: the number of evaluations of each context (0 for
                a definition without evaluations)
        """
        with closing(self._connect()) as connection:
            rows = connection.execute(
                "SELECT c.context, COALESCE(r.evaluations, 0) FROM "
                + "(SELECT context FROM results WHERE namespace = ? UNION "
                + "SELECT context FROM definitions WHERE namespace = ?) AS c "
                + "LEFT JOIN (SELECT context, COUNT(*) AS evaluations FROM results "
                + "WHERE namespace = ? GROUP BY context) AS r ON r.context = c.context "
                + "LEFT JOIN definitions AS d ON d.namespace = ? "
                + "AND d.context = c.context "
                + "ORDER BY d.defined, d.rowid, c.context",
                [self.namespace] * 4,
            ).fetchall()
        return dict(rows)

    def __len__(self) -> int:
        return sum(self.contexts().values())

    def remove(self, contexts: str | Iterable[str]) -> int:
        """
        Removes the stored records (and definitions) of contexts of the
        namespace, such as those of earlier definitions of a study (see
        `Study.prune`, which removes those of the study with other code,
        also of definitions without records, and not those of other studies
        that share the namespace).

        Args:
            contexts (str | Iterable[str]): The contexts.

        Returns:
            int: the number of evaluations removed
        """
        contexts = [contexts] if isinstance(contexts, str) else list(contexts)
        removed = 0
        with closing(self._connect()) as connection:
            with connection:
                for context in contexts:
                    removed += connection.execute(
                        "DELETE FROM results WHERE namespace = ? AND context = ?",
                        [self.namespace, context],
                    ).rowcount
                    connection.execute(
                        "DELETE FROM definitions WHERE namespace = ? AND context = ?",
                        [self.namespace, context],
                    )
        return removed

    def clear(self) -> None:
        """
        Removes the stored records (and definitions) of the namespace (of
        every context).
        """
        with closing(self._connect()) as connection:
            with connection:
                for table in ("results", "definitions"):
                    connection.execute(
                        f"DELETE FROM {table} WHERE namespace = ?", [self.namespace]
                    )

    def __repr__(self) -> str:
        return f"ResultStore({str(self.path)!r}, {self.namespace!r}: {len(self)} evaluations)"
