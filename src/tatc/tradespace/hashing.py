"""
Canonical text of plain values, and short hashes of it, which identify
designs (see `Design.key`), definitions (see `get_fingerprint`), and the
contents of caches and stores.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import hashlib
import json
from typing import Any


def _json_default(value: Any) -> Any:
    """
    Gets a value that JSON holds of one it does not: a numpy scalar as
    Python's, or another value as its string.
    """
    item = getattr(value, "item", None)
    if callable(item) and type(value).__module__ == "numpy":
        return item()
    return str(value)


def _canonical(value: Any) -> str:
    """
    Gets the canonical text of a plain value: its JSON, with the keys of
    mappings sorted (numpy scalars as Python's, and other values as their
    strings).
    """
    return json.dumps(value, sort_keys=True, default=_json_default)


def _digest(data: str | bytes, length: int | None = None) -> str:
    """Gets a hash of text (or bytes), as hexadecimal digits (the first `length`)."""
    raw = data.encode() if isinstance(data, str) else data
    text = hashlib.sha1(raw).hexdigest()
    return text if length is None else text[:length]
