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


def _canonical(value: Any) -> str:
    """
    Gets the canonical text of a plain value: its JSON, with the keys of
    mappings sorted (and other values as their strings).
    """
    return json.dumps(value, sort_keys=True, default=str)


def _digest(data: str | bytes, length: int | None = None) -> str:
    """Gets a hash of text (or bytes), as hexadecimal digits (the first `length`)."""
    raw = data.encode() if isinstance(data, str) else data
    text = hashlib.sha1(raw).hexdigest()
    return text if length is None else text[:length]
