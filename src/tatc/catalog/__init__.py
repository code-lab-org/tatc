"""
Catalog of instruments of past and present missions, in configurations compared
with their data products.

An entry's `validation` states how far each configuration has been validated:
not every entry is fully `validated`, and only `validated` configurations have
passed every test of the standard (see `tatc.catalog.standard`). The schemas of
entries and their records are in `tatc.catalog.schemas`.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .loader import build_payload, get_entry, get_entry_ids, list_entries
from .schemas import (
    ConfigurationValidation,
    EntryValidation,
    GateResult,
    InstrumentEntry,
    InstrumentPayload,
)
from .standard import STANDARD_VERSION, Grade, InstrumentClass, Tier

__all__ = [
    "ConfigurationValidation",
    "EntryValidation",
    "GateResult",
    "Grade",
    "InstrumentClass",
    "InstrumentEntry",
    "InstrumentPayload",
    "STANDARD_VERSION",
    "Tier",
    "build_payload",
    "get_entry",
    "get_entry_ids",
    "list_entries",
]
