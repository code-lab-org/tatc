"""
Catalog of instruments of past and present missions, in configurations compared
with their data products.

An entry's `validation` states how far each configuration has been validated:
not every entry is fully `validated`, and only `validated` configurations have
passed every test of the standard (see `tatc.catalog.standard`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .loader import get_entry, get_entry_ids, list_entries
from .schemas import (
    Configuration,
    ConfigurationRecord,
    ConfigurationValidation,
    DataProduct,
    EntryValidation,
    GateResult,
    InstrumentEntry,
    InstrumentPayload,
    InstrumentSpecification,
    Metric,
    Platform,
    Provenance,
    ValidationRecord,
)
from .standard import (
    MIN_HOLDOUT_GAP_DAYS,
    STANDARD_VERSION,
    Gate,
    Grade,
    InstrumentClass,
    Tier,
    get_gates,
    get_required_grades,
)

__all__ = [
    "Configuration",
    "ConfigurationRecord",
    "ConfigurationValidation",
    "DataProduct",
    "EntryValidation",
    "Gate",
    "GateResult",
    "Grade",
    "InstrumentClass",
    "InstrumentEntry",
    "InstrumentPayload",
    "InstrumentSpecification",
    "MIN_HOLDOUT_GAP_DAYS",
    "Metric",
    "Platform",
    "Provenance",
    "STANDARD_VERSION",
    "Tier",
    "ValidationRecord",
    "get_entry",
    "get_entry_ids",
    "get_gates",
    "get_required_grades",
    "list_entries",
]
