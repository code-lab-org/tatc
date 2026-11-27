"""
Catalog of instruments of past and present missions, in configurations compared
with their data products, and of systems (satellite programs, radar networks,
surface networks, and networks of ground stations), in configurations compared
with reference data such as element sets and station metadata, from which to
build baseline systems.

An entry's `validation` states how far each configuration has been validated:
not every entry is fully `validated`, and only `validated` configurations have
passed every test of the standard (see `tatc.catalog.standard`, which also
evaluates the records, and `tatc.catalog.system_standard` for systems). The
schemas of entries are in `tatc.catalog.schemas` and `tatc.catalog.system_schemas`,
and of their records in `tatc.catalog.records`. Instrument entries list the
satellites that carry them by the identifiers of the systems' spacecraft
(`get_spacecraft`).

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from .loader import (
    build_payload,
    clear_cache,
    get_entry,
    get_entry_ids,
    get_spacecraft,
    get_spacecraft_system,
    get_system,
    get_system_ids,
    get_transmitters,
    list_entries,
    list_systems,
)
from .schemas import InstrumentEntry, InstrumentPayload
from .standard import (
    STANDARD_VERSION,
    ConfigurationValidation,
    EntryValidation,
    GateResult,
    Grade,
    InstrumentClass,
    Tier,
)
from .system_schemas import (
    GroundNetwork,
    RadarNetwork,
    SatelliteSystem,
    Spacecraft,
    SurfaceNetwork,
    SystemEntry,
)
from .system_standard import SYSTEM_STANDARD_VERSION, SystemClass, SystemGrade

__all__ = [
    "ConfigurationValidation",
    "EntryValidation",
    "GateResult",
    "GroundNetwork",
    "Grade",
    "InstrumentClass",
    "InstrumentEntry",
    "InstrumentPayload",
    "RadarNetwork",
    "STANDARD_VERSION",
    "SYSTEM_STANDARD_VERSION",
    "SatelliteSystem",
    "Spacecraft",
    "SurfaceNetwork",
    "SystemClass",
    "SystemEntry",
    "SystemGrade",
    "Tier",
    "build_payload",
    "clear_cache",
    "get_entry",
    "get_entry_ids",
    "get_spacecraft",
    "get_spacecraft_system",
    "get_system",
    "get_transmitters",
    "get_system_ids",
    "list_entries",
    "list_systems",
]
