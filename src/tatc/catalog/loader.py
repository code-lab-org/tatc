"""
Methods to load catalog entries from packaged resources.

Files are parsed on first use, and only for the entries requested, to avoid
slowing the import of the package.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import functools
import json
import pathlib

import yaml

from .schemas import InstrumentEntry
from .standard import InstrumentClass, Tier

_RESOURCES = pathlib.Path(__file__).parents[1] / "resources"
_ENTRIES = _RESOURCES / "instruments"
_RECORDS = _RESOURCES / "validation"


def get_entry_ids() -> list[str]:
    """
    Gets the identifiers of the catalog entries.

    Returns:
        list[str]: The sorted identifiers.
    """
    return sorted(path.stem for path in _ENTRIES.glob("*.yml"))


@functools.lru_cache(maxsize=None)
def _load_entry(entry_id: str) -> InstrumentEntry:
    """
    Loads an entry and its validation record, if any, from files.

    Args:
        entry_id (str): The identifier of the entry.

    Returns:
        InstrumentEntry: The entry, which callers must not modify.
    """
    with open(_ENTRIES / f"{entry_id}.yml", "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if data.get("id") != entry_id:
        raise ValueError(f"the file of entry {entry_id!r} has id {data.get('id')!r}")
    record_path = _RECORDS / f"{entry_id}.json"
    if record_path.exists():
        with open(record_path, "r", encoding="utf-8") as f:
            data["record"] = json.load(f)
    return InstrumentEntry(**data)


def get_entry(entry_id: str) -> InstrumentEntry:
    """
    Gets a catalog entry by its identifier or one of its aliases (ignoring case).

    Args:
        entry_id (str): The identifier or alias.

    Returns:
        InstrumentEntry: A copy of the entry.

    Raises:
        KeyError: If there is no such entry.
    """
    wanted = entry_id.casefold()
    for candidate in get_entry_ids():
        if wanted == candidate:
            return _load_entry(candidate).model_copy(deep=True)
    for candidate in get_entry_ids():
        entry = _load_entry(candidate)
        names = (entry.name, *entry.aliases)
        if wanted in (name.casefold() for name in names):
            return entry.model_copy(deep=True)
    raise KeyError(f"no catalog entry {entry_id!r}: {get_entry_ids()}")


def list_entries(
    instrument_class: InstrumentClass | str | None = None,
    tier: Tier | str | None = None,
) -> list[InstrumentEntry]:
    """
    Lists catalog entries, optionally of one class or validation tier.

    Args:
        instrument_class (InstrumentClass | str | None): The class, if any.
        tier (Tier | str | None): The validation tier, if any.

    Returns:
        list[InstrumentEntry]: Copies of the entries, sorted by identifier.
    """
    wanted_class = (
        None if instrument_class is None else InstrumentClass(instrument_class)
    )
    wanted_tier = None if tier is None else Tier(tier)
    return [
        entry
        for entry in map(get_entry, get_entry_ids())
        if (wanted_class is None or entry.instrument_class == wanted_class)
        and (wanted_tier is None or entry.validation.tier == wanted_tier)
    ]
