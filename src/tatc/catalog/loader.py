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
from collections.abc import Sequence
from datetime import datetime

import yaml

from ..schemas import AttitudeSchedule
from ..schemas.orbit import AllOrbits
from ..schemas.space.base import AllAttitudes
from .schemas import InstrumentEntry, InstrumentPayload
from .standard import Grade, InstrumentClass, Tier
from .system_schemas import (
    SatelliteSystem,
    Spacecraft,
    SystemEntry,
    _read_table,
    parse_system,
)
from .system_standard import SystemClass

_RESOURCES = pathlib.Path(__file__).parents[1] / "resources"
_ENTRIES = _RESOURCES / "instruments"
_RECORDS = _RESOURCES / "validation"
_SYSTEMS = _RESOURCES / "systems"
_SYSTEM_RECORDS = _RECORDS / "systems"


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
    *,
    platform: str | None = None,
    grade: Grade | str | None = None,
) -> list[InstrumentEntry]:
    """
    Lists catalog entries, optionally of one class or validation tier, that list
    a platform, or that have a configuration for a use (a grade). With a grade,
    the tier is the entry's for that use (see `InstrumentEntry.get_tier`): the
    lowest tier of its configurations that serve the grade, so that, for
    example, `list_entries(tier="validated", grade="coverage")` gives the
    entries validated for coverage, whatever their tier at other grades.

    Args:
        instrument_class (InstrumentClass | str | None): The class, if any.
        tier (Tier | str | None): The validation tier, if any: the entry's, or,
            with a grade, the entry's for the grade.
        platform (str | None): The identifier of a spacecraft of the catalog's
            systems (see `get_spacecraft`) that the entries list, if any.
        grade (Grade | str | None): The grade of a use that the entries have a
            configuration for, if any.

    Returns:
        list[InstrumentEntry]: Copies of the entries, sorted by identifier.
    """
    wanted_class = (
        None if instrument_class is None else InstrumentClass(instrument_class)
    )
    wanted_tier = None if tier is None else Tier(tier)
    wanted_grade = None if grade is None else Grade(grade)
    selected = []
    for entry in map(_load_entry, get_entry_ids()):
        if wanted_class is not None and entry.instrument_class != wanted_class:
            continue
        if platform is not None and platform not in entry.platforms:
            continue
        if wanted_grade is not None or wanted_tier is not None:
            entry_tier = entry.get_tier(wanted_grade)
            if entry_tier is None or (
                wanted_tier is not None and entry_tier != wanted_tier
            ):
                continue
        selected.append(entry.model_copy(deep=True))
    return selected


def build_payload(
    references: str | Sequence[str],
    *,
    orbit: AllOrbits,
    grade: Grade | str | None = None,
    time_step: float | None = None,
    attitude: AllAttitudes | AttitudeSchedule | None = None,
    platform: str | None = None,
    start: datetime | None = None,
    end: datetime | None = None,
) -> InstrumentPayload:
    """
    Builds the instruments of one or more catalog configurations for a
    satellite that carries them together, with the attitude they share, to be
    analyzed together. Each reference is the identifier of an entry (or of one
    of its aliases), optionally with a configuration after a colon (such as
    `"atms:swath_noaa_20"`, or `"atms:swath"` for the variant tested on the
    platform), whose configuration is otherwise selected by grade and platform
    (see `InstrumentEntry.select_configuration`). A configuration of
    `alternatives` gives only the instrument it analyzes by default. Analyze a
    satellite of the payload (`InstrumentPayload.to_satellite`) with its
    `instrument_index`.

    Args:
        references (str | Sequence[str]): The references.
        orbit (AllOrbits): The orbit.
        grade (Grade | str | None): The grade of the use, if any, which every
            configuration must serve.
        time_step (float | None): The analysis time step (seconds), for
            configurations that tailor a field of view to it.
        attitude (AllAttitudes | AttitudeSchedule | None): An attitude to use
            instead of the validated one, in the velocity frame of every
            configuration.
        platform (str | None): The identifier of the platform modeled, if any.
        start (datetime.datetime | None): The start of the analysis, for
            configurations that schedule their tilt by day.
        end (datetime.datetime | None): The end of the analysis, for
            configurations that schedule their tilt by day.

    Returns:
        InstrumentPayload: The instruments, with an `instrument_index` of 0 for
        one instrument and `None` (analyzed together) for several.

    Raises:
        KeyError: If an entry is unknown.
        ValueError: If there are no references, a configuration is unknown or
            does not serve the grade, or the configurations have different
            attitudes.
    """
    references = [references] if isinstance(references, str) else list(references)
    if not references:
        raise ValueError("a payload needs at least one catalog reference")
    payloads = []
    for reference in references:
        entry_id, _, configuration = reference.partition(":")
        payloads.append(
            get_entry(entry_id).build(
                configuration=configuration or None,
                grade=grade,
                orbit=orbit,
                time_step=time_step,
                attitude=attitude,
                platform=platform,
                start=start,
                end=end,
            )
        )
    if any(payload.sounding is not None for payload in payloads):
        if len(payloads) > 1:
            raise ValueError(
                f"a sounding is analyzed alone, not with other references: {references}"
            )
        return payloads[0]
    attitude = payloads[0].attitude
    for payload in payloads[1:]:
        if payload.attitude != attitude:
            raise ValueError(
                "the configurations have different attitudes of the body, so they "
                f"cannot be carried together: {references} (give an attitude to "
                "replace theirs)"
            )
    instruments = [
        instrument
        for payload in payloads
        for instrument in (
            payload.instruments
            if payload.instrument_index is None
            else payload.instruments[payload.instrument_index :][:1]
        )
    ]
    return InstrumentPayload(
        instruments=instruments,
        attitude=attitude,
        instrument_index=0 if len(instruments) == 1 else None,
        references=tuple(r for payload in payloads for r in payload.references),
        transferred=any(payload.transferred for payload in payloads),
        limitations=tuple(
            dict.fromkeys(text for payload in payloads for text in payload.limitations)
        ),
    )


def get_system_ids() -> list[str]:
    """
    Gets the identifiers of the catalog's systems.

    Returns:
        list[str]: The sorted identifiers.
    """
    return sorted(path.stem for path in _SYSTEMS.glob("*.yml"))


@functools.lru_cache(maxsize=None)
def _load_system(system_id: str) -> SystemEntry:
    """
    Loads a system and its validation record, if any, from files.

    Args:
        system_id (str): The identifier of the system.

    Returns:
        SystemEntry: The system, which callers must not modify.
    """
    with open(_SYSTEMS / f"{system_id}.yml", "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if data.get("id") != system_id:
        raise ValueError(f"the file of system {system_id!r} has id {data.get('id')!r}")
    record_path = _SYSTEM_RECORDS / f"{system_id}.json"
    if record_path.exists():
        with open(record_path, "r", encoding="utf-8") as f:
            data["record"] = json.load(f)
    return parse_system(data)


def get_system(system_id: str) -> SystemEntry:
    """
    Gets a catalog system by its identifier or one of its aliases (ignoring
    case).

    Args:
        system_id (str): The identifier or alias.

    Returns:
        SystemEntry: A copy of the system.

    Raises:
        KeyError: If there is no such system.
    """
    wanted = system_id.casefold()
    for candidate in get_system_ids():
        if wanted == candidate:
            return _load_system(candidate).model_copy(deep=True)
    for candidate in get_system_ids():
        entry = _load_system(candidate)
        if wanted in (name.casefold() for name in (entry.name, *entry.aliases)):
            return entry.model_copy(deep=True)
    raise KeyError(f"no catalog system {system_id!r}: {get_system_ids()}")


def list_systems(
    system_class: SystemClass | str | None = None,
    tier: Tier | str | None = None,
    *,
    satellite: bool | None = None,
) -> list[SystemEntry]:
    """
    Lists catalog systems, optionally of one class or validation tier, or only
    satellite systems (or only networks).

    Args:
        system_class (SystemClass | str | None): The class, if any.
        tier (Tier | str | None): The validation tier, if any.
        satellite (bool | None): `True` for satellite systems only, `False` for
            networks only, or `None` (the default) for both.

    Returns:
        list[SystemEntry]: Copies of the systems, sorted by identifier.
    """
    wanted_class = None if system_class is None else SystemClass(system_class)
    wanted_tier = None if tier is None else Tier(tier)
    return [
        entry
        for entry in map(get_system, get_system_ids())
        if (wanted_class is None or entry.system_class == wanted_class)
        and (satellite is None or entry.system_class.is_satellite() == satellite)
        and (wanted_tier is None or entry.validation.tier == wanted_tier)
    ]


@functools.lru_cache(maxsize=None)
def _get_spacecraft_systems() -> dict[str, str]:
    """
    Gets the system of each spacecraft, by the spacecraft's identifier.

    Returns:
        dict[str, str]: The identifiers of the systems.

    Raises:
        ValueError: If a spacecraft belongs to several systems.
    """
    systems: dict[str, str] = {}
    for system_id in get_system_ids():
        system = _load_system(system_id)
        if not isinstance(system, SatelliteSystem):
            continue
        for spacecraft in system.spacecraft:
            if spacecraft.id in systems:
                raise ValueError(
                    f"spacecraft {spacecraft.id!r} belongs to systems "
                    f"{systems[spacecraft.id]!r} and {system_id!r}"
                )
            systems[spacecraft.id] = system_id
    return systems


def _find_spacecraft(spacecraft: str) -> str:
    """
    Gets the identifier of a spacecraft given by its identifier or name
    (ignoring case).
    """
    systems = _get_spacecraft_systems()
    if spacecraft in systems:
        return spacecraft
    wanted = spacecraft.casefold()
    for spacecraft_id, system_id in systems.items():
        found = _load_system(system_id).get_spacecraft(spacecraft_id)
        if wanted in (found.id.casefold(), found.name.casefold()):
            return spacecraft_id
    raise KeyError(f"no catalog spacecraft {spacecraft!r}")


def get_spacecraft(spacecraft: str) -> Spacecraft:
    """
    Gets a spacecraft of a catalog system by its identifier (which instrument
    entries list) or name, ignoring case.

    Args:
        spacecraft (str): The identifier or name of the spacecraft (such as
            `noaa-20` or `NOAA-20`).

    Returns:
        Spacecraft: The spacecraft.

    Raises:
        KeyError: If no system has such a spacecraft.
    """
    spacecraft_id = _find_spacecraft(spacecraft)
    return _load_system(_get_spacecraft_systems()[spacecraft_id]).get_spacecraft(
        spacecraft_id
    )


def get_spacecraft_system(spacecraft: str) -> SystemEntry:
    """
    Gets the system of a spacecraft.

    Args:
        spacecraft (str): The identifier or name of the spacecraft.

    Returns:
        SystemEntry: A copy of its system.

    Raises:
        KeyError: If no system has such a spacecraft.
    """
    return get_system(_get_spacecraft_systems()[_find_spacecraft(spacecraft)])


def clear_cache() -> None:
    """
    Clears the catalog's cached entries, systems, records, and tables of
    stations, so that files rewritten during a session (as by the validation
    scripts) are read again.
    """
    _load_entry.cache_clear()
    _load_system.cache_clear()
    _get_spacecraft_systems.cache_clear()
    _read_table.cache_clear()
