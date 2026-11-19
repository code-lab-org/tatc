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
    attitude = payloads[0].attitude
    for payload in payloads[1:]:
        if payload.attitude != attitude:
            raise ValueError(
                "the configurations have different attitudes, so they cannot be "
                f"carried together: {references}"
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
