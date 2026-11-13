"""
Object schemas for instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from collections.abc import Mapping
from typing import Annotated, Any

from pydantic import Discriminator, Tag

from .simple import Instrument
from .pointed import PointedInstrument
from .conical import ConicalInstrument

_BASE_FIELDS = set(Instrument.model_fields)
_POINTED_FIELDS = set(PointedInstrument.model_fields) - _BASE_FIELDS
_CONICAL_FIELDS = set(ConicalInstrument.model_fields) - _BASE_FIELDS - _POINTED_FIELDS


def _get_instrument_type(value: Any) -> str | None:
    """
    Gets the discriminator of an instrument: its `type` or, for data that omits
    it (as written before the field existed), the type its fields imply.

    Args:
        value (Any): The instrument object or data to validate or serialize.

    Returns:
        str | None: The type, or `None` if it cannot be determined.
    """
    if not isinstance(value, Mapping):
        return getattr(value, "type", None)
    if isinstance(value.get("type"), str):
        return value["type"]
    if _CONICAL_FIELDS.intersection(value):
        return "conical"
    if _POINTED_FIELDS.intersection(value):
        return "pointed"
    return "instrument"


AllInstruments = Annotated[
    Annotated[Instrument, Tag("instrument")]
    | Annotated[PointedInstrument, Tag("pointed")]
    | Annotated[ConicalInstrument, Tag("conical")],
    Discriminator(_get_instrument_type),
]
"""Instruments, identified by their `type` or, if omitted, by their fields."""

__all__ = ["AllInstruments", "ConicalInstrument", "Instrument", "PointedInstrument"]
