"""
Object schemas for instruments.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from collections.abc import Mapping
from typing import Annotated, Any

from pydantic import Discriminator, Tag

from ...utils.projection import VelocityFrame

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


def _check_platform_instruments(
    instruments: list[Instrument],
    platform: str,
    direction: str,
    allow_conical: bool = True,
) -> None:
    """
    Checks that the views of instruments are defined on a platform other than
    a satellite (an aircraft or a surface platform), raising a `ValueError`
    otherwise: without roll or pitch angle profiles (which vary with an
    orbit's argument of latitude) or an inertial velocity frame (aligned with
    an orbit plane), and, unless `allow_conical`, not a `ConicalInstrument`.

    Args:
        instruments (list[Instrument]): The instruments.
        platform (str): The platform, as named in an error (for example,
            "an aircraft").
        direction (str): What defines the platform's along-track direction,
            as named in an error (for example, "track").
        allow_conical (bool): `True` to accept a `ConicalInstrument`.
    """
    for instrument in instruments:
        if isinstance(instrument, ConicalInstrument) and not allow_conical:
            raise ValueError(
                f"instrument {instrument.name!r}: a ConicalInstrument is not "
                f"supported on {platform}"
            )
        if isinstance(instrument, PointedInstrument):
            if (
                instrument.roll_angle_profile is not None
                or instrument.pitch_angle_profile is not None
            ):
                raise ValueError(
                    f"instrument {instrument.name!r}: roll and pitch angle "
                    f"profiles are not supported on {platform}"
                )
            if instrument.velocity_frame != VelocityFrame.EARTH_FIXED:
                raise ValueError(
                    f"instrument {instrument.name!r}: the {direction} of "
                    f"{platform} defines the along-track direction, which "
                    "requires velocity_frame='earth_fixed'"
                )


__all__ = ["AllInstruments", "ConicalInstrument", "Instrument", "PointedInstrument"]
