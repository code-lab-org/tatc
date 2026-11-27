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
from .in_situ import InSituInstrument
from .sounder import LimbSounder, OccultationReceiver, ScanDirection, ScanStart

_BASE_FIELDS = set(Instrument.model_fields)
_POINTED_FIELDS = set(PointedInstrument.model_fields) - _BASE_FIELDS
_CONICAL_FIELDS = set(ConicalInstrument.model_fields) - _BASE_FIELDS - _POINTED_FIELDS
_IN_SITU_FIELDS = set(InSituInstrument.model_fields) - _BASE_FIELDS
# the fields of each kind of sounder that the other does not have (both have a
# sample_elevation)
_LIMB_FIELDS = set(LimbSounder.model_fields) - set(OccultationReceiver.model_fields)
_OCCULTATION_FIELDS = set(OccultationReceiver.model_fields) - set(
    LimbSounder.model_fields
)


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
    if _IN_SITU_FIELDS.intersection(value):
        return "in_situ"
    if _LIMB_FIELDS.intersection(value):
        return "limb"
    if _OCCULTATION_FIELDS.intersection(value):
        return "occultation"
    if _CONICAL_FIELDS.intersection(value):
        return "conical"
    if _POINTED_FIELDS.intersection(value):
        return "pointed"
    return "instrument"


AllInstruments = Annotated[
    Annotated[Instrument, Tag("instrument")]
    | Annotated[PointedInstrument, Tag("pointed")]
    | Annotated[ConicalInstrument, Tag("conical")]
    | Annotated[InSituInstrument, Tag("in_situ")]
    | Annotated[LimbSounder, Tag("limb")]
    | Annotated[OccultationReceiver, Tag("occultation")],
    Discriminator(_get_instrument_type),
]
"""Instruments, identified by their `type` or, if omitted, by their fields."""


def _check_platform_instruments(
    instruments: list[Instrument],
    platform: str,
    direction: str,
    allow_conical: bool = True,
    allow_in_situ: bool = False,
) -> None:
    """
    Checks that the views of instruments are defined on a platform other than
    a satellite (an aircraft or a surface platform), raising a `ValueError`
    otherwise: without roll or pitch angle profiles (which vary with an
    orbit's argument of latitude) or an inertial velocity frame (aligned with
    an orbit plane), not a sounder (a `LimbSounder` or `OccultationReceiver`,
    which only a satellite carries), unless `allow_conical`, not a
    `ConicalInstrument`, and, unless `allow_in_situ`, not an `InSituInstrument`.

    Args:
        instruments (list[Instrument]): The instruments.
        platform (str): The platform, as named in an error (for example,
            "an aircraft").
        direction (str): What defines the platform's along-track direction,
            as named in an error (for example, "track").
        allow_conical (bool): `True` to accept a `ConicalInstrument`.
        allow_in_situ (bool): `True` to accept an `InSituInstrument`.
    """
    for instrument in instruments:
        if isinstance(instrument, (LimbSounder, OccultationReceiver)):
            raise ValueError(
                f"instrument {instrument.name!r}: a {type(instrument).__name__} is "
                f"only supported on a satellite, not on {platform}"
            )
        if isinstance(instrument, InSituInstrument) and not allow_in_situ:
            raise ValueError(
                f"instrument {instrument.name!r}: an InSituInstrument is only "
                f"supported on a surface platform, not on {platform}"
            )
        if isinstance(instrument, ConicalInstrument) and not allow_conical:
            raise ValueError(
                f"instrument {instrument.name!r}: a ConicalInstrument is not "
                f"supported on {platform}"
            )
        if isinstance(instrument, PointedInstrument) and (
            instrument.roll_angle_profile is not None
            or instrument.pitch_angle_profile is not None
        ):
            raise ValueError(
                f"instrument {instrument.name!r}: roll and pitch angle "
                f"profiles are not supported on {platform}"
            )
        # of a PointedInstrument or a ConicalInstrument
        velocity_frame = getattr(instrument, "velocity_frame", None)
        if velocity_frame not in (None, VelocityFrame.EARTH_FIXED):
            raise ValueError(
                f"instrument {instrument.name!r}: the {direction} of "
                f"{platform} defines the along-track direction, which "
                "requires velocity_frame='earth_fixed'"
            )


__all__ = [
    "AllInstruments",
    "ConicalInstrument",
    "InSituInstrument",
    "Instrument",
    "LimbSounder",
    "OccultationReceiver",
    "PointedInstrument",
    "ScanDirection",
    "ScanStart",
]
