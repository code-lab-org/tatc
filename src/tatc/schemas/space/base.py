"""
Base object schemas for space systems.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Annotated

from pydantic import BaseModel, Field, model_validator

from ..attitude import Attitude, AttitudeSchedule, GlintAttitude, TargetAttitude
from ..instrument import AllInstruments, Instrument

AllAttitudes = Annotated[
    Attitude | TargetAttitude | GlintAttitude | AttitudeSchedule,
    Field(discriminator="type"),
]
"""Attitudes of a space system: about nadir, pointed at targets or the sun glint, or scheduled."""


class SpaceSystem(BaseModel):
    """
    Base class for space systems.
    """

    name: str = Field(
        ...,
        description="Space system name.",
        examples=["International Space Station"],
    )
    instruments: list[AllInstruments] = Field(
        default=[Instrument()],
        min_length=1,
        description="List of assigned instruments.",
    )
    attitude: AllAttitudes | None = Field(
        default=None,
        description="Attitude of the spacecraft's body (or a schedule of "
        + "attitudes), from which the views of its instruments are rotated, or "
        + "`None` for a body aligned with the nadir frame (each instrument's "
        + "`velocity_frame` then orients it).",
    )

    @model_validator(mode="after")
    def _check_velocity_frames(self) -> "SpaceSystem":
        """
        Checks that instruments do not set a velocity frame if the attitude
        sets the body's.
        """
        if self.attitude is not None:
            for instrument in self.instruments:
                if "velocity_frame" in instrument.model_fields_set:
                    raise ValueError(
                        f"Instrument {instrument.name!r} sets a velocity frame, "
                        + "but the attitude sets the body's (see Attitude.velocity_frame)."
                    )
        return self

    def get_instrument(self, index: int = 0) -> Instrument:
        """
        Gets an instrument carried by this space system, with its views
        rotated from the axes of the body in its attitude (see
        `Instrument.bind_attitude`).

        Args:
            index (int): The index of the instrument.

        Returns:
            Instrument: the instrument
        """
        instrument = self.instruments[index]
        if self.attitude is None:
            return instrument
        return instrument.bind_attitude(self.attitude)
