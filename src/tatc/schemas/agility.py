"""
Object schemas for the agility of spacecraft: how quickly they turn to point
their instruments at targets.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, ConfigDict, Field


class Agility(BaseModel):
    """
    The agility of a spacecraft: the time it takes to turn (slew) its body
    through an angle, about the axis between its attitudes (eigenaxis), and
    settle before it observes. With a maximum acceleration, a slew
    accelerates to the maximum rate, coasts, and decelerates (a trapezoidal
    profile of rate), or, if too short to reach the maximum rate, accelerates
    and decelerates (a triangular profile): it takes `2 * sqrt(angle /
    acceleration)` for angles up to `rate ** 2 / acceleration`, and otherwise
    `angle / rate + rate / acceleration`. Without one, a slew turns at the
    maximum rate throughout (`angle / rate`). The settle time is added to
    both.

    Agile small satellites of Earth observation turn at about 1 to 4 degrees
    per second; large agile imagers, with control moment gyroscopes, at up
    to 4 or more, with accelerations of about 0.5 to 2 degrees per second
    squared and settle times of a few seconds.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    max_rate: float = Field(
        ..., gt=0, description="The maximum rate (degrees per second) of a slew."
    )
    max_acceleration: float | None = Field(
        None,
        gt=0,
        description="The maximum acceleration (degrees per second squared) of a "
        + "slew, or `None` to turn at the maximum rate throughout.",
    )
    settle_time: float = Field(
        0.0,
        ge=0,
        description="The time (seconds) to settle after a slew before observing.",
    )

    def get_slew_time(self, angle: float | npt.ArrayLike) -> float | npt.NDArray:
        """
        Gets the time to slew through an angle and settle (see `Agility`).
        No slew (an angle of zero) takes no time, not even to settle.

        Args:
            angle (float | numpy.typing.ArrayLike): The angle(s) (degrees).

        Returns:
            float | numpy.typing.NDArray: The time(s) (seconds).

        Raises:
            ValueError: If an angle is negative.
        """
        angle = np.asarray(angle, dtype=float)
        if np.any(angle < 0):
            raise ValueError("slew angles must not be negative")
        rate = self.max_rate
        if self.max_acceleration is None:
            turn = angle / rate
        else:
            acceleration = self.max_acceleration
            turn = np.where(
                angle <= rate**2 / acceleration,
                2 * np.sqrt(angle / acceleration),
                angle / rate + rate / acceleration,
            )
        time = np.where(angle > 0, turn + self.settle_time, 0.0)
        return float(time) if np.ndim(time) == 0 else time
