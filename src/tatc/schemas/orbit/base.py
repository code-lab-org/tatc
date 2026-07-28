"""
Object schemas for orbits.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime, timezone

from pydantic import BaseModel, Field

from ... import utils


class OrbitBase(BaseModel):
    """
    Base class for orbits.
    """

    true_anomaly: float = Field(0, description="True anomaly (degrees).", ge=0, lt=360)
    epoch: datetime = Field(
        datetime.now(tz=timezone.utc),
        description="Timestamp (epoch) of the initial orbital state.",
    )

    def get_mean_anomaly(self) -> float:
        """
        Gets the mean anomaly (decimal degrees).

        Returns:
            float: the mean anomaly
        """
        return utils.orbital.true_anomaly_to_mean_anomaly(self.true_anomaly)
