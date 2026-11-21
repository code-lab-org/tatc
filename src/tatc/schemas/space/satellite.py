"""
Object schema for satellites.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime
from typing import Literal

from pydantic import Field
from skyfield.positionlib import Geocentric
from skyfield.timelib import Time

from ...utils.time import _to_time
from ..orbit import AllOrbits
from .base import SpaceSystem


class Satellite(SpaceSystem):
    """
    Single satellite.
    """

    type: Literal["satellite"] = Field(
        default="satellite", description="Space system type discriminator."
    )
    orbit: AllOrbits = Field(..., description="Orbit specification.")

    def get_track(self, times: datetime | list[datetime] | Time) -> Geocentric:
        """
        Gets the track of the satellite at times (see
        `GeneralPerturbationsOrbit.get_orbit_track_at_time`), from which its
        instruments project their views, as `tatc.schemas.Aircraft.get_track`
        and `tatc.schemas.SurfacePlatform.get_track` get those of an aircraft
        and a fixed surface platform.

        Args:
            times (datetime.datetime | list[datetime.datetime] | skyfield.timelib.Time):
                The time(s).

        Returns:
            skyfield.positionlib.Geocentric: The position/velocity.
        """
        return self.orbit.to_gp_orbit().get_orbit_track_at_time(_to_time(times))
