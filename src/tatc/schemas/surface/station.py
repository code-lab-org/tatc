"""
Object schemas for ground stations.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import timedelta

from pydantic import Field

from .point import Point


class GroundStation(Point):
    """
    Ground station in the WGS 84 coordinate system.
    """

    name: str = Field(..., description="Ground station name", examples=["station 1"])
    min_elevation_angle: float = Field(
        default=0,
        description="The minimum elevation angle (decimal degrees) required "
        + "for satellite communication.",
        ge=0,
        le=90,
    )
    min_access_time: timedelta = Field(
        default=timedelta(0),
        description="Minimum access (integration) time required for satellite communication.",
        examples=[timedelta(seconds=10)],
    )
