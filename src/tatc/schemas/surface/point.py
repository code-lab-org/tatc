"""
Object schemas for sampling points.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from pydantic import BaseModel, Field, NonNegativeInt


class Point(BaseModel):
    """
    Geodetic point in the WGS 84 coordinate system.
    """

    id: NonNegativeInt = Field(..., description="Unique point identifier.")
    latitude: float = Field(
        ...,
        description="Latitude (decimal degrees) in the WGS 84 coordinate system.",
        ge=-90,
        le=90,
        examples=[40.74259],
    )
    longitude: float = Field(
        ...,
        description="Longitude (decimal degrees) in the WGS 84 coordinate system.",
        ge=-180,
        le=180,
        examples=[-74.02686],
    )
    elevation: float = Field(
        0,
        description="Elevation (meters) above datum in the WGS 84 coordinate system.",
    )
