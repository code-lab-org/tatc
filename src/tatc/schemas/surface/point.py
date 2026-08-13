"""
Base classes for surface objects.

@author Paul T. Grogan <paul.grogan@asu.edu>
"""

from pydantic import BaseModel, Field, NonNegativeInt


class Point(BaseModel):
    """
    Surface point in the WGS 84 coordinate system.
    """

    id: NonNegativeInt = Field(default=0, description="Unique point identifier.")
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
        default=0,
        description="Elevation (meters) above datum in the WGS 84 coordinate system.",
    )
