# -*- coding: utf-8 -*-
"""
Object schemas for sampling points.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import timedelta

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


class GroundStation(BaseModel):
    """
    Ground station in the WGS 84 coordinate system.
    """

    name: str = Field(..., description="Ground station name", example="station 1")
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
    min_elevation_angle: float = Field(
        0,
        description="The minimum elevation angle (decimal degrees) required "
        + "for satellite communication.",
        ge=0,
        le=90,
    )
    min_access_time: timedelta = Field(
        timedelta(0),
        description="Minimum access (integration) time required for satellite communication.",
        examples=[timedelta(seconds=10)],
    )
