# -*- coding: utf-8 -*-
"""
Object schemas for sampling points.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

from pydantic import BaseModel, Field, NonNegativeInt
from datetime import timedelta


class Point(BaseModel):
    """
    Representation of a geodetic point in the WGS 84 coordinate system.
    """

    id: NonNegativeInt = Field(..., description="Unique point identifier.")
    latitude: float = Field(
        ..., description="Latitude (decimal degrees).", ge=-90, le=90, example=40.74259
    )
    longitude: float = Field(
        ...,
        description="Longitude (decimal degrees).",
        ge=-180,
        le=180,
        example=-74.02686,
    )


class GroundStation(BaseModel):
    """
    Representation of a ground station in the WGS 84 coordinate system.
    """

    name: str = Field(
        ..., description="Ground station name", example="station 1"
    )
    latitude: float = Field(
        ..., description="Latitude (decimal degrees).", ge=-90, le=90, example=40.74259
    )
    longitude: float = Field(
        ...,
        description="Longitude (decimal degrees).",
        ge=-180,
        le=180,
        example=-74.02686,
    )
    min_elevation_angle: float = Field(
        0,
        description="The minimum elevation angle (decimal degrees) required for satellite communication.",
        ge=0,
        le=90,
    )
    min_access_time: timedelta = Field(
        timedelta(0),
        description="Minimum access (integration) time required for satellite communication.",
        example=timedelta(seconds=10),
    )
