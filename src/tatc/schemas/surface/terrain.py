"""
Object schemas for terrain masks of fixed surface objects.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, model_validator


class TerrainMask(BaseModel):
    """
    Azimuthal terrain blockage mask for a fixed surface object (a
    `RadarStation` or a `SurfacePlatform`): the minimum usable elevation
    angle at each azimuth, due to terrain blocking the line of sight at low
    grazing angles. Typically derived by ray-tracing a digital elevation
    model (DEM) outward from the antenna or instruments to find the highest
    terrain obstruction angle along each azimuthal direction (see
    `tatc.preprocess.compute_terrain_mask`).

    This represents one blocking angle per azimuth (the common case for a
    DEM-derived terrain horizon, where the nearest/highest ridge along a
    ray dominates): any elevation angle at or above the masked value is
    assumed clear, and any below it is blocked. It does not model
    multiple, separately-blocked elevation bands along the same azimuth
    (e.g. a nearby building blocking only a narrow low-elevation slice
    while a distant mountain blocks a separate, higher slice).
    """

    azimuth: list[float] = Field(
        ...,
        description="Azimuth samples (decimal degrees, clockwise from "
        + "north, in [0, 360)), strictly increasing, at which the minimum "
        + "elevation angle is specified. Values are linearly interpolated "
        + "between samples, wrapping around 0/360 degrees.",
        min_length=2,
    )
    min_elevation_angle: list[float] = Field(
        ...,
        description="Minimum usable elevation angle (decimal degrees), "
        + "due to terrain blockage, at each corresponding azimuth.",
        min_length=2,
    )

    @model_validator(mode="after")
    def validate_profile(self) -> TerrainMask:
        """
        Validates that `azimuth` and `min_elevation_angle` are the same
        length, and that `azimuth` values fall in [0, 360) and are
        strictly increasing.
        """
        if len(self.azimuth) != len(self.min_elevation_angle):
            raise ValueError(
                "azimuth and min_elevation_angle must have the same length"
            )
        if any(a < 0 or a >= 360 for a in self.azimuth):
            raise ValueError("azimuth values must be in [0, 360)")
        if any(b <= a for a, b in zip(self.azimuth, self.azimuth[1:])):
            raise ValueError("azimuth values must be strictly increasing")
        return self

    def get_min_elevation_angle(self, azimuth: float) -> float:
        """
        Interpolates (periodically, wrapping at 0/360 degrees) the minimum
        usable elevation angle at a specified azimuth.

        Args:
            azimuth (float): Azimuth (decimal degrees, clockwise from north).

        Returns:
            float: The interpolated minimum usable elevation angle (degrees).
        """
        return float(self.get_min_elevation_angles(azimuth))

    def get_min_elevation_angles(self, azimuths: npt.ArrayLike) -> npt.NDArray:
        """
        Interpolates (periodically, wrapping at 0/360 degrees) the minimum
        usable elevation angles at specified azimuths, vectorized across
        azimuths (see `get_min_elevation_angle`).

        Args:
            azimuths (numpy.typing.ArrayLike): Azimuths (decimal degrees,
                clockwise from north).

        Returns:
            numpy.typing.NDArray: The interpolated minimum usable elevation
            angles (degrees).
        """
        samples = np.concatenate(
            ([self.azimuth[-1] - 360], self.azimuth, [self.azimuth[0] + 360])
        )
        angles = np.concatenate(
            (
                [self.min_elevation_angle[-1]],
                self.min_elevation_angle,
                [self.min_elevation_angle[0]],
            )
        )
        return np.interp(np.asarray(azimuths) % 360, samples, angles)
