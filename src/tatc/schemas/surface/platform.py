"""
Object schemas for fixed surface platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime

import numpy as np
from pydantic import Field, model_validator
from shapely.geometry import MultiPolygon, Polygon
from skyfield.positionlib import Geocentric

from ...utils.projection import (
    VelocityFrame,
    _compute_horizon_region,
    compute_static_track,
)
from ...utils.time import _to_time
from ..instrument import (
    AllInstruments,
    ConicalInstrument,
    Instrument,
    PointedInstrument,
)
from .point import Point
from .terrain import TerrainMask


class SurfacePlatform(Point):
    """
    Platform fixed above the Earth's surface (e.g. a tower or a tethered
    balloon) in the WGS 84 coordinate system, hosting one or more
    instruments. Unlike a satellite, the platform does not move in the
    Earth-fixed frame, so the views of its instruments, and the footprints
    they project, are static.

    Instruments are oriented as on a satellite whose direction of motion is
    the platform's `heading`: a view at zero roll and pitch looks straight
    down (nadir), a positive roll angle looks to the left of the heading,
    a positive pitch angle looks toward the heading, and a roll angle of 180
    degrees looks straight up (zenith). A generic `Instrument` observes a
    cone about nadir; a `PointedInstrument` observes a view rotated from
    nadir by fixed roll and pitch angles, for example looking up at a
    target above the instruments (such as clouds).
    """

    name: str = Field(..., description="Platform name.", examples=["tower 1"])
    elevation: float = Field(
        default=0,
        description="Elevation (meters) of the ground below the platform above "
        + "the WGS 84 ellipsoid (which differs from a height above mean sea "
        + "level by the local geoid undulation).",
    )
    height: float = Field(
        ...,
        description="Height (meters) of the instruments above the ground, for "
        + "example the height of a tower or the altitude of a tethered balloon "
        + "above its tether point.",
        gt=0,
        examples=[30],
    )
    heading: float = Field(
        default=0,
        description="Azimuth (decimal degrees, clockwise from north) of the "
        + "platform's forward (along-track) direction, which orients the views "
        + "of its instruments.",
        ge=0,
        lt=360,
    )
    instruments: list[AllInstruments] = Field(
        default=[Instrument()],
        min_length=1,
        description="List of hosted instruments (`Instrument` or "
        + "`PointedInstrument`).",
    )
    terrain_mask: TerrainMask | None = Field(
        default=None,
        description="Optional azimuthal terrain blockage mask as seen from the "
        + "instruments (e.g. from "
        + "`tatc.preprocess.compute_terrain_mask_for_platform`), which, with "
        + "the horizon of the ground, blocks views toward a target above the "
        + "instruments (see `compute_visible_region`). It does not apply to "
        + "targets below the instruments.",
    )

    @model_validator(mode="after")
    def check_instruments(self) -> SurfacePlatform:
        """
        Checks that each instrument's view is defined for a fixed platform:
        not a `ConicalInstrument` (whose footprint is swept by motion along
        track), and without roll or pitch angle profiles (which vary with an
        orbit's argument of latitude) or an inertial velocity frame (aligned
        with an orbit plane).
        """
        for instrument in self.instruments:
            if isinstance(instrument, ConicalInstrument):
                raise ValueError(
                    f"instrument {instrument.name!r}: a ConicalInstrument is not "
                    "supported on a surface platform"
                )
            if isinstance(instrument, PointedInstrument):
                if (
                    instrument.roll_angle_profile is not None
                    or instrument.pitch_angle_profile is not None
                ):
                    raise ValueError(
                        f"instrument {instrument.name!r}: roll and pitch angle "
                        "profiles are not supported on a surface platform"
                    )
                if instrument.velocity_frame != VelocityFrame.EARTH_FIXED:
                    raise ValueError(
                        f"instrument {instrument.name!r}: a surface platform's "
                        "heading defines the along-track direction, which "
                        "requires velocity_frame='earth_fixed'"
                    )
        return self

    def get_altitude(self) -> float:
        """
        Gets the elevation of the instruments: the ground elevation plus the
        platform height.

        Returns:
            float: the elevation (meters) above the WGS 84 ellipsoid
        """
        return self.elevation + self.height

    def compute_visible_region(
        self, elevation: float, number_azimuths: int = 360
    ) -> Polygon | MultiPolygon:
        """
        Computes the region at a target elevation above the instruments that
        is visible over the horizon: bounded at each azimuth by the lowest
        line of sight that clears both the terrain mask (if any) and the
        horizon of the ground, projected to the target elevation. Footprints
        toward a target above the instruments are clipped to this region
        (see `tatc.analysis.collect_platform_track`).

        Args:
            elevation (float): The target elevation (meters) above the WGS 84
                ellipsoid, above the instruments.
            number_azimuths (int): The number of evenly spaced azimuths at
                which to sample the region (in addition to those of the
                terrain mask).

        Returns:
            shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The
                visible region.
        """
        if elevation <= self.get_altitude():
            raise ValueError(
                f"target elevation ({elevation} m) must be above the "
                f"instruments of platform {self.name!r} ({self.get_altitude()} m)"
            )
        azimuths = np.linspace(0, 360, number_azimuths, endpoint=False)
        if self.terrain_mask is None:
            # below the horizon of the ground at every azimuth
            angles = np.full(len(azimuths), -89.0)
        else:
            azimuths = np.union1d(azimuths, self.terrain_mask.azimuth)
            angles = self.terrain_mask.get_min_elevation_angles(azimuths)
        return _compute_horizon_region(
            self.longitude,
            self.latitude,
            self.get_altitude(),
            elevation,
            self.elevation,
            azimuths,
            angles,
        )

    def get_track(self, times: datetime | list[datetime] | None = None) -> Geocentric:
        """
        Gets the static track of the platform's instruments (see
        `tatc.utils.compute_static_track`), from which they project their
        views as from a satellite's orbit track.

        Args:
            times (datetime.datetime | list[datetime.datetime] | None): The
                time(s). The views and footprints do not depend on time; times
                are needed only to evaluate time-varying conditions (such as
                the sunlit conditions of `Instrument.is_valid_observation`).
                Defaults to 2020-01-01 00:00:00 UTC.

        Returns:
            skyfield.positionlib.Geocentric: The position/velocity.
        """
        return compute_static_track(
            self.longitude,
            self.latitude,
            self.get_altitude(),
            self.heading,
            None if times is None else _to_time(times),
        )
