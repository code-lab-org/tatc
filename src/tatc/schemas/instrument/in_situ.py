"""
Object schemas for in-situ instruments of fixed surface platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from pydantic import Field
from shapely import MultiPolygon, Polygon
from skyfield.api import wgs84
from skyfield.positionlib import Geocentric
from skyfield.toposlib import GeographicPosition

from ...utils.geometry import (
    _WGS84_GEOD,
    project_polygon_to_elevation,
    split_polygon,
)
from .simple import Instrument

# the number of vertices of a footprint, by default
_NUMBER_POINTS = 72


class InSituInstrument(Instrument):
    """
    In-situ instrument of a fixed surface platform (such as the thermometer or
    rain gauge of a weather station), which samples the air or ground at the
    platform rather than viewing a scene. Its footprint is the area its samples
    represent: a geodesic disk of `radius` about the platform, at the target
    elevation, rather than the projection of a view (which, from the height of a
    weather station's sensors, could not reach beyond a horizon a few
    kilometers away). The radius is an assumption of representativeness (for
    example, of the scale over which a variable is spatially correlated), not a
    property of the sensor. Only a `SurfacePlatform` carries it.
    """

    type: Literal["in_situ"] = Field(  # type: ignore[assignment]
        default="in_situ", description="Instrument type discriminator."
    )
    radius: float = Field(
        ...,
        gt=0,
        description="Radius (meters) of the area the instrument's samples "
        + "represent, about its platform.",
        examples=[10000],
    )

    def get_swath_width(self, height: float) -> float:
        """
        Gets the width of the footprint: its diameter.

        Args:
            height (float): Not used.

        Returns:
            float: the width (meters)
        """
        del height
        return 2 * self.radius

    def compute_footprint(  # pylint: disable=too-many-locals
        self,
        orbit_track: Geocentric,
        number_points: int | None = None,
        elevation: float = 0,
        ground_elevation: float = 0,
    ) -> list[Polygon | MultiPolygon]:
        """
        Computes the footprint at each position of a track: a geodesic disk of
        the instrument's `radius` about the position's subpoint, at the
        elevation.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The platform's positions.
            number_points (int | None): The number of vertices of each disk
                (by default, 72).
            elevation (float): The elevation (meters) of the footprint.
            ground_elevation (float): Not used (the samples are not blocked).

        Returns:
            list[shapely.geometry.Polygon | shapely.geometry.MultiPolygon]: The footprints.
        """
        del ground_elevation
        number = _NUMBER_POINTS if number_points is None else number_points
        position = wgs84.geographic_position_of(orbit_track)
        latitudes = np.atleast_1d(position.latitude.degrees)
        longitudes = np.atleast_1d(position.longitude.degrees)
        azimuths = np.linspace(0, 360, number, endpoint=False)
        footprints = []
        for latitude, longitude in zip(latitudes, longitudes):
            ring_longitudes, ring_latitudes, _ = _WGS84_GEOD.fwd(
                np.full(number, longitude),
                np.full(number, latitude),
                azimuths,
                np.full(number, self.radius),
            )
            # unwrap longitudes about the center, so that a disk across the
            # anti-meridian is split rather than spanning the globe
            ring_longitudes = longitude + (
                (np.asarray(ring_longitudes) - longitude + 180) % 360 - 180
            )
            disk = split_polygon(Polygon(zip(ring_longitudes, ring_latitudes)))
            footprints.append(project_polygon_to_elevation(disk, elevation))
        return footprints

    def compute_footprint_center(
        self,
        orbit_track: Geocentric,
        elevation: float = 0,
        ground_elevation: float = 0,
    ) -> GeographicPosition:
        """
        Computes the center of the footprint: the platform's subpoint, at the
        elevation.

        Args:
            orbit_track (skyfield.positionlib.Geocentric): The platform's positions.
            elevation (float): The elevation (meters) of the footprint.
            ground_elevation (float): Not used.

        Returns:
            skyfield.toposlib.GeographicPosition: The center.
        """
        del ground_elevation
        position = wgs84.geographic_position_of(orbit_track)
        return wgs84.latlon(
            position.latitude.degrees,
            position.longitude.degrees,
            elevation_m=elevation,
        )
