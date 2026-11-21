"""
Methods to collect the static tracks of instruments on fixed surface platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime

import geopandas as gpd
import pandas as pd
import shapely
from shapely.geometry import MultiPolygon, Polygon

from .._warnings import warn
from ..schemas import SurfacePlatform
from ._tracks import (
    Mask,
    _assemble,
    _collect_tracks,
    _FootprintRequest,
)
from .check import SURFACE_PLATFORM, _dissolve_footprints


def _get_target_elevation(
    platform: SurfacePlatform, request: _FootprintRequest
) -> float:
    """
    Gets the elevation at which to project a platform's footprints (see
    `_FootprintRequest.get_elevation`: by default, the platform's ground
    elevation), which must differ from the elevation of its instruments
    (below them, looking down, or above them, looking up).

    Args:
        platform (SurfacePlatform): The platform.
        request (_FootprintRequest): The request.

    Returns:
        float: the target elevation (meters) above the WGS 84 ellipsoid
    """
    target = request.get_elevation(platform)
    if target == platform.get_altitude():
        raise ValueError(
            f"target elevation ({target} m) must differ from the elevation of "
            f"the instruments of platform {platform.name!r} "
            f"({platform.get_altitude()} m)"
        )
    return target


def _clip_to_region(
    footprint: Polygon | MultiPolygon, region: Polygon | MultiPolygon, elevation: float
) -> Polygon | MultiPolygon:
    """
    Clips a footprint to a region (such as the region visible over a terrain
    mask), keeping only its polygonal parts.

    Args:
        footprint (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The footprint.
        region (shapely.geometry.Polygon | shapely.geometry.MultiPolygon): The region.
        elevation (float): The elevation (meters) of the footprint.

    Returns:
        shapely.geometry.Polygon | shapely.geometry.MultiPolygon: The clipped footprint.
    """
    clipped = shapely.intersection(shapely.force_2d(footprint), region)
    parts = [
        part
        for part in shapely.get_parts(shapely.get_parts(clipped))
        if isinstance(part, Polygon) and not part.is_empty
    ]
    if len(parts) == 0:
        return Polygon()
    return shapely.force_3d(
        parts[0] if len(parts) == 1 else MultiPolygon(parts), elevation
    )


def collect_surface_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    times: list[datetime] | None = None,
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Mask | None = None,
    platform_altaz: bool = False,
    solar_altaz: bool = False,
) -> gpd.GeoDataFrame:
    """
    Collect the static footprint of each instrument on fixed surface
    platforms (e.g. towers or tethered balloons), projected to a target
    elevation by computing the exact ray/WGS 84 ellipsoid intersection, as
    for a satellite's ground track (see `collect_ground_track`). Rays that
    miss a target below the instruments (looking above the horizon) are
    projected to the horizon. A target above the instruments (as for clouds
    observed by an instrument looking up, e.g. with a roll angle of 180
    degrees) is observed only above the horizon: views are clipped to the
    lines of sight tangent to the ground, which blocks the rest (a view
    entirely below the horizon has an empty footprint), and to the
    platform's terrain mask, if any (see
    `SurfacePlatform.compute_visible_region`). An empty footprint (as of a
    view entirely blocked by the ground or terrain) is not a valid
    observation (`valid_obs`), as for satellites and aircraft.

    A platform does not move in the Earth-fixed frame, so its footprints do
    not change with time. Without `times`, this returns one row per platform
    and instrument, without a time, and an instrument's time-varying
    requirements (such as a sunlit target, see
    `Instrument.is_valid_observation`) are not evaluated. With `times`, it
    returns one row per time, platform, and instrument, as
    `collect_ground_track` does for satellites, with the validity of each
    observation evaluated at its time (at the footprint center) and, if
    `solar_altaz`, the solar angles.

    Args:
        platforms (SurfacePlatform | list[SurfacePlatform]): The platform(s).
        times (list[datetime.datetime] | None): The times to sample, or
                `None` for a static footprint of each instrument.
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` for every instrument of each
                platform.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid of the observed target for which to compute
                footprints, below (looking down) or above (looking up) each
                platform's instruments. Defaults to each platform's ground
                elevation.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.
        platform_altaz (bool): `True` to include the altitude/azimuth
                angles (`platform_alt` and `platform_az`) of the instruments
                seen from the center of their view at the target elevation
                (as for a satellite),
                `NaN` for an empty footprint. Below the instruments, the
                altitude is the elevation angle of the view; above them, it
                is negative.
        solar_altaz (bool): `True` to include the solar altitude/azimuth
                angles (`solar_alt` and `solar_az`) and the glint angle
                (`glint_angle`, see `tatc.utils.compute_glint_angle`) at the
                footprint center, `NaN` for an empty footprint. Requires
                `times`.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected surface track
            results, with columns `surface_platform`, `instrument`,
            `valid_obs`, and `geometry` (whose z-coordinate is the target
            elevation), and the angles requested: without `times`, in the order of the platforms and
            their instruments; with `times`, also `time`, for a single
            platform and instrument (an integer `instrument_index`) sorted
            by time, otherwise those of each platform and instrument,
            concatenated and sorted by time.
    """
    return _collect_tracks(
        platforms,
        SURFACE_PLATFORM,
        _FootprintRequest.create(
            times, instrument_index, elevation, mask, platform_altaz, solar_altaz
        ),
        _collect_surface_platform,
    )


def _collect_surface_platform(
    platform: SurfacePlatform, instrument_index: int, request: _FootprintRequest
) -> gpd.GeoDataFrame | None:
    """
    Collect a surface platform instrument's footprints (see
    `collect_surface_track`) or, if the request is of pixels, its pixels, or
    `None` if there are none.

    Args:
        platform (SurfacePlatform): The platform.
        instrument_index (int): The index of the instrument.
        request (_FootprintRequest): The request.

    Returns:
        geopandas.GeoDataFrame | None: the footprints or pixels
    """
    times = request.times
    if request.solar_altaz and times is None:
        raise ValueError(
            "solar_altaz requires times: the solar angles of a fixed "
            "platform's footprint change with time"
        )
    target = _get_target_elevation(platform, request)
    # the static track, from which the footprints are projected
    static = platform.get_track()
    # the track at each time (or at the time of the static track, without
    # times)
    track = platform.get_track(
        [static.t.utc_datetime()] if times is None else times  # type: ignore
    )
    timestamps = None if times is None else pd.to_datetime(times, utc=True)
    # the region visible over the terrain mask, toward a target above
    region = (
        platform.compute_visible_region(target)
        if platform.terrain_mask is not None and target > platform.get_altitude()
        else None
    )
    instrument = platform.get_instrument(instrument_index)
    geometries = None
    if not request.pixels:
        footprint = instrument.compute_footprint(
            static, elevation=target, ground_elevation=platform.elevation
        )[0]
        if region is not None:
            footprint = _clip_to_region(footprint, region, target)
        geometries = [footprint] * (1 if times is None else len(times))
    return _assemble(
        platform,
        instrument,
        timestamps,
        track,
        target,
        request,
        geometries,
        ground_elevation=platform.elevation,
        region=region,
    )


def compute_surface_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    times: list[datetime] | None = None,
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Mask | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated (dissolved) coverage of instruments on fixed
    surface platforms into a single composite geometry (e.g. a network
    coverage map): unlike `collect_surface_track` (one footprint polygon per
    platform and instrument, and time), this dissolves the footprints of
    valid observations into one geometry (see
    `tatc.analysis.footprints.dissolve_footprints`).

    Args:
        platforms (SurfacePlatform | list[SurfacePlatform]): The platform(s).
        times (list[datetime.datetime] | None): The times to sample (see
                `collect_surface_track`), or `None` for static footprints
                (whose time-varying requirements are not evaluated).
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` for every instrument of each
                platform.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid of the observed target for which to compute
                footprints, below (looking down) or above (looking up) each
                platform's instruments. Defaults to each platform's ground
                elevation.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.

    Returns:
        geopandas.GeoDataFrame: The data frame of the aggregated surface
            track: one row with the observed area, or none if no observation
            is valid.
    """
    return _dissolve_footprints(
        collect_surface_track(
            platforms,
            times,
            instrument_index=instrument_index,
            elevation=elevation,
            mask=mask,
        )
    )


def _warn_renamed(old: str, new: str) -> None:
    """
    Warns that a function was renamed.

    Args:
        old (str): The deprecated name.
        new (str): The new name.
    """
    warn(
        f"{old} is deprecated and will be removed in a future release: use {new} "
        "(whose results name each platform in a surface_platform column)",
        DeprecationWarning,
    )


def collect_platform_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Mask | None = None,
    **kwargs,
) -> gpd.GeoDataFrame:
    """
    Deprecated: use `collect_surface_track`, of which this is an alias (with
    the former order of arguments, without times). Its results name each
    platform in a `surface_platform` column (formerly `platform`).
    """
    _warn_renamed("collect_platform_track", "collect_surface_track")
    return collect_surface_track(
        platforms,
        instrument_index=instrument_index,
        elevation=elevation,
        mask=mask,
        **kwargs,
    )


def compute_platform_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Mask | None = None,
) -> gpd.GeoDataFrame:
    """
    Deprecated: use `compute_surface_track`, of which this is an alias (with
    the former order of arguments, without times).
    """
    _warn_renamed("compute_platform_track", "compute_surface_track")
    return compute_surface_track(
        platforms, instrument_index=instrument_index, elevation=elevation, mask=mask
    )
