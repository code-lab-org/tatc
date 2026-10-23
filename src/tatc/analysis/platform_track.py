"""
Methods to collect the static tracks of instruments on fixed surface platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon

from ..schemas import SurfacePlatform
from .check import _get_instrument_indices
from .ground_track import _get_mask_geometry


def _get_empty_platform_track() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for platform track results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "platform": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "elevation": pd.Series([], dtype="float"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def _get_target_elevation(platform: SurfacePlatform, elevation: float | None) -> float:
    """
    Gets the elevation at which to project a platform's footprints: the
    given elevation, or the platform's ground elevation if None, which must
    be below its instruments.

    Args:
        platform (SurfacePlatform): The platform.
        elevation (float | None): The target elevation (meters) above the
            WGS 84 ellipsoid, or None for the platform's ground elevation.

    Returns:
        float: the target elevation (meters).
    """
    target = platform.elevation if elevation is None else elevation
    if target >= platform.get_altitude():
        raise ValueError(
            f"target elevation ({target} m) must be below the instruments of "
            f"platform {platform.name!r} ({platform.get_altitude()} m)"
        )
    return target


def collect_platform_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect the static footprint of each instrument on fixed surface
    platforms (e.g. towers or tethered balloons), projected to a target
    elevation by computing the exact ray/WGS 84 ellipsoid intersection, as
    for a satellite's ground track (see `collect_ground_track`). Rays that
    miss the target surface (looking above the horizon) are projected to the
    horizon. Unlike `collect_ground_track`, there is no time dimension: a
    platform does not move in the Earth-fixed frame, so this returns one row
    per platform and instrument rather than one row per time step.

    Args:
        platforms (SurfacePlatform | list[SurfacePlatform]): The platform(s).
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` for every instrument of each
                platform.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid of the observed target for which to compute
                footprints, which must be below each platform's instruments.
                Defaults to each platform's ground elevation.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected platform track
            results, in the order of the platforms and their instruments.
    """
    platforms = [platforms] if isinstance(platforms, SurfacePlatform) else platforms
    records = []
    for platform in platforms:
        target = _get_target_elevation(platform, elevation)
        track = platform.get_track()
        for index in _get_instrument_indices(platform, instrument_index, "platform"):
            instrument = platform.instruments[index]
            records.append(
                {
                    "platform": platform.name,
                    "instrument": instrument.name,
                    "elevation": target,
                    "geometry": instrument.compute_footprint(track, elevation=target)[
                        0
                    ],
                }
            )
    if len(records) == 0:
        return _get_empty_platform_track()
    track = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if mask is not None:
        # clipping does not preserve the order of rows
        track = gpd.clip(track, _get_mask_geometry(mask)).sort_index()
        track = track.reset_index(drop=True)
    return track


def compute_platform_track(
    platforms: SurfacePlatform | list[SurfacePlatform],
    instrument_index: int | None = 0,
    elevation: float | None = None,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated (dissolved) coverage of instruments on fixed
    surface platforms into a single composite geometry (e.g. a network
    coverage map): unlike `collect_platform_track` (one footprint polygon per
    platform and instrument), this dissolves every footprint into one
    geometry.

    Args:
        platforms (SurfacePlatform | list[SurfacePlatform]): The platform(s).
        instrument_index (int | None): The index of the observing instrument
                on each platform, or `None` for every instrument of each
                platform.
        elevation (float | None): The elevation (meters) above the WGS 84
                ellipsoid of the observed target for which to compute
                footprints, which must be below each platform's instruments.
                Defaults to each platform's ground elevation.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.

    Returns:
        geopandas.GeoDataFrame: The data frame of aggregated platform track results.
    """
    track = collect_platform_track(platforms, instrument_index, elevation, mask)
    if track.empty:
        return track
    return track[["geometry"]].dissolve()
