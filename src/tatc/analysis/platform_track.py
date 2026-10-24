"""
Methods to collect the static tracks of instruments on fixed surface platforms.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
import shapely
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
    differ from the elevation of its instruments (below them, looking down,
    or above them, looking up).

    Args:
        platform (SurfacePlatform): The platform.
        elevation (float | None): The target elevation (meters) above the
            WGS 84 ellipsoid, or None for the platform's ground elevation.

    Returns:
        float: the target elevation (meters).
    """
    target = platform.elevation if elevation is None else elevation
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
    miss a target below the instruments (looking above the horizon) are
    projected to the horizon. A target above the instruments (as for clouds
    observed by an instrument looking up, e.g. with a roll angle of 180
    degrees) is observed only above the horizon: views are clipped to the
    lines of sight tangent to the ground, which blocks the rest (a view
    entirely below the horizon has an empty footprint), and to the
    platform's terrain mask, if any (see
    `SurfacePlatform.compute_visible_region`). Unlike
    `collect_ground_track`, there is no time dimension: a platform does not
    move in the Earth-fixed frame, so this returns one row per platform and
    instrument rather than one row per time step.

    Args:
        platforms (SurfacePlatform | list[SurfacePlatform]): The platform(s).
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
        geopandas.GeoDataFrame: The data frame of collected platform track
            results, in the order of the platforms and their instruments.
    """
    platforms = [platforms] if isinstance(platforms, SurfacePlatform) else platforms
    records = []
    for platform in platforms:
        target = _get_target_elevation(platform, elevation)
        track = platform.get_track()
        # the region visible over the terrain mask, toward a target above
        region = (
            platform.compute_visible_region(target)
            if platform.terrain_mask is not None and target > platform.get_altitude()
            else None
        )
        for index in _get_instrument_indices(platform, instrument_index, "platform"):
            instrument = platform.instruments[index]
            footprint = instrument.compute_footprint(
                track, elevation=target, ground_elevation=platform.elevation
            )[0]
            if region is not None:
                footprint = _clip_to_region(footprint, region, target)
            records.append(
                {
                    "platform": platform.name,
                    "instrument": instrument.name,
                    "elevation": target,
                    "geometry": footprint,
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
                footprints, below (looking down) or above (looking up) each
                platform's instruments. Defaults to each platform's ground
                elevation.
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
