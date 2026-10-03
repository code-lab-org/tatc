"""
Methods to collect ground-based radar coverage.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon

from ..schemas import RadarStation


def _get_empty_radar_track() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for radar track results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "station": pd.Series([], dtype="str"),
        "band": pd.Series([], dtype="str"),
        "elevation": pd.Series([], dtype="float"),
        "inner_ground_range": pd.Series([], dtype="float"),
        "outer_ground_range": pd.Series([], dtype="float"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_radar_track(
    stations: RadarStation | list[RadarStation],
    elevation: float,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect each radar station's static coverage footprint, projected to a
    specified target elevation. Unlike `collect_ground_track`, there is no
    time dimension: a ground-based radar station does not move, so this
    returns one row per station rather than one row per time step.

    Args:
        stations (RadarStation | list[RadarStation]): The radar station(s).
        elevation (float): The elevation (meters) above the WGS 84 datum
                of the observed target for which to compute coverage (an
                absolute elevation, common to all stations, not a height
                relative to each station).
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected radar track results.
    """
    stations_list = [stations] if isinstance(stations, RadarStation) else stations
    if len(stations_list) == 0:
        return _get_empty_radar_track()
    records = [
        {
            "station": station.name,
            "band": station.band.value if station.band is not None else None,
            "elevation": elevation,
            "inner_ground_range": ranges[0] if ranges is not None else None,
            "outer_ground_range": ranges[1] if ranges is not None else None,
            "geometry": station.compute_footprint(elevation),
        }
        for station in stations_list
        for ranges in [station.compute_ground_ranges(elevation)]
    ]
    track = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if mask is not None:
        track = gpd.clip(track, mask).reset_index(drop=True)
    return track


def compute_radar_track(
    stations: RadarStation | list[RadarStation],
    elevation: float,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated (dissolved) radar coverage across all stations
    into a single composite geometry (e.g. a network coverage map): unlike
    `collect_radar_track` (one footprint polygon per station), this
    dissolves every station's footprint into one geometry.

    Args:
        stations (RadarStation | list[RadarStation]): The radar station(s).
        elevation (float): The elevation (meters) above the WGS 84 datum
                of the observed target for which to compute coverage (an
                absolute elevation, common to all stations, not a height
                relative to each station).
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.

    Returns:
        geopandas.GeoDataFrame: The data frame of aggregated radar track results.
    """
    track = collect_radar_track(stations, elevation, mask)
    if track.empty:
        return track
    return track[["geometry"]].dissolve()
