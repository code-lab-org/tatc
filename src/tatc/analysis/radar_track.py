"""
Methods to collect ground-based radar coverage.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

import geopandas as gpd
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon

from ..schemas import RadarStation
from .check import _check_platforms, _dissolve_footprints
from .ground_track import _clip_track, _get_mask_geometry


def _get_empty_radar_track() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for radar track results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "radar_station": pd.Series([], dtype="str"),
        "band": pd.Series([], dtype="str"),
        "elevation": pd.Series([], dtype="float"),
        "inner_ground_range": pd.Series([], dtype="float"),
        "outer_ground_range": pd.Series([], dtype="float"),
        "valid_obs": pd.Series([], dtype="bool"),
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
    returns one row per station rather than one row per time step. An empty
    footprint (as of a target that no beam reaches) is not a valid
    observation (`valid_obs`), as for other platforms.

    Args:
        stations (RadarStation | list[RadarStation]): The radar station(s).
        elevation (float): The elevation (meters) of the observed target
                for which to compute coverage, above the WGS 84
                ellipsoid (as the stations' elevations): an absolute
                elevation, common to all stations, not a height relative
                to each station.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected radar track
            results, with columns `radar_station`, `band`, `elevation`,
            `inner_ground_range`, `outer_ground_range`, `valid_obs`, and
            `geometry`, in the order of the stations.
    """
    stations = _check_platforms(stations, RadarStation, "stations")
    if len(stations) == 0:
        return _get_empty_radar_track()
    records = []
    for station in stations:
        ranges = station.compute_ground_ranges(elevation)
        footprint = station.compute_footprint(elevation)
        records.append(
            {
                "radar_station": station.name,
                "band": station.band.value if station.band is not None else None,
                "elevation": elevation,
                "inner_ground_range": ranges[0] if ranges is not None else None,
                "outer_ground_range": ranges[1] if ranges is not None else None,
                "valid_obs": not footprint.is_empty,
                "geometry": footprint,
            }
        )
    track = gpd.GeoDataFrame(records, crs="EPSG:4326")
    return _clip_track(track, _get_mask_geometry(mask) if mask is not None else None)


def compute_radar_track(
    stations: RadarStation | list[RadarStation],
    elevation: float,
    mask: Polygon | MultiPolygon | gpd.GeoDataFrame | gpd.GeoSeries | None = None,
) -> gpd.GeoDataFrame:
    """
    Compute the aggregated (dissolved) radar coverage across all stations
    into a single composite geometry (e.g. a network coverage map): unlike
    `collect_radar_track` (one footprint polygon per station), this
    dissolves the footprints of valid observations into one geometry (see
    `tatc.analysis.footprints.dissolve_footprints`).

    Args:
        stations (RadarStation | list[RadarStation]): The radar station(s).
        elevation (float): The elevation (meters) of the observed target
                for which to compute coverage, above the WGS 84
                ellipsoid (as the stations' elevations): an absolute
                elevation, common to all stations, not a height relative
                to each station.
        mask (shapely.geometry.Polygon | shapely.geometry.MultiPolygon | geopandas.GeoDataFrame | geopandas.GeoSeries | None):
                An optional mask, always interpreted in WGS84 (lon/lat)
                coordinates, to constrain results.

    Returns:
        geopandas.GeoDataFrame: The data frame of the aggregated radar
            track: one row with the observed area, or none if no
            observation is valid.
    """
    return _dissolve_footprints(collect_radar_track(stations, elevation, mask))
