"""
Methods to collect ground-based radar coverage.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime

import geopandas as gpd
import numpy as np
import pandas as pd

from ..schemas import RadarStation
from ._tracks import Mask, _clip_track, _get_mask_geometry
from .check import (
    RADAR_STATION,
    _check_distinct_names,
    _check_platforms,
    _check_times,
    _dissolve_footprints,
)


def _get_empty_radar_track(timed: bool = False) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for radar track results.

    Args:
        timed (bool): `True` to include a `time` column.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        **({"time": pd.Series([], dtype="datetime64[ns, utc]")} if timed else {}),
        "radar_station": pd.Series([], dtype="str"),
        "band": pd.Series([], dtype="str"),
        "inner_ground_range": pd.Series([], dtype="float"),
        "outer_ground_range": pd.Series([], dtype="float"),
        "valid_obs": pd.Series([], dtype="bool"),
        "geometry": pd.Series([], dtype="object"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_radar_track(
    stations: RadarStation | list[RadarStation],
    elevation: float,
    mask: Mask | None = None,
    *,
    times: list[datetime] | None = None,
) -> gpd.GeoDataFrame:
    """
    Collect each radar station's static coverage footprint, projected to a
    specified target elevation. A ground-based radar station does not move,
    so its footprint does not change with time: without `times`, this
    returns one row per station, without a time. With `times`, it returns
    the same footprint at each time (one row per time and station), as
    `collect_surface_track` does for fixed surface platforms, so that the
    coverage is sampled at the times of other platforms' footprints (see
    `tatc.analysis.footprints.collect_footprint_observations`). An empty
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
        times (list[datetime.datetime] | None): The times at which to sample
                the coverage, or `None` for a static footprint of each
                station (keyword-only).

    Returns:
        geopandas.GeoDataFrame: The data frame of collected radar track
            results, with columns `radar_station`, `band`,
            `inner_ground_range`, `outer_ground_range`, `valid_obs`, and
            `geometry` (whose z-coordinate is the target elevation): without
            `times`, in the order of the stations; with
            `times`, also `time` (first), sorted by time and then in the
            order of the stations.
    """
    stations = _check_platforms(stations, RADAR_STATION)
    _check_distinct_names(stations, RADAR_STATION, RADAR_STATION.argument)
    _check_times(times)
    timed = times is not None
    if len(stations) == 0 or (timed and len(times) == 0):  # type: ignore
        return _get_empty_radar_track(timed)
    records = []
    for station in stations:
        ranges = station.compute_ground_ranges(elevation)
        footprint = station.compute_footprint(elevation)
        records.append(
            {
                "radar_station": station.name,
                "band": station.band.value if station.band is not None else None,
                "inner_ground_range": ranges[0] if ranges is not None else None,
                "outer_ground_range": ranges[1] if ranges is not None else None,
                "valid_obs": not footprint.is_empty,
                "geometry": footprint,
            }
        )
    track = gpd.GeoDataFrame(records, crs="EPSG:4326")
    if timed:
        # the footprint of each station at each time
        track = track.iloc[np.tile(np.arange(len(stations)), len(times))]  # type: ignore
        track.insert(
            0,
            "time",
            np.repeat(pd.to_datetime(times, utc=True).as_unit("ns"), len(stations)),
        )
        track = track.reset_index(drop=True)
    return _clip_track(track, _get_mask_geometry(mask) if mask is not None else None)


def compute_radar_track(
    stations: RadarStation | list[RadarStation],
    elevation: float,
    mask: Mask | None = None,
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
