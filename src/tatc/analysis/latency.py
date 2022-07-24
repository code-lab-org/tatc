# -*- coding: utf-8 -*-
"""
Methods to perform latency analysis.

@author: Isaac Feldman, Paul T. Grogan <pgrogan@stevens.edu>
"""

import numpy as np
import pandas as pd
import geopandas as gpd
from typing import List, Union
from shapely import geometry as geo
from numba import njit
from datetime import datetime, timedelta
from skyfield.api import load, wgs84, EarthSatellite

from ..schemas.point import Point, GroundStation
from ..schemas.satellite import Satellite
from ..schemas.instrument import Instrument

from .coverage import _get_visible_interval_series
from ..utils import (
    compute_min_elevation_angle,
    swath_width_to_field_of_regard,
    compute_max_access_time,
)
from ..constants import de421, timescale


def _get_empty_downlinks_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for downlink results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "station": pd.Series([], dtype="str"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_downlinks(
    stations: Union[GroundStation, List[GroundStation]],
    satellite: Satellite,
    start: datetime,
    end: datetime,
) -> gpd.GeoDataFrame:
    """
    Collect satellite downlink opportunities to ground station(s) of interest.

    Args:
        stations (GroundStation or typing.List[GroundStation]): The ground stations.
        satellite (Satellite): The observing satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected downlink results.
    """
    # convert orbit to tle
    orbit = satellite.orbit.to_tle()
    # construct a satellite for propagation
    sat = EarthSatellite(orbit.tle[0], orbit.tle[1], satellite.name)
    # collect the records of ground station overpasses
    records = [
        {
            "station": station.name,
            "geometry": geo.Point(
                station.longitude, station.latitude, station.elevation
            ),
            "satellite": satellite.name,
            "start": period.left,
            "end": period.right,
            "epoch": period.mid,
        }
        for station in (stations if type(stations) == list else [stations])
        for period in _get_visible_interval_series(
            wgs84.latlon(station.latitude, station.longitude, station.elevation),
            sat,
            station.min_elevation_angle,
            start,
            end,
        )
        if (station.min_access_time <= period.right - period.left)
    ]
    # build the dataframe
    if len(records) > 0:
        gdf = (
            gpd.GeoDataFrame(records, crs="EPSG:4326")
            .sort_values("start")
            .reset_index(drop=True)
        )
    else:
        gdf = _get_empty_downlinks_frame()
    return gdf


def _get_empty_latency_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for downlink results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "observed": pd.Series([], dtype="datetime64[ns, utc]"),
        "station": pd.Series([], dtype="str"),
        "downlinked": pd.Series([], dtype="datetime64[ns, utc]"),
        "latency": pd.Series([], dtype="timedelta64[ns]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def compute_latencies(
    observations: gpd.GeoDataFrame, downlinks: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Collect latencies between an observation and the first downlink opportunity.

    Args:
        observations (geopandas.GeoDataFrame): The data frame of observations to downlink.
        downlinks (geopandas.GeoDataFrame): The data frame of downlink opportunities.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected latency results.
    """
    if observations.empty or downlinks.empty:
        return _get_empty_latency_frame()

    def _align_downlinks(r):
        # filter downlinks after observation occurs
        dls = downlinks[
            np.logical_and(r.satellite == downlinks.satellite, r.end < downlinks.start)
        ]
        # append latency-specific columns
        if dls.empty:
            r["station"] = pd.NA
            r["downlinked"] = pd.NA
            r["latency"] = pd.NA
        else:
            r["station"] = dls.iloc[0].station
            r["downlinked"] = dls.iloc[0].epoch
            r["latency"] = dls.iloc[0].epoch - r.epoch
        return r

    # append the latency-specific columns
    observations = observations.apply(_align_downlinks, axis=1)
    # add observed column
    observations["observed"] = observations["epoch"]
    # drop start, epoch, and end columns
    observations = observations.drop(["start", "epoch", "end"], axis=1)
    return observations


def _get_empty_reduce_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for reduced latency results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "latency": pd.Series([], dtype="timedelta64[ns]"),
        "samples": pd.Series([], dtype="int"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def reduce_latencies(latency_observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reduce observation latencies. Computes descriptive statistics for each
    pair of observation and first downlink opportunities.

    Args:
        latency_observations (geopandas.GeoDataFrame): The latency observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced latencies.
    """
    if latency_observations.notna().empty:
        return _get_empty_reduce_frame()
    # operate on a copy of the dataframe
    gdf = latency_observations.copy()
    # convert latency to a numeric value before aggregation
    gdf.loc[gdf.latency.notna(), "latency"] = gdf.loc[
        gdf.latency.notna(), "latency"
    ] / timedelta(seconds=1)
    # assign each record to one observation
    gdf["samples"] = 1
    # perform the aggregation operation
    gdf = gdf.dissolve(
        "point_id",
        aggfunc={
            "latency": "mean",
            "samples": "sum",
        },
    )
    # convert latency from numeric values after aggregation
    gdf["latency"] = gdf["latency"].apply(
        lambda t: pd.NaT if pd.isna(t) else timedelta(seconds=t)
    )
    return gdf
