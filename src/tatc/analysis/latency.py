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
    Get an empty downlinks data frame.
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
    Collect satellite downlink opportunities to a ground station of interest.

    :param stations: The ground station(s) of interest.
    :type stations: :class:`tatc.schemas.GroundStation`
    :param satellite: The satellite performing the downlink
    :type satellite: :class:`tatc.schemas.Satellite`
    :param start: The start of the mission window
    :type start: :class:`datetime.datetime`
    :param end: The end of the mission window
    :type end: :class:`datetime.datetime`
    :return: An instance of  :class:`geopandas.GeoDataFrame`
        with all recorded downlink opportunities.
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    # convert orbit to tle
    orbit = satellite.orbit.to_tle()
    # construct a satellite for propagation
    sat = EarthSatellite(orbit.tle[0], orbit.tle[1], satellite.name)
    # collect the records of ground station overpasses
    records = [
        {
            "station": station.name,
            "geometry": geo.Point(station.longitude, station.latitude),
            "satellite": satellite.name,
            "start": period.left,
            "end": period.right,
            "epoch": period.mid,
        }
        for station in (stations if type(stations) == list else [stations])
        for period in _get_visible_interval_series(
            wgs84.latlon(station.latitude, station.longitude),
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
    Get an empty latency data frame.
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
    Collect data latencies between an observation and the first downlink opportunity.

    :param observation: An instance of  :class:`geopandas.GeoDataFrame` containing
        a single observation
    :type observation: :class:`geopandas.GeoDataFrame`
    :param downlinks: The full collection of downlink opportunities
    :type downlinks: :class:`geopandas.GeoDataFrame`
    :return: An instance of  :class:`geopandas.GeoDataFrame` with data latency information.
    :rtype: :class:`geopandas.GeoDataFrame`
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
    Get an empty reduce data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "latency": pd.Series([], dtype="timedelta64[ns]"),
        "samples": pd.Series([], dtype="int"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def reduce_latencies(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reduce observation latencies between observations and first downlink opportunities.

    :param gdf: The aggregated latencies
    :type gdf: :class:`geopandas.GeodataFrame`
    :return: An instance of :class:`geopandas.GeoDataFrame`: The data frame
        with reduced latencies.
    :rtype: :class:`geopanadas.GeodataFrame`
    """
    if gdf.notna().empty:
        return _get_empty_reduce_frame()
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
