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
    station: GroundStation, satellite: Satellite, start: datetime, end: datetime
) -> gpd.GeoDataFrame:
    """
    Collect satellite downlink opportunities to a ground station of interest.

    :param station: The ground station of interest.
    :type station: :class:`tatc.schemas.GroundStation`
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
    # build a topocentric point at the designated ground station
    topos = wgs84.latlon(station.latitude, station.longitude)
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
        for period in _get_visible_interval_series(
            topos, sat, station.min_elevation_angle, start, end
        )
        if (station.min_access_time <= period.right - period.left)
    ]
    # build the dataframe
    if len(records) > 0:
        gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
    else:
        gdf = _get_empty_downlinks_frame()
    return gdf


def collect_multi_downlinks(
    stations: List[GroundStation],
    satellite: Satellite,
    start: datetime,
    end: datetime,
) -> gpd.GeoDataFrame:
    """
    Collect downlink opportunities to multiple ground stations.

    :param stations: The ground stations of interest
    :type stations: :class:`tatc.schemas.GroundStation`
    :param satellites: The observing satellites
    :type satellites: list of :class:`tatc.schemas.SpaceSystem`
    :param start: The start of the mission window
    :type start: :`datetime.datetime`
    :param end: The end of the mission window
    :type end: :class:`datetime.datetime`
    :return: an instance of :class:`geopandas.GeoDataFrame` containing all
        recorded downlinks opportunities
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    # collect downlinks for each ground station
    gdfs = [collect_downlinks(station, satellite, start, end) for station in stations]
    # merge the observations into one data frame
    df = pd.concat(gdfs, ignore_index=True)
    # sort the values by start datetime
    df = df.sort_values("start").reset_index(drop=True)
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")


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
    if gdf.empty:
        return _get_empty_reduce_frame()
    # convert latency to a numeric value before aggregation
    gdf["latency"] = gdf["latency"] / timedelta(seconds=1)
    # assign each record to one observation
    gdf["samples"] = 1
    # perform the aggregation operation
    gdf = gpd.GeoDataFrame(
        gdf.groupby("point_id").agg(
            {
                "point_id": "first",
                "geometry": "first",
                "latency": "mean",
                "samples": "sum",
            }
        ),
        crs="EPSG:4326",
    )
    # convert latency from numeric values after aggregation
    gdf["latency"] = gdf["latency"].apply(
        lambda t: pd.NaT if pd.isna(t) else timedelta(seconds=t)
    )
    return gdf
