# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 21:31:32 2021

@author: Isaac
"""
# IMPORTED MODULES-------------------------------------------------------------
import numpy as np
import pandas as pd
import geopandas as gpd
from typing import List
from shapely import geometry as geo
from numba import njit
from datetime import datetime, timedelta
from skyfield.api import load, wgs84, EarthSatellite
# from tatc.visualization import plot_points_latency, plot_cells_latency
# from tatc_core.generation.cells import generate_cubed_sphere_cells
from ..schemas.point import Point, GroundStation
from ..schemas.satellite import Satellite
from ..schemas.instrument import Instrument

from ..utils import compute_min_altitude, swath_width_to_field_of_regard, compute_max_access_time

# FUNCTIONS---------------------------------------------------------------------
def collect_downlinks(
        station: GroundStation,
        satellite: Satellite,
        start: datetime,
        end: datetime
    ):
    """
    Collect satellite downlik opportunities to a groundstation of interest.

    :param station: The ground station of interest.
    :type station: class:`tatc.schemas.point.GroundStation`
    :param satellite: The satellite performing the downlink
    :type satellite: class:`tatc.schemas.satellite.Satellite`
    :param start: The start of the mission window
    :type start: class:`datetime.datetime`
    :param end: The end of the mission window
    :type end: class:`datetime.datetime`
    :return: An instance of :class:`geopandas.GeoDataFrame`
        with all recorded downlink opportunities.
    :rtype: class:`geopandas.GeoDataFrame`
    """
    # build a topocentric point at the designated geodetic point (the ground station)
    topos = wgs84.latlon(station.latitude, station.longitude)
    # load the timescale and define starting and ending points
    ts = load.timescale()
    t0 = ts.utc(start)
    t1 = ts.utc(end)
    # load the ephemerides
    eph = load('de421.bsp')
    # convert orbit to tle
    orbit = satellite.orbit.to_tle()
    # construct a satellite for propagation
    sat = EarthSatellite(orbit.tle[0], orbit.tle[1], satellite.name)
    # compute the initial satellite height (altitude)
    satellite_height = wgs84.subpoint(sat.at(t0)).elevation.m
    # compute the maximum access time to filter bad data
    max_access_time = timedelta(seconds=compute_max_access_time(satellite_height,
                                                                 station.min_elevation_angle))
    # find the set of observation events
    t, events = sat.find_events(topos, t0, t1, altitude_degrees=station.min_elevation_angle)
    # initialize dataframe
    df = pd.DataFrame({
            'station_name': pd.Series([], dtype='str'),
            'geometry': pd.Series([]),
            'satellite': pd.Series([], dtype='str'),
            'start': pd.Series([], dtype='datetime64[ns, utc]'),
            'epoch': pd.Series([], dtype='datetime64[ns, utc]'),
            'end': pd.Series([], dtype='datetime64[ns, utc]')
            })
    # check for geocentricity
    if numpy.all(events == 1):
        df = df.append({
            'station_name': station.name,
            'geometry': geo.Point(station.longitude, station.latitude),
            'satellite': satellite.name,
            'start': start,
            'epoch': start + (end - start) / 2,
            'end': end
        }, ignore_index=True)
    # define variables for stepping through the events list
    t_rise = None
    t_culminate = None
    for j in range(len(events)):
        if events[j] == 0:
            # record the rise time
            t_rise = t[j].utc_datetime()
        elif events[j] == 1:
            # record the culmination time
            t_culminate = t[j].utc_datetime()
        elif events[j] == 2:
            # record the set time
            t_set = t[j].utc_datetime()
            # only record an observation if a previous rise and culminate
            # events were recorded (sometimes they are out-of-order)
            if t_rise is not None and t_culminate is not None:
                # check if the observation meets minimum access duration,
                # ground sunlit conditions, and satellite sunlit conditions
                if (station.min_access_time <= t_set - t_rise <= max_access_time*2):
                    df = df.append({
                        'station_name': station.name,
                        'geometry': geo.Point(station.longitude, station.latitude),
                        'satellite': satellite.name,
                        'start': pd.Timestamp(t_rise),
                        'epoch': pd.Timestamp(t_culminate),
                        'end': pd.Timestamp(t_set)
                    }, ignore_index=True)
            # reset the variables for stepping through the event list
            t_rise = None
            t_culminate = None
            sat_alt = None
            sat_az = None

    # compute the access time for each downlink(end - start)
    df['access'] = df['end'] - df['start']
    # compute the revisit time for each downlink (previous end - start)
    df['revisit'] = df['end'] - df['start'].shift()
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")

def aggregate_downlinks(gdfs: List[gpd.GeoDataFrame]):
    """
    Aggregate constellation downlink opportunities for a ground station of interest.

    :param gdfs: The individual downlink opportunities.
    :type gdfs: list of class:`geopandas.GeoDataFrame`
    :return: An instance of :class:`geopandas.GeoDataFrame` with aggregated observations.
    :rtype: class:`geopandas.GeoDataFrame`
    """
    if all(len(gdf.index)==0 for gdf in gdfs):
        empty_df = pd.DataFrame({
            'station_name': pd.Series([]),
            'geometry': pd.Series([]),
            'satellite': pd.Series([]),
            'start': pd.Series([]),
            'epoch': pd.Series([]),
            'end': pd.Series([]),
            'access': pd.Series([]),
            'revisit': pd.Series([])
        })
        return gpd.GeoDataFrame(
            empty_df,
            geometry=empty_df.geometry,
            crs="EPSG:4326"
        )
    # merge the downlink oportunitiess into one data frame
    df = pd.concat(gdfs, ignore_index=True)
    # sort the values by start datetime
    df = df.sort_values('epoch')
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")

def compute_latency(
        observation:gpd.GeoDataFrame,
        downlinks: pd.DataFrame):
    """
    Compute data latency between an observation and the first downlink opportunity.

    :param observation: An instance of :class:`geopandas.GeoDataFrame` containing
        a single observation
    :type observation: class:`geopandas.GeoDataFrame`
    :param: downlink_ops_full: The full collection of dowlink opportunities
    :type: downlink_ops_full: class:`geopandas.GeoDataFrame`
    :return: An instance of :class:`geopandas.GeoDataFrame` with data latency information.
    :rtype: class:`geopandas.GeoDataFrame`
    """
    # Create a database to save latencies
    df = pd.DataFrame({
        'point_id': pd.Series([], dtype='int'),
        'geometry': pd.Series([]),
        'satellite': pd.Series([], dtype='str'),
        'observing_instrument':pd.Series([], dtype='str'),
        'station_name':pd.Series([], dtype='str'),
        'observed': pd.Series([], dtype='datetime64[ns, utc]'),
        'downlinked': pd.Series([], dtype='datetime64[ns, utc]'),
        'latency': pd.Series([], dtype='timedelta64[ns]')
    })

    # Pull out the epoch time of the observation
    obs_time = observation['epoch']
    # return an empty dataframe if no downlink opportunities input
    if downlinks.empty:
        df = df.append({'point_id': observation.point_id,
                        'geometry': observation.geometry,
                        'satellite': observation.satellite,
                        'observing_instrument': observation.instrument,
                        'station_name': None,
                        'observed': obs_time,
                        'downlinked': None,
                        'latency': None
                        }, ignore_index=True)

        return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")
    # filter the downlink opportunities to be only those after the
    # observation
    downlinks = downlinks.loc[downlinks['epoch'] >= obs_time]
    # return an empty data frame if there are no downlink opportunities
    # after the observation
    if downlinksd.empty:
        df = df.append({'point_id': observation.point_id,
                        'geometry': observation.geometry,
                        'satellite': observation.satellite,
                        'observing_instrument': observation.instrument,
                        'station_name': None,
                        'observed': obs_time,
                        'downlinked': None,
                        'latency': None
                        }, ignore_index=True)

        return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")
    # reindex the downlink opportunities starting at zero
    downlinks.sort_values(by='epoch')
    downlinks = downlinks.reset_index()
    # extract the first downlink opportunity
    first_dl_op = downlinks.loc[0]
    # pull out the time of the first downlink opportunity
    first_dl_time = first_dl_op.epoch
    # Calculate the latency as the difference between the first downlink and
    # the first observation
    latency = pd.Timestamp(first_dl_time) - pd.Timestamp(obs_time)

    df = df.append({'point_id': observation.point_id,
                    'geometry': observation.geometry,
                    'satellite': observation.satellite,
                    'observing_instrument': observation.instrument,
                    'station_name': first_dl_op.station_name,
                    'observed': obs_time,
                    'downlinked': first_dl_time,
                    'latency': latency
                    }, ignore_index=True)

    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")

def aggregate_latencies(dfs: List[pd.DataFrame]):
    """
    Aggregate constellation latency information for a ground station of interest

    :param dfs: The individual downlink opportunities.
    :type dfs: list class:`pandas.DataFrame`
    :return: An instance of :class:`geopandas.GeoDataFrame` with
        aggregated latency information.
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    if all(len(df.index)==0 for df in dfs):
        empty_df = pd.DataFrame({
                'point_id': pd.Series([]),
                'geometry': pd.Series([]),
                'satellite': pd.Series([]),
                'observing_instrument': pd.Series([]),
                'station_name': pd.Series([]),
                'observed': pd.Series([]),
                'downlinked': pd.Series([]),
                'latency': pd.Series([]),
                }, ignore_index=True)
        return gpd.GeoDataFrame(
            empty_df,
            geometry=empty_df.geometry,
            crs="EPSG:4326"
        )
    # merge the downlink oportunitiess into one data frame
    df = pd.concat(dfs, ignore_index=True)
    # sort the values by start datetime
    df = df.sort_values('latency')
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")
