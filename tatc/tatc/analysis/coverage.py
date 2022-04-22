# -*- coding: utf-8 -*-
"""
Methods to perform coverage analysis.

@author: Paul T. Grogan <pgrogan@stevens.edu>
"""

import pandas as pd
import numpy as np
import geopandas as gpd
from typing import List, Optional
from shapely import geometry as geo
from datetime import datetime, timedelta
from skyfield.api import load, wgs84, EarthSatellite

from ..schemas.point import Point
from ..schemas.satellite import Satellite
from ..schemas.instrument import Instrument, DutyCycleScheme

from ..utils import (
    compute_min_altitude,
    swath_width_to_field_of_regard,
    compute_max_access_time,
    compute_orbit_period,
)


def collect_observations(
    point: Point,
    satellite: Satellite,
    instrument: Instrument,
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
    sample_distance: Optional[float] = None,
) -> gpd.GeoDataFrame:
    """
    Collect single satellite observations of a geodetic point of interest.

    :param point: The ground point of interest
    :type point: :class:`tatc.schemas.point.Point`
    :param satellite: The observing satellite
    :type satellite: :class:`tatc.schemas.satellite.Satellite`
    :param instrument: The instrument used to make observations
    :type instrument::`tatc.schemas.instrument.instrument`
    :param start: The start of the mission window
    :type start::`datetime.datetime`
    :param end: The end of the mission window
    :type end::`datetime.datetime`
    :param omit_solar: True, if solar angles should be omitted
        to improve computational efficiency, defaults to True
    :type omit_solar: bool, optional
    :param sample_distance: Ground sample distance (m) to override
        instrument field of regard, defaults to None
    :type sample_distance: int, optional
    :return: An instance of :class:`geopandas.GeoDataFrame` containing all
        recorded reduce_observations
    :rtype::`geopandas.GeoDataFrame`
    """

    # build a topocentric point at the designated geodetic point
    topos = wgs84.latlon(point.latitude, point.longitude)
    # load the timescale and define starting and ending points
    ts = load.timescale()
    t0 = ts.from_datetime(start)
    t1 = ts.from_datetime(end)
    # load the ephemerides
    eph = load("de421.bsp")
    # convert orbit to tle
    orbit = satellite.orbit.to_tle()
    # construct a satellite for propagation
    sat = EarthSatellite(orbit.tle[0], orbit.tle[1], satellite.name)
    # compute the initial satellite height (altitude)
    satellite_height = wgs84.subpoint(sat.at(t0)).elevation.m
    # compute the minimum altitude angle required for observation
    min_altitude = compute_min_altitude(
        satellite_height,
        instrument.field_of_regard
        if sample_distance is None
        else swath_width_to_field_of_regard(satellite_height, sample_distance),
    )
    # compute the maximum access time to filter bad data
    max_access_time = timedelta(
        seconds=compute_max_access_time(satellite_height, min_altitude)
    )
    # TODO: consider instrument operational intervals
    ops_intervals = pd.Series(
        [pd.Interval(pd.Timestamp(start), pd.Timestamp(end), "both")]
    )

    # find the set of observation events
    t, events = sat.find_events(topos, t0, t1, altitude_degrees=min_altitude)

    if omit_solar:
        # basic dataframe without solar angles
        df = pd.DataFrame(
            {
                "point_id": pd.Series([], dtype="int"),
                "geometry": pd.Series([], dtype="object"),
                "satellite": pd.Series([], dtype="str"),
                "instrument": pd.Series([], dtype="str"),
                "start": pd.Series([], dtype="datetime64[ns, utc]"),
                "end": pd.Series([], dtype="datetime64[ns, utc]"),
                "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
                "sat_alt": pd.Series([], dtype="float64"),
                "sat_az": pd.Series([], dtype="float64"),
            }
        )
    else:
        # extended dataframe including solar angles
        df = pd.DataFrame(
            {
                "point_id": pd.Series([], dtype="int"),
                "geometry": pd.Series([], dtype="object"),
                "satellite": pd.Series([], dtype="str"),
                "instrument": pd.Series([], dtype="str"),
                "start": pd.Series([], dtype="datetime64[ns, utc]"),
                "end": pd.Series([], dtype="datetime64[ns, utc]"),
                "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
                "sat_alt": pd.Series([], dtype="float64"),
                "sat_az": pd.Series([], dtype="float64"),
                "sat_sunlit": pd.Series([], dtype="bool"),
                "solar_alt": pd.Series([], dtype="float64"),
                "solar_az": pd.Series([], dtype="float64"),
                "solar_time": pd.Series([], dtype="float64"),
            }
        )
    # define variables for stepping through the events list
    t_rise = None
    t_culminate = None
    sat_sunlit = None
    solar_time = None
    sat_alt = None
    sat_az = None
    solar_alt = None
    solar_az = None
    # check for geocentricity
    if np.all(events == 1) and events:
        # find the satellite altitude, azimuth, and distance at t0
        sat_alt, sat_az, sat_dist = (sat - topos).at(t[0]).altaz()
        # if ommiting solar angles
        if omit_solar:
            df = pd.concat([
                df,
                pd.DataFrame.from_records(
                    {
                        "point_id": point.id,
                        "geometry": geo.Point(point.longitude, point.latitude),
                        "satellite": satellite.name,
                        "instrument": instrument.name,
                        "start": start,
                        "epoch": start + (end - start) / 2,
                        "end": end,
                        "sat_alt": sat_alt.degrees,
                        "sat_az": sat_az.degrees,
                    }, index=[0]
                )
            ], ignore_index=True)
        # otherwise if solar angles are included
        else:
            df = pd.concat([
                df,
                pd.DataFrame.from_records(
                    {
                        "point_id": point.id,
                        "geometry": geo.Point(point.longitude, point.latitude),
                        "satellite": satellite.name,
                        "instrument": instrument.name,
                        "start": start,
                        "epoch": start + (end - start) / 2,
                        "end": end,
                        "sat_alt": sat_alt.degrees,
                        "sat_az": sat_az.degrees,
                        "sat_sunlit": None,
                        "solar_alt": None,
                        "solar_az": None,
                        "solar_time": None
                    }, index=[0]
                )
            ], ignore_index=True)
        # compute the access time for the observation (end - start)
        df["access"] = df["end"] - df["start"]
        # compute the revisit time for each observation (previous end - start)
        df["revisit"] = df["end"] - df["start"].shift()
        return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")

    for j in range(len(events)):
        if events[j] == 0:
            # record the rise time
            t_rise = t[j].utc_datetime()
        elif events[j] == 1:
            # record the culmination time
            t_culminate = t[j].utc_datetime()
            # find the satellite altitude, azimuth, and distance
            sat_alt, sat_az, sat_dist = (sat - topos).at(t[j]).altaz()
            if not omit_solar or instrument.req_target_sunlit is not None:
                # find the solar altitude, azimuth, and distance
                solar_obs = (
                    (eph["earth"] + topos).at(t[j]).observe(eph["sun"]).apparent()
                )
                solar_alt, solar_az, solar_dist = solar_obs.altaz()
                # find the local solar time
                solar_time = solar_obs.hadec()[0].hours + 12
            if not omit_solar or instrument.req_self_sunlit is not None:
                # find whether the satellite is sunlit
                sat_sunlit = sat.at(t[j]).is_sunlit(eph)
        elif events[j] == 2:
            # record the set time
            t_set = t[j].utc_datetime()
            # only record an observation if a previous rise and culminate
            # events were recorded (sometimes they are out-of-order)
            if t_rise is not None and t_culminate is not None:
                # check if the observation meets minimum access duration,
                # ground sunlit conditions, and satellite sunlit conditions
                if (
                    instrument.min_access_time <= t_set - t_rise <= max_access_time * 2
                    and instrument.is_valid_observation(
                        eph,
                        ts.from_datetime(t_culminate),
                        sat.at(ts.from_datetime(t_culminate)),
                    )
                    and (
                        instrument.duty_cycle >= 1
                        or any(ops_intervals.apply(lambda i: t_culminate in i))
                    )
                ):
                    # if omitting solar angles
                    if omit_solar:
                        df = pd.concat([
                            df,
                            pd.DataFrame.from_records(
                                {
                                    "point_id": point.id,
                                    "geometry": geo.Point(point.longitude, point.latitude),
                                    "satellite": satellite.name,
                                    "instrument": instrument.name,
                                    "start": pd.Timestamp(t_rise),
                                    "epoch": pd.Timestamp(t_culminate),
                                    "end": pd.Timestamp(t_set),
                                    "sat_alt": sat_alt.degrees,
                                    "sat_az": sat_az.degrees,
                                }, index=[0]
                            )
                        ], ignore_index=True)
                    # otherwise if solar angles are included
                    else:
                        df = pd.concat([
                            df,
                            pd.DataFrame.from_records(
                                {
                                    "point_id": point.id,
                                    "geometry": geo.Point(point.longitude, point.latitude),
                                    "satellite": satellite.name,
                                    "instrument": instrument.name,
                                    "start": pd.Timestamp(t_rise),
                                    "epoch": pd.Timestamp(t_culminate),
                                    "end": pd.Timestamp(t_set),
                                    "sat_alt": sat_alt.degrees,
                                    "sat_az": sat_az.degrees,
                                    "sat_sunlit": sat_sunlit,
                                    "solar_alt": solar_alt.degrees,
                                    "solar_az": solar_az.degrees,
                                    "solar_time": solar_time,
                                }, index=[0]
                            )
                        ], ignore_index=True)
            # reset the variables for stepping through the event list
            t_rise = None
            t_culminate = None
            sat_sunlit = None
            solar_time = None
            sat_alt = None
            sat_az = None
            solar_alt = None
            solar_az = None

    # compute the access time for each observation (end - start)
    df["access"] = df["end"] - df["start"]
    # compute the revisit time for each observation (previous end - start)
    df["revisit"] = df["end"] - df["start"].shift()
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")


def collect_multi_observations(
    point: Point,
    satellites: List[Satellite],
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
    sample_distance: Optional[float] = None,
) -> gpd.GeoDataFrame:
    """
    Collect multiple satellite observations of a geodetic point of interest.

    :param point: The ground point of interest
    :type point: :class:`tatc.schemas.point.Point`
    :param satellites: The observing satellites
    :type satellites: list of :class:`tatc.schemas.satellite.Satellite`
    :param start: The start of the mission window
    :type start: :`datetime.datetime`
    :param end: The end of the mission window
    :type end: :class:`datetime.datetime`
    :param omit_solar: True, if solar angles should be omitted
        to improve computational efficiency, defaults to True
    :type omit_solar: bool, optional
    :param sample_distance: Ground sample distance (m) to override
        instrument field of regard, defaults to None
    :type sample_distance: int, optional
    :return: an instance of :class:`geopandas.GeoDataFrame` containing all
        recorded observations
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    gdfs = [
        collect_observations(
            point, satellite, instrument, start, end, omit_solar, sample_distance
        )
        for constellation in satellites
        for satellite in (constellation.generate_members())
        for instrument in satellite.instruments
    ]
    # merge the observations into one data frame
    df = pd.concat(gdfs, ignore_index=True)
    # sort the values by start datetime
    df = df.sort_values("start")
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")


def aggregate_observations(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Aggregate constellation observations for a geodetic point of interest.

    :param gdf: The individual observations.
    :type gdf: :class:`geopandas.GeoDataFrame`
    :return: An instance of :class:`geopandas.GeoDataFrame` The data frame
        with aggregated observations.
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    if len(gdf.index) == 0:
        empty_df = pd.DataFrame(
            {
                "point_id": pd.Series([], dtype="int"),
                "geometry": pd.Series([], dtype="object"),
                "start": pd.Series([], dtype="datetime64[ns, utc]"),
                "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
                "end": pd.Series([], dtype="datetime64[ns, utc]"),
                "satellite": pd.Series([], dtype="str"),
                "instrument": pd.Series([], dtype="str"),
                "access": pd.Series([], dtype="timedelta64[ns]"),
                "revisit": pd.Series([], dtype="timedelta64[ns]")
            }
        )
        return gpd.GeoDataFrame(empty_df, geometry=empty_df.geometry, crs="EPSG:4326")

    # sort the values by start datetime
    df = gdf.sort_values("start")
    # assign the observation group number based on overlapping start/end times
    df["obs"] = (df["start"] > df["end"].shift().cummax()).cumsum()
    if all(key in gdf.columns for key in ["solar_alt", "solar_az", "solar_time"]):
        # reduce solar angles
        df = df.groupby("obs").agg(
            {
                "point_id": "first",
                "geometry": "first",
                "start": "min",
                "epoch": "first",
                "end": "max",
                "solar_alt": "mean",
                "solar_az": "mean",
                "solar_time": "mean",
                "satellite": ", ".join,
                "instrument": ", ".join,
            }
        )
    else:
        # reduce only core attributes
        df = df.groupby("obs").agg(
            {
                "point_id": "first",
                "geometry": "first",
                "start": "min",
                "epoch": "first",
                "end": "max",
                "satellite": ", ".join,
                "instrument": ", ".join,
            }
        )
    # compute the access time for each observation (end - start)
    df["access"] = df["end"] - df["start"]
    # compute the revisit time for each observation (previous end - start)
    df["revisit"] = df["end"] - df["start"].shift()
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")


def reduce_observations(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reduce constellation observations for a geodetic point of interest.

    :param gdf: The aggregated observations
    :type gdf: :class:`geopandas.GeodataFrame`
    :return: An instance of :class:`geopandas.GeoDataFrame`: The data frame
        with reduced observations.
    :rtype: :class:`geopanadas.GeodataFrame`
    """
    if len(gdf.index) == 0:
        empty_df = pd.DataFrame(
            {
                "point_id": pd.Series([], dtype="int"),
                "geometry": pd.Series([], dtype="object"),
                "access": pd.Series([], dtype="timedelta64[ns]"),
                "revisit": pd.Series([], dtype="timedetla64[ns]"),
                "num_obs": pd.Series([], dtype="int"),
            }
        )
        return gpd.GeoDataFrame(empty_df, geometry=empty_df.geometry, crs="EPSG:4326")
    gdf["access"] = gdf["access"] / timedelta(seconds=1)
    gdf["revisit"] = gdf["revisit"] / timedelta(seconds=1)
    gdf["num_obs"] = 1
    df = gdf.groupby("point_id").agg(
        {
            "point_id": "first",
            "geometry": "first",
            "access": "mean",
            "revisit": "mean",
            "num_obs": "sum",
        }
    )
    df["access"] = df["access"].apply(lambda t: timedelta(seconds=t))
    df["revisit"] = df["revisit"].apply(
        lambda t: pd.NaT if pd.isna(t) else timedelta(seconds=t)
    )
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")
