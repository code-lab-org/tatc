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
from ..schemas.satellite import Satellite, SpaceSystem
from ..schemas.instrument import Instrument, DutyCycleScheme

from ..utils import (
    compute_min_altitude,
    swath_width_to_field_of_regard,
    compute_max_access_time,
)
from ..constants import de421, timescale


def _get_visible_interval_series(
    point: Point,
    satellite: Satellite,
    instrument: Instrument,
    start: datetime,
    end: datetime,
) -> pd.Series:
    """
    Get the series of visible intervals based on observation angle constraints.
    """
    # build a topocentric point at the designated geodetic point
    topos = wgs84.latlon(point.latitude, point.longitude)
    # define starting and ending points
    t0 = timescale.from_datetime(start)
    t1 = timescale.from_datetime(end)
    # convert orbit to tle
    orbit = satellite.orbit.to_tle()
    # construct a satellite for propagation
    sat = EarthSatellite(orbit.tle[0], orbit.tle[1], satellite.name)
    # compute the initial satellite height (altitude)
    satellite_height = wgs84.geographic_position_of(sat.at(t0)).elevation.m
    # compute the minimum altitude angle required for observation
    min_altitude = compute_min_altitude(
        satellite_height,
        instrument.field_of_regard,
    )
    # compute the maximum access time to filter bad data
    max_access_time = timedelta(
        seconds=compute_max_access_time(satellite_height, min_altitude)
    )
    # find the set of observation events
    t, events = sat.find_events(topos, t0, t1, altitude_degrees=min_altitude)

    # build the observation periods
    obs_periods = []
    if len(events) > 0 and np.all(events == 1):
        # if all events are type 1 (culminate), create a period from start to end
        obs_periods += [pd.Interval(left=pd.Timestamp(start), right=pd.Timestamp(end))]
    elif len(events) > 0:
        # otherwise, match rise/set events
        rises = t[events == 0]
        sets = t[events == 2]
        if len(sets) > 0 and (
            len(rises) == 0 or sets[0].utc_datetime() < rises[0].utc_datetime()
        ):
            # if first event is a set, create a period from the start
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(start), right=pd.Timestamp(sets[0].utc_datetime())
                )
            ]
        # create an observation period to match with each rise event if
        # there is a following set event within twice the maximum access time
        obs_periods += [
            pd.Interval(
                left=pd.Timestamp(rise.utc_datetime()),
                right=pd.Timestamp(
                    sets[
                        np.logical_and(
                            rise.utc_datetime() < sets.utc_datetime(),
                            sets.utc_datetime()
                            < rise.utc_datetime() + 2 * max_access_time,
                        )
                    ][0].utc_datetime()
                ),
            )
            for rise in rises
            if np.any(
                np.logical_and(
                    rise.utc_datetime() < sets.utc_datetime(),
                    sets.utc_datetime() < rise.utc_datetime() + 2 * max_access_time,
                )
            )
        ]
        if len(rises) > 0 and (
            len(sets) == 0 or rises[-1].utc_datetime() > sets[-1].utc_datetime()
        ):
            # if last event is a rise, create a period to the end
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(rises[-1].utc_datetime()), right=pd.Timestamp(end)
                )
            ]
    return pd.Series(obs_periods, dtype="interval")


def _get_satellite_altaz_series(
    observations: gpd.GeoDataFrame, sat: EarthSatellite
) -> pd.Series:
    """
    Get a series with the satellite altitude/azimuth for each observation.
    """
    sat_altaz = observations.apply(
        lambda r: (sat - wgs84.latlon(r.geometry.y, r.geometry.x))
        .at(timescale.from_datetime(r.epoch))
        .altaz(),
        axis=1,
    )
    return pd.Series(sat_altaz, dtype="object")


def _get_satellite_sunlit_series(
    observations: gpd.GeoDataFrame, sat: EarthSatellite
) -> pd.Series:
    """
    Get a series with the satellite sunlit condition for each observation.
    """
    sat_sunlit = observations.apply(
        lambda r: sat.at(timescale.from_datetime(r.epoch)).is_sunlit(de421),
        axis=1,
    )
    return pd.Series(sat_sunlit, dtype="bool")


def _get_solar_altaz_series(observations: gpd.GeoDataFrame) -> pd.Series:
    """
    Get a series with the solar altitude/azimuth for each observation.
    """
    sun_altaz = observations.apply(
        lambda r: (de421["earth"] + wgs84.latlon(r.geometry.y, r.geometry.x))
        .at(timescale.from_datetime(r.epoch))
        .observe(de421["sun"])
        .apparent()
        .altaz(),
        axis=1,
    )
    return pd.Series(sun_altaz, dtype="object")


def _get_solar_time_series(observations: gpd.GeoDataFrame) -> pd.Series:
    """
    Get a series with the local solar time for each observation.
    """
    solar_time = observations.apply(
        lambda r: (de421["earth"] + wgs84.latlon(r.geometry.y, r.geometry.x))
        .at(timescale.from_datetime(r.epoch))
        .observe(de421["sun"])
        .apparent()
        .hadec()[0]
        .hours
        + 12,
        axis=1,
    )
    return pd.Series(solar_time, dtype="float")


def _get_access_series(observations: gpd.GeoDataFrame) -> pd.Series:
    """
    Get a series with access times for each observation.
    """
    # compute the access time for the observation (end - start)
    return observations["end"] - observations["start"]


def _get_revisit_series(observations: gpd.GeoDataFrame) -> pd.Series:
    """
    Get a series with revisit times for each observation.
    """
    # compute the revisit time for each observation (previous end - start)
    return observations["end"] - observations["start"].shift()


def _get_empty_coverage_frame(omit_solar: bool) -> gpd.GeoDataFrame:
    """
    Get an empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "sat_alt": pd.Series(dtype="float"),
        "sat_az": pd.Series(dtype="float"),
    }
    if not omit_solar:
        columns = {
            **columns,
            **{
                "sat_sunlit": pd.Series(dtype="bool"),
                "solar_alt": pd.Series(dtype="float"),
                "solar_az": pd.Series(dtype="float"),
                "solar_time": pd.Series(dtype="float"),
            },
        }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def collect_observations(
    point: Point,
    satellite: Satellite,
    instrument: Instrument,
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
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
    :return: An instance of :class:`geopandas.GeoDataFrame` containing all
        recorded reduce_observations
    :rtype::`geopandas.GeoDataFrame`
    """
    # build a topocentric point at the designated geodetic point
    topos = wgs84.latlon(point.latitude, point.longitude)
    # convert orbit to tle
    orbit = satellite.orbit.to_tle()
    # construct a satellite for propagation
    sat = EarthSatellite(orbit.tle[0], orbit.tle[1], satellite.name)
    records = [
        {
            "point_id": point.id,
            "geometry": geo.Point(point.longitude, point.latitude),
            "satellite": satellite.name,
            "instrument": instrument.name,
            "start": period.left,
            "end": period.right,
            "epoch": period.mid,
        }
        for period in _get_visible_interval_series(
            point, satellite, instrument, start, end
        )
        if (
            instrument.min_access_time <= period.right - period.left
            and instrument.is_valid_observation(
                sat, timescale.from_datetime(period.mid)
            )
        )
    ]

    # build the dataframe
    if len(records) > 0:
        gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
        # append satellite altitude/azimuth columns
        sat_altaz = _get_satellite_altaz_series(gdf, sat)
        gdf["sat_alt"] = sat_altaz.apply(lambda r: r[0].degrees)
        gdf["sat_az"] = sat_altaz.apply(lambda r: r[1].degrees)
        if not omit_solar:
            # append satellite sunlit column
            gdf["sat_sunlit"] = _get_satellite_sunlit_series(gdf, sat)
            # append solar altitude/azimuth columns
            sun_altaz = _get_solar_altaz_series(gdf)
            gdf["sun_alt"] = sun_altaz.apply(lambda r: r[0].degrees)
            gdf["sun_az"] = sun_altaz.apply(lambda r: r[1].degrees)
            # append local solar time column
            gdf["solar_time"] = _get_solar_time_series(gdf)
    else:
        gdf = _get_empty_coverage_frame(omit_solar)
    # compute access and revisit metrics
    gdf["access"] = _get_access_series(gdf)
    gdf["revisit"] = _get_revisit_series(gdf)
    return gdf


def collect_multi_observations(
    point: Point,
    satellites: List[SpaceSystem],
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
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
    :return: an instance of :class:`geopandas.GeoDataFrame` containing all
        recorded observations
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    gdfs = [
        collect_observations(point, satellite, instrument, start, end, omit_solar)
        for constellation in satellites
        for satellite in (constellation.generate_members())
        for instrument in satellite.instruments
    ]
    # merge the observations into one data frame
    df = pd.concat(gdfs, ignore_index=True)
    # sort the values by start datetime
    df = df.sort_values("start")
    return gpd.GeoDataFrame(df, geometry=df.geometry, crs="EPSG:4326")


def _get_empty_aggregate_frame() -> gpd.GeoDataFrame:
    """
    Get an empty aggregate data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "access": pd.Series([], dtype="timedelta64[ns]"),
        "revisit": pd.Series([], dtype="timedelta64[ns]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def aggregate_observations(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Aggregate constellation observations for a geodetic point of interest.

    :param gdf: The individual observations.
    :type gdf: :class:`geopandas.GeoDataFrame`
    :return: An instance of :class:`geopandas.GeoDataFrame` The data frame
        with aggregated observations.
    :rtype: :class:`geopandas.GeoDataFrame`
    """
    if gdf.empty:
        return _get_empty_aggregate_frame()

    # sort the values by start datetime
    gdf = gdf.sort_values("start")
    # assign the observation group number based on overlapping start/end times
    gdf["obs"] = (gdf["start"] > gdf["end"].shift().cummax()).cumsum()
    # perform the aggregation to group overlapping observations
    gdf = gpd.GeoDataFrame(
        gdf.groupby("obs").agg(
            {
                "point_id": "first",
                "geometry": "first",
                "start": "min",
                "epoch": "first",
                "end": "max",
                "satellite": ", ".join,
                "instrument": ", ".join,
            }
        ),
        crs="EPSG:4326",
    )
    # compute access and revisit metrics
    gdf["access"] = _get_access_series(gdf)
    gdf["revisit"] = _get_revisit_series(gdf)
    return gdf


def _get_empty_reduce_frame() -> gpd.GeoDataFrame:
    """
    Get an empty reduce data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "access": pd.Series([], dtype="timedelta64[ns]"),
        "revisit": pd.Series([], dtype="timedelta64[ns]"),
        "num_obs": pd.Series([], dtype="int"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def reduce_observations(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reduce constellation observations for a geodetic point of interest.

    :param gdf: The aggregated observations
    :type gdf: :class:`geopandas.GeodataFrame`
    :return: An instance of :class:`geopandas.GeoDataFrame`: The data frame
        with reduced observations.
    :rtype: :class:`geopanadas.GeodataFrame`
    """
    if gdf.empty:
        return _get_empty_reduce_frame()

    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"] / timedelta(seconds=1)
    gdf["revisit"] = gdf["revisit"] / timedelta(seconds=1)
    # assign each record to one observation
    gdf["num_obs"] = 1
    # perform the aggregation operation
    gdf = gpd.GeoDataFrame(
        gdf.groupby("point_id").agg(
            {
                "point_id": "first",
                "geometry": "first",
                "access": "mean",
                "revisit": "mean",
                "num_obs": "sum",
            }
        ),
        crs="EPSG:4326",
    )
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = gdf["access"].apply(lambda t: timedelta(seconds=t))
    gdf["revisit"] = gdf["revisit"].apply(
        lambda t: pd.NaT if pd.isna(t) else timedelta(seconds=t)
    )
    return gdf
