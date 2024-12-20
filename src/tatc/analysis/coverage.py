# -*- coding: utf-8 -*-
"""
Methods to perform coverage analysis.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import List, Union, Tuple
from datetime import datetime, timedelta

import pandas as pd
import numpy as np
import numpy.typing as npt
import geopandas as gpd
from shapely import geometry as geo
from skyfield.api import Angle, Distance, wgs84

from ..schemas.point import Point
from ..schemas.satellite import Satellite

from ..utils import (
    compute_min_elevation_angle,
    compute_max_access_time,
)
from ..constants import de421, timescale


def _get_visible_interval_series(
    point: Point,
    satellite: Satellite,
    min_elevation_angle: float,
    start: datetime,
    end: datetime,
) -> pd.Series:
    """
    Get the series of visible intervals based on altitude angle constraints.

    Args:
        point (Point): Point to observe.
        satellite (Satellite): Satellite doing the observation.
        min_elevation_angle (float): Minimum elevation angle (degrees) for valid observation.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.

    Returns:
        pandas.Series: Series of observation intervals.
    """
    # define starting and ending points
    t_0 = timescale.from_datetime(start)
    t_1 = timescale.from_datetime(end)
    # build skyfield objects
    topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
    sat = satellite.as_skyfield()
    # compute the initial satellite altitude
    satellite_altitude = wgs84.geographic_position_of(sat.at(t_0)).elevation.m
    # compute the maximum access time to filter bad data
    max_access_time = timedelta(
        seconds=compute_max_access_time(satellite_altitude, min_elevation_angle)
    )
    # find the set of observation events
    times, events = sat.find_events(
        topos, t_0, t_1, altitude_degrees=min_elevation_angle
    )

    # build the observation periods
    obs_periods = []
    if len(events) > 0 and np.all(events == 1):
        # if all events are type 1 (culminate), create a period from start to end
        obs_periods += [pd.Interval(left=pd.Timestamp(start), right=pd.Timestamp(end))]
    elif len(events) > 0:
        # otherwise, match rise/set events
        rises = times[events == 0]
        sets = times[events == 2]
        if (
            len(sets) > 0
            and (len(rises) == 0 or sets[0].utc_datetime() < rises[0].utc_datetime())
            and start < sets[0].utc_datetime()
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
        if (
            len(rises) > 0
            and (len(sets) == 0 or rises[-1].utc_datetime() > sets[-1].utc_datetime())
            and rises[-1].utc_datetime() < end
        ):
            # if last event is a rise, create a period to the end
            obs_periods += [
                pd.Interval(
                    left=pd.Timestamp(rises[-1].utc_datetime()), right=pd.Timestamp(end)
                )
            ]
    return pd.Series(obs_periods, dtype="interval")


def _get_satellite_altaz_series(
    point: Point,
    satellite: Satellite,
    times: List[datetime],
) -> Tuple[Angle, Angle, Distance]:
    """
    Get a series with the satellite altitude/azimuth for each observation.

    Args:
        point (skyfield.positionlib.GeographicPosition): Observed location.
        satellite (skyfield.api.EarthSatellite): Satellite doing the observation.
        times (List[datetime]): Times of observation.

    Returns:
        Tuple[skyfield.api.Angle, skyfield.api.Angle, skyfield.api.Distance]: Altitude-azimuth objects associated with observations.
    """
    sat = satellite.as_skyfield()
    topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
    ts = timescale.from_datetimes(times)
    return (sat - topos).at(ts).altaz()


def _get_satellite_sunlit_series(
    satellite: Satellite,
    times: List[datetime],
) -> npt.NDArray:
    """
    Get a series with the satellite sunlit condition for each observation.

    Args:
        satellite (skyfield.api.EarthSatellite): Satellite doing the observation.
        times (List[datetime]): Times of observation.

    Returns:
        numpy.typing.NDArray: Array of indicators whether the satellite is sunlit.
    """
    return satellite.get_orbit_track(times).is_sunlit(de421)


def _get_solar_altaz_series(
    point: Point,
    times: List[datetime],
) -> Tuple[Angle, Angle, Distance]:
    """
    Get a series with the solar altitude/azimuth for each observation.

    Args:
        point (skyfield.positionlib.GeographicPosition): Observed location.
        times (List[datetime]): Times of observation.

    Returns:
        Tuple[skyfield.api.Angle, skyfield.api.Angle, skyfield.api.Distance]: Altitude-azimuth objects associated with observations.
    """
    topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
    ts = timescale.from_datetimes(times)
    return (de421["earth"] + topos).at(ts).observe(de421["sun"]).apparent().altaz()


def _get_solar_time_series(
    point: Point, times: List[datetime]
) -> npt.NDArray:
    """
    Get a series with the local solar time for each observation.

    Args:
        observations (geopandas.GeoDataFrame): Data frame of observation records.

    Returns:
        numpy.typing.NDArray: Array of floats (local solar time in hours) associated with observations.
    """
    topos = wgs84.latlon(point.latitude, point.longitude, point.elevation)
    ts = timescale.from_datetimes(times)
    return (de421["earth"] + topos).at(ts).observe(de421["sun"]).apparent().hadec()[
        0
    ].hours + 12


def _get_empty_coverage_frame(omit_solar: bool) -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
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
    start: datetime,
    end: datetime,
    instrument_index: int = 0,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect single satellite observations of a geodetic point of interest.

    Args:
        point (Point): The ground point of interest.
        satellite (Satellite): The observing satellite.
        instrument (Instrument): The observing instrument.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        instrument_index (int): The index of the observing instrument in satellite.
        omit_solar (bool): `True`, to omit solar angles to improve
            computational efficiency; otherwise `False`.

    Returns:
        geopandas.GeoDataFrame: The data frame with recorded observations.
    """
    instrument = satellite.instruments[instrument_index]
    # compute the initial satellite altitude
    satellite_altitude = wgs84.geographic_position_of(
        satellite.get_orbit_track(start)
    ).elevation.m
    # compute the minimum altitude angle required for observation
    min_elevation_angle = compute_min_elevation_angle(
        satellite_altitude,
        instrument.field_of_regard,
    )
    records = [
        {
            "point_id": point.id,
            "geometry": geo.Point(point.longitude, point.latitude, point.elevation),
            "satellite": satellite.name,
            "instrument": instrument.name,
            "start": (
                period.left
                if not instrument.access_time_fixed
                else period.mid - instrument.min_access_time / 2
            ),
            "end": (
                period.right
                if not instrument.access_time_fixed
                else period.mid + instrument.min_access_time / 2
            ),
            "epoch": period.mid,
        }
        for period in _get_visible_interval_series(
            point, satellite, min_elevation_angle, start, end
        )
        if (
            instrument.min_access_time <= period.right - period.left
            and instrument.is_valid_observation(satellite.get_orbit_track(period.mid))
        )
    ]

    # build the dataframe
    if len(records) > 0:
        gdf = gpd.GeoDataFrame(records, crs="EPSG:4326")
        # append satellite altitude/azimuth columns
        sat_altaz = _get_satellite_altaz_series(point, satellite, gdf.epoch)
        gdf["sat_alt"] = sat_altaz[0].degrees
        gdf["sat_az"] = sat_altaz[1].degrees
        if not omit_solar:
            # append satellite sunlit column
            gdf["sat_sunlit"] = _get_satellite_sunlit_series(satellite, gdf.epoch)
            # append solar altitude/azimuth columns
            sun_altaz = _get_solar_altaz_series(point, gdf.epoch)
            gdf["solar_alt"] = sun_altaz[0].degrees
            gdf["solar_az"] = sun_altaz[1].degrees
            # append local solar time column
            gdf["solar_time"] = _get_solar_time_series(point, gdf.epoch)
    else:
        gdf = _get_empty_coverage_frame(omit_solar)
    return gdf


def collect_multi_observations(
    point: Point,
    satellites: Union[Satellite, List[Satellite]],
    start: datetime,
    end: datetime,
    omit_solar: bool = True,
) -> gpd.GeoDataFrame:
    """
    Collect multiple satellite observations of a geodetic point of interest.

    Args:
        point (Point): The ground point of interest.
        satellites (Satellite or List[Satellite]): The observing satellite(s).
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.
        omit_solar (bool): `True`, to omit solar angles to improve
            computational efficiency; otherwise `False`.

    Returns:
        geopandas.GeoDataFrame: The data frame with all recorded observations.
    """
    gdfs = [
        collect_observations(point, satellite, start, end, instrument_index, omit_solar)
        for constellation in (
            satellites if isinstance(satellites, list) else [satellites]
        )
        for satellite in (constellation.generate_members())
        for instrument_index in range(len(satellite.instruments))
    ]
    # concatenate into one data frame, sort by start time, and re-index
    return pd.concat(gdfs).sort_values("start").reset_index(drop=True)


def _get_empty_aggregate_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for aggregated coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "start": pd.Series([], dtype="datetime64[ns, utc]"),
        "epoch": pd.Series([], dtype="datetime64[ns, utc]"),
        "end": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def aggregate_observations(observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Aggregate constellation observations. Interleaves observations by multiple
    satellites to compute aggregate performance metrics including access
    (observation duration) and revisit (duration between observations).

    Args:
        observations (geopandas.GeoDataFrame): The collected observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with aggregated observations.
    """
    if observations.empty:
        return _get_empty_aggregate_frame()
    gdfs = []
    # split into constituent data frames based on point_id
    for _, gdf in observations.groupby("point_id"):
        # sort the values by start datetime
        gdf = gdf.sort_values("start")
        # assign the observation group number based on overlapping start/end times
        gdf["obs"] = (gdf["start"] > gdf["end"].shift().cummax()).cumsum()
        # perform the aggregation to group overlapping observations
        gdf = gdf.dissolve(
            "obs",
            aggfunc={
                "point_id": "first",
                "satellite": ", ".join,
                "instrument": ", ".join,
                "start": "min",
                "epoch": "mean",
                "end": "max",
            },
        )
        # compute access and revisit metrics
        gdf["access"] = gdf["end"] - gdf["start"]
        gdf["revisit"] = gdf["start"] - gdf["end"].shift()
        # append to the list of data frames
        gdfs.append(gdf)
    # return a concatenated data frame and re-index
    return pd.concat(gdfs).reset_index(drop=True)


def _get_empty_reduce_frame() -> gpd.GeoDataFrame:
    """
    Gets an empty data frame for reduced coverage analysis results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "access": pd.Series([], dtype="timedelta64[ns]"),
        "revisit": pd.Series([], dtype="timedelta64[ns]"),
        "samples": pd.Series([], dtype="int"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def reduce_observations(aggregated_observations: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Reduce constellation observations. Computes descriptive statistics for each
    geodetic point of interest contained in aggregated observations.

    Args:
        aggregated_observations (geopandas.GeoDataFrame): The aggregated observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced observations.
    """
    if aggregated_observations.empty:
        return _get_empty_reduce_frame()
    # operate on a copy of the data frame
    gdf = aggregated_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"] / timedelta(seconds=1)
    gdf["revisit"] = gdf["revisit"] / timedelta(seconds=1)
    # assign each record to one observation
    gdf["samples"] = 1
    # perform the aggregation operation
    gdf = gdf.dissolve(
        "point_id",
        aggfunc={
            "access": "mean",
            "revisit": "mean",
            "samples": "sum",
        },
    ).reset_index()
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = gdf["access"].apply(lambda t: timedelta(seconds=t))
    gdf["revisit"] = gdf["revisit"].apply(
        lambda t: pd.NaT if pd.isna(t) else timedelta(seconds=t)
    )
    return gdf


def grid_observations(
    reduced_observations: gpd.GeoDataFrame, cells: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Grid reduced observations to cells.

    Args:
        reduced_observations (geopandas.GeoDataFrame): The reduced observations.
        cells (geopandas.GeoDataFrame): The cell specification.

    Returns:
        geopandas.GeoDataFrame: The data frame with gridded observations.
    """
    if reduced_observations.empty:
        gdf = cells.copy()
        gdf["samples"] = 0
        gdf["access"] = None
        gdf["revisit"] = None
        return gdf
    # operate on a copy of the data frame
    gdf = reduced_observations.copy()
    # convert access and revisit to numeric values before aggregation
    gdf["access"] = gdf["access"] / timedelta(seconds=1)
    gdf["revisit"] = gdf["revisit"] / timedelta(seconds=1)
    gdf = (
        cells.sjoin(gdf, how="inner", predicate="contains")
        .dissolve(
            by="cell_id",
            aggfunc={
                "samples": "sum",
                "access": lambda r: np.average(r, weights=gdf.loc[r.index, "samples"]),
                "revisit": lambda r: np.average(r, weights=gdf.loc[r.index, "samples"]),
            },
        )
        .reset_index()
    )
    # convert access and revisit from numeric values after aggregation
    gdf["access"] = gdf["access"].apply(lambda t: timedelta(seconds=t))
    gdf["revisit"] = gdf["revisit"].apply(
        lambda t: pd.NaT if pd.isna(t) else timedelta(seconds=t)
    )
    return gdf
