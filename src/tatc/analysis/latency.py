"""
Methods to perform latency analysis.

@author: Isaac Feldman
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime

import geopandas as gpd
import pandas as pd
from shapely import geometry as geo

from ..constants import EARTH_MEAN_RADIUS
from ..schemas import GroundStation, Satellite
from ..utils.orbital import compute_apoapsis_radius
from .coverage import _get_visible_interval_series


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
    stations: GroundStation | list[GroundStation],
    satellite: Satellite,
    start: datetime,
    end: datetime,
) -> gpd.GeoDataFrame:
    """
    Collect satellite downlink opportunities to ground station(s) of interest.

    Args:
        stations (GroundStation | list[GroundStation]): The ground stations.
        satellite (Satellite): The observing satellite.
        start (datetime.datetime): Start of analysis period.
        end (datetime.datetime): End of analysis period.

    Returns:
        geopandas.GeoDataFrame: The data frame of collected downlink results.
    """
    # use the orbit's apogee altitude as a conservative upper bound
    max_altitude = (
        compute_apoapsis_radius(
            satellite.orbit.get_semimajor_axis(), satellite.orbit.get_eccentricity()
        )
        - EARTH_MEAN_RADIUS
    )
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
        for station in ([stations] if isinstance(stations, GroundStation) else stations)
        for period in _get_visible_interval_series(
            station,
            satellite,
            station.min_elevation_angle,
            max_altitude,
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
        geopandas.GeoDataFrame: The data frame of collected latency results, sorted by the 'observed'
        column in ascending order. It includes the following columns:
        - 'point_id' (int64): Identifier for the observation point.
        - 'geometry' (geometry): Geometry representing the observation point.
        - 'satellite' (object): Name or identifier of the satellite.
        - 'instrument' (object): Name or identifier of the instrument.
        - 'sat_alt' (float64): Altitude of the satellite at the time of observation.
        - 'sat_az' (float64): Azimuth of the satellite at the time of observation.
        - 'station' (object): Name or identifier of the ground station for downlink.
        - 'downlinked' (datetime64[ns, UTC]): Timestamp when the observation data was downlinked.
        - 'latency' (timedelta64[ns]): Latency between observation and downlink.
        - 'observed' (datetime64[ns, UTC]): Timestamp when the observation was made.
    """
    if observations.empty or downlinks.empty:
        return _get_empty_latency_frame()

    # merge observations with downlinks to find matching satellite downlinks
    obs = pd.merge_asof(
        observations.sort_values(by="end"),
        downlinks.sort_values(by="start"),
        by="satellite",
        left_on="end",
        right_on="start",
        direction="forward",
    )

    # compute latency
    obs["latency"] = obs["epoch_y"] - obs["epoch_x"]

    # rename and select relevant columns. Only "epoch" and "geometry" exist
    # in both `observations` and `downlinks`, so merge_asof only suffixes
    # those two with "_x"/"_y"; "station", "sat_alt", and "sat_az" exist in
    # just one of the two frames each and so are never suffixed at all.
    obs.rename(
        columns={
            "epoch_y": "downlinked",
            "epoch_x": "observed",
            "geometry_x": "geometry",
        },
        inplace=True,
    )

    # reorder columns
    obs = obs[
        [
            "point_id",
            "geometry",
            "satellite",
            "instrument",
            "sat_alt",
            "sat_az",
            "station",
            "downlinked",
            "latency",
            "observed",
        ]
    ].copy()

    # handle rows without matching downlinks (if any)
    no_downlink_rows = obs["downlinked"].isna()
    if no_downlink_rows.any():
        obs.loc[no_downlink_rows, ["station", "downlinked", "latency"]] = [
            None,
            pd.NaT,
            pd.NaT,
        ]

    # ensure result_df is a GeoDataFrame with geometry set
    obs = gpd.GeoDataFrame(obs, geometry="geometry")

    # set CRS if observations is a GeoDataFrame and has a defined CRS
    if isinstance(observations, gpd.GeoDataFrame) and observations.crs:
        obs.set_crs(observations.crs)

    # sort observations by observed time
    obs.sort_values(by="observed", inplace=True)

    obs.reset_index(drop=True, inplace=True)
    return obs


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
    Reduce observation latencies: for each unique point_id in
    `latency_observations`, computes the mean latency and the total number
    of samples (observation/downlink pairs). An observation with no
    matching downlink has an undefined (NaT) latency (see
    `compute_latencies`), which is excluded from the mean rather than
    counted as zero, while still counting toward `samples`.

    Args:
        latency_observations (geopandas.GeoDataFrame): The latency observations.

    Returns:
        geopandas.GeoDataFrame: The data frame with reduced latencies.
    """
    if latency_observations.empty:
        return _get_empty_reduce_frame()
    # operate on a copy of the dataframe
    gdf = latency_observations.copy()
    # convert latency to a numeric value before aggregation
    gdf["latency"] = gdf["latency"].dt.total_seconds()
    # assign each record to one observation
    gdf["samples"] = 1
    # perform the aggregation operation
    gdf = gdf.dissolve(
        "point_id",
        aggfunc={
            "latency": "mean",
            "samples": "sum",
        },
    ).reset_index()
    # convert latency from numeric values after aggregation
    gdf["latency"] = pd.to_timedelta(gdf["latency"], unit="s")
    return gdf


def grid_latencies(
    reduced_latencies: gpd.GeoDataFrame, cells: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Grid reduced latencies to cells: for every cell, sums the number of
    samples across every point it contains, and combines those points'
    latency into a single sample-weighted arithmetic mean per cell. Latency
    (a per-observation duration, like `access` in `grid_observations`, not
    a time-between-events/rate quantity like `revisit`) does not need a
    harmonic-mean treatment.

    Args:
        reduced_latencies (geopandas.GeoDataFrame): The reduced latencies.
        cells (geopandas.GeoDataFrame): The cell specification.

    Returns:
        geopandas.GeoDataFrame: The data frame with gridded latencies.
    """
    if reduced_latencies.empty:
        gdf = cells.copy()
        gdf["samples"] = 0
        gdf["latency"] = None
        return gdf
    # operate on a copy of the data frame
    gdf = reduced_latencies.copy()
    # convert latency to numeric values before aggregation
    gdf["latency"] = gdf["latency"].dt.total_seconds()
    # pre-multiply so the sample-weighted mean below reduces to a plain sum:
    # groupby().agg() with a dict of {column: function} only ever hands a
    # custom callable its own column's Series, never a sibling column like
    # "samples" needed to compute a weighted statistic within the callable
    gdf["latency_x_samples"] = gdf["latency"] * gdf["samples"]
    gdf = (
        cells.sjoin(gdf, how="inner", predicate="contains")
        .dissolve(
            by="cell_id",
            aggfunc={
                "samples": "sum",
                "latency_x_samples": "sum",
            },
        )
        .reset_index()
    )
    # finish the weighted mean
    gdf["latency"] = gdf["latency_x_samples"] / gdf["samples"]
    gdf = gdf.drop(columns=["latency_x_samples"])
    # convert latency from numeric values after aggregation
    gdf["latency"] = pd.to_timedelta(gdf["latency"], unit="s")
    return gdf
