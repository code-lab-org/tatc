"""
Methods to perform latency analysis.

@author: Isaac Feldman
@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from __future__ import annotations

from datetime import datetime
from typing import Literal

import geopandas as gpd
import pandas as pd
from shapely import geometry as geo

from ..constants import EARTH_MEAN_RADIUS
from ..schemas import GroundStation, Satellite
from ..utils.orbital import compute_apoapsis_radius
from .coverage import _get_visible_interval_series
from .validation import _check_satellite


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
    _check_satellite(satellite)
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
    Gets an empty data frame for latency results.

    Returns:
        geopandas.GeoDataFrame: Empty data frame.
    """
    columns = {
        "point_id": pd.Series([], dtype="int"),
        "geometry": pd.Series([], dtype="object"),
        "satellite": pd.Series([], dtype="str"),
        "instrument": pd.Series([], dtype="str"),
        "sat_alt": pd.Series([], dtype="float"),
        "sat_az": pd.Series([], dtype="float"),
        "station": pd.Series([], dtype="str"),
        "downlinked": pd.Series([], dtype="datetime64[ns, utc]"),
        "latency": pd.Series([], dtype="timedelta64[ns]"),
        "observed": pd.Series([], dtype="datetime64[ns, utc]"),
    }
    return gpd.GeoDataFrame(columns, crs="EPSG:4326")


def compute_latencies(
    observations: gpd.GeoDataFrame,
    downlinks: gpd.GeoDataFrame,
    during_contact: Literal["next", "end", "immediate"] = "end",
) -> gpd.GeoDataFrame:
    """
    Collect latencies between an observation and the first downlink opportunity.

    An observation that ends before a downlink starts is downlinked at that
    downlink's epoch (midpoint). The `during_contact` option sets how an
    observation that ends while a downlink is in progress is downlinked:
    `"end"` (default) downlinks it at the end of the downlink in progress
    (stored data is played back after the data recorded before the contact);
    `"next"` waits for the next downlink to start (stored data is only
    played back from the start of a contact); `"immediate"` downlinks it as it is
    observed (real-time downlink), at the later of the observation epoch and
    the start of the downlink in progress. Each observation is assigned the
    earliest of these downlink times.

    Args:
        observations (geopandas.GeoDataFrame): The data frame of observations to downlink.
        downlinks (geopandas.GeoDataFrame): The data frame of downlink opportunities.
        during_contact (str): Downlink of observations that end during a
            downlink opportunity: `"end"` (default), `"next"`, or `"immediate"`.

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
    if during_contact not in ("next", "end", "immediate"):
        raise ValueError(
            f"during_contact must be 'next', 'end', or 'immediate', not {during_contact!r}"
        )
    if observations.empty or downlinks.empty:
        return _get_empty_latency_frame()

    obs = observations.sort_values(by="end").reset_index(drop=True)
    contacts = downlinks[["satellite", "station", "start", "epoch", "end"]].rename(
        columns={
            "start": "downlink_start",
            "epoch": "downlinked",
            "end": "downlink_end",
        }
    )
    # pair each observation with the first downlink that starts after it ends
    obs = pd.merge_asof(
        obs,
        contacts.sort_values(by="downlink_start"),
        by="satellite",
        left_on="end",
        right_on="downlink_start",
        direction="forward",
    )
    if during_contact != "next":
        # find the first downlink that ends after the observation ends, which
        # is in progress if it started before the observation ends
        current = pd.merge_asof(
            obs[["satellite", "end"]],
            contacts.sort_values(by="downlink_end"),
            by="satellite",
            left_on="end",
            right_on="downlink_end",
            direction="forward",
        )
        in_progress = current["downlink_start"] <= current["end"]
        if during_contact == "end":
            downlinked = current["downlink_end"]
        else:
            downlinked = current["downlink_start"].where(
                current["downlink_start"] > obs["epoch"], obs["epoch"]
            )
        earlier = in_progress & (
            obs["downlinked"].isna() | (downlinked < obs["downlinked"])
        )
        obs.loc[earlier, "station"] = current.loc[earlier, "station"]
        obs.loc[earlier, "downlinked"] = downlinked[earlier]

    # compute latency
    obs["latency"] = obs["downlinked"] - obs["epoch"]
    obs.rename(columns={"epoch": "observed"}, inplace=True)

    # select relevant columns
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
    ]

    # ensure the result is a GeoDataFrame with the observations' CRS
    obs = gpd.GeoDataFrame(
        obs,
        geometry="geometry",
        crs=(
            observations.crs
            if isinstance(observations, gpd.GeoDataFrame) and observations.crs
            else "EPSG:4326"
        ),
    )

    # sort observations by observed time
    obs = obs.sort_values(by="observed").reset_index(drop=True)
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
    latency into a single sample-weighted arithmetic mean per cell.

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
