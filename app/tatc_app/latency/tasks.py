from datetime import datetime, timedelta
import geopandas as gpd
# -*- coding: utf-8 -*-
"""
Task specifications for latency analysis endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pandas as pd

from tatc.schemas.point import Point, GroundStation
from tatc.schemas.satellite import Satellite
from tatc.analysis.latency import (
    collect_downlinks,
    compute_latencies,
    reduce_latencies,
    grid_latencies,
)
from tatc.analysis.coverage import collect_multi_observations

from .schemas import LatencyAnalysisResult
from ..worker import app


@app.task
def collect_downlinks_task(stations: list, satellite: str, start: str, end: str) -> str:
    """
    Task to collect downlinks for latency analysis.

    Args:
        stations (str): List of JSON serialized :class:`tatc.schemas.GroundStation` objects.
        satellite (str): JSON serialized :class:`tatc.schemas.Satellite` object.
        start (str): ISO 8601 serialized start time.
        end (str): ISO 8601 serialized end time.
    """
    # call analysis function, parsing the serialized arguments
    results = collect_downlinks(
        [GroundStation.parse_raw(station) for station in stations],
        Satellite.parse_raw(satellite),
        datetime.fromisoformat(start),
        datetime.fromisoformat(end),
    )
    # re-serialize constituent data
    results["start"] = results["start"].apply(lambda t: t.isoformat())
    results["epoch"] = results["epoch"].apply(lambda t: t.isoformat())
    results["end"] = results["end"].apply(lambda t: t.isoformat())
    return results.to_json(show_bbox=False, drop_id=True)


@app.task
def run_latency_analysis_task(
    downlinks: str, point: str, satellites: list, start: str, end: str
) -> str:
    """
    Task to run latency analysis.

    Args:
        downlinks (str): GeoJSON serialized `FeatureCollection` containing downlink results.
        point (str): JSON serialized :class:`tatc.schemas.Point` object.
        satellites (list): List of JSON serialized :class:`tatc.schemas.Satellite` objects.
        start (str): ISO 8601 serialized start time.
        end (str): ISO 8601 serialized end time.

    Returns:
        str: GeoJSON serialized `FeatureCollection` containing latency analysis results.
    """
    # parse downlinks
    downlinks = gpd.GeoDataFrame.from_features(json.loads(downlinks), crs="EPSG:4326")
    # coerce types for constituent data
    downlinks["start"] = downlinks["start"].astype("datetime64[ns, utc]")
    downlinks["epoch"] = downlinks["epoch"].astype("datetime64[ns, utc]")
    downlinks["end"] = downlinks["end"].astype("datetime64[ns, utc]")
    # call analysis function, parsing the serialized arguments
    observations = collect_multi_observations(
        Point.parse_raw(point),
        [Satellite.parse_raw(satellite) for satellite in satellites],
        datetime.fromisoformat(start),
        datetime.fromisoformat(end),
    )
    # compute latencies
    results = reduce_latencies(compute_latencies(observations, downlinks))
    # re-serialize constituent data
    results["latency"] = results["latency"].apply(lambda t: t.value)
    return results.to_json(show_bbox=False, drop_id=True)


@app.task
def grid_latency_analysis_task(latency_results: str, cells: str) -> str:
    """
    Task to grid latency analysis over specified cells.

    Args:
        latency_results (str): GeoJSON serialized `FeatureCollection` containing latency analysis results.
        cells (str): GeoJSON serialized `FeatureCollection` containing cell polygons.

    Returns:
        str: JSON serialized `LatencyAnalysisResult` containing latency analysis results.
    """
    # deserialize the coverage statistics
    gdf = gpd.GeoDataFrame.from_features(json.loads(latency_results), crs="EPSG:4326")
    # de-serialize constituent data
    gdf["latency"] = gdf["latency"].astype("timedelta64[ns]")
    # deserialize the cells
    grid_cells = gpd.GeoDataFrame.from_features(json.loads(cells), crs="EPSG:4326")
    # grid the results
    grid_data = grid_latencies(gdf, grid_cells)
    # re-serialize the constituent data
    grid_data["latency"] = grid_data["latency"].apply(lambda t: t.value)
    # return the results object in json format
    return LatencyAnalysisResult(
        points=json.loads(latency_results),
        cells=json.loads(grid_data.to_json(show_bbox=False, drop_id=True)),
    ).json()
