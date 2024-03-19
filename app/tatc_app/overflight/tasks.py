from datetime import datetime, timedelta
import geopandas as gpd
# -*- coding: utf-8 -*-
"""
Task specifications for overflight analysis endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
import pandas as pd
from tatc.schemas.instrument import Instrument
from tatc.schemas.point import Point
from tatc.schemas.satellite import Satellite
from tatc.analysis.coverage import collect_observations, aggregate_observations

from ..worker import app


@app.task
def collect_observations_task(
    point: str, satellite: str, instrument: str, start: str, end: str, omit_solar: bool
) -> str:
    """
    Task to collect observations.

    Args:
        point (str): JSON serialized :class:`tatc.schemas.Point` object.
        satellites (str): JSON serialized :class:`tatc.schemas.Satellite` object.
        instrument (str): JSON serialized :class:`tatc.schemas.Instrument` object.
        start (str): ISO 8601 serialized start time.
        end (str): ISO 8601 serialized end time.
        omit_solar (bool): `True`, if solar angles can be omitted to improve computation speed.

    Returns:
        FeatureCollection: GeoJSON serialized observations.
    """
    # call analysis function
    results = collect_observations(
        Point.parse_raw(point),
        Satellite.parse_raw(satellite),
        Instrument.parse_raw(instrument),
        datetime.fromisoformat(start),
        datetime.fromisoformat(end),
        omit_solar,
    )
    # serialize constituent data
    results["start"] = results["start"].apply(lambda t: t.isoformat())
    results["epoch"] = results["epoch"].apply(lambda t: t.isoformat())
    results["end"] = results["end"].apply(lambda t: t.isoformat())
    return results.to_json(show_bbox=False, drop_id=True)


@app.task
def aggregate_observations_task(observations_results):
    """
    Task to aggregate observations.

    Args:
        observations_results (str): GeoJSON serialized observations.

    Returns:
        FeatureCollection: GeoJSON serialized aggregated observations.
    """
    gdf = gpd.GeoDataFrame.from_features(
        json.loads(observations_results), crs="EPSG:4326"
    )
    # de-serialize constituent data
    gdf["start"] = gdf["start"].astype("datetime64[ns, utc]")
    gdf["epoch"] = gdf["epoch"].astype("datetime64[ns, utc]")
    gdf["end"] = gdf["end"].astype("datetime64[ns, utc]")
    # call analysis function
    results = aggregate_observations(gdf)
    # re-serialize constituent data
    results["start"] = results["start"].apply(lambda t: t.isoformat())
    results["epoch"] = results["epoch"].apply(lambda t: t.isoformat())
    results["end"] = results["end"].apply(lambda t: t.isoformat())
    results["access"] = results["access"].apply(lambda t: t.value)
    results["revisit"] = results["revisit"].apply(lambda t: t.value)
    return results.to_json(show_bbox=False, drop_id=True)
