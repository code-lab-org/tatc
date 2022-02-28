from datetime import datetime, timedelta
import geopandas as gpd
import json
import pandas as pd
from tatc.schemas.instrument import Instrument
from tatc.schemas.point import Point
from tatc.schemas.satellite import Satellite
from tatc.analysis.coverage import collect_observations, aggregate_observations

from ..worker import app


@app.task
def collect_observations_task(
    point, satellite, instrument, start, end, omit_solar, sample_distance
):
    # call analysis function
    results = collect_observations(
        Point.parse_raw(point),
        Satellite.parse_raw(satellite),
        Instrument.parse_raw(instrument),
        datetime.fromisoformat(start),
        datetime.fromisoformat(end),
        omit_solar,
        sample_distance,
    )
    # serialize constituent data
    results["start"] = results["start"].apply(lambda t: t.isoformat())
    results["epoch"] = results["epoch"].apply(lambda t: t.isoformat())
    results["end"] = results["end"].apply(lambda t: t.isoformat())
    results["access"] = results["access"].apply(lambda t: t / timedelta(seconds=1))
    results["revisit"] = results["revisit"].apply(lambda t: t / timedelta(seconds=1))
    return results.to_json(show_bbox=False, drop_id=True)


@app.task
def aggregate_observations_task(feature_collection):
    collection = json.loads(feature_collection)
    if len(collection.get("features")) == 0:
        return feature_collection
    gdf = gpd.GeoDataFrame.from_features(collection, crs="EPSG:4326")
    # de-serialize constituent data
    gdf["start"] = gdf["start"].apply(lambda t: datetime.fromisoformat(t))
    gdf["epoch"] = gdf["epoch"].apply(lambda t: datetime.fromisoformat(t))
    gdf["end"] = gdf["end"].apply(lambda t: datetime.fromisoformat(t))
    gdf["access"] = gdf["access"].apply(lambda t: timedelta(seconds=t))
    gdf["revisit"] = gdf["revisit"].apply(
        lambda t: pd.NaT if pd.isna(t) else timedelta(seconds=t)
    )
    # call analysis function
    results = aggregate_observations(gdf)
    # re-serialize constituent data
    results["start"] = results["start"].apply(lambda t: t.isoformat())
    results["epoch"] = results["epoch"].apply(lambda t: t.isoformat())
    results["end"] = results["end"].apply(lambda t: t.isoformat())
    results["access"] = results["access"].apply(lambda t: t / timedelta(seconds=1))
    results["revisit"] = results["revisit"].apply(lambda t: t / timedelta(seconds=1))
    return results.to_json(show_bbox=False, drop_id=True)
