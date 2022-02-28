from datetime import datetime, timedelta
import geopandas as gpd
import json
import pandas as pd

from tatc.schemas.instrument import Instrument
from tatc.schemas.point import Point, GroundStation
from tatc.schemas.satellite import Satellite
from tatc.analysis.latency import (
    collect_downlinks,
    compute_latency,
    aggregate_downlinks,
    aggregate_latencies
)
from tatc.analysis.coverage import collect_observations

from ..worker import app

@app.task
def collect_downlinks_task(station, satellite, start, end):
    results = collect_downlinks(
        GroundStation.parse_raw(station),
        Satellite.parse_raw(satellite),
        datetime.fromisoformat(start),
        datetime.fromisoformat(end)
    )
    results['start'] = results['start'].apply(lambda t: t.isoformat())
    results['epoch'] = results['epoch'].apply(lambda t: t.isoformat())
    results['end'] = results['end'].apply(lambda t: t.isoformat())
    results['access'] = results['access'].apply(lambda t: t/timedelta(seconds=1))
    results['revisit'] = results['revisit'].apply(lambda t: t/timedelta(seconds=1))

    return results.to_json(
        show_bbox=False,
        drop_id=True
    )

@app.task
def aggregate_downlinks_task(downlinks):
    collections = [
        json.loads(dl)
        for dl in downlinks
    ]
    gdfs = [
        gpd.GeoDataFrame.from_features(
            collection,
            crs="EPSG:4326"
        )
        for collection in collections
        if len(collection.get("features", [])) > 0
    ]
    for gdf in gdfs:
        gdf['start'] = gdf['start'].apply(lambda t: datetime.fromisoformat(t))
        gdf['epoch'] = gdf['epoch'].apply(lambda t: datetime.fromisoformat(t))
        gdf['end'] = gdf['end'].apply(lambda t: datetime.fromisoformat(t))
        gdf['access'] = gdf['access'].apply(lambda t: timedelta(seconds=t))
        gdf['revisit'] = gdf['revisit'].apply(
            lambda t: pd.NaT if pd.isna(t) else timedelta(seconds=t)
        )
    results = aggregate_downlinks(gdfs)
    results['start'] = results['start'].apply(lambda t: t.isoformat())
    results['epoch'] = results['epoch'].apply(lambda t: t.isoformat())
    results['end'] = results['end'].apply(lambda t: t.isoformat())
    results['access'] = results['access'].apply(lambda t: t/timedelta(seconds=1))
    results['revisit'] = results['revisit'].apply(lambda t: t/timedelta(seconds=1))
    return results.to_json(
        show_bbox=False,
        drop_id=True
    )

@app.task
def compute_latency_task(downlinks, point, satellite, instrument, start, end, omit_solar):
    observations = collect_observations(
        Point.parse_raw(point),
        Satellite.parse_raw(satellite),
        Instrument.parse_raw(instrument),
        datetime.fromisoformat(start),
        datetime.fromisoformat(end),
        omit_solar)
    observations['epoch'] = observations['epoch'].apply(lambda t: t.isoformat())

    results = gpd.GeoDataFrame()
    downlinks = gpd.GeoDataFrame.from_features(
        json.loads(downlinks),
        crs="EPSG:4326"
    )

    for _, observation in observations.iterrows():
        results = results.append(compute_latency(
                                    observation,
                                    downlinks
                                    )
                                )

    results['latency'] = results['latency'].apply(lambda t: t/timedelta(seconds=1))

    return results.to_json(
        show_bbox=False,
        drop_id=True
    )

@app.task
def aggregate_latencies_task(latencies):
    collections = [
        json.loads(latency)
        for latency in latencies
    ]

    gdfs = [
        gpd.GeoDataFrame.from_features(
            collection,
            crs="EPSG:4326"
        )
        for collection in collections
        if len(collection.get("features", [])) > 0
    ]

    for gdf in gdfs:
        gdf['observed'] = gdf['observed'].apply(lambda t: datetime.fromisoformat(t))
        gdf['downlinked'] = gdf['downlinked'].apply(lambda t: datetime.fromisoformat(t))
        gdf['latency'] =gdf['latency'].apply(lambda t: timedelta(seconds=t))

    results = aggregate_latencies(gdfs)

    results['observed'] = results['observed'].apply(lambda t: datetime.isoformat(t))
    results['downlinked'] = results['downlinked'].apply(lambda t: datetime.isoformat(t))
    results['latency'] = results['latency'].apply(lambda t: t/timedelta(seconds=1))

    return results.to_json(
        show_bbox=False,
        drop_id=True
    )
