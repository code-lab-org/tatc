from datetime import datetime, timedelta
import geopandas as gpd
from geojson_pydantic import FeatureCollection
import json
from itertools import chain
import pandas as pd
from tatc.schemas.instrument import Instrument
from tatc.schemas.point import Point
from tatc.schemas.satellite import Satellite
from tatc.analysis.coverage import (
    collect_multi_observations,
    reduce_observations,
    aggregate_observations,
)

from .schemas import CoverageStatisticsAnalysisResult
from ..worker import app


@app.task
def compute_coverage_statistics_task(
    point, satellites, start, end, omit_solar, sample_distance
):
    # call analysis function
    results = reduce_observations(
        aggregate_observations(
            collect_multi_observations(
                Point.parse_raw(point),
                [Satellite.parse_raw(satellite) for satellite in satellites],
                datetime.fromisoformat(start),
                datetime.fromisoformat(end),
                omit_solar,
                sample_distance,
            )
        )
    )
    # re-serialize constituent data
    results["access"] = results["access"].apply(lambda t: t / timedelta(seconds=1))
    results["revisit"] = results["revisit"].apply(
        lambda t: None if pd.isna(t) else t / timedelta(seconds=1)
    )
    return results.to_json(show_bbox=False, drop_id=True)


@app.task
def aggregate_coverage_statistics_task(feature_collection, cells):
    collection = json.loads(feature_collection)
    if len(collection.get("features")) == 0:
        return feature_collection
    gdf = gpd.GeoDataFrame.from_features(collection, crs="EPSG:4326")
    grid_cells = gpd.GeoDataFrame.from_features(json.loads(cells), crs="EPSG:4326")
    grid_data = (
        gpd.sjoin(grid_cells, gdf, how="inner", op="contains")
        .dissolve(by="cell_id", aggfunc="mean")
        .reset_index()
        .drop(columns=["index_right", "point_id"])
    )

    return CoverageStatisticsAnalysisResult(
        points=collection,
        cells=json.loads(grid_data.to_json(show_bbox=False, drop_id=True)),
    ).json()
