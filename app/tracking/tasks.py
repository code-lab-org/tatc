from datetime import datetime
from shapely.geometry import shape
import time
import json
from geojson_pydantic import FeatureCollection

from tatc.analysis.track import collect_ground_track, collect_ground_track_swath
from tatc import schemas

from ..worker import app


@app.task
def collect_ground_track_task(satellite, instrument, times, mask):
    results = collect_ground_track(
        schemas.satellite.Satellite.parse_raw(satellite),
        schemas.instrument.Instrument.parse_raw(instrument),
        [datetime.fromisoformat(time) for time in times],
        shape(mask) if mask is not None else None,
    )
    # serialize Timestamp
    results["time"] = results["time"].apply(lambda t: t.isoformat())
    return results.to_json(show_bbox=False, drop_id=True)


@app.task
def collect_ground_track_swath_task(
    satellite, instrument, times, mask, fast, resolution
):
    results = collect_ground_track_swath(
        schemas.satellite.Satellite.parse_raw(satellite),
        schemas.instrument.Instrument.parse_raw(instrument),
        [datetime.fromisoformat(time) for time in times],
        shape(mask) if mask is not None else None,
        fast,
        resolution,
    )
    # serialize Timestamp
    results["time"] = results["time"].apply(lambda t: t.isoformat())
    return results.to_json(show_bbox=False, drop_id=True)


@app.task
def aggregate_ground_track_swath_task(results):
    return FeatureCollection(
        features=[
            feature
            for result in results
            for feature in json.loads(result).get("features")
        ]
    ).json()
