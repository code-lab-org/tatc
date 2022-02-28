from geojson_pydantic import FeatureCollection
import json
from itertools import chain

from ..worker import app


@app.task
def merge_feature_collections_task(collections):
    if isinstance(collections, list):
        return FeatureCollection(
            features=list(
                chain(
                    *list(
                        json.loads(collection).get("features")
                        for collection in collections
                    )
                )
            )
        ).json()
    return collections
