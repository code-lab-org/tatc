from geojson_pydantic import FeatureCollection
import json
# -*- coding: utf-8 -*-
"""
Task specifications for general utilities.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from itertools import chain

from ..worker import app


@app.task
def merge_feature_collections_task(collections: str) -> str:
    """
    Task to merge a list of feature collections into a single feature collection.

    Args:
        collections (str): GeoJSON serialized list of feature collections.

    Results:
        str: GeoJSON serialized feature collection.
    """
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
