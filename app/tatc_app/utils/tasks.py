# -*- coding: utf-8 -*-
"""
Task specifications for general utilities.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from typing import Union
from itertools import chain

from geojson_pydantic import FeatureCollection

from ..worker import app


@app.task
def merge_feature_collections_task(collections: Union[list, str]) -> str:
    """
    Task to merge a list of feature collections into a single feature collection.

    Args:
        collections (Union[list,str]): GeoJSON serialized list of feature collections.

    Results:
        str: GeoJSON serialized feature collection.
    """
    if isinstance(collections, list):
        return FeatureCollection(
            type="FeatureCollection",
            features=list(
                chain(
                    *list(
                        FeatureCollection.model_validate_json(collection).features
                        for collection in collections
                    )
                )
            ),
        ).model_dump_json()
    return collections
