# -*- coding: utf-8 -*-
"""
Schema specifications for tracking analysis endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime, timedelta
from fastapi_utils.api_model import APIModel
from geojson_pydantic import Polygon
from pydantic import Field, validator
from shapely.geometry import box, mapping
from tatc.schemas.satellite import Satellite, TrainConstellation, WalkerConstellation
from tatc.schemas.instrument import Instrument
from typing import List, Optional, Union

from ..generation.schemas import TimeGenerator


class OrbitTrackAnalysisRequest(APIModel):
    """
    User request to perform orbit track analysis.
    """

    satellite: Union[Satellite, TrainConstellation, WalkerConstellation] = Field(
        ..., description="Satellite from which to observe."
    )
    instrument: Optional[Union[int, Instrument]] = Field(
        None, description="Instrument (or index thereof) performing observation."
    )
    times: Union[List[datetime], TimeGenerator] = Field(
        ...,
        description="Start date time of the observation period.",
        example=TimeGenerator(
            start=datetime.fromisoformat("2021-01-01T00:00:00+00:00"),
            end=datetime.fromisoformat("2021-01-01T01:00:00+00:00"),
            delta=timedelta(seconds=30),
        ),
    )
    elevation: float = Field(
        0,
        description="Elevation (meters) above datum in the WGS 84 coordinate system.",
    )
    mask: Optional[Polygon] = Field(
        None,
        description="Mask to limit the extent of generated points.",
        example=Polygon.parse_obj(mapping(box(-180, -90, 180, 90))),
    )

    # TODO add validator for mask


class GroundTrackAnalysisRequest(OrbitTrackAnalysisRequest):
    """
    User request to perform ground track analysis.
    """

    crs: str = Field(
        "EPSG:4087",
        description="Coordinate reference system for which to compute ground track. (Note: `utm` uses Universal Transverse Mercator (UTM) zones).",
    )

    # TODO add validator for crs
