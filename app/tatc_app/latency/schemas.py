# -*- coding: utf-8 -*-
"""
Schema specifications for latency analysis endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from datetime import datetime
from typing import List, Union

from geojson_pydantic import FeatureCollection
from pydantic import BaseModel, ConfigDict, Field
from pydantic.alias_generators import to_camel

from tatc.schemas import GroundStation, Point, Satellite, TrainConstellation, WalkerConstellation
from ..generation.schemas import (
    PointGenerator,
    Cell,
    CellGenerator,
)


class LatencyAnalysisRequest(BaseModel):
    """
    User request to perform latency analysis.
    """

    model_config = ConfigDict(
        from_attributes=True, populate_by_name=True, alias_generator=to_camel
    )
    points: Union[List[Point], PointGenerator] = Field(
        ..., description="Points from which to collect observations."
    )
    cells: Union[List[Cell], CellGenerator] = Field(
        ..., description="Cells from which to aggregate observations."
    )
    stations: List[GroundStation] = Field(
        ..., description="Stations to which observations can be downlinked"
    )
    satellites: List[Union[Satellite, TrainConstellation, WalkerConstellation]] = Field(
        ..., description="Satellite from which to observe."
    )
    start: datetime = Field(
        ...,
        description="Start date time of the observation period.",
        examples=[datetime.fromisoformat("2021-01-01T00:00:00+00:00")],
    )
    end: datetime = Field(
        ...,
        description="End date time of the observation period.",
        examples=[datetime.fromisoformat("2021-01-02T00:00:00+00:00")],
    )


class LatencyAnalysisResult(BaseModel):
    """
    Response to provide latency analysis results.
    """

    model_config = ConfigDict(
        from_attributes=True, populate_by_name=True, alias_generator=to_camel
    )
    points: FeatureCollection = Field(..., description="Point-level latency analysis.")
    cells: FeatureCollection = Field(..., description="Cell-level latency analysis.")
