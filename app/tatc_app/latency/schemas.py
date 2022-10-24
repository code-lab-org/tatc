from datetime import datetime, timedelta
from fastapi_utils.api_model import APIModel
from geojson_pydantic import FeatureCollection
from typing import List, Optional, Union
from pydantic import Field, conlist, validator
import pandas as pd

from tatc.schemas.satellite import Satellite, TrainConstellation, WalkerConstellation
from tatc.schemas.instrument import Instrument
from ..generation.schemas import (
    Point,
    PointGenerator,
    GroundStation,
    Cell,
    CellGenerator,
)
from typing import List, Optional, Union


class LatencyAnalysisRequest(APIModel):
    """
    User request to perform latency analysis.
    """

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
        example=datetime.fromisoformat("2021-01-01T00:00:00+00:00"),
    )
    end: datetime = Field(
        ...,
        description="End date time of the observation period.",
        example=datetime.fromisoformat("2021-01-02T00:00:00+00:00"),
    )


class LatencyAnalysisResult(APIModel):
    """
    Response to provide latency analysis results.
    """

    points: FeatureCollection = Field(..., description="Point-level latency analysis.")
    cells: FeatureCollection = Field(..., description="Cell-level latency analysis.")
