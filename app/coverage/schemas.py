from fastapi_utils.api_model import APIModel
from geojson_pydantic import FeatureCollection
from typing import List, Optional, Union
from pydantic import Field
from datetime import datetime
from tatc.schemas.satellite import Satellite, TrainConstellation, WalkerConstellation

from ..generation.schemas import Point, PointGenerator, Cell, CellGenerator


class CollectCoverageStatisticsRequest(APIModel):
    satellites: List[Union[Satellite, TrainConstellation, WalkerConstellation]] = Field(
        ..., description="Satellites from which to observe."
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
    points: Union[List[Point], PointGenerator] = Field(
        ..., description="Points from which to collect observations."
    )
    cells: Union[List[Cell], CellGenerator] = Field(
        ..., description="Cells from which to aggregate observations."
    )
    sample_distance: Optional[float] = Field(
        None,
        description="Sample distance to override instrument field of regard (optional, default: None).",
    )


class CoverageStatisticsAnalysisResult(APIModel):
    points: FeatureCollection = Field(
        ..., description="Point-level coverage statistics."
    )
    cells: FeatureCollection = Field(..., description="Cell-level coverage statistics.")
