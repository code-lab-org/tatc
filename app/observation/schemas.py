from datetime import datetime
from fastapi_utils.api_model import APIModel
from pydantic import Field, conlist, validator
from tatc.schemas.satellite import Satellite, TrainConstellation, WalkerConstellation
from tatc.schemas.instrument import Instrument
from typing import List, Optional, Union

from ..generation.schemas import Point, PointGenerator, CellGenerator


class CollectObservationsRequest(APIModel):
    points: Union[List[Point], PointGenerator] = Field(
        ..., description="Points from which to collect observations."
    )
    satellite: Union[Satellite, TrainConstellation, WalkerConstellation] = Field(
        ..., description="Satellite from which to observe."
    )
    instrument: Union[int, Instrument] = Field(
        0, description="Instrument (or index thereof) performing the observation."
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
    omit_solar: bool = Field(
        True,
        description="True, if solar angles can be omitted to improve computation speed.",
    )
    sample_distance: Optional[float] = Field(
        None,
        description="Sample distance to override instrument field of regard (optional, default: None).",
    )

    @validator("instrument")
    def valid_instrument_index(cls, v, values):
        if (
            isinstance(v, int)
            and "satellite" in values
            and v >= len(values["satellite"].instruments)
        ):
            raise ValueError("Instrument index out of bounds.")
        return v
