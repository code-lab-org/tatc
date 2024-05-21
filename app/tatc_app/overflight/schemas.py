# -*- coding: utf-8 -*-
"""
Schema specifications for overflight analysis endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""


from datetime import datetime
from typing import List, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator
from pydantic.alias_generators import to_camel
from tatc.schemas import (
    Instrument,
    Point,
    Satellite,
    TrainConstellation,
    WalkerConstellation,
)

from ..generation.schemas import PointGenerator


class OverflightAnalysisRequest(BaseModel):
    """
    User request to collect overflights.
    """

    model_config = ConfigDict(
        from_attributes=True, populate_by_name=True, alias_generator=to_camel
    )
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
        examples=[datetime.fromisoformat("2021-01-01T00:00:00+00:00")],
    )
    end: datetime = Field(
        ...,
        description="End date time of the observation period.",
        examples=[datetime.fromisoformat("2021-01-02T00:00:00+00:00")],
    )
    omit_solar: bool = Field(
        True,
        description="`True`, if solar angles can be omitted to improve computation speed.",
    )

    @model_validator(mode="after")
    def valid_instrument_index(self) -> "OverflightAnalysisRequest":
        """
        Validates the instrument index is in bounds.
        """
        # pylint: disable=E1101
        if (
            self.satellite is not None
            and self.instrument is not None
            and isinstance(self.instrument, int)
            and self.instrument >= len(self.satellite.instruments)
        ):
            raise ValueError("Instrument index out of bounds.")
        return self
