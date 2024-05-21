# -*- coding: utf-8 -*-
"""
Schema specifications for tracking analysis endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""


from datetime import datetime, timedelta
from typing import List, Optional, Union

from geojson_pydantic import Polygon
from pydantic import BaseModel, ConfigDict, Field, model_validator
from pydantic.alias_generators import to_camel
from shapely.geometry import box, mapping
from tatc.schemas.satellite import Satellite, TrainConstellation, WalkerConstellation
from tatc.schemas.instrument import Instrument

from ..generation.schemas import TimeGenerator


class OrbitTrackAnalysisRequest(BaseModel):
    """
    User request to perform orbit track analysis.
    """

    model_config = ConfigDict(
        from_attributes=True, populate_by_name=True, alias_generator=to_camel
    )
    satellite: Union[Satellite, TrainConstellation, WalkerConstellation] = Field(
        ..., description="Satellite from which to observe."
    )
    instrument: Optional[Union[int, Instrument]] = Field(
        None, description="Instrument (or index thereof) performing observation."
    )
    times: Union[List[datetime], TimeGenerator] = Field(
        ...,
        description="Start date time of the observation period.",
        examples=[
            TimeGenerator(
                start=datetime.fromisoformat("2021-01-01T00:00:00+00:00"),
                end=datetime.fromisoformat("2021-01-01T01:00:00+00:00"),
                delta=timedelta(seconds=30),
            )
        ],
    )
    elevation: float = Field(
        0,
        description="Elevation (meters) above datum in the WGS 84 coordinate system.",
    )
    # TODO add validator for mask
    mask: Optional[Polygon] = Field(
        None,
        description="Mask to limit the extent of generated points.",
        examples=[Polygon(**mapping(box(-180, -90, 180, 90)))],
    )

    @model_validator(mode="after")
    def valid_instrument_index(self) -> "OrbitTrackAnalysisRequest":
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


class GroundTrackAnalysisRequest(OrbitTrackAnalysisRequest):
    """
    User request to perform ground track analysis.
    """

    # TODO add validator for crs
    crs: str = Field(
        "EPSG:4087",
        description="Coordinate reference system for which to compute ground track. (Note: `utm` uses Universal Transverse Mercator (UTM) zones).",
    )
