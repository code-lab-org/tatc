from datetime import datetime, timedelta
from fastapi_utils.api_model import APIModel
from geojson_pydantic import Polygon
from pydantic import Field, validator
from shapely.geometry import box, mapping
from tatc.schemas.satellite import Satellite, TrainConstellation, WalkerConstellation
from tatc.schemas.instrument import Instrument
from typing import List, Optional, Union

from ..generation.schemas import TimeGenerator


class GroundTrackRequest(APIModel):
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
    mask: Optional[Polygon] = Field(
        None,
        description="Mask to limit the extent of generated points.",
        example=Polygon.parse_obj(mapping(box(-180, -90, 180, 90))),
    )


class GroundTrackSwathRequest(GroundTrackRequest):
    fast: bool = Field(
        True,
        description="Toggles fast (True) or accurate (False) projection. \
            Fast uses EPSG:4087 (World Equidistant Cylindrical) \
            to project distances, accurate finds the appropriate Unified \
            Transverse Mercator (UTM) zone for each point.",
    )
    resolution: int = Field(
        4,
        ge=1,
        description="Resolution of the projected swath (1=square; 16=99.8% accurate).",
    )
