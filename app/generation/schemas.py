from fastapi_utils.api_model import APIModel
from geojson_pydantic import Polygon
from enum import Enum
from tatc.schemas.point import Point as TatcPoint
from tatc.schemas.point import GroundStation as TatcGroundStation
from typing import Optional, Union
from pydantic import Field, validator, constr
from shapely.geometry import box, mapping
from datetime import datetime, timedelta


class PointGeneratorMethod(str, Enum):
    fibonacci_lattice = "fibonacci_lattice"
    cubed_square = "cubed_square"


class KnownShape(str, Enum):
    conus = "conus"


class PointGenerator(APIModel):
    method: PointGeneratorMethod = Field(
        PointGeneratorMethod.cubed_square, description="Point generation method."
    )
    distance: float = Field(
        1e6, description="Characteristic distance (meters) between points."
    )
    mask: Optional[
        Union[Polygon, KnownShape, constr(min_length=3, max_length=3)]
    ] = Field(
        None,
        description="Mask to limit the extent of generated points. Allows ISO 3166-1 alpha-3 country codes.",
        example=Polygon.parse_obj(mapping(box(-180, -90, 180, 90))),
    )


class Point(TatcPoint):
    pass


class GroundStation(TatcGroundStation):
    pass


class Cell(Polygon):
    pass


class CellGeneratorMethod(str, Enum):
    cubed_square = "cubed_square"


class CellStrips(str, Enum):
    latitude = "lat"
    longitude = "lon"


class CellGenerator(APIModel):
    method: CellGeneratorMethod = Field(
        CellGeneratorMethod.cubed_square, description="Cell generation method."
    )
    distance: float = Field(
        1e6, description="Characteristic distance (meters) between cell centroids."
    )
    mask: Optional[
        Union[Polygon, KnownShape, constr(min_length=3, max_length=3)]
    ] = Field(
        None,
        description="Mask to limit the extent of generated cells. Allows ISO 3166-1 alpha-3 country codes.",
        example=Polygon.parse_obj(mapping(box(-180, -90, 180, 90))),
    )
    strips: Optional[CellStrips] = Field(
        None, description="Option for latitude or longitude strips for cells."
    )


class TimeGenerator(APIModel):
    start: datetime = Field(
        ...,
        description="Start date time of the observation period.",
        example=datetime.fromisoformat("2021-01-01T00:00:00+00:00"),
    )
    end: datetime = Field(
        ...,
        description="End date time of the observation period.",
        example=datetime.fromisoformat("2021-01-01T01:00:00+00:00"),
    )
    delta: timedelta = Field(
        ..., description="Time step duration.", example=timedelta(seconds=30)
    )

    @validator("end")
    def validate_start_end(cls, v, values, **kwargs):
        if "start" in values and v < values["start"]:
            raise ValueError("Invalid end: must precede start.")
        return v

    @validator("delta")
    def validate_delta(cls, v, values, **kwargs):
        if (
            "start" in values
            and "end" in values
            and v > values["end"] - values["start"]
        ):
            raise ValueError("Invalid delta: must be smaller than (end-start).")
        return v
