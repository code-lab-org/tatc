from fastapi_utils.api_model import APIModel
from geojson_pydantic import Polygon
# -*- coding: utf-8 -*-
"""
Schema specifications for generation endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

from enum import Enum
from tatc.schemas import Point as TatcPoint
from tatc.schemas import GroundStation as TatcGroundStation
from typing import Optional, Union
from pydantic import Field, validator, constr
from shapely.geometry import box, mapping
from datetime import datetime, timedelta


class PointGeneratorMethod(str, Enum):
    """
    Method for point generation.
    """

    fibonacci_lattice = "fibonacci_lattice"
    cubed_square = "cubed_square"


class KnownShape(str, Enum):
    """
    Known shapes to use as masks.
    """

    conus = "conus"


class PointGenerator(APIModel):
    """
    Specification to generate points.
    """

    method: PointGeneratorMethod = Field(
        PointGeneratorMethod.cubed_square, description="Point generation method."
    )
    distance: float = Field(
        1e6, description="Characteristic distance (meters) between points."
    )
    elevation: float = Field(
        0,
        description="Elevation (meters) above datum in the WGS 84 coordinate system.",
    )
    mask: Optional[
        Union[Polygon, KnownShape, constr(min_length=3, max_length=3)]
    ] = Field(
        None,
        description="Mask to limit the extent of generated points. Allows ISO 3166-1 alpha-3 country codes.",
        example=Polygon.parse_obj(mapping(box(-180, -90, 180, 90))),
    )


class Point(TatcPoint):
    """
    Point.
    """

    pass


class GroundStation(TatcGroundStation):
    """
    Ground station.
    """

    pass


class Cell(Polygon):
    """
    Cell represented as a polygon.
    """

    pass


class CellGeneratorMethod(str, Enum):
    """
    Method for cell generation.
    """

    cubed_square = "cubed_square"


class CellStrips(str, Enum):
    """
    Options for cells in strips.
    """

    latitude = "lat"
    longitude = "lon"


class CellGenerator(APIModel):
    """
    Specification to generate cells.
    """

    method: CellGeneratorMethod = Field(
        CellGeneratorMethod.cubed_square, description="Cell generation method."
    )
    distance: float = Field(
        1e6, description="Characteristic distance (meters) between cell centroids."
    )
    elevation: float = Field(
        0,
        description="Elevation (meters) above datum in the WGS 84 coordinate system.",
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
    """
    Specification to generate time series.
    """

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
