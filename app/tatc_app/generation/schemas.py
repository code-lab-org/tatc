# -*- coding: utf-8 -*-
"""
Schema specifications for generation endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""


from datetime import datetime, timedelta
from enum import Enum
from typing import Annotated, Optional, Union

from geojson_pydantic import Polygon
from pydantic import BaseModel, ConfigDict, Field, model_validator, StringConstraints
from pydantic.alias_generators import to_camel
from shapely.geometry import box, mapping


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


class PointGenerator(BaseModel):
    """
    Specification to generate points.
    """

    model_config = ConfigDict(
        from_attributes=True, populate_by_name=True, alias_generator=to_camel
    )
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
        Union[
            Polygon,
            KnownShape,
            Annotated[str, StringConstraints(min_length=3, max_length=3)],
        ]
    ] = Field(
        None,
        description="Mask to limit the extent of generated points. Allows ISO 3166-1 alpha-3 country codes.",
        examples=[Polygon(**mapping(box(-180, -90, 180, 90)))],
    )


class Cell(Polygon):
    """
    Cell represented as a polygon.
    """


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


class CellGenerator(BaseModel):
    """
    Specification to generate cells.
    """

    model_config = ConfigDict(
        from_attributes=True, populate_by_name=True, alias_generator=to_camel
    )
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
        Union[
            Polygon,
            KnownShape,
            Annotated[str, StringConstraints(min_length=3, max_length=3)],
        ]
    ] = Field(
        None,
        description="Mask to limit the extent of generated cells. Allows ISO 3166-1 alpha-3 country codes.",
        examples=[Polygon(**mapping(box(-180, -90, 180, 90)))],
    )
    strips: Optional[CellStrips] = Field(
        None, description="Option for latitude or longitude strips for cells."
    )


class TimeGenerator(BaseModel):
    """
    Specification to generate time series.
    """

    model_config = ConfigDict(
        from_attributes=True, populate_by_name=True, alias_generator=to_camel
    )
    start: datetime = Field(
        ...,
        description="Start date time of the observation period.",
        examples=[datetime.fromisoformat("2021-01-01T00:00:00+00:00")],
    )
    end: datetime = Field(
        ...,
        description="End date time of the observation period.",
        examples=[datetime.fromisoformat("2021-01-01T01:00:00+00:00")],
    )
    delta: timedelta = Field(
        ..., description="Time step duration.", examples=[timedelta(seconds=30)]
    )

    @model_validator(mode="after")
    def validate_start_end(self) -> "TimeGenerator":
        """
        Validates the end time precedes the start time.
        """
        if self.start is not None and self.end is not None and self.end < self.start:
            raise ValueError("Invalid end: start must precede end.")
        return self

    @model_validator(mode="after")
    def validate_delta(self) -> "TimeGenerator":
        """
        Validates the delta is not smaller than end - start.
        """
        if (
            self.start is not None
            and self.end is not None
            and self.delta is not None
            and self.delta > self.end - self.start
        ):
            raise ValueError("Invalid delta: must be smaller than (end-start).")
        return self
