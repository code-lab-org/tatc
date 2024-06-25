# -*- coding: utf-8 -*-
"""
Utility functions for generation endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json
from typing import Union, List

import geopandas as gpd
from tatc.schemas import Point

from .schemas import (
    PointGenerator,
    PointGeneratorMethod,
    Cell,
    CellGenerator,
    CellGeneratorMethod,
)
from .tasks import (
    generate_equally_spaced_points_task,
    generate_fibonacci_lattice_points_task,
    generate_equally_spaced_cells_task,
)


def generate_points(points: Union[List[Point], PointGenerator]) -> List[Point]:
    """
    Generate points from either a generator specification or an existing list.

    Args:
        points (Union[PointGenerator, List[Point]]): points to generate

    Returns:
        List[Point]: generated points
    """
    if isinstance(points, PointGenerator):
        if points.method == PointGeneratorMethod.equally_spaced:
            return gpd.GeoDataFrame.from_features(
                json.loads(
                    generate_equally_spaced_points_task.delay(
                        points.distance,
                        points.elevation,
                        points.model_dump().get("mask", None),
                    ).get()
                )
            ).apply(
                lambda r: Point(
                    id=r.point_id, latitude=r.geometry.y, longitude=r.geometry.x
                ),
                axis=1,
            )
        elif points.method == PointGeneratorMethod.fibonacci_lattice:
            return gpd.GeoDataFrame.from_features(
                json.loads(
                    generate_fibonacci_lattice_points_task.delay(
                        points.distance,
                        points.elevation,
                        points.model_dump().get("mask", None),
                    ).get()
                )
            ).apply(
                lambda r: Point(
                    id=r.point_id, latitude=r.geometry.y, longitude=r.geometry.x
                ),
                axis=1,
            )
    return points


def generate_cells(cells: Union[CellGenerator, List[Cell]]) -> List[Cell]:
    """
    Generate cells from either a generator specification or an existing list.

    Args:
        points (Union[CellGenerator, List[Cell]]): cells to generate

    Returns:
        List[Cell]: generated cells
    """
    if isinstance(cells, CellGenerator):
        if cells.method == CellGeneratorMethod.equally_spaced:
            return gpd.GeoDataFrame.from_features(
                json.loads(
                    generate_equally_spaced_cells_task.delay(
                        cells.distance,
                        cells.elevation,
                        cells.model_dump().get("mask", None),
                        cells.model_dump().get("strips", None),
                    ).get()
                )
            )
    return cells
