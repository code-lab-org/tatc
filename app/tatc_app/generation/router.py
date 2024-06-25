# -*- coding: utf-8 -*-
"""
Router specifications for generation endpoints.

@author: Paul T. Grogan <paul.grogan@asu.edu>
"""

import json

from fastapi import APIRouter
from geojson_pydantic import FeatureCollection

from .schemas import (
    PointGenerator,
    PointGeneratorMethod,
    CellGenerator,
    CellGeneratorMethod,
)
from .tasks import (
    generate_equally_spaced_points_task,
    generate_fibonacci_lattice_points_task,
    generate_equally_spaced_cells_task,
)

router = APIRouter()


@router.post(
    "/points",
    response_model=FeatureCollection,
    response_model_exclude_none=True,
    responses={"200": {"content": {"application/x-msgpack": {}}}},
)
async def generate_points(generator: PointGenerator):
    """
    Endpoint to generate a set of points.

    Args:
        generator (PointGenerator): specification for point generation

    Returns:
        FeatureCollection: generated points
    """
    if generator.method == PointGeneratorMethod.equally_spaced:
        task = generate_equally_spaced_points_task.delay(
            generator.distance,
            generator.elevation,
            generator.model_dump().get("mask", None),
        )
    elif generator.method == PointGeneratorMethod.fibonacci_lattice:
        task = generate_fibonacci_lattice_points_task.delay(
            generator.distance,
            generator.elevation,
            generator.model_dump().get("mask", None),
        )
    return json.loads(task.get())


@router.post(
    "/cells",
    response_model=FeatureCollection,
    response_model_exclude_none=True,
    responses={"200": {"content": {"application/x-msgpack": {}}}},
)
async def generate_cells(generator: CellGenerator):
    """
    Endpoint to generate a set of cells.

    Args:
        generator (CellGenerator): specification for cell generation

    Returns:
        FeatureCollection: generated cells
    """
    if generator.method == CellGeneratorMethod.equally_spaced:
        task = generate_equally_spaced_cells_task.delay(
            generator.distance,
            generator.elevation,
            generator.model_dump().get("mask", None),
            generator.strips,
        )
    return json.loads(task.get())
