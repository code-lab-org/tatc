from fastapi import APIRouter
from geojson_pydantic import FeatureCollection
import json

from .schemas import (
    PointGenerator,
    PointGeneratorMethod,
    CellGenerator,
    CellGeneratorMethod,
)
from .tasks import (
    generate_cubed_sphere_points_task,
    generate_fibonacci_lattice_points_task,
    generate_cubed_sphere_cells_task,
)
from ..worker import app as celery_app

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
    if generator.method == PointGeneratorMethod.cubed_square:
        task = generate_cubed_sphere_points_task.delay(
            generator.distance, generator.elevation, generator.dict().get("mask", None)
        )
    elif generator.method == PointGeneratorMethod.fibonacci_lattice:
        task = generate_fibonacci_lattice_points_task.delay(
            generator.distance, generator.elevation, generator.dict().get("mask", None)
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
    if generator.method == CellGeneratorMethod.cubed_square:
        task = generate_cubed_sphere_cells_task.delay(
            generator.distance,
            generator.elevation,
            generator.dict().get("mask", None),
            generator.strips,
        )
    return json.loads(task.get())
